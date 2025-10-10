import scanpy as sc
import pandas as pd
import numpy as np
import argparse
import os
import scipy.sparse
import scipy.io

def clean_barcode_suffix(barcode):
    """Remove '-' and following numbers from cell barcode."""
    if '-' in barcode:
        return barcode.split('-')[0]
    return barcode

def load_data(matrix_path, barcodes_path, features_path):
    import scipy.io

    # Load the matrix
    matrix = scipy.io.mmread(matrix_path).tocsr()

    # Load barcodes (cell IDs)
    barcodes_raw = pd.read_csv(barcodes_path, header=None)[0].tolist()
    barcodes = [clean_barcode_suffix(b) for b in barcodes_raw]

    # Load features (gene IDs and gene symbols)
    features = pd.read_csv(features_path, sep='\t', header=None)
    features.columns = ['gene_id', 'gene_symbol', 'feature_type']

    # Debugging output
    print(f"Matrix shape: {matrix.shape}")
    print(f"Features shape: {features.shape}")
    print(f"Barcodes count: {len(barcodes)}")

    # Clean up hidden characters
    features['gene_id'] = features['gene_id'].astype(str).str.strip()
    features['gene_symbol'] = features['gene_symbol'].astype(str).str.strip()

    # Check for duplicate gene IDs
    if not features['gene_id'].is_unique:
        print("⚠️ Found duplicate gene IDs. Removing duplicates.")
        features = features.drop_duplicates(subset='gene_id', keep='first')
        matrix = matrix[:features.shape[0], :]
    else:
        print("✅ Gene IDs are unique.")

    # Check for duplicate barcodes
    barcode_duplicates = pd.Series(barcodes)[pd.Series(barcodes).duplicated()]
    if not barcode_duplicates.empty:
        print(f"⚠️ Found {barcode_duplicates.nunique()} duplicate barcodes.")
        barcodes = pd.Series(barcodes).drop_duplicates().tolist()

    # ✅ FINAL FIX: Transpose the matrix
    if matrix.shape[0] == len(features) and matrix.shape[1] == len(barcodes):
        print("🔄 Transposing matrix to match AnnData format (cells × genes).")
        matrix = matrix.transpose()
    else:
        print("✅ Matrix orientation matches AnnData expectations.")

    # Final assignment
    adata = sc.AnnData(matrix)
    adata.var_names = pd.Index(features['gene_id'].values, name='gene_id')
    adata.obs_names = pd.Index(barcodes, name='cell_id')
    adata.var['gene_symbol'] = features['gene_symbol'].values

    return adata

def calculate_metrics(adata):
    """ Calculate basic QC metrics """
    sc.pp.calculate_qc_metrics(adata, percent_top=[20], log1p=False, inplace=True)

def mitochondrial_percentage(adata, mito_prefix='MT-'):
    """ Calculate mitochondrial gene expression percentage """
    adata.var['mt'] = adata.var['gene_symbol'].str.startswith(mito_prefix, na=False)
    adata.obs['pct_counts_mt'] = (
        np.sum(adata[:, adata.var['mt']].X, axis=1).A1 / np.sum(adata.X, axis=1).A1 * 100
    )

def compute_complexity(adata):
    """ Compute complexity score for each cell """
    adata.obs['complexity'] = adata.obs['n_genes_by_counts'] / adata.obs['total_counts']

def compute_mad_outliers(data, factor=5):
    """ Identify outliers using Median Absolute Deviation (MAD) """
    median = np.median(data)
    mad = np.median(np.abs(data - median))
    threshold_upper = median + factor * mad
    threshold_lower = median - factor * mad
    return (data > threshold_upper) | (data < threshold_lower)

def classify_outliers(adata, factors):
    """ Classify cells as outliers based on QC metrics """
    outliers_total_counts = compute_mad_outliers(adata.obs['total_counts'], factor=factors['total_counts'])
    outliers_n_genes = compute_mad_outliers(adata.obs['n_genes_by_counts'], factor=factors['n_genes'])
    outliers_pct_top_20 = compute_mad_outliers(adata.obs['pct_counts_in_top_20_genes'], factor=factors['top_20_genes'])
    outliers_pct_mt = compute_mad_outliers(adata.obs['pct_counts_mt'], factor=factors['pct_mt'])

    adata.obs['outlier'] = outliers_total_counts & outliers_n_genes & outliers_pct_top_20 & outliers_pct_mt

def export_metrics(adata, output_dir, filename="cell_metrics.csv"):
    """ Export QC metrics to a CSV file """
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, filename)
    metrics_df = adata.obs[['total_counts', 'n_genes_by_counts', 'pct_counts_mt', 'complexity', 'pct_counts_in_top_20_genes', 'outlier']]
    metrics_df.index.name = 'CellId'
    metrics_df.reset_index().to_csv(output_path, index=False)

def get_mito_prefix(species):
    """ Return mitochondrial gene prefix based on species """
    prefix_dict = {
        'saccharomyces-cerevisiae': 'MT-',
        'homo-sapiens': 'MT-',
        'mus-musculus': 'mt-',
        'rattus-norvegicus': 'MT-',
        'danio-rerio': 'mt-',
        'drosophila-melanogaster': 'mt-',
        'arabidopsis-thaliana': 'ATMT-',
        'caenorhabditis-elegans': 'mt-',
        'gallus-gallus': 'MT-',
        'bos-taurus': 'MT-',
        'sus-scrofa': 'MT-'
    }
    return prefix_dict.get(species.lower(), 'MT-')

def main():
    parser = argparse.ArgumentParser(description='Calculate QC metrics for scRNA-seq data.')
    parser.add_argument('--matrix', type=str, required=True, help='Path to matrix.mtx.gz file')
    parser.add_argument('--barcodes', type=str, required=True, help='Path to barcodes.tsv.gz file')
    parser.add_argument('--features', type=str, required=True, help='Path to features.tsv.gz file')
    parser.add_argument('--species', type=str, required=True, help='Species (e.g., homo-sapiens)')
    parser.add_argument('--output', type=str, required=True, help='Output directory')

    args = parser.parse_args()

    # Load data
    adata = load_data(args.matrix, args.barcodes, args.features)

    # Calculate QC metrics
    calculate_metrics(adata)

    # Calculate mitochondrial gene expression percentage
    mito_prefix = get_mito_prefix(args.species)
    mitochondrial_percentage(adata, mito_prefix)

    # Compute complexity
    compute_complexity(adata)

    # Classify outliers
    factors = {'total_counts': 5, 'n_genes': 5, 'top_20_genes': 5, 'pct_mt': 3}
    classify_outliers(adata, factors)

    # Export QC metrics
    export_metrics(adata, args.output)

if __name__ == "__main__":
    main()
