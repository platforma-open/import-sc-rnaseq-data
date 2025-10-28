import scanpy as sc
import pandas as pd
import numpy as np
import argparse
import os
import scipy.sparse
import scipy.io

def clean_barcode_suffix(barcode):
    """Remove '-' and following numbers from cell barcode (for 10X Genomics format).
    
    Only cleans if:
    - There's exactly one dash
    - The suffix after the dash is purely numeric (e.g., "-1", "-2")
    - The prefix contains only nucleotide letters (A, T, G, C, N)
    
    This distinguishes 10X barcodes from barcodes that use dashes as part of the identifier.
    """
    if '-' in barcode:
        parts = barcode.split('-')
        # Only clean if there's exactly one dash and the suffix is purely numeric
        if len(parts) == 2 and parts[1].isdigit():
            prefix = parts[0].upper()
            # Check if prefix contains only nucleotide letters (A, T, G, C, N)
            if all(c in 'ATCGN' for c in prefix):
                return parts[0]
    return barcode

def load_data_mtx(matrix_path, barcodes_path, features_path):
    """Load data from 10X Genomics format files."""
    import scipy.io

    # Load the matrix
    matrix = scipy.io.mmread(matrix_path).tocsr()

    # Load barcodes (cell IDs)
    barcodes_raw = pd.read_csv(barcodes_path, header=None)[0].tolist()
    barcodes = [clean_barcode_suffix(b) for b in barcodes_raw]

    # Load features (gene IDs and gene symbols)
    features = pd.read_csv(features_path, sep='\t', header=None)
    
    # Handle both 2-column and 3-column formats
    if features.shape[1] == 3:
        features.columns = ['gene_id', 'gene_symbol', 'feature_type']
    elif features.shape[1] == 2:
        features.columns = ['gene_id', 'gene_symbol']
        features['feature_type'] = 'Gene Expression'  # Add default feature type
    else:
        raise ValueError(f"Features file must have 2 or 3 columns, found {features.shape[1]}")

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
    barcode_duplicates = pd.Series(barcodes).duplicated()
    if barcode_duplicates.any():
        print(f"⚠️ Found {barcode_duplicates.sum()} duplicate barcodes.")
        # Keep first occurrence of each unique barcode
        unique_indices = ~barcode_duplicates
        barcodes = pd.Series(barcodes)[unique_indices].tolist()
        # Filter matrix columns to match unique barcodes
        if matrix.shape[1] == len(unique_indices):
            matrix = matrix[:, unique_indices]
            print(f"Filtered matrix to {unique_indices.sum()} unique cells")

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

def load_data_csv(csv_path):
    """Load data from long format CSV file (output from generate_counts_csv.py)."""
    print(f"Loading long format CSV: {csv_path}")
    
    # Read the CSV file
    df = pd.read_csv(csv_path)
    
    # Validate required columns
    required_columns = ['CellId', 'GeneId', 'Count']
    if not all(col in df.columns for col in required_columns):
        raise ValueError(f"CSV file must contain columns: {required_columns}")
    
    print(f"Found {len(df)} entries")
    print(f"Unique cells: {df['CellId'].nunique()}")
    print(f"Unique genes: {df['GeneId'].nunique()}")
    
    # Check for duplicate CellId-GeneId combinations
    duplicates = df.duplicated(subset=['CellId', 'GeneId'], keep=False)
    if duplicates.any():
        print(f"⚠️ Found {duplicates.sum()} duplicate CellId-GeneId combinations. Keeping maximum counts.")
        df = df.loc[df.groupby(['CellId', 'GeneId'])['Count'].idxmax()]
        print(f"After deduplication: {len(df)} entries")
    
    # Pivot to create matrix format (cells × genes)
    print("Converting to matrix format...")
    matrix_df = df.pivot(index='CellId', columns='GeneId', values='Count').fillna(0)
    
    # Convert to sparse matrix for efficiency
    matrix = scipy.sparse.csr_matrix(matrix_df.values)
    
    # Create AnnData object
    adata = sc.AnnData(matrix)
    adata.obs_names = pd.Index(matrix_df.index, name='cell_id')
    adata.var_names = pd.Index(matrix_df.columns, name='gene_id')
    
    # For gene symbols, we'll use the gene IDs as symbols (since we don't have separate symbols)
    adata.var['gene_symbol'] = matrix_df.columns.values
    
    print(f"Created AnnData object: {adata.shape[0]} cells × {adata.shape[1]} genes")
    
    return adata

def load_data_h5ad(h5ad_path):
    """Load data from AnnData h5ad file."""
    print(f"Loading h5ad file: {h5ad_path}")
    
    # Load the h5ad file
    adata = sc.read_h5ad(h5ad_path)
    
    print(f"AnnData shape: {adata.shape[0]} cells × {adata.shape[1]} genes")
    
    # Clean cell barcodes if needed
    original_barcodes = list(adata.obs_names)
    cleaned_barcodes = [clean_barcode_suffix(b) for b in original_barcodes]
    
    if cleaned_barcodes != original_barcodes:
        print("Cleaning barcode suffixes...")
        adata.obs_names = cleaned_barcodes
        print(f"Example: {original_barcodes[0]} -> {cleaned_barcodes[0]}")
    
    # Check for gene_symbol in var
    if 'gene_symbol' not in adata.var.columns:
        # Use var_names as gene_symbol if not present
        adata.var['gene_symbol'] = adata.var_names
    
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
    parser.add_argument('--format', required=True, choices=['mtx', 'csv', 'h5ad'], 
                       help="Input format: 'mtx' for 10X Genomics format, 'csv' for long format CSV, 'h5ad' for AnnData format")
    parser.add_argument('--matrix', help='Path to matrix.mtx.gz file (required for mtx format)')
    parser.add_argument('--barcodes', help='Path to barcodes.tsv.gz file (required for mtx format)')
    parser.add_argument('--features', help='Path to features.tsv.gz file (required for mtx format)')
    parser.add_argument('--csv', help='Path to long format CSV file (required for csv format)')
    parser.add_argument('--h5ad', help='Path to h5ad file (required for h5ad format)')
    parser.add_argument('--species', type=str, required=True, help='Species (e.g., homo-sapiens)')
    parser.add_argument('--output', type=str, required=True, help='Output directory')

    args = parser.parse_args()
    
    # Load data based on format
    if args.format == 'mtx':
        if not all([args.matrix, args.barcodes, args.features]):
            parser.error("For mtx format, --matrix, --barcodes, and --features are required")
        adata = load_data_mtx(args.matrix, args.barcodes, args.features)
    elif args.format == 'csv':
        if not args.csv:
            parser.error("For csv format, --csv is required")
        adata = load_data_csv(args.csv)
    elif args.format == 'h5ad':
        if not args.h5ad:
            parser.error("For h5ad format, --h5ad is required")
        adata = load_data_h5ad(args.h5ad)

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
