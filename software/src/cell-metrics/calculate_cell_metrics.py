import scanpy as sc
import pandas as pd
import numpy as np
import argparse
import os
import scipy.sparse
import scipy.io
from datetime import datetime

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
        print("Found duplicate gene IDs. Removing duplicates.")
        features = features.drop_duplicates(subset='gene_id', keep='first')
        matrix = matrix[:features.shape[0], :]
    else:
        print("Gene IDs are unique.")

    # Check for duplicate barcodes
    barcode_duplicates = pd.Series(barcodes).duplicated()
    if barcode_duplicates.any():
        print(f"Found {barcode_duplicates.sum()} duplicate barcodes.")
        # Keep first occurrence of each unique barcode
        unique_indices = ~barcode_duplicates
        barcodes = pd.Series(barcodes)[unique_indices].tolist()
        # Filter matrix columns to match unique barcodes
        if matrix.shape[1] == len(unique_indices):
            matrix = matrix[:, unique_indices]
            print(f"Filtered matrix to {unique_indices.sum()} unique cells")

    # FINAL FIX: Transpose the matrix
    if matrix.shape[0] == len(features) and matrix.shape[1] == len(barcodes):
        print("Transposing matrix to match AnnData format (cells x genes).")
        matrix = matrix.transpose()
    else:
        print("Matrix orientation matches AnnData expectations.")

    # Final assignment
    adata = sc.AnnData(matrix)
    adata.var_names = pd.Index(features['gene_id'].values, name='gene_id')
    adata.obs_names = pd.Index(barcodes, name='cell_id')
    adata.var['gene_symbol'] = features['gene_symbol'].values

    return adata

def load_data_csv(csv_path):
    """Load data from long format Parquet file (output from generate_counts_csv.py).
    Supports both Parquet and CSV files for backward compatibility.
    """
    print(f"Loading long format file: {csv_path}")
    
    # Read Parquet file if .parquet extension, otherwise read CSV (backward compatibility)
    if csv_path.endswith('.parquet'):
        df = pd.read_parquet(csv_path)
    else:
        # Read the CSV file (pandas auto-detects gzip compression for .gz files)
        compression = 'gzip' if csv_path.endswith('.gz') else None
        df = pd.read_csv(csv_path, compression=compression)
    
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
        print(f"Found {duplicates.sum()} duplicate CellId-GeneId combinations. Keeping maximum counts.")
        df = df.loc[df.groupby(['CellId', 'GeneId'])['Count'].idxmax()]
        print(f"After deduplication: {len(df)} entries")
    
    # Create sparse matrix directly from long format (much more efficient than pivot)
    print("Converting to matrix format...")
    
    # Get unique cells and genes
    unique_cells = df['CellId'].unique()
    unique_genes = df['GeneId'].unique()
    
    # Create mappings for indices
    cell_to_idx = {cell: idx for idx, cell in enumerate(unique_cells)}
    gene_to_idx = {gene: idx for idx, gene in enumerate(unique_genes)}
    
    # Build sparse matrix using COO format
    row_indices = [cell_to_idx[cell] for cell in df['CellId']]
    col_indices = [gene_to_idx[gene] for gene in df['GeneId']]
    values = df['Count'].values
    
    # Create sparse matrix (cells × genes)
    matrix = scipy.sparse.csr_matrix(
        (values, (row_indices, col_indices)),
        shape=(len(unique_cells), len(unique_genes))
    )
    
    # Create AnnData object
    adata = sc.AnnData(matrix)
    adata.obs_names = pd.Index(unique_cells, name='cell_id')
    adata.var_names = pd.Index(unique_genes, name='gene_id')
    
    # For gene symbols, we'll use the gene IDs as symbols (since we don't have separate symbols)
    adata.var['gene_symbol'] = unique_genes
    
    print(f"Created AnnData object: {adata.shape[0]} cells × {adata.shape[1]} genes")
    
    return adata

def find_sample_column(obs_df, target_sample):
    """
    Finds the column in the .obs DataFrame that contains the target sample name.
    Prioritizes columns with 'sample' in the name.
    """
    # Prioritize columns with 'sample' in the name, then check all others
    candidate_columns = sorted(obs_df.columns, key=lambda col: 'sample' not in col.lower())

    for col_name in candidate_columns:
        # Check only string or categorical columns
        if obs_df[col_name].dtype.name not in ['category', 'object']:
            continue
        
        unique_values = set(obs_df[col_name].unique())
        
        # Check if the target sample is present in the column's unique values
        if target_sample in unique_values:
            print(f"Found sample column: '{col_name}' containing sample '{target_sample}'")
            return col_name
            
    return None

def load_data_h5ad(h5ad_path, sample_name=None, sample_column_name=None):
    """Load data from AnnData h5ad file.
    
    Args:
        h5ad_path: Path to the h5ad file
        sample_name: Optional sample name to filter cells by
        sample_column_name: Optional column name to use for sample filtering
    """
    print(f"Loading h5ad file: {h5ad_path}")
    
    # Load the h5ad file
    adata = sc.read_h5ad(h5ad_path)
    
    # Filter by sample if sample_name is provided
    if sample_name:
        print(f"Filtering for sample: {sample_name}")
        
        if sample_column_name:
            # Use provided column name
            print(f"Using provided sample column name: '{sample_column_name}'")
            if sample_column_name not in adata.obs.columns:
                raise ValueError(f"Provided sample column '{sample_column_name}' not found in the h5ad file's .obs dataframe. Available columns: {list(adata.obs.columns)}")
            sample_column = sample_column_name
        else:
            # Auto-detect column name
            sample_column = find_sample_column(adata.obs, sample_name)
            
            if sample_column is None:
                raise ValueError(f"Could not automatically identify the sample column containing '{sample_name}' in the h5ad file's .obs dataframe.")
        
        # Filter AnnData for the target sample
        adata = adata[adata.obs[sample_column] == sample_name].copy()
        
        if adata.n_obs == 0:
            raise ValueError(f"No cells found for sample '{sample_name}'.")
        
        print(f"Filtered to {adata.n_obs} cells for sample '{sample_name}'")
    
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
    """ Calculate mitochondrial gene expression percentage (optimized for sparse matrices) """
    adata.var['mt'] = adata.var['gene_symbol'].str.startswith(mito_prefix, na=False)
    
    # Use sparse-aware operations to avoid unnecessary dense conversions
    if scipy.sparse.issparse(adata.X):
        # For sparse matrices, use scanpy's built-in calculation which is optimized
        mito_genes = adata.var['mt'].values
        mito_counts = np.array(adata.X[:, mito_genes].sum(axis=1)).flatten()
        total_counts = np.array(adata.X.sum(axis=1)).flatten()
    else:
        # For dense matrices
        mito_counts = adata[:, adata.var['mt']].X.sum(axis=1)
        total_counts = adata.X.sum(axis=1)
    
    # Calculate percentage, avoiding division by zero
    adata.obs['pct_counts_mt'] = np.divide(
        mito_counts, 
        total_counts, 
        out=np.zeros_like(mito_counts, dtype=float),
        where=total_counts != 0
    ) * 100

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

def export_metrics(adata, output_dir, filename="cell_metrics.csv", sample_id=None):
    """ Export QC metrics to a CSV file """
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, filename)
    metrics_df = adata.obs[['total_counts', 'n_genes_by_counts', 'pct_counts_mt', 'complexity', 'pct_counts_in_top_20_genes', 'outlier']]
    metrics_df.index.name = 'CellId'
    metrics_df = metrics_df.reset_index()
    
    if sample_id:
        # Add SampleId as first column
        metrics_df.insert(0, 'SampleId', sample_id)
    
    metrics_df.to_csv(output_path, index=False)

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
    # Setup timing logging
    script_name = os.path.splitext(os.path.basename(__file__))[0]
    log_file = f"{script_name}.time.log"
    start_time = datetime.now()
    
    with open(log_file, 'w') as f:
        f.write(f"Script: {script_name}\n")
        f.write(f"Start time: {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    parser = argparse.ArgumentParser(description='Calculate QC metrics for scRNA-seq data.')
    parser.add_argument('--format', required=True, choices=['mtx', 'parquet', 'h5ad'], 
                       help="Input format: 'mtx' for 10X Genomics format, 'parquet' for long format Parquet/CSV, 'h5ad' for AnnData format")
    parser.add_argument('--matrix', help='Path to matrix.mtx.gz file (required for mtx format)')
    parser.add_argument('--barcodes', help='Path to barcodes.tsv.gz file (required for mtx format)')
    parser.add_argument('--features', help='Path to features.tsv.gz file (required for mtx format)')
    parser.add_argument('--input-file', help='Path to long format Parquet(required for parquet format)')
    parser.add_argument('--h5ad', help='Path to h5ad file (required for h5ad format)')
    parser.add_argument('--sample-name', help="Sample name to filter cells by (optional, for h5ad format only)")
    parser.add_argument('--sample-column-name', help="Column name in h5ad .obs that contains sample identifiers (optional, for h5ad format only)")
    parser.add_argument('--sample-id', help="Sample ID to add as first column in output CSV")
    parser.add_argument('--species', type=str, required=True, help='Species (e.g., homo-sapiens)')
    parser.add_argument('--output', type=str, required=True, help='Output directory')

    args = parser.parse_args()
    
    try:
        # Load data based on format
        if args.format == 'mtx':
            if not all([args.matrix, args.barcodes, args.features]):
                parser.error("For mtx format, --matrix, --barcodes, and --features are required")
            adata = load_data_mtx(args.matrix, args.barcodes, args.features)
        elif args.format == 'parquet':
            if not args.input_file:
                parser.error("For parquet format, --input-file is required")
            adata = load_data_csv(args.input_file)
        elif args.format == 'h5ad':
            if not args.h5ad:
                parser.error("For h5ad format, --h5ad is required")
            adata = load_data_h5ad(args.h5ad, args.sample_name, args.sample_column_name)

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
        export_metrics(adata, args.output, sample_id=args.sample_id)
    except Exception as e:
        # Log end time even on error
        end_time = datetime.now()
        duration = end_time - start_time
        with open(log_file, 'a') as f:
            f.write(f"End time: {end_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Total duration: {duration.total_seconds():.2f} seconds ({duration})\n")
            f.write(f"Status: ERROR - {str(e)}\n")
        raise
    
    # Log end time and duration
    end_time = datetime.now()
    duration = end_time - start_time
    
    with open(log_file, 'a') as f:
        f.write(f"End time: {end_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Total duration: {duration.total_seconds():.2f} seconds ({duration})\n")

if __name__ == "__main__":
    main()
