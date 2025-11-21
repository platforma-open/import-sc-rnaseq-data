import argparse
import pandas as pd
import polars as pl
import scipy.io
import gzip
import scanpy as sc
import anndata
import os
import sys
from pathlib import Path
from datetime import datetime
from io import BytesIO 

def write_parquet(df, output_path):
    """Write DataFrame to Parquet format using Polars for optimized performance.
    
    Args:
        df: Pandas DataFrame to write
        output_path: Path to output Parquet file
    """
    start_time = datetime.now()
    total_rows = len(df)
    
    # Ensure output path ends with .parquet
    if not output_path.endswith('.parquet'):
        output_path = output_path + '.parquet'
    
    # Convert pandas DataFrame to Polars DataFrame
    conversion_start = datetime.now()
    pl_df = pl.from_pandas(df)
    conversion_time = (datetime.now() - conversion_start).total_seconds()
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Converted to Polars DataFrame in {conversion_time:.1f}s")
    
    # Write to Parquet format
    pl_df.write_parquet(output_path, compression='zstd')
    
    end_time = datetime.now()
    elapsed_time = (end_time - start_time).total_seconds()
    print(f"[{end_time.strftime('%Y-%m-%d %H:%M:%S')}] Completed writing {total_rows:,} rows to {output_path} - Total time: {elapsed_time:.1f}s")

def read_gzip_tsv(file_path):
    """Reads a gzipped TSV file into a pandas DataFrame."""
    with gzip.open(file_path, 'rt') as f:
        return pd.read_csv(f, sep='\t', header=None)

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

def clean_gene_name(gene_name):
    """Clean gene name by removing special characters like newlines, tabs, etc.
    
    Removes:
    - Newlines (\n, \r)
    - Tabs (\t)
    - Other control characters
    - Leading/trailing whitespace
    """
    if not isinstance(gene_name, str):
        gene_name = str(gene_name)
    
    # Remove newlines, carriage returns, tabs
    gene_name = gene_name.replace('\n', '').replace('\r', '').replace('\t', '')
    
    # Remove other control characters (ASCII 0-31 except space)
    gene_name = ''.join(char for char in gene_name if ord(char) >= 32)
    
    # Strip leading/trailing whitespace
    gene_name = gene_name.strip()
    
    return gene_name

def process_input_files(matrix_path, barcodes_path, features_path, output_csv_path, sample_id=None):
    # Load the input files
    print("Loading matrix.mtx.gz...")
    matrix = scipy.io.mmread(matrix_path).tocoo()  # Sparse COO format

    print("Loading barcodes.tsv.gz...")
    barcodes = read_gzip_tsv(barcodes_path)
    
    print("Loading features.tsv.gz...")
    features = read_gzip_tsv(features_path)
    
    barcodes_list = [clean_barcode_suffix(barcode) for barcode in barcodes[0].tolist()]
    print(f"Cleaned barcode suffixes. Example: {barcodes[0].tolist()[0]} -> {barcodes_list[0]}")
    features_list = [clean_gene_name(gene) for gene in features[0].tolist()]

    data = []

    print(f"Processing {matrix.nnz} nonzero entries...")
    count_debug = 0

    for i, j, value in zip(matrix.row, matrix.col, matrix.data):
        if count_debug < 10:
            print(f"Row: {i}, Column: {j}, Count: {value}")
        cell_id = barcodes_list[j]
        gene_id = features_list[i]
        count = value
        if sample_id:
            data.append([sample_id, cell_id, gene_id, count])
        else:
            data.append([cell_id, gene_id, count])
        count_debug += 1

    print(f"Total rows written: {count_debug}")
    if sample_id:
        df = pd.DataFrame(data, columns=["SampleId", "CellId", "GeneId", "Count"])
    else:
        df = pd.DataFrame(data, columns=["CellId", "GeneId", "Count"])

    print(f"Writing raw count matrix to {output_csv_path}...")
    write_parquet(df, output_csv_path)

    # Normalize counts
    print("Normalizing counts...")

    # Reconstruct sparse matrix in CSC format (Scanpy prefers this)
    matrix = matrix.tocsc()

    adata = anndata.AnnData(X=matrix.transpose())
    adata.obs_names = barcodes_list
    adata.var_names = features_list

    sc.pp.normalize_total(adata, target_sum=1e4)
    normalized_matrix = adata.X.tocoo()

    # Extract normalized counts
    norm_data = []
    for i, j, value in zip(normalized_matrix.row, normalized_matrix.col, normalized_matrix.data):
        cell_id = adata.obs_names[i]
        gene_id = adata.var_names[j]
        if sample_id:
            norm_data.append([sample_id, cell_id, gene_id, value])
        else:
            norm_data.append([cell_id, gene_id, value])

    if sample_id:
        norm_df = pd.DataFrame(norm_data, columns=["SampleId", "CellId", "GeneId", "NormalizedCount"])
    else:
        norm_df = pd.DataFrame(norm_data, columns=["CellId", "GeneId", "NormalizedCount"])
    # Generate normalized output path with .parquet extension
    if output_csv_path.endswith(".parquet"):
        normalized_output_path = output_csv_path.replace(".parquet", "_normalized.parquet")
    elif output_csv_path.endswith(".csv.gz"):
        normalized_output_path = output_csv_path.replace(".csv.gz", "_normalized.parquet")
    elif output_csv_path.endswith(".csv"):
        normalized_output_path = output_csv_path.replace(".csv", "_normalized.parquet")
    else:
        normalized_output_path = output_csv_path + "_normalized.parquet"

    print(f"Writing normalized count matrix to {normalized_output_path}...")
    write_parquet(norm_df, normalized_output_path)

    print("Done!")

def detect_orientation(df):
    """Detect if genes are in rows (index) or columns."""
    # Get first row and first column to analyze
    first_row = df.iloc[0].tolist()
    first_col = df.iloc[:, 0].tolist()
    
    # Count numeric vs non-numeric values
    def count_numeric(values):
        numeric_count = 0
        for val in values:
            try:
                float(val)
                numeric_count += 1
            except (ValueError, TypeError):
                pass
        return numeric_count
    
    row_numeric = count_numeric(first_row)
    col_numeric = count_numeric(first_col)
    
    # If more numeric values in columns, genes are in rows
    # If more numeric values in rows, genes are in columns (transposed)
    if row_numeric > col_numeric:
        return "transposed"  # genes in columns, cells in rows
    else:
        return "normal"      # genes in rows, cells in columns

def load_gene_annotation(annotation_path):
    """Load gene annotation file and create symbol to Ensembl ID mapping."""
    print(f"Loading gene annotation file: {annotation_path}")
    
    # Read the annotation file
    annotation_df = pd.read_csv(annotation_path)
    
    # Create mapping from gene symbol to Ensembl ID
    symbol_to_ensembl = {}
    for _, row in annotation_df.iterrows():
        symbol = row['Gene symbol']
        ensembl_id = row['Ensembl Id']
        if pd.notna(symbol) and pd.notna(ensembl_id):
            # Clean both symbol and Ensembl ID
            clean_symbol = clean_gene_name(str(symbol))
            clean_ensembl_id = clean_gene_name(str(ensembl_id))
            symbol_to_ensembl[clean_symbol] = clean_ensembl_id
    
    print(f"Loaded {len(symbol_to_ensembl)} gene symbol to Ensembl ID mappings")
    return symbol_to_ensembl

def convert_gene_symbols_to_ensembl(df, symbol_to_ensembl, gene_names):
    """Convert gene symbols to Ensembl IDs in the dataframe."""
    print("Converting gene symbols to Ensembl IDs...")
    
    # Create mapping for gene names
    gene_mapping = {}
    converted_count = 0
    not_found_genes = []
    
    for gene in gene_names:
        if gene in symbol_to_ensembl:
            gene_mapping[gene] = symbol_to_ensembl[gene]
            converted_count += 1
        else:
            gene_mapping[gene] = gene  # Keep original if not found
            not_found_genes.append(gene)
    
    print(f"Converted {converted_count} gene symbols to Ensembl IDs")
    if not_found_genes:
        print(f"Warning: {len(not_found_genes)} genes not found in annotation file (keeping original names)")
        if len(not_found_genes) <= 10:
            print(f"Not found genes: {', '.join(not_found_genes)}")
        else:
            print(f"Not found genes (first 10): {', '.join(not_found_genes[:10])}...")
    
    # Update the dataframe index with Ensembl IDs
    df.index = [gene_mapping[gene] for gene in df.index]
    
    return df, list(df.index)

def handle_duplicate_combinations(df_long):
    """Handle duplicate CellId-GeneId combinations by keeping the entry with maximum count."""
    original_count = len(df_long)
    
    # Check for duplicates
    duplicates = df_long.duplicated(subset=['CellId', 'GeneId'], keep=False)
    duplicate_count = duplicates.sum()
    
    if duplicate_count > 0:
        print(f"Warning: Found {duplicate_count} duplicate CellId-GeneId combinations")
        
        # Group by CellId and GeneId, keep the entry with maximum count
        df_long = df_long.loc[df_long.groupby(['CellId', 'GeneId'])['Count'].idxmax()]
        
        final_count = len(df_long)
        removed_count = original_count - final_count
        
        print(f"Removed {removed_count} duplicate entries (kept entries with maximum count)")
        print(f"Final unique combinations: {final_count}")
    else:
        print("No duplicate CellId-GeneId combinations found")
    
    return df_long

def handle_duplicate_combinations_normalized(df_long):
    """Handle duplicate CellId-GeneId combinations in normalized data by keeping the entry with maximum normalized count."""
    original_count = len(df_long)
    
    # Check for duplicates
    duplicates = df_long.duplicated(subset=['CellId', 'GeneId'], keep=False)
    duplicate_count = duplicates.sum()
    
    if duplicate_count > 0:
        print(f"Warning: Found {duplicate_count} duplicate CellId-GeneId combinations in normalized data")
        
        # Group by CellId and GeneId, keep the entry with maximum normalized count
        df_long = df_long.loc[df_long.groupby(['CellId', 'GeneId'])['NormalizedCount'].idxmax()]
        
        final_count = len(df_long)
        removed_count = original_count - final_count
        
        print(f"Removed {removed_count} duplicate entries from normalized data (kept entries with maximum normalized count)")
        print(f"Final unique combinations in normalized data: {final_count}")
    else:
        print("No duplicate CellId-GeneId combinations found in normalized data")
    
    return df_long

def process_csv_file(csv_path, output_csv_path, gene_format=None, annotation_path=None, sample_id=None):
    """Process a CSV/TSV file with automatic format detection, orientation detection and optional gene format conversion."""
    print(f"Loading file: {csv_path}")
    
    # Detect file format and read accordingly
    file_path = Path(csv_path)
    if file_path.suffix.lower() == '.csv' or (file_path.suffix.lower() == '.gz' and file_path.stem.endswith('.csv')):
        # CSV format
        df = pd.read_csv(csv_path, index_col=0)
    elif file_path.suffix.lower() == '.tsv' or (file_path.suffix.lower() == '.gz' and file_path.stem.endswith('.tsv')):
        # TSV format
        df = pd.read_csv(csv_path, sep='\t', index_col=0)
    else:
        # Try to auto-detect format
        try:
            df = pd.read_csv(csv_path, index_col=0)
        except Exception:
            try:
                df = pd.read_csv(csv_path, sep='\t', index_col=0)
            except Exception as e:
                raise ValueError(f"Could not read file as CSV or TSV: {e}")
    
    # Detect orientation
    orientation = detect_orientation(df)
    print(f"Detected orientation: {orientation}")
    
    if orientation == "transposed":
        # Genes in columns, cells in rows - transpose the dataframe
        print("Transposing matrix: genes in columns -> genes in rows")
        df = df.T
        gene_names = [clean_gene_name(gene) for gene in df.index.tolist()]
        cell_names = df.columns.tolist()
    else:
        # Genes in rows, cells in columns - use as is
        gene_names = [clean_gene_name(gene) for gene in df.index.tolist()]
        cell_names = df.columns.tolist()
    
    print(f"Found {len(gene_names)} genes and {len(cell_names)} cells")
    
    # Convert gene symbols to Ensembl IDs if requested
    if gene_format == 'gene symbol' and annotation_path:
        symbol_to_ensembl = load_gene_annotation(annotation_path)
        df, gene_names = convert_gene_symbols_to_ensembl(df, symbol_to_ensembl, gene_names)
    elif gene_format == 'gene symbol' and not annotation_path:
        print("Warning: Gene format is 'gene symbol' but no annotation file provided. Using original gene names.")
    
    # Clean cell barcodes (remove -1, -2, etc. suffixes)
    cleaned_cell_names = [clean_barcode_suffix(cell) for cell in cell_names]
    print(f"Cleaned barcode suffixes. Example: {cell_names[0]} -> {cleaned_cell_names[0]}")
    
    # Convert to long format using pandas stack (much faster than nested loops)
    print(f"Processing {df.shape[0]} genes x {df.shape[1]} cells...")
    
    # Update dataframe with cleaned names
    df.index = gene_names
    df.columns = cleaned_cell_names
    
    # Stack to convert wide to long format: creates a Series with MultiIndex (gene, cell)
    df_stacked = df.stack()
    
    # Filter out zeros (much faster than doing it in a loop)
    df_stacked = df_stacked[df_stacked > 0]
    
    print(f"Total non-zero entries: {len(df_stacked)}")
    
    # Convert to DataFrame with proper column names
    df_long = df_stacked.reset_index()
    df_long.columns = ["GeneId", "CellId", "Count"]
    
    # Reorder columns and add sample_id if needed
    if sample_id:
        df_long.insert(0, "SampleId", sample_id)
        df_long = df_long[["SampleId", "CellId", "GeneId", "Count"]]
    else:
        df_long = df_long[["CellId", "GeneId", "Count"]]
    
    # Check for and handle duplicate CellId-GeneId combinations
    df_long = handle_duplicate_combinations(df_long)
    
    print(f"Writing raw count matrix to {output_csv_path}...")
    write_parquet(df_long, output_csv_path)
    
    # Normalize counts
    print("Normalizing counts...")
    
    # Create AnnData object for normalization
    adata = anndata.AnnData(X=df.values.T)  # Transpose: cells x genes
    adata.obs_names = cleaned_cell_names
    adata.var_names = gene_names
    
    sc.pp.normalize_total(adata, target_sum=1e4)
    normalized_matrix = adata.X
    
    # Extract normalized counts using pandas (much faster than nested loops)
    # Convert normalized matrix back to DataFrame
    norm_df_wide = pd.DataFrame(normalized_matrix, index=cleaned_cell_names, columns=gene_names)
    
    # Stack to convert wide to long format
    norm_stacked = norm_df_wide.stack()
    
    # Filter out zeros
    norm_stacked = norm_stacked[norm_stacked > 0]
    
    # Convert to DataFrame with proper column names
    norm_df = norm_stacked.reset_index()
    norm_df.columns = ["CellId", "GeneId", "NormalizedCount"]
    
    # Add sample_id if needed
    if sample_id:
        norm_df.insert(0, "SampleId", sample_id)
        norm_df = norm_df[["SampleId", "CellId", "GeneId", "NormalizedCount"]]
    else:
        norm_df = norm_df[["CellId", "GeneId", "NormalizedCount"]]
    
    # Check for and handle duplicate CellId-GeneId combinations in normalized data
    norm_df = handle_duplicate_combinations_normalized(norm_df)
    
    # Generate normalized output path with .parquet extension
    if output_csv_path.endswith(".parquet"):
        normalized_output_path = output_csv_path.replace(".parquet", "_normalized.parquet")
    elif output_csv_path.endswith(".csv.gz"):
        normalized_output_path = output_csv_path.replace(".csv.gz", "_normalized.parquet")
    elif output_csv_path.endswith(".csv"):
        normalized_output_path = output_csv_path.replace(".csv", "_normalized.parquet")
    else:
        normalized_output_path = output_csv_path + "_normalized.parquet"
    
    print(f"Writing normalized count matrix to {normalized_output_path}...")
    write_parquet(norm_df, normalized_output_path)
    
    print("Done!")

def find_sample_column(obs_df, target_sample):
    """
    Finds the column in the .obs DataFrame that contains the target sample name.
    Prioritizes columns with 'sample' in the name.
    """
    # Prioritize columns with 'sample' in the name, then check all others
    candidate_columns = [col for col in obs_df.columns if 'sample' in col.lower()]
    candidate_columns.extend([col for col in obs_df.columns if col not in candidate_columns])

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

def process_h5ad_file(h5ad_path, output_csv_path, sample_name=None, sample_id=None, sample_column_name=None, gene_format=None, annotation_path=None):
    """Process an AnnData h5ad file and convert to long format CSV.
    
    Args:
        h5ad_path: Path to the h5ad file
        output_csv_path: Path to save the output CSV
        sample_name: Optional sample name to filter cells by
        sample_id: Optional sample ID to add as first column
        sample_column_name: Optional column name to use for sample filtering
        gene_format: Optional gene identifier format ('gene symbol' or 'Ensembl Id')
        annotation_path: Optional path to gene annotation CSV file (required when gene_format is 'gene symbol')
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
    
    print(f"AnnData shape: {adata.shape} (cells × genes)")
    
    # Extract gene and cell identifiers
    gene_names = [clean_gene_name(gene) for gene in adata.var_names]
    cell_names = list(adata.obs_names)
    
    # Convert gene symbols to Ensembl IDs if requested
    if gene_format == 'gene symbol' and annotation_path:
        symbol_to_ensembl = load_gene_annotation(annotation_path)
        # Create a mapping for gene names
        gene_mapping = {}
        for gene in gene_names:
            if gene in symbol_to_ensembl:
                gene_mapping[gene] = symbol_to_ensembl[gene]
            else:
                gene_mapping[gene] = gene  # Keep original if not found
        # Update gene names with Ensembl IDs
        gene_names = [gene_mapping[gene] for gene in gene_names]
    elif gene_format == 'gene symbol' and not annotation_path:
        print("Warning: Gene format is 'gene symbol' but no annotation file provided. Using original gene names.")
    
    # Clean cell barcodes if needed
    cleaned_cell_names = [clean_barcode_suffix(cell) for cell in cell_names]
    if any(c != o for c, o in zip(cleaned_cell_names, cell_names)):
        print(f"Cleaned barcode suffixes. Example: {cell_names[0]} -> {cleaned_cell_names[0]}")
    
    # Get the count matrix
    X = adata.X
    
    # Convert to long format for raw counts (optimized)
    print(f"Processing {X.shape[0]} cells × {X.shape[1]} genes...")
    
    if hasattr(X, 'tocoo'):
        # Sparse matrix - create DataFrame directly from COO format
        X_coo = X.tocoo()
        df_raw = pd.DataFrame({
            'CellId': [cleaned_cell_names[i] for i in X_coo.row],
            'GeneId': [gene_names[j] for j in X_coo.col],
            'Count': X_coo.data
        })
        # Filter out zeros
        df_raw = df_raw[df_raw['Count'] > 0]
    else:
        # Dense matrix - use pandas stack for efficiency
        df_wide = pd.DataFrame(X, index=cleaned_cell_names, columns=gene_names)
        df_stacked = df_wide.stack()
        df_stacked = df_stacked[df_stacked > 0]
        df_raw = df_stacked.reset_index()
        df_raw.columns = ["CellId", "GeneId", "Count"]
    
    # Add sample_id if needed
    if sample_id:
        df_raw.insert(0, "SampleId", sample_id)
        df_raw = df_raw[["SampleId", "CellId", "GeneId", "Count"]]
    else:
        df_raw = df_raw[["CellId", "GeneId", "Count"]]
    
    print(f"Total non-zero entries: {len(df_raw)}")
    
    # Check for and handle duplicate CellId-GeneId combinations
    df_raw = handle_duplicate_combinations(df_raw)
    
    print(f"Writing raw count matrix to {output_csv_path}...")
    write_parquet(df_raw, output_csv_path)
    
    # Normalize counts
    print("Normalizing counts...")
    
    # Create a copy of adata for normalization
    adata_norm = adata.copy()
    sc.pp.normalize_total(adata_norm, target_sum=1e4)
    
    # Extract normalized counts (optimized)
    X_norm = adata_norm.X
    
    if hasattr(X_norm, 'tocoo'):
        # Sparse matrix - create DataFrame directly from COO format
        X_norm_coo = X_norm.tocoo()
        norm_df = pd.DataFrame({
            'CellId': [cleaned_cell_names[i] for i in X_norm_coo.row],
            'GeneId': [gene_names[j] for j in X_norm_coo.col],
            'NormalizedCount': X_norm_coo.data
        })
        # Filter out zeros
        norm_df = norm_df[norm_df['NormalizedCount'] > 0]
    else:
        # Dense matrix - use pandas stack for efficiency
        norm_df_wide = pd.DataFrame(X_norm, index=cleaned_cell_names, columns=gene_names)
        norm_stacked = norm_df_wide.stack()
        norm_stacked = norm_stacked[norm_stacked > 0]
        norm_df = norm_stacked.reset_index()
        norm_df.columns = ["CellId", "GeneId", "NormalizedCount"]
    
    # Add sample_id if needed
    if sample_id:
        norm_df.insert(0, "SampleId", sample_id)
        norm_df = norm_df[["SampleId", "CellId", "GeneId", "NormalizedCount"]]
    else:
        norm_df = norm_df[["CellId", "GeneId", "NormalizedCount"]]
    
    # Check for and handle duplicate CellId-GeneId combinations in normalized data
    norm_df = handle_duplicate_combinations_normalized(norm_df)
    
    # Generate normalized output path with .parquet extension
    if output_csv_path.endswith(".parquet"):
        normalized_output_path = output_csv_path.replace(".parquet", "_normalized.parquet")
    elif output_csv_path.endswith(".csv.gz"):
        normalized_output_path = output_csv_path.replace(".csv.gz", "_normalized.parquet")
    elif output_csv_path.endswith(".csv"):
        normalized_output_path = output_csv_path.replace(".csv", "_normalized.parquet")
    else:
        normalized_output_path = output_csv_path + "_normalized.parquet"
    
    print(f"Writing normalized count matrix to {normalized_output_path}...")
    write_parquet(norm_df, normalized_output_path)
    
    print("Done!")

def main():
    # Setup timing logging
    script_name = os.path.splitext(os.path.basename(__file__))[0]
    log_file = f"{script_name}.time.log"
    start_time = datetime.now()
    
    with open(log_file, 'w') as f:
        f.write(f"Script: {script_name}\n")
        f.write(f"Start time: {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    parser = argparse.ArgumentParser(description="Convert .mtx.gz, .tsv.gz files or CSV files into a count matrix Parquet file.")
    parser.add_argument('--format', required=True, choices=['mtx', 'xsv', 'h5ad'], 
                       help="Input format: 'mtx' for 10X Genomics format, 'xsv' for CSV/TSV format, 'h5ad' for AnnData format")
    parser.add_argument('--matrix', help="Path to the matrix.mtx.gz file (required for mtx format)")
    parser.add_argument('--barcodes', help="Path to the barcodes.tsv.gz file (required for mtx format)")
    parser.add_argument('--features', help="Path to the features.tsv.gz file (required for mtx format)")
    parser.add_argument('--xsv', help="Path to the XSV file (required for xsv format)")
    parser.add_argument('--h5ad', help="Path to the h5ad file (required for h5ad format)")
    parser.add_argument('--gene-format', choices=['gene symbol', 'Ensembl Id'], 
                       help="Gene identifier format: 'gene symbol' or 'Ensembl Id'")
    parser.add_argument('--annotation', help="Path to gene annotation CSV file (required when --gene-format is 'gene symbol')")
    parser.add_argument('--sample-name', help="Sample name to filter cells by (optional, for h5ad format only)")
    parser.add_argument('--sample-column-name', help="Column name in h5ad .obs that contains sample identifiers (optional, for h5ad format only)")
    parser.add_argument('--sample-id', help="Sample ID to add as first column in output CSV")
    parser.add_argument('--output', required=True, help="Path to output the raw CSV file")

    args = parser.parse_args()
    
    try:
        if args.format == 'mtx':
            if not all([args.matrix, args.barcodes, args.features]):
                parser.error("For mtx format, --matrix, --barcodes, and --features are required")
            process_input_files(args.matrix, args.barcodes, args.features, args.output, args.sample_id)
        elif args.format == 'xsv':
            if not args.xsv:
                parser.error("For xsv format, --xsv is required")
            if args.gene_format == 'gene symbol' and not args.annotation:
                parser.error("For gene format 'gene symbol', --annotation is required")
            process_csv_file(args.xsv, args.output, args.gene_format, args.annotation, args.sample_id)
        elif args.format == 'h5ad':
            if not args.h5ad:
                parser.error("For h5ad format, --h5ad is required")
            if args.gene_format == 'gene symbol' and not args.annotation:
                parser.error("For gene format 'gene symbol', --annotation is required")
            process_h5ad_file(args.h5ad, args.output, args.sample_name, args.sample_id, args.sample_column_name, args.gene_format, args.annotation)
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
