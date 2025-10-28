import argparse
import pandas as pd
import scipy.io
import gzip
import scanpy as sc
import anndata
from pathlib import Path 

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

def process_input_files(matrix_path, barcodes_path, features_path, output_csv_path):
    # Load the input files
    print("Loading matrix.mtx.gz...")
    matrix = scipy.io.mmread(matrix_path).tocoo()  # Sparse COO format

    print("Loading barcodes.tsv.gz...")
    barcodes = read_gzip_tsv(barcodes_path)
    
    print("Loading features.tsv.gz...")
    features = read_gzip_tsv(features_path)
    
    barcodes_list = [clean_barcode_suffix(barcode) for barcode in barcodes[0].tolist()]
    print(f"Cleaned barcode suffixes. Example: {barcodes[0].tolist()[0]} -> {barcodes_list[0]}")
    features_list = features[0].tolist()

    data = []

    print(f"Processing {matrix.nnz} nonzero entries...")
    count_debug = 0

    for i, j, value in zip(matrix.row, matrix.col, matrix.data):
        if count_debug < 10:
            print(f"Row: {i}, Column: {j}, Count: {value}")
        cell_id = barcodes_list[j]
        gene_id = features_list[i]
        count = value
        data.append([cell_id, gene_id, count])
        count_debug += 1

    print(f"Total rows written: {count_debug}")
    df = pd.DataFrame(data, columns=["CellId", "GeneId", "Count"])

    print(f"Writing raw count matrix to {output_csv_path}...")
    df.to_csv(output_csv_path, index=False)

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
        norm_data.append([cell_id, gene_id, value])

    norm_df = pd.DataFrame(norm_data, columns=["CellId", "GeneId", "NormalizedCount"])
    normalized_output_csv_path = output_csv_path.replace(".csv", "_normalized.csv")

    print(f"Writing normalized count matrix to {normalized_output_csv_path}...")
    norm_df.to_csv(normalized_output_csv_path, index=False)

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
            symbol_to_ensembl[symbol] = ensembl_id
    
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

def process_csv_file(csv_path, output_csv_path, gene_format=None, annotation_path=None):
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
        gene_names = df.index.tolist()
        cell_names = df.columns.tolist()
    else:
        # Genes in rows, cells in columns - use as is
        gene_names = df.index.tolist()
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
    
    # Convert to long format
    data = []
    print(f"Processing {df.shape[0]} genes x {df.shape[1]} cells...")
    
    for i, gene in enumerate(gene_names):
        for j, cell in enumerate(cleaned_cell_names):
            count = df.iloc[i, j]
            if count > 0:  # Only store non-zero counts
                data.append([cell, gene, count])
    
    print(f"Total non-zero entries: {len(data)}")
    df_long = pd.DataFrame(data, columns=["CellId", "GeneId", "Count"])
    
    # Check for and handle duplicate CellId-GeneId combinations
    df_long = handle_duplicate_combinations(df_long)
    
    print(f"Writing raw count matrix to {output_csv_path}...")
    df_long.to_csv(output_csv_path, index=False)
    
    # Normalize counts
    print("Normalizing counts...")
    
    # Create AnnData object for normalization
    adata = anndata.AnnData(X=df.values.T)  # Transpose: cells x genes
    adata.obs_names = cleaned_cell_names
    adata.var_names = gene_names
    
    sc.pp.normalize_total(adata, target_sum=1e4)
    normalized_matrix = adata.X
    
    # Extract normalized counts
    norm_data = []
    for i, cell in enumerate(cleaned_cell_names):
        for j, gene in enumerate(gene_names):
            value = normalized_matrix[i, j]
            if value > 0:  # Only store non-zero counts
                norm_data.append([cell, gene, value])
    
    norm_df = pd.DataFrame(norm_data, columns=["CellId", "GeneId", "NormalizedCount"])
    
    # Check for and handle duplicate CellId-GeneId combinations in normalized data
    norm_df = handle_duplicate_combinations_normalized(norm_df)
    
    normalized_output_csv_path = output_csv_path.replace(".csv", "_normalized.csv")
    
    print(f"Writing normalized count matrix to {normalized_output_csv_path}...")
    norm_df.to_csv(normalized_output_csv_path, index=False)
    
    print("Done!")

def process_h5ad_file(h5ad_path, output_csv_path):
    """Process an AnnData h5ad file and convert to long format CSV."""
    print(f"Loading h5ad file: {h5ad_path}")
    
    # Load the h5ad file
    adata = sc.read_h5ad(h5ad_path)
    
    print(f"AnnData shape: {adata.shape} (cells × genes)")
    
    # Extract gene and cell identifiers
    gene_names = list(adata.var_names)
    cell_names = list(adata.obs_names)
    
    # Clean cell barcodes if needed
    cleaned_cell_names = [clean_barcode_suffix(cell) for cell in cell_names]
    if any(c != o for c, o in zip(cleaned_cell_names, cell_names)):
        print(f"Cleaned barcode suffixes. Example: {cell_names[0]} -> {cleaned_cell_names[0]}")
    
    # Get the count matrix
    X = adata.X
    
    # Convert to long format for raw counts
    print(f"Processing {X.shape[0]} cells × {X.shape[1]} genes...")
    raw_data = []
    
    if hasattr(X, 'tocoo'):
        # Sparse matrix
        for i, j, value in zip(X.tocoo().row, X.tocoo().col, X.tocoo().data):
            if value > 0:
                raw_data.append([cleaned_cell_names[i], gene_names[j], value])
    else:
        # Dense matrix
        for i, cell in enumerate(cleaned_cell_names):
            for j, gene in enumerate(gene_names):
                value = X[i, j]
                if value > 0:
                    raw_data.append([cell, gene, value])
    
    print(f"Total non-zero entries: {len(raw_data)}")
    df_raw = pd.DataFrame(raw_data, columns=["CellId", "GeneId", "Count"])
    
    # Check for and handle duplicate CellId-GeneId combinations
    df_raw = handle_duplicate_combinations(df_raw)
    
    print(f"Writing raw count matrix to {output_csv_path}...")
    df_raw.to_csv(output_csv_path, index=False)
    
    # Normalize counts
    print("Normalizing counts...")
    
    # Create a copy of adata for normalization
    adata_norm = adata.copy()
    sc.pp.normalize_total(adata_norm, target_sum=1e4)
    
    # Extract normalized counts
    X_norm = adata_norm.X
    norm_data = []
    
    if hasattr(X_norm, 'tocoo'):
        # Sparse matrix
        for i, j, value in zip(X_norm.tocoo().row, X_norm.tocoo().col, X_norm.tocoo().data):
            if value > 0:
                norm_data.append([cleaned_cell_names[i], gene_names[j], value])
    else:
        # Dense matrix
        for i, cell in enumerate(cleaned_cell_names):
            for j, gene in enumerate(gene_names):
                value = X_norm[i, j]
                if value > 0:
                    norm_data.append([cell, gene, value])
    
    norm_df = pd.DataFrame(norm_data, columns=["CellId", "GeneId", "NormalizedCount"])
    
    # Check for and handle duplicate CellId-GeneId combinations in normalized data
    norm_df = handle_duplicate_combinations_normalized(norm_df)
    
    normalized_output_csv_path = output_csv_path.replace(".csv", "_normalized.csv")
    
    print(f"Writing normalized count matrix to {normalized_output_csv_path}...")
    norm_df.to_csv(normalized_output_csv_path, index=False)
    
    print("Done!")

def main():
    parser = argparse.ArgumentParser(description="Convert .mtx.gz, .tsv.gz files or CSV files into a count matrix CSV.")
    parser.add_argument('--format', required=True, choices=['mtx', 'xsv', 'h5ad'], 
                       help="Input format: 'mtx' for 10X Genomics format, 'xsv' for CSV/TSV format, 'h5ad' for AnnData format")
    parser.add_argument('--matrix', help="Path to the matrix.mtx.gz file (required for mtx format)")
    parser.add_argument('--barcodes', help="Path to the barcodes.tsv.gz file (required for mtx format)")
    parser.add_argument('--features', help="Path to the features.tsv.gz file (required for mtx format)")
    parser.add_argument('--xsv', help="Path to the XSV file (required for xsv format)")
    parser.add_argument('--h5ad', help="Path to the h5ad file (required for h5ad format)")
    parser.add_argument('--gene-format', choices=['gene symbol', 'Ensembl Id'], 
                       help="Gene identifier format: 'gene symbol' or 'Ensembl Id' (for xsv format only)")
    parser.add_argument('--annotation', help="Path to gene annotation CSV file (required when --gene-format is 'gene symbol')")
    parser.add_argument('--output', required=True, help="Path to output the raw CSV file")

    args = parser.parse_args()
    
    if args.format == 'mtx':
        if not all([args.matrix, args.barcodes, args.features]):
            parser.error("For mtx format, --matrix, --barcodes, and --features are required")
        process_input_files(args.matrix, args.barcodes, args.features, args.output)
    elif args.format == 'xsv':
        if not args.xsv:
            parser.error("For xsv format, --xsv is required")
        if args.gene_format == 'gene symbol' and not args.annotation:
            parser.error("For gene format 'gene symbol', --annotation is required")
        process_csv_file(args.xsv, args.output, args.gene_format, args.annotation)
    elif args.format == 'h5ad':
        if not args.h5ad:
            parser.error("For h5ad format, --h5ad is required")
        process_h5ad_file(args.h5ad, args.output)

if __name__ == "__main__":
    main()
