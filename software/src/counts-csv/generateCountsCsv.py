import argparse
import pandas as pd
import scipy.io
import gzip
import scanpy as sc
import anndata 

def read_gzip_tsv(file_path):
    """Reads a gzipped TSV file into a pandas DataFrame."""
    with gzip.open(file_path, 'rt') as f:
        return pd.read_csv(f, sep='\t', header=None)

def clean_barcode_suffix(barcode):
    """Remove '-' and following numbers from cell barcode."""
    if '-' in barcode:
        return barcode.split('-')[0]
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

def main():
    parser = argparse.ArgumentParser(description="Convert .mtx.gz, .tsv.gz files into a count matrix CSV.")
    parser.add_argument('--matrix', required=True, help="Path to the matrix.mtx.gz file")
    parser.add_argument('--barcodes', required=True, help="Path to the barcodes.tsv.gz file")
    parser.add_argument('--features', required=True, help="Path to the features.tsv.gz file")
    parser.add_argument('--output', required=True, help="Path to output the raw CSV file")

    args = parser.parse_args()
    process_input_files(args.matrix, args.barcodes, args.features, args.output)

if __name__ == "__main__":
    main()
