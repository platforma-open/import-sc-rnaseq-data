import polars as pl
import argparse
import os

def filter_genes(raw_counts_path, normalized_counts_path, output_raw_path, output_normalized_path):
    """
    Filters genes with zero counts across all cells and samples from count matrices.
    Gene filtering is based on the raw counts file and applied to both raw and normalized counts.

    Args:
        raw_counts_path (str): Path to the raw counts CSV file.
        normalized_counts_path (str): Path to the normalized counts CSV file.
        output_raw_path (str): Path to save the filtered raw counts CSV.
        output_normalized_path (str): Path to save the filtered normalized counts CSV.
    """
    # Load data
    raw_counts_df = pl.read_csv(raw_counts_path)
    normalized_counts_df = pl.read_csv(normalized_counts_path)

    # Log initial counts
    initial_gene_count = raw_counts_df['Ensembl Id'].n_unique()
    initial_cell_count = raw_counts_df.select(pl.concat_str([pl.col('Sample'), pl.col('Cell ID')], separator='_')).n_unique()
    print(f"Total number of unique cells: {initial_cell_count}")
    print(f"Total number of genes: {initial_gene_count}")

    # Identify genes with zero counts in all cells across all samples, based on raw counts
    gene_sums = raw_counts_df.group_by('Ensembl Id').agg(pl.sum('Raw gene expression').alias('TotalCount'))
    genes_to_drop = gene_sums.filter(pl.col('TotalCount') == 0).select('Ensembl Id')
    
    num_genes_to_drop = len(genes_to_drop)
    print(f"Number of genes with zero counts across all cells (based on raw counts): {num_genes_to_drop}")

    if num_genes_to_drop > 0:
        # Filter out genes with zero counts from both dataframes
        filtered_raw_counts_df = raw_counts_df.join(genes_to_drop, on='Ensembl Id', how='anti')
        filtered_normalized_counts_df = normalized_counts_df.join(genes_to_drop, on='Ensembl Id', how='anti')
    else:
        filtered_raw_counts_df = raw_counts_df
        filtered_normalized_counts_df = normalized_counts_df

    # Log final gene counts
    final_gene_count = filtered_raw_counts_df['Ensembl Id'].n_unique()
    print(f"Number of genes after filtering out zero-count genes: {final_gene_count}")

    # Save filtered data
    if os.path.isdir(output_raw_path):
        output_raw_path = os.path.join(output_raw_path, 'filtered_raw_counts.csv')
    if os.path.isdir(output_normalized_path):
        output_normalized_path = os.path.join(output_normalized_path, 'filtered_normalized_counts.csv')
    
    filtered_raw_counts_df.write_csv(output_raw_path)
    filtered_normalized_counts_df.write_csv(output_normalized_path)

def main():
    parser = argparse.ArgumentParser(description='Filter genes with zero counts from raw and normalized count matrices.')
    parser.add_argument('--raw_counts', type=str, required=True, help='Path to raw counts CSV file.')
    parser.add_argument('--normalized_counts', type=str, required=True, help='Path to normalized counts CSV file.')
    parser.add_argument('--output_raw', type=str, required=True, help='Path to save filtered raw counts CSV.')
    parser.add_argument('--output_normalized', type=str, required=True, help='Path to save filtered normalized counts CSV.')

    args = parser.parse_args()

    filter_genes(args.raw_counts, args.normalized_counts, args.output_raw, args.output_normalized)

if __name__ == "__main__":
    main()