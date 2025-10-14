import polars as pl
import argparse
import time
import os
import polars.selectors as cs


def filter_outliers(raw_counts_path, normalized_counts_path, metrics_path, output_raw_path, output_normalized_path):
    """
    Filters outlier cells from count matrices based on a metrics file.

    Args:
        raw_counts_path (str): Path to the raw counts CSV file.
        normalized_counts_path (str): Path to the normalized counts CSV file.
        metrics_path (str): Path to the cell metrics CSV file containing outlier flags.
        output_raw_path (str): Path to save the filtered raw counts CSV.
        output_normalized_path (str): Path to save the filtereThed normalized counts CSV.
    """
    # Load data
    metrics_df = pl.read_csv(metrics_path)
    raw_counts_long = pl.read_csv(raw_counts_path)
    normalized_counts_long = pl.read_csv(normalized_counts_path)

    # Identify outliers
    outlier_cells_df = metrics_df.filter(pl.col('outlier') == True).select('CellId')

    # Log initial counts
    initial_cell_count = raw_counts_long['CellId'].n_unique()
    outlier_cell_count = len(outlier_cells_df)
    print(f"Total number of cells: {initial_cell_count}")
    print(f"Number of cells flagged as outliers: {outlier_cell_count}")

    # Filter out outliers
    final_raw_long = raw_counts_long.join(outlier_cells_df, on='CellId', how='anti')
    final_normalized_long = normalized_counts_long.join(outlier_cells_df, on='CellId', how='anti')

    # Log final counts
    final_cell_count = final_raw_long['CellId'].n_unique()
    print(f"Number of cells after filtering: {final_cell_count}")

    # Save filtered data
    if os.path.isdir(output_raw_path):
        output_raw_path = os.path.join(output_raw_path, 'filtered_raw_counts.csv')
    if os.path.isdir(output_normalized_path):
        output_normalized_path = os.path.join(output_normalized_path, 'filtered_normalized_counts.csv')
    final_raw_long.write_csv(output_raw_path)
    final_normalized_long.write_csv(output_normalized_path)

def main():
    parser = argparse.ArgumentParser(description='Filter outlier cells from count matrices.')
    parser.add_argument('--raw_counts', type=str, required=True, help='Path to raw counts CSV file.')
    parser.add_argument('--normalized_counts', type=str, required=True, help='Path to normalized counts CSV file.')
    parser.add_argument('--metrics', type=str, required=True, help='Path to cell metrics CSV file.')
    parser.add_argument('--output_raw', type=str, required=True, help='Path to save filtered raw counts CSV.')
    parser.add_argument('--output_normalized', type=str, required=True, help='Path to save filtered normalized counts CSV.')

    args = parser.parse_args()

    filter_outliers(args.raw_counts, args.normalized_counts, args.metrics, args.output_raw, args.output_normalized)

if __name__ == "__main__":
    main()