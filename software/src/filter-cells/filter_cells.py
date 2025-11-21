import polars as pl
import argparse
import time
import os
import sys
import polars.selectors as cs
from datetime import datetime


def write_csv_chunked(pl_df, output_path, chunk_size=500000):
    """Write Polars DataFrame to CSV using optimized chunked approach.
    
    Args:
        pl_df: Polars DataFrame to write
        output_path: Path to output CSV file (uncompressed)
        chunk_size: Number of rows to write per chunk (default: 500000)
    """
    start_time = datetime.now()
    total_rows = len(pl_df)
    
    if total_rows == 0:
        # Write empty DataFrame with headers only
        pl_df.write_csv(output_path)
        print(f"[{start_time.strftime('%Y-%m-%d %H:%M:%S')}] Written empty DataFrame with headers to {output_path}")
        return
    
    print(f"[{start_time.strftime('%Y-%m-%d %H:%M:%S')}] Writing CSV in chunks of {chunk_size:,} rows...")
    
    # Write in chunks using file handle
    if total_rows <= chunk_size:
        # Small enough to write in one go
        pl_df.write_csv(output_path)
    else:
        # Write in chunks using Polars directly with file handle
        with open(output_path, 'w') as f:
            # Write header
            f.write(",".join(pl_df.columns) + "\n")
            for i in range(0, len(pl_df), chunk_size):
                batch = pl_df.slice(i, chunk_size)
                batch.write_csv(f, include_header=False)
                written_rows = min(i + chunk_size, total_rows)
                if (i // chunk_size) % 10 == 0 or written_rows == total_rows:
                    current_time = datetime.now()
                    elapsed = (current_time - start_time).total_seconds()
                    print(f"[{current_time.strftime('%Y-%m-%d %H:%M:%S')}] Written {written_rows:,} / {total_rows:,} rows ({100 * written_rows / total_rows:.1f}%) - Elapsed: {elapsed:.1f}s")
        
    end_time = datetime.now()
    elapsed_time = (end_time - start_time).total_seconds()
    print(f"[{end_time.strftime('%Y-%m-%d %H:%M:%S')}] Completed writing {total_rows:,} rows to {output_path} - Total time: {elapsed_time:.1f}s")


def filter_outliers(raw_counts_path, normalized_counts_path, metrics_path, output_raw_path, output_normalized_path, sample_id=None):
    """
    Filters outlier cells from count matrices based on a metrics file.
    Supports Parquet and CSV files as input.

    Args:
        raw_counts_path (str): Path to the raw counts Parquet or CSV file.
        normalized_counts_path (str): Path to the normalized counts Parquet or CSV file.
        metrics_path (str): Path to the cell metrics CSV file containing outlier flags.
        output_raw_path (str): Path to save the filtered raw counts CSV.
        output_normalized_path (str): Path to save the filtered normalized counts CSV.
        sample_id (str): Optional sample ID to add as first column in output CSV.
    """
    # Load data - use Parquet if available, otherwise CSV
    metrics_df = pl.read_csv(metrics_path)
    if raw_counts_path.endswith('.parquet'):
        raw_counts_long = pl.read_parquet(raw_counts_path)
    else:
        raw_counts_long = pl.read_csv(raw_counts_path)
    if normalized_counts_path.endswith('.parquet'):
        normalized_counts_long = pl.read_parquet(normalized_counts_path)
    else:
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

    # Add SampleId as first column if provided
    if sample_id:
        if 'SampleId' not in final_raw_long.columns:
            final_raw_long = final_raw_long.with_columns(pl.lit(sample_id).alias('SampleId'))
            final_raw_long = final_raw_long.select(['SampleId'] + [col for col in final_raw_long.columns if col != 'SampleId'])
        
        if 'SampleId' not in final_normalized_long.columns:
            final_normalized_long = final_normalized_long.with_columns(pl.lit(sample_id).alias('SampleId'))
            final_normalized_long = final_normalized_long.select(['SampleId'] + [col for col in final_normalized_long.columns if col != 'SampleId'])

    # Save filtered data using optimized chunked approach
    if os.path.isdir(output_raw_path):
        output_raw_path = os.path.join(output_raw_path, 'filtered_raw_counts.csv')
    if os.path.isdir(output_normalized_path):
        output_normalized_path = os.path.join(output_normalized_path, 'filtered_normalized_counts.csv')
    
    print(f"Writing filtered raw counts to {output_raw_path}...")
    write_csv_chunked(final_raw_long, output_raw_path)
    
    print(f"Writing filtered normalized counts to {output_normalized_path}...")
    write_csv_chunked(final_normalized_long, output_normalized_path)

def main():
    # Setup timing logging
    script_name = os.path.splitext(os.path.basename(__file__))[0]
    log_file = f"{script_name}.time.log"
    start_time = datetime.now()
    
    with open(log_file, 'w') as f:
        f.write(f"Script: {script_name}\n")
        f.write(f"Start time: {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    parser = argparse.ArgumentParser(description='Filter outlier cells from count matrices.')
    parser.add_argument('--raw_counts', type=str, required=True, help='Path to raw counts Parquet or CSV file.')
    parser.add_argument('--normalized_counts', type=str, required=True, help='Path to normalized counts Parquet or CSV file.')
    parser.add_argument('--metrics', type=str, required=True, help='Path to cell metrics CSV file.')
    parser.add_argument('--sample-id', type=str, help='Sample ID to add as first column in output CSV.')
    parser.add_argument('--output_raw', type=str, required=True, help='Path to save filtered raw counts CSV.')
    parser.add_argument('--output_normalized', type=str, required=True, help='Path to save filtered normalized counts CSV.')

    args = parser.parse_args()

    try:
        filter_outliers(args.raw_counts, args.normalized_counts, args.metrics, args.output_raw, args.output_normalized, args.sample_id)
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