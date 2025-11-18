#!/usr/bin/env python3
"""
Calculate basic metrics per sample from long-format count data.

This script processes a long-format CSV file with columns:
- Sample: Sample identifier
- Cell ID: Cell barcode identifier  
- Ensembl Id: Gene identifier
- Raw gene expression: Count value

And calculates various metrics per sample including:
- Number of cells
- Mean/median reads per cell
- Number of genes detected
- Mean/median genes per cell
- Complexity (log10(genes)/log10(UMIs))
- Total reads
- Mean/median reads per gene
- Fraction of cells with >0 reads
- Fraction of genes with >0 reads
"""

import polars as pl
import numpy as np
import argparse
from pathlib import Path
import sys
from typing import Dict, Any, List

SAMPLE_COLUMN_NAME = "SampleId"
CELL_COLUMN_NAME = "CellId"
GENE_COLUMN_NAME = "GeneId"
COUNT_COLUMN_NAME = "Count"

def calculate_sample_metrics(df: pl.DataFrame) -> Dict[str, Any]:
    """
    Calculate basic metrics for a single sample using Polars.
    
    Args:
        df: Polars DataFrame with columns ['SampleId', 'CellId', 'GeneId', 'Count']
        
    Returns:
        Dictionary with calculated metrics
    """
    # Remove header row if present
    df = df.filter(pl.col(SAMPLE_COLUMN_NAME) != SAMPLE_COLUMN_NAME)
    
    # Basic counts
    n_cells = df.select(pl.col(CELL_COLUMN_NAME).n_unique()).item()
    n_genes = df.select(pl.col(GENE_COLUMN_NAME).n_unique()).item()
    total_umis = df.select(pl.col(COUNT_COLUMN_NAME).sum()).item()
    
    # Calculate UMIs per cell
    umis_per_cell = (df.group_by(CELL_COLUMN_NAME)
                    .agg(pl.col(COUNT_COLUMN_NAME).sum())
                    .select(COUNT_COLUMN_NAME))
    
    mean_umis_per_cell = umis_per_cell.select(pl.col(COUNT_COLUMN_NAME).mean()).item()
    median_umis_per_cell = umis_per_cell.select(pl.col(COUNT_COLUMN_NAME).median()).item()
    
    # Calculate genes per cell (count non-zero values per cell)
    # Use a more robust approach: count distinct genes per cell
    genes_per_cell = (df.filter(pl.col(COUNT_COLUMN_NAME) > 0)
                     .group_by(CELL_COLUMN_NAME)
                     .agg(pl.col(GENE_COLUMN_NAME).n_unique().alias(COUNT_COLUMN_NAME))
                     .select(COUNT_COLUMN_NAME))
    
    mean_genes_per_cell = genes_per_cell.select(pl.col(COUNT_COLUMN_NAME).mean()).item()
    median_genes_per_cell = genes_per_cell.select(pl.col(COUNT_COLUMN_NAME).median()).item()
    
    return {
        'n_cells': n_cells,
        'n_genes': n_genes,
        'total_umis': total_umis,
        'mean_umis_per_cell': mean_umis_per_cell,
        'median_umis_per_cell': median_umis_per_cell,
        'mean_genes_per_cell': mean_genes_per_cell,
        'median_genes_per_cell': median_genes_per_cell
    }


def process_files(input_files: List[str], output_file: str) -> None:
    """
    Process multiple input files and calculate metrics for each sample using Polars.
    
    Args:
        input_files: List of paths to input CSV files
        output_file: Path to output CSV file
    """
    print(f"Loading data from {len(input_files)} file(s)...")
    
    # Read all CSV files and concatenate them
    dataframes = []
    for input_file in input_files:
        input_path = Path(input_file)
        if not input_path.exists():
            print(f"Warning: Input file {input_file} does not exist, skipping...")
            continue
        
        print(f"  Reading {input_file}...")
        try:
            df = pl.read_csv(input_file)
            dataframes.append(df)
            print(f"    Loaded {len(df)} rows")
        except Exception as e:
            print(f"  Error reading {input_file}: {e}")
            continue
    
    if not dataframes:
        print("Error: No valid CSV files could be read!")
        sys.exit(1)
    
    # Concatenate all dataframes
    print("Concatenating dataframes...")
    df = pl.concat(dataframes, how="vertical")
    print(f"Total rows: {len(df)}")
    
    # Get all unique samples
    print("Identifying samples...")
    samples = df.select(pl.col(SAMPLE_COLUMN_NAME).unique()).to_series().to_list()
    
    # Remove 'Sample' header if present
    if SAMPLE_COLUMN_NAME in samples:
        samples.remove(SAMPLE_COLUMN_NAME)
    
    print(f"Found {len(samples)} samples: {samples}")
    
    all_metrics = []
    
    # Process each sample
    for sample in samples:
        print(f"Processing sample: {sample}")
        
        # Filter data for this sample
        sample_df = df.filter(pl.col(SAMPLE_COLUMN_NAME) == sample)
        
        if not sample_df.is_empty():
            metrics = calculate_sample_metrics(sample_df)
            metrics['sample'] = sample
            all_metrics.append(metrics)
            
            print(f"  - Cells: {metrics['n_cells']}")
            print(f"  - Genes: {metrics['n_genes']}")
            print(f"  - Total UMIs: {metrics['total_umis']:,.0f}")
            print(f"  - Mean UMIs/cell: {metrics['mean_umis_per_cell']:.1f}")
            print(f"  - Mean genes/cell: {metrics['mean_genes_per_cell']:.1f}")
    
    # Convert to DataFrame and save
    if all_metrics:
        metrics_df = pl.DataFrame(all_metrics)
        
        # Reorder columns to match requested order
        column_order = ['sample', 'n_cells', 'n_genes', 'total_umis', 
                       'mean_umis_per_cell', 'median_umis_per_cell',
                       'mean_genes_per_cell', 'median_genes_per_cell']
        
        metrics_df = metrics_df.select(column_order)
        
        # Round numeric columns
        metrics_df = metrics_df.with_columns([
            pl.col(col).round(3) for col in metrics_df.columns 
            if col != 'sample' and metrics_df[col].dtype in [pl.Float64, pl.Float32]
        ])
        
        print(f"\nSaving metrics to {output_file}...")
        metrics_df.write_csv(output_file)
        print(f"Successfully calculated metrics for {len(all_metrics)} samples")
        
        # Print summary
        print("\nSummary:")
        summary_cols = ['sample', 'n_cells', 'total_umis', 'mean_umis_per_cell', 'mean_genes_per_cell']
        print(metrics_df.select(summary_cols))
        
    else:
        print("No valid samples found!")


def main():
    parser = argparse.ArgumentParser(
        description="Calculate basic metrics per sample from long-format count data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python calculate_sample_metrics.py -i rawCounts1.csv -i rawCounts2.csv -o sample_metrics.csv
  
The input file(s) should have columns:
- SampleId: Sample identifier
- CellId: Cell barcode identifier  
- GeneId: Gene identifier
- Count: Count value

Output metrics include:
- n_cells: Number of cells
- n_genes: Number of genes detected
- total_umis: Total UMI counts across all cells (NOT raw reads)
- mean_umis_per_cell: Average UMI counts per cell
- median_umis_per_cell: Median UMI counts per cell
- mean_genes_per_cell: Average genes detected per cell
- median_genes_per_cell: Median genes detected per cell
        """
    )
    
    parser.add_argument('-i', '--input', required=True, action='append', dest='inputs',
                       help='Input CSV file(s) with long-format count data (can be specified multiple times)')
    parser.add_argument('-o', '--output', required=True,
                       help='Output CSV file with sample metrics')
    
    args = parser.parse_args()
    
    # Check input files exist
    for input_file in args.inputs:
        input_path = Path(input_file)
        if not input_path.exists():
            print(f"Error: Input file {input_file} does not exist!")
            sys.exit(1)
    
    # Create output directory if needed
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    try:
        process_files(args.inputs, args.output)
    except Exception as e:
        print(f"Error processing file: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
