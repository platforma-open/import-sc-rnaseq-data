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
from typing import Dict, Any


def calculate_sample_metrics(df: pl.DataFrame) -> Dict[str, Any]:
    """
    Calculate basic metrics for a single sample using Polars.
    
    Args:
        df: Polars DataFrame with columns ['Sample', 'Cell ID', 'Ensembl Id', 'Raw gene expression']
        
    Returns:
        Dictionary with calculated metrics
    """
    # Remove header row if present
    df = df.filter(pl.col("Sample") != "Sample")
    
    # Basic counts
    n_cells = df.select(pl.col("Cell ID").n_unique()).item()
    n_genes = df.select(pl.col("Ensembl Id").n_unique()).item()
    total_umis = df.select(pl.col("Raw gene expression").sum()).item()
    
    # Calculate UMIs per cell
    umis_per_cell = (df.group_by("Cell ID")
                    .agg(pl.col("Raw gene expression").sum())
                    .select("Raw gene expression"))
    
    mean_umis_per_cell = umis_per_cell.select(pl.col("Raw gene expression").mean()).item()
    median_umis_per_cell = umis_per_cell.select(pl.col("Raw gene expression").median()).item()
    
    # Calculate genes per cell (count non-zero values per cell)
    # Use a more robust approach: count distinct genes per cell
    genes_per_cell = (df.filter(pl.col("Raw gene expression") > 0)
                     .group_by("Cell ID")
                     .agg(pl.col("Ensembl Id").n_unique().alias("Raw gene expression"))
                     .select("Raw gene expression"))
    
    mean_genes_per_cell = genes_per_cell.select(pl.col("Raw gene expression").mean()).item()
    median_genes_per_cell = genes_per_cell.select(pl.col("Raw gene expression").median()).item()
    
    return {
        'n_cells': n_cells,
        'n_genes': n_genes,
        'total_umis': total_umis,
        'mean_umis_per_cell': mean_umis_per_cell,
        'median_umis_per_cell': median_umis_per_cell,
        'mean_genes_per_cell': mean_genes_per_cell,
        'median_genes_per_cell': median_genes_per_cell
    }


def process_file(input_file: str, output_file: str) -> None:
    """
    Process the input file and calculate metrics for each sample using Polars.
    
    Args:
        input_file: Path to input CSV file
        output_file: Path to output CSV file
    """
    print(f"Loading data from {input_file}...")
    
    # Read the entire file with Polars (more efficient for this use case)
    df = pl.read_csv(input_file)
    
    # Get all unique samples
    print("Identifying samples...")
    samples = df.select(pl.col("Sample").unique()).to_series().to_list()
    
    # Remove 'Sample' header if present
    if 'Sample' in samples:
        samples.remove('Sample')
    
    print(f"Found {len(samples)} samples: {samples}")
    
    all_metrics = []
    
    # Process each sample
    for sample in samples:
        print(f"Processing sample: {sample}")
        
        # Filter data for this sample
        sample_df = df.filter(pl.col("Sample") == sample)
        
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
  python calculate_sample_metrics.py -i rawCounts.csv -o sample_metrics.csv
  
The input file should have columns:
- Sample: Sample identifier
- Cell ID: Cell barcode identifier  
- Ensembl Id: Gene identifier
- Raw gene expression: Count value

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
    
    parser.add_argument('-i', '--input', required=True, 
                       help='Input CSV file with long-format count data')
    parser.add_argument('-o', '--output', required=True,
                       help='Output CSV file with sample metrics')
    
    args = parser.parse_args()
    
    # Check input file exists
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: Input file {args.input} does not exist!")
        sys.exit(1)
    
    # Create output directory if needed
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    try:
        process_file(args.input, args.output)
    except Exception as e:
        print(f"Error processing file: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
