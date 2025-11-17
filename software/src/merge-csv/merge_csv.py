#!/usr/bin/env python3
"""
Merge multiple CSV files into a single CSV file.

This script takes multiple CSV files with the same schema and concatenates them
into a single output file using Polars for efficient processing.
"""

import polars as pl
import argparse
from pathlib import Path
import sys
from typing import List


def merge_csv_files(input_files: List[str], output_file: str) -> None:
    """
    Merge multiple CSV files into a single output file.
    
    Args:
        input_files: List of paths to input CSV files
        output_file: Path to output CSV file
    """
    if not input_files:
        print("Error: No input files provided!")
        sys.exit(1)
    
    print(f"Merging {len(input_files)} CSV files...")
    
    # Read all CSV files
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
    merged_df = pl.concat(dataframes, how="vertical")
    
    print(f"Merged dataframe contains {len(merged_df)} total rows")
    
    # Rename columns to match expected format
    column_mapping = {
        "SampleId": "Sample",
        "CellId": "Cell ID",
        "GeneId": "Ensembl Id",
        "Count": "Raw gene expression"
    }
    
    # Only rename columns that exist
    rename_dict = {old: new for old, new in column_mapping.items() if old in merged_df.columns}
    if rename_dict:
        print(f"Renaming columns: {rename_dict}")
        merged_df = merged_df.rename(rename_dict)
    
    # Create output directory if needed
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Write merged CSV
    print(f"Writing merged data to {output_file}...")
    merged_df.write_csv(output_file)
    
    print(f"Successfully merged {len(dataframes)} files into {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Merge multiple CSV files into a single CSV file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python merge_csv.py sample1.csv sample2.csv sample3.csv -o merged.csv
  
The input files should have the same schema. All files will be concatenated
vertically into a single output file.
        """
    )
    
    parser.add_argument('input_files', nargs='+', 
                       help='Input CSV files to merge')
    parser.add_argument('-o', '--output', required=True,
                       help='Output CSV file')
    
    args = parser.parse_args()
    
    try:
        merge_csv_files(args.input_files, args.output)
    except Exception as e:
        print(f"Error merging files: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

