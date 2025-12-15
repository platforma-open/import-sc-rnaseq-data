#!/usr/bin/env python3
"""
File Format Check Script for Single-Cell RNA-seq Data

This script analyzes count matrices from single-cell RNA-seq data to ensure
they contain recognizable gene and cell identifiers.

Usage:
    python check_format.py <input_file>
    
Output:
    Prints a success message if the format is valid, otherwise an error.
"""

import argparse
import sys
import re
from pathlib import Path
from typing import List
import polars as pl


class FileFormatError(Exception):
    """Custom exception for file format errors."""
    pass


class CountMatrixAnalyzer:
    """Analyzes count matrices to determine structure and extract gene identifiers."""
    
    def __init__(self, file_path: str):
        self.file_path = Path(file_path)
        self.df = None
        self.gene_orientation = None
        self.gene_identifiers = []
    
    def analyze(self) -> List[str]:
        """
        Analyze the count matrix and return gene identifiers.
        
        Returns:
            List of gene identifiers
        """
        self._load_data()
        self._check_column_types()
        self._determine_gene_orientation()
        return self.gene_identifiers
    
    def _load_data(self):
        """Load the data file using polars."""
        if not self.file_path.exists():
            raise FileFormatError(f"1\tFile not found: {self.file_path}")
        
        # Determine file format
        if self.file_path.suffix.lower() == '.csv':
            self.df = pl.read_csv(self.file_path, has_header=True)
        elif self.file_path.suffix.lower() in ['.tsv', '.txt']:
            self.df = pl.read_csv(self.file_path, separator='\t', has_header=True)
        else:
            # Try CSV first, then TSV
            try:
                self.df = pl.read_csv(self.file_path, has_header=True)
            except Exception:
                try:
                    self.df = pl.read_csv(self.file_path, separator='\t', has_header=True)
                except Exception as e:
                    raise FileFormatError(f"2\tCould not read one of the provided files: {e}. File format different from CSV or TSV")
    
    def _check_column_types(self):
        """
        Check value types of loaded columns.
        For non-numeric columns, check if values can be handled by scanpy.
        Scanpy can handle: numeric strings, "NaN"/"NA"/"inf" (missing values).
        Raise error if non-numeric column contains strings that scanpy cannot convert.
        """
        if self.df is None:
            raise FileFormatError("3\tData not loaded")
        
        # Values that scanpy can handle even if they're strings
        # These are missing value representations that pandas/scanpy recognizes
        scanpy_handled_missing = {
            "NaN", "nan", "NAN", "Na", "na",
            "NA", "N/A", "n/a",
            "inf", "Inf", "INF", "-inf", "-Inf", "-INF",
            "infinity", "Infinity", "INFINITY",
            "",  # Empty string is handled as missing
        }

        # Use tuple for faster membership testing (slightly more efficient than list)
        numeric_types = (
            pl.Int8, pl.Int16, pl.Int32, pl.Int64,
            pl.UInt8, pl.UInt16, pl.UInt32, pl.UInt64,
            pl.Float32, pl.Float64
        )
        
        # Skip first column (typically gene/cell identifiers)
        for col_name in self.df.columns[1:]:
            col = self.df[col_name]
            col_dtype = col.dtype
            
            # Check if column is numeric
            if col_dtype not in numeric_types:
                # Get first 10 distinct non-null values in one chain
                distinct_values = (
                    col.drop_nulls()
                    .unique()
                    .head(10)
                    .to_list()
                )
                
                # Check each value to see if scanpy can handle it
                for val in distinct_values:
                    # Convert to string for checking
                    val_str = str(val).strip()
                    
                    # check scanpy-handled missing values first
                    if val_str in scanpy_handled_missing:
                        continue
                    
                    # Check if it's a numeric string (scanpy can convert)
                    try:
                        float(val_str)
                        # If conversion succeeds, it's numeric - scanpy can handle it
                        continue
                    except (ValueError, TypeError):
                        # Not numeric and not in missing set - this will break scanpy
                        raise FileFormatError(
                            f"6\tColumn '{col_name}' contains a non-numeric string value '{val_str}' "
                            f"that scanpy cannot handle. The input file does not appear to be a correctly formatted "
                            f"count matrix. Please check the file format and try again."
                        )
    
    def _determine_gene_orientation(self):
        """Determine if genes are in rows or columns."""
        if self.df is None:
            raise FileFormatError("3\tData not loaded")
                
        # Get first column name and first few values
        first_col = self.df.columns[0]
        first_col_values = self.df.select(pl.col(first_col)).head(10).to_series().to_list()
        
        # Check if first column contains gene identifiers vs cell barcodes
        # Gene identifiers: short names like "ACT1", "TP53", "YDL247W-A"
        # Cell barcodes: long strings like "AAACCTGTCGGAAACG-2"
        
        gene_like_count = 0
        cell_barcode_like_count = 0
        
        for val in first_col_values:
            if isinstance(val, str):
                # Check for cell barcode patterns (long strings with dashes and numbers)
                if len(val) > 15 and '-' in val and any(c.isdigit() for c in val):
                    cell_barcode_like_count += 1
                # Check for generic cell identifiers (cell1, cell2, etc.)
                elif val.startswith('cell') and val[4:].isdigit():
                    cell_barcode_like_count += 1
                # Check for gene identifier patterns (shorter, more gene-like)
                elif len(val) <= 15 and not val.replace('.', '').replace('-', '').isdigit():
                    gene_like_count += 1
        
        # If we have more cell barcode-like patterns, genes are likely in columns
        if cell_barcode_like_count > gene_like_count:
            self.gene_orientation = 'columns'
        elif cell_barcode_like_count + gene_like_count == 0:
            raise FileFormatError("4\tNo gene identifiers were found in the first column or row, which indicates that the format of at least one input file is not suitable for this block.")
        else:
            # Check if first row contains gene identifiers
            first_row = self.df.head(1).to_pandas().iloc[0].tolist()
            gene_like_count_row = sum(1 for val in first_row[1:]  # Skip first element
                                    if isinstance(val, str) and len(val) <= 15 and not val.replace('.', '').replace('-', '').isdigit())
            
            if gene_like_count_row >= len(first_row) * 0.7:
                self.gene_orientation = 'columns'
            else:
                # Default to rows if uncertain
                self.gene_orientation = 'rows'    

def main():
    """Main function to handle command line arguments and execute the analysis."""
    parser = argparse.ArgumentParser(
        description="Check the format of a single-cell RNA-seq count matrix.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python check_format.py data.csv
    python check_format.py data.tsv
        """
    )
    
    parser.add_argument(
        'input_file',
        help='Path to the input CSV or TSV file containing the count matrix'
    )
    parser.add_argument(
        '--output', '-o',
        default='check_error.tsv',
        help='Output file name (default: check_error.tsv)'
    )

    args = parser.parse_args()
    
    try:
        # Analyze the count matrix to check format
        analyzer = CountMatrixAnalyzer(args.input_file)
        _ = analyzer.analyze()
        
        # A successful analysis implies a valid format based on the heuristics.
        print(f"File format check passed for: {args.input_file}")
        print(f"Detected gene orientation: genes in {analyzer.gene_orientation}")

        # Store empty log file
        with open(args.output, 'w') as f:
            f.write("Index\tLog\n")
        
    except FileFormatError as e:
        print(f"Error: {e}", file=sys.stderr)
        # Store error message in output file
        with open(args.output, 'w') as f:
            f.write(f"Index\tLog\n{e}")
    except Exception as e:
        print(f"Data structure not correct: {e}", file=sys.stderr)
        # Store error message in output file
        with open(args.output, 'w') as f:
            f.write(f"Index\tLog\n5\tData structure not correct\n")


if __name__ == '__main__':
    main()
