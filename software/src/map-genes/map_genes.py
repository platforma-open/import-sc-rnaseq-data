import polars as pl
import argparse

SAMPLE_COLUMN_NAME = "SampleId"
CELL_COLUMN_NAME = "CellId"
GENE_COLUMN_NAME = "GeneId"
COUNT_COLUMN_NAME = "Count"

def map_ensembl_to_gene_symbol(raw_counts_paths, annotation_path, output_path):
    # Load only the Ensembl Id column from raw count data (much more efficient)
    print(f"Loading unique Ensembl IDs from {len(raw_counts_paths)} raw count file(s)...")
    
    # Collect unique Ensembl IDs from all files
    all_ensembl_ids = []
    
    for i, raw_counts_path in enumerate(raw_counts_paths):
        print(f"  Processing file {i+1}/{len(raw_counts_paths)}: {raw_counts_path}")
        
        # First check if file has headers we need
        header_df = pl.read_csv(raw_counts_path, n_rows=1)
        
        # Validate required columns in raw count data (support both legacy and new headers)
        base_required = {SAMPLE_COLUMN_NAME, GENE_COLUMN_NAME, COUNT_COLUMN_NAME}
        cell_headers = {"Cell Barcode", "Cell ID", CELL_COLUMN_NAME}
        missing_base = base_required - set(header_df.columns)
        has_cell_header = any(h in header_df.columns for h in cell_headers)
        if missing_base or not has_cell_header:
            expected_desc = f"{sorted(base_required)} and one of {sorted(cell_headers)}"
            raise ValueError(f"Raw count data must contain columns: {expected_desc}")
        
        # Extract unique Ensembl IDs - only read the column we need
        ensembl_ids = (pl.read_csv(raw_counts_path, columns=[GENE_COLUMN_NAME])
                      .select(pl.col(GENE_COLUMN_NAME).unique())
                      .to_series()
                      .to_list())
        
        all_ensembl_ids.extend(ensembl_ids)
        print(f"    Found {len(ensembl_ids)} unique Ensembl IDs in this file")
    
    # Get overall unique Ensembl IDs
    unique_ensembl_ids = list(set(all_ensembl_ids))
    print(f"Found {len(unique_ensembl_ids)} total unique Ensembl IDs across all files.")

    # Load the annotation data
    print("Loading annotation data...")
    annotation_df = pl.read_csv(annotation_path)

    # Validate required columns in annotation data
    if not {"Ensembl Id", "Gene symbol"}.issubset(annotation_df.columns):
        raise ValueError("Annotation file must contain 'Ensembl Id' and 'Gene symbol' columns")

    # Create DataFrame with unique Ensembl IDs
    mapped_df = pl.DataFrame({"Ensembl Id": unique_ensembl_ids})
    
    # Use join instead of map for better performance
    print("Mapping Ensembl IDs to gene symbols...")
    annotation_unique = annotation_df.unique(subset=["Ensembl Id"], keep="first")
    
    # Left join to get gene symbols
    mapped_df = mapped_df.join(
        annotation_unique.select(["Ensembl Id", "Gene symbol"]),
        on="Ensembl Id",
        how="left"
    )

    # Report missing mappings
    missing = mapped_df.filter(pl.col("Gene symbol").is_null()).height
    print(f"Mapping complete. {missing} Ensembl IDs could not be mapped to gene symbols.")

    # Add Ensembl ID as label when gene symbol is missing
    mapped_df = mapped_df.with_columns(
        pl.when(pl.col("Gene symbol").is_null())
        .then(pl.col("Ensembl Id"))
        .otherwise(pl.col("Gene symbol"))
        .alias("Gene symbol")
    )

    # Save output
    print(f"Saving output to {output_path}...")
    mapped_df.write_csv(output_path)
    print("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Map Ensembl IDs in single-cell RNA-seq count data to gene symbols.")
    parser.add_argument("--raw_counts", required=True, action='append', dest='raw_counts_files',
                       help="Path to raw count CSV file (can be specified multiple times)")
    parser.add_argument("--annotation", required=True, help="Path to annotation CSV file (e.g. mus_musculus_gene_annotations.csv)")
    parser.add_argument("--output", required=True, help="Path to output CSV file with Ensembl Id and Gene symbol")

    args = parser.parse_args()
    map_ensembl_to_gene_symbol(args.raw_counts_files, args.annotation, args.output)
