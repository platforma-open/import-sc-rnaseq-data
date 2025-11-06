import anndata as ad
import pandas as pd
import sys
import json
from pathlib import Path
import os


def find_sample_column(obs_df, expected_samples):
    """
    Finds the column in the .obs DataFrame that contains the sample names.
    """
    expected_samples_set = set(expected_samples)
    
    # Prioritize columns with 'sample' in the name, then check all others
    candidate_columns = [col for col in obs_df.columns if 'sample' in col.lower()]
    candidate_columns.extend([col for col in obs_df.columns if col not in candidate_columns])

    for col_name in candidate_columns:
        # Check only string or categorical columns
        if obs_df[col_name].dtype.name not in ['category', 'object']:
            continue
        
        unique_values = set(obs_df[col_name].unique())
        
        # Check if all expected samples are present in the column's unique values
        if expected_samples_set.issubset(unique_values):
            print(f"Found sample column: '{col_name}'")
            return col_name
            
    return None


def main():
    if len(sys.argv) < 2:
        print("Usage: python main.py <input.h5ad>", file=sys.stderr)
        sys.exit(1)

    input_h5ad = sys.argv[1]
    
    try:
        with open('sampleGroups.json', 'r') as f:
            sample_groups = json.load(f)
    except FileNotFoundError:
        print("Error: sampleGroups.json not found.", file=sys.stderr)
        sys.exit(1)
    except json.JSONDecodeError:
        print("Error: Could not decode sampleGroups.json.", file=sys.stderr)
        sys.exit(1)

    adata = ad.read_h5ad(input_h5ad)

    # Assuming a single sample group from the input file
    group_id = list(sample_groups.keys())[0]
    sample_mapping = sample_groups[group_id]
    expected_sample_names = list(sample_mapping.values())

    sample_column = find_sample_column(adata.obs, expected_sample_names)

    if sample_column is None:
        print("Error: Could not automatically identify the sample column in the h5ad file's .obs dataframe.", file=sys.stderr)
        print(f"Please ensure one of the columns contains the following sample names: {expected_sample_names}", file=sys.stderr)
        sys.exit(1)

    output_dir = Path("output_h5ads")
    output_dir.mkdir(exist_ok=True)

    linker_data = {"samples": []}

    for sample_id, sample_name in sample_mapping.items():
        # Filter AnnData for the current sample
        sample_adata = adata[adata.obs[sample_column] == sample_name].copy()

        if sample_adata.n_obs == 0:
            print(f"Warning: No cells found for sample '{sample_name}'. Skipping.", file=sys.stderr)
            continue

        # Define output path and save the new h5ad file
        output_path = output_dir / f"{sample_id}.h5ad"
        sample_adata.write_h5ad(output_path)

        # Add entry to linker data
        linker_data["samples"].append({
            "sampleId": sample_id,
            "filePath": str(output_path),
            "sampleName": sample_name
        })

    with open("linker.json", "w") as f:
        json.dump(linker_data, f)

    print("Successfully split h5ad file.")


if __name__ == "__main__":
    main()
