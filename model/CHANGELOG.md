# @platforma-open/milaboratories.import-sc-rnaseq-data.model

## 1.5.0

### Minor Changes

- 34f4c2e: H5 file format support

## 1.4.0

### Minor Changes

- 81fa688: Add Seurat (rds) datasets support, optimize performance

  - Added support for Seurat RDS file format with automatic import mode detection
  - Implemented multi-sample Seurat file processing with automatic sample detection and filtering
  - Created new R-based software components (extract_counts.R, infer_species.R) for Seurat object processing
  - Switched output format from CSV to Parquet with zstd compression for improved I/O performance and reduced storage footprint
  - Optimized sparse matrix handling and data.table operations for improved memory efficiency
  - Enhanced memory calculation utilities to dynamically allocate resources based on data size

## 1.3.0

### Minor Changes

- c7bd46d: h5ad format support

  - Added support for h5ad format input files (AnnData format)
  - Support for multisample h5ad files with automatic sample detection
  - Automatic import mode inference from selected dataset
  - Performance improvements: replaced xsv.exportFrame with custom merge-csv template to improve performance
  - Removed gene filtering functionality from the workflow

## 1.2.0

### Minor Changes

- febbe4e: Added support for cell ranger output files.

## 1.1.2

### Patch Changes

- ce5778a: Add file check in prerun

## 1.1.1

### Patch Changes

- 2c32685: Technical release

## 1.1.0

### Minor Changes

- 778be21: First block version supporting xsv counts import
