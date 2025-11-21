# @platforma-open/milaboratories.import-sc-rnaseq-data

## 1.3.0

### Minor Changes

- 81fa688: Add Seurat (rds) datasets support, optimize performance

  - Added support for Seurat RDS file format with automatic import mode detection
  - Implemented multi-sample Seurat file processing with automatic sample detection and filtering
  - Created new R-based software components (extract_counts.R, infer_species.R) for Seurat object processing
  - Switched output format from CSV to Parquet with zstd compression for improved I/O performance and reduced storage footprint
  - Optimized sparse matrix handling and data.table operations for improved memory efficiency
  - Enhanced memory calculation utilities to dynamically allocate resources based on data size

### Patch Changes

- Updated dependencies [81fa688]
  - @platforma-open/milaboratories.import-sc-rnaseq-data.workflow@1.4.0
  - @platforma-open/milaboratories.import-sc-rnaseq-data.model@1.4.0
  - @platforma-open/milaboratories.import-sc-rnaseq-data.ui@1.3.1

## 1.2.0

### Minor Changes

- c7bd46d: h5ad format support

  - Added support for h5ad format input files (AnnData format)
  - Support for multisample h5ad files with automatic sample detection
  - Automatic import mode inference from selected dataset
  - Performance improvements: replaced xsv.exportFrame with custom merge-csv template to improve performance
  - Removed gene filtering functionality from the workflow

### Patch Changes

- Updated dependencies [c7bd46d]
  - @platforma-open/milaboratories.import-sc-rnaseq-data.workflow@1.3.0
  - @platforma-open/milaboratories.import-sc-rnaseq-data.model@1.3.0
  - @platforma-open/milaboratories.import-sc-rnaseq-data.ui@1.3.0

## 1.1.6

### Patch Changes

- Updated dependencies [4623894]
  - @platforma-open/milaboratories.import-sc-rnaseq-data.ui@1.2.1

## 1.1.5

### Patch Changes

- Updated dependencies [fbcd123]
  - @platforma-open/milaboratories.import-sc-rnaseq-data.workflow@1.2.1

## 1.1.4

### Patch Changes

- Updated dependencies [febbe4e]
  - @platforma-open/milaboratories.import-sc-rnaseq-data.workflow@1.2.0
  - @platforma-open/milaboratories.import-sc-rnaseq-data.model@1.2.0
  - @platforma-open/milaboratories.import-sc-rnaseq-data.ui@1.2.0

## 1.1.3

### Patch Changes

- Updated dependencies [2bcce25]
  - @platforma-open/milaboratories.import-sc-rnaseq-data.workflow@1.1.3

## 1.1.2

### Patch Changes

- Updated dependencies [ce5778a]
  - @platforma-open/milaboratories.import-sc-rnaseq-data.workflow@1.1.2
  - @platforma-open/milaboratories.import-sc-rnaseq-data.model@1.1.2
  - @platforma-open/milaboratories.import-sc-rnaseq-data.ui@1.1.2

## 1.1.1

### Patch Changes

- 2c32685: Technical release
- Updated dependencies [2c32685]
  - @platforma-open/milaboratories.import-sc-rnaseq-data.model@1.1.1
  - @platforma-open/milaboratories.import-sc-rnaseq-data.ui@1.1.1
  - @platforma-open/milaboratories.import-sc-rnaseq-data.workflow@1.1.1

## 1.1.0

### Minor Changes

- 778be21: First block version supporting xsv counts import

### Patch Changes

- Updated dependencies [778be21]
  - @platforma-open/milaboratories.import-sc-rnaseq-data.workflow@1.1.0
  - @platforma-open/milaboratories.import-sc-rnaseq-data.model@1.1.0
  - @platforma-open/milaboratories.import-sc-rnaseq-data.ui@1.1.0
