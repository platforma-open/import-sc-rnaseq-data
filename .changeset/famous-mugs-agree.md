---
'@platforma-open/milaboratories.import-sc-rnaseq-data.software': minor
'@platforma-open/milaboratories.import-sc-rnaseq-data.workflow': minor
'@platforma-open/milaboratories.import-sc-rnaseq-data.model': minor
'@platforma-open/milaboratories.import-sc-rnaseq-data': minor
---

Add Seurat (rds) datasets support, optimize performance

- Added support for Seurat RDS file format with automatic import mode detection
- Implemented multi-sample Seurat file processing with automatic sample detection and filtering
- Created new R-based software components (extract_counts.R, infer_species.R) for Seurat object processing
- Switched output format from CSV to Parquet with zstd compression for improved I/O performance and reduced storage footprint
- Optimized sparse matrix handling and data.table operations for improved memory efficiency
- Enhanced memory calculation utilities to dynamically allocate resources based on data size
