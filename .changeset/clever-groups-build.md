---
'@platforma-open/milaboratories.import-sc-rnaseq-data.software': minor
'@platforma-open/milaboratories.import-sc-rnaseq-data.workflow': minor
'@platforma-open/milaboratories.import-sc-rnaseq-data': minor
'@platforma-open/milaboratories.import-sc-rnaseq-data.model': minor
'@platforma-open/milaboratories.import-sc-rnaseq-data.ui': minor
---

h5ad format support

- Added support for h5ad format input files (AnnData format)
- Support for multisample h5ad files with automatic sample detection
- Automatic import mode inference from selected dataset
- Performance improvements: replaced xsv.exportFrame with custom merge-csv template for better memory efficiency
- Removed gene filtering functionality from the workflow
- Updated Python dependencies to support h5ad file processing
- Enhanced species inference logic for h5ad files
- Added extract-single-file template for handling single file extraction
- Updated cell metrics calculation to work with h5ad format
- UI improvements for dataset selection and import mode handling
