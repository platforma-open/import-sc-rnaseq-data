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
