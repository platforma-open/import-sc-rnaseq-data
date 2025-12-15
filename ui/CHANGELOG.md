# @platforma-open/milaboratories.import-sc-rnaseq-data.ui

## 1.3.3

### Patch Changes

- 8d8bf74: Add validation of input seurat files

## 1.3.2

### Patch Changes

- Updated dependencies [34f4c2e]
  - @platforma-open/milaboratories.import-sc-rnaseq-data.model@1.5.0

## 1.3.1

### Patch Changes

- Updated dependencies [81fa688]
  - @platforma-open/milaboratories.import-sc-rnaseq-data.model@1.4.0

## 1.3.0

### Minor Changes

- c7bd46d: h5ad format support

  - Added support for h5ad format input files (AnnData format)
  - Support for multisample h5ad files with automatic sample detection
  - Automatic import mode inference from selected dataset
  - Performance improvements: replaced xsv.exportFrame with custom merge-csv template to improve performance
  - Removed gene filtering functionality from the workflow

### Patch Changes

- Updated dependencies [c7bd46d]
  - @platforma-open/milaboratories.import-sc-rnaseq-data.model@1.3.0

## 1.2.1

### Patch Changes

- 4623894: Do not show h5ad option until fully supported

## 1.2.0

### Minor Changes

- febbe4e: Added support for cell ranger output files.

### Patch Changes

- Updated dependencies [febbe4e]
  - @platforma-open/milaboratories.import-sc-rnaseq-data.model@1.2.0

## 1.1.2

### Patch Changes

- ce5778a: Add file check in prerun
- Updated dependencies [ce5778a]
  - @platforma-open/milaboratories.import-sc-rnaseq-data.model@1.1.2

## 1.1.1

### Patch Changes

- 2c32685: Technical release
- Updated dependencies [2c32685]
  - @platforma-open/milaboratories.import-sc-rnaseq-data.model@1.1.1

## 1.1.0

### Minor Changes

- 778be21: First block version supporting xsv counts import

### Patch Changes

- Updated dependencies [778be21]
  - @platforma-open/milaboratories.import-sc-rnaseq-data.model@1.1.0
