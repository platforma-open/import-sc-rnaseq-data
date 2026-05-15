# @platforma-open/milaboratories.import-sc-rnaseq-data.software

## 1.6.3

### Patch Changes

- 6800989: Make `seurat_operations.R` executable in git (mode 0755) so the docker `cmd` `/app/seurat_operations.R` actually runs. The script has a `#!/usr/bin/env Rscript` shebang, but `COPY *.R` preserves the source mode — and the file was committed as 0644, so direct execution failed with `Permission denied`. Matches the existing pattern in other R-using blocks (`differential-expression`, `functional-analysis`, `star-read-mapping` all keep their entry scripts at 0755).

  The block patch bump is needed so the new software image actually ships — the block re-packs and the published manifest references the rebuilt image.

## 1.6.2

### Patch Changes

- 98057d3: Fix runtime `Permission denied` when the container runs as a non-root UID (MILAB-6263). The entrypoint re-invoked `renv::restore()` on every start, which tries to reconcile the system R library at `/usr/local/lib/R/site-library/` with the project lockfile. When the `r-base:4.4.2` base image preinstalled a version of a locked package (e.g. `rlang`) that differs from `renv.lock`, renv attempted to back up the system-library copy before replacing it — failing on hosts that run the container unprivileged. renv now installs into a project-local library at `/app/renv/library` and `R_LIBS_USER` points R at the same path, so the obsolete `/app/run.sh` wrapper and runtime restore are removed.

## 1.6.1

### Patch Changes

- 8d8bf74: Add validation of input seurat files

## 1.6.0

### Minor Changes

- 34f4c2e: H5 file format support

## 1.5.0

### Minor Changes

- 81fa688: Add Seurat (rds) datasets support, optimize performance

  - Added support for Seurat RDS file format with automatic import mode detection
  - Implemented multi-sample Seurat file processing with automatic sample detection and filtering
  - Created new R-based software components (extract_counts.R, infer_species.R) for Seurat object processing
  - Switched output format from CSV to Parquet with zstd compression for improved I/O performance and reduced storage footprint
  - Optimized sparse matrix handling and data.table operations for improved memory efficiency
  - Enhanced memory calculation utilities to dynamically allocate resources based on data size

## 1.4.0

### Minor Changes

- c7bd46d: h5ad format support

  - Added support for h5ad format input files (AnnData format)
  - Support for multisample h5ad files with automatic sample detection
  - Automatic import mode inference from selected dataset
  - Performance improvements: replaced xsv.exportFrame with custom merge-csv template to improve performance
  - Removed gene filtering functionality from the workflow

## 1.3.0

### Minor Changes

- c281e5a: Fixed script compatibility with other OS

## 1.2.0

### Minor Changes

- febbe4e: Added support for cell ranger output files.

## 1.1.3

### Patch Changes

- 2bcce25: Allow tSV inputs

## 1.1.2

### Patch Changes

- ce5778a: Add file check in prerun

## 1.1.1

### Patch Changes

- 2c32685: Technical release

## 1.1.0

### Minor Changes

- 778be21: First block version supporting xsv counts import
