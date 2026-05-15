---
"@platforma-open/milaboratories.import-sc-rnaseq-data": patch
"@platforma-open/milaboratories.import-sc-rnaseq-data.software": patch
---

Make `seurat_operations.R` executable in git (mode 0755) so the docker `cmd` `/app/seurat_operations.R` actually runs. The script has a `#!/usr/bin/env Rscript` shebang, but `COPY *.R` preserves the source mode — and the file was committed as 0644, so direct execution failed with `Permission denied`. Matches the existing pattern in other R-using blocks (`differential-expression`, `functional-analysis`, `star-read-mapping` all keep their entry scripts at 0755).

The block patch bump is needed so the new software image actually ships — the block re-packs and the published manifest references the rebuilt image.
