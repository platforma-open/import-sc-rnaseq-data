---
"@platforma-open/milaboratories.import-sc-rnaseq-data": patch
---

Re-pack the block so it picks up the latest `…import-sc-rnaseq-data.software` image (MILAB-6263 follow-up). The previous MILAB-6263 PR bumped only the `.software` package — because the workflow declares `.software` as a `devDependency`, changesets' `updateInternalDependencies` didn't propagate the bump up to the block, the block stayed at 1.4.3, and the published block manifest kept pointing at the pre-fix seurat-utils image (which still ran `/app/run.sh` and hit `Permission denied`).
