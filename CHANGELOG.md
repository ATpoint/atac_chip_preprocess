# atac_chip_preprocess

# v2.1
- input is now a samplesheet indicating sample name and the R1/R2 files plus a group information.
- separated trimming and alignment/sort into separate processes
- allow multiple fastq files per sample, will be merged prior to alignment by a dedicated process
- called peaks are now filtered against provided blacklists, currently mm10 and hg38 are supported
- now publish software versions and command lines (the latter separated by process and sample) into a dedicated pipeline directory
- major code cleanup
- new container versions
- better documentation

## 2.0
- added error strategy `finish` to all modules
- all modules use the same publishing mode which by default is `copy`
- some eye-candy in the startup message
- only run with <= 21.04.3 as in later versions include statements must be outside the workflow definitions, which we have not yet done
- remove the conda CI test, only run with Docker

## < 2.0
(...)