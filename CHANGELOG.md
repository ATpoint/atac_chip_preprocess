# atac_chip_preprocess

## 2.0
- added error strategy `finish` to all modules
- all modules use the same publishing mode which by default is `copy`
- some eye-candy in the startup message
- only run with <= 21.04.3 as in later versions include statements must be outside the workflow definitions, which we have not yet done
- remove the conda CI test, only run with Docker
