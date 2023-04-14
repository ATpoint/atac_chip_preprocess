# atac_chip_preprocess

  
![CI](https://github.com/ATpoint/atac_chip_preprocess/actions/workflows/CI.yml/badge.svg)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with apptainer/singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000&logo=data%3Aimage%2Fjpeg%3Bbase64%2C%2F9j%2F4AAQSkZJRgABAQAAkACQAAD%2F2wBDABwcHBwcHDAcHDBEMDAwRFxEREREXHRcXFxcXHSMdHR0dHR0jIyMjIyMjIyoqKioqKjExMTExNzc3Nzc3Nzc3Nz%2F2wBDASIkJDg0OGA0NGDmnICc5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ub%2FwAARCACAAHgDASIAAhEBAxEB%2F8QAGgAAAgMBAQAAAAAAAAAAAAAABAUAAgMBBv%2FEADMQAAIBAwIDBgYBAwUAAAAAAAECAAMEESExEpHRBRVBUVJxEyIyYYGxQhQzwUNTgpLh%2F8QAGAEAAwEBAAAAAAAAAAAAAAAAAAIDAQT%2FxAAfEQEBAQACAgMBAQAAAAAAAAAAAQIDESExEkFRIjL%2F2gAMAwEAAhEDEQA%2FADe5rX1PzHSTua19T8x0jYbTsAUdzWvqfmOknc1r6n5jpG8yq1kpDLn8eMAW9zWvqfmOko3ZNmmrOw9yOk0qXdR9F%2BUfbeCHU5OspOP9SvJ%2BIbHs4f6jn2x0lk7OsnVmDVMLvt0lYZb%2FANqr7CG8yS0Z3bQw7OsD%2FNx746TVeybNvpdz%2BR0kk%2B4nLOSn7W7mtfU%2FMdJO5rX1PzHSbJcOujfMIalRagyplM6lb2WdzWvqfmOknc1r6n5jpG8kZpR3Na%2Bp%2BY6SRvJAODadnBtMq9UUU4vHwEILelLi4FEYGrGKWZnPExyTOMxYlmOSZtRoNWbTQDcy8kzHPbdVkqM54UGTDqdid6h%2FA6w6nTSmvCgxLxLu%2FSk459h1taA%2Fjn31mop01BCqADvLznEvmIndp%2BoyNCkf449pg9r4ofwYZkHYzsS4lHRQyspwwwZFYqeJTgxo6K4wwzF1Wk1I66jwMjrHXllg2jWFQYOjTeJwSpBG4jOlUFRc%2BPjKY334rZWskkko1wbRPc1fi1DjYaCMq7%2FDolhvsPzEsrxz7S5L9NKdM1XCL4x0iLTUIuwgtlT4U%2BId2%2FUNi7vd6bjPU7Ud0poXc4UakmILjtaq5K244F8zqeXhO9r1yai242UcR9ztEsMw9rZ69aocu7H8n%2FGJjgeIEk0p06lVxTpqWY%2BAjMZj5dV09tP1N0ubin9FRh%2Bc%2FvM2fs%2B8QZamSPsQYHAHFn2jd1KyUWw%2FEcZIwceJ0noGUOpVtjEHY9HiqPXP8Rwj3OpnoZPTYUuhRipl6L8DjyOhhNymU4xuP1AJzWfGsOZJlRbjpgyS8vZgt63yIvmcxbvp5w69yXT2MFpqTUUHzE6M%2BnPvzo7VQqhR4CWkkkHQ8l2ln%2BtqZ%2B2PbEBjzti3IIul20Vv8RJKy%2BC1yF2d0bSr8Th4gRgjx%2FEEkgHsaF7bXGlNhnyOh5TC9sEuV40wKg2Pn9jPKw%2B37QuaGnFxr5N13i%2FH8b29BYUDb2qIww27e5hkGtbqndoXp5GDgg%2BBhMStcYcSlT4xRtpHEU1Bh2H3MlyRlF2p0ZfIySlr9Te0kbHoRtVdkxjxEzFZiQDiXrj5QfKCx2mckqp4lB85aAUcIylKgBDaYPjEdz2QQS1qdPS3%2BDL9sViPh0VODniOPttNLPtRKgFO5PC%2B3F4HoY079sIKlKpRbhqqUP36zOe3dadVOFwGU%2BeonjrhaaV6iUjlAcCNL2yxjJJJGAuzuGt7hXB0JCsPsek9hPCHY4nuV%2Bke0TTYtF7XNQMQMYz5Q52CqWPgIonNy6666V4537H29V6hIbGkk5aD5S3mZI%2BP8%2BS79%2BBDrxIRF8ZDaB1k4WyNjHK1oNkcB8IRFqsVIYeEPVhUXIgHkr6r8a6dxsDwj2H%2FALmCR1V7GqDWjUDfZtDzHSBP2fepvTJ9iDKSwoMEgcIJA8gTicmzW9wv1UnH%2FEyvwq3%2B2%2F8A1MYM5IStpdP9NJvyMfuH0eyKznNdgg8hqekzsA7G3NzcKP4qQzH22H5nr5jRoUrdPh0hgfv3l6jimvEZPWvs0ga6fAFMeOpgMszF2LNuZtb0%2BN8nZZyW%2FLTon8wdSTgphZJpJOmTrw564NpV1Drgyw2nZoLmUqcGdRyhyIZUphxrv5wJlKHDQA5KiuNOUvFoJGomy12H1awAySYCuh3yJf4qFS2dBvANJIMbqmNsmYPdO2i%2FL%2B4l5Mw8xaMqVVpjXfyi2pUao2W5ShJJydZZEZzwqJDW7rwrnMy4ql2CruY1poKahRK0qS0h5k7may3HjrzUt67SSSSUI4Np2J%2B%2Bbb0vyHWd75tvS%2FIdYA3lWUMMMMxV3zbel%2BQ6yd823pfkOsAMegRqmswKlfqGJl3zbel%2BQ6znfNqf4PyHWAazVf7VT2gR7VsjvTbkOsnetngrwPg77dZlnhsWnVVm%2BkZmY7TsRtTbkOs075tRsj8h1kZw%2FtVvJ%2BCUtWOtQ4%2BwhqoqDCjAirvm29L8h1k75tvS%2FIdZXOJPSd1abyRR3zbel%2BQ6yd823pfkOsYpvJE%2FfNt6X5DrJAP%2F2Q%3D%3D)](https://sylabs.io/docs/)  

## Introduction

**atac_chip_preprocess** is a containerized Nextflow pipeline for preprocessing of ATAC-seq and ChIP-seq data. 

The workflow consists of these steps:  

- validation of the provided samplesheet
- initial QC with `fastqc`
- merging of lane replicates per sample into one fastq file per R1/R2
- adapter and quality trimming with `fastp`
- mapping with `bowtie2` using the `--very-sensitive` (and `-X 2000` for paired-end data) flags  
- duplicate marking with `samblaster`
- removal of MAPQ < 20, non-primary or supplementary, reads mapped to non-primary (random/unplaced) chromosomes (regex `chr[1-9,X,Y]`), mitochondrial alignments and duplicate reads with `samtools`
- for paired-end data fetching of insert size metrics with `picard`
- for ATAC-seq data extraction of transposome insertion events (cutsites) using custom GNU tool combinations
- peak calling with `macs2` and filtering of peaks against NGS blacklists (ENCODE+mitochondrial homologs in the nuclear genome, the latter for ATAC-seq only) using `bedtools`
- calculation of Fractions Of Reads in Peaks (FRiPs) as a QC metric with `featureCounts`
- creation of raw bigwig tracks for visual inspection of data quality with `bedtools`
- summary report with `MultiQC`
- output of all used software versions and the exact command lines per process step and sample using custom scripts

Run the following test profile to see all possible outputs that the pipeline produces. Default output directory is `./atac_chip_preprocess_results/`). 
Downloading the [Docker image](https://hub.docker.com/r/atpoint/atac_chip_preprocess) may take a minute or two (automated).

```bash
NXF_VER=21.10.6 nextflow run atpoint/atac_chip_preprocess -r main -profile docker,test --keep_merge --keep_trim
```

An overview of current software versions and exact command lines when using default settings of the pipeline can be found in the [misc directory](misc/).  

## Usage

The minimal parameters the user has to provide are the following ones:

- `--samplesheet`: path to a [samplesheet csv file](test/samplesheet.csv) with three columns, being `sample` (the sample name), `r1` (path to R1) and `r2` (path to R2), where r2 can be empty. If empty, then the sample is considered single-end.  
- `--index`: path to a folder containing a `bowtie2` index with the typical `*.bt2` files. Note, it is the path to the folder, not the path to the index basename, as the pipeline will find the bt2 files automatically.  
- `--species`: either of `mm` or `hs` to let the peak caller know whether mouse or human data are provided, so it gets the effective genome length right.  

Note that the `bowtie2` index must be produced beforehand, we did not include that into the pipeline as it is trivially just `bowtie2-build genome.fa idx`.

On our HPC we typically use:  

```bash
# Example for mouse ATAC-seq data
NXF_VER=21.10.6 nextflow run atpoint/atac_chip_preprocess -r main -profile apptainer,test --samplesheet path/to/samplesheet.csv --index path/to/index_folder --species mm

# Example for mouse ChIP-seq data
NXF_VER=21.10.6 nextflow run atpoint/atac_chip_preprocess -r main -profile apptainer,test --samplesheet path/to/samplesheet.csv --index path/to/index_folder --species mm --atacseq false
```

Use either of `-profile docker/singularity/apptainer` to use any of these container engines.

## Options

We used reasonable defaults for all processing steps that should be used without modifications. Still, the following options exist for customization:  

**General options**  

- `--outdir`, path to desired output folder collecting all results  
- `--atacseq`, a logical, set to `false` if processing something like ChIP-seq data, by default `true` for ATAC-seq data  

**Filtering options**  

- `--blacklist`: path to a BED file to filter peaks against. By default when `--species` is `mm` then the provided mm10 blacklist is used, for `hs` the hg38 one is used.  
- `--filter_blacklist`: logical, set to `false` to turn off any blacklist filtering, default `true`.  
- `--flag_remove`: a numeric flag to be used with `samtools view -F`, so indicating which alignments to remove. Default is 3332, so discard unmapped, not primary, supplementary and duplicates . See [here](https://broadinstitute.github.io/picard/explain-flags.html) for details.  
- `--chr_regex`: a groovy-compatible regex to indicate which chromosomes to keep in the BAM alignments. Default is `chr[1-9,X,Y]` which means keep everything starting with `chr` and then a number of X/Y. That in turn removes decoys (`chrEBV`) and unplaced/random contigs such as `chrU...`, therefore keeping only the primary autosomes and sex chromosomes.  
- `--min_mapq`: an integer, keep only alignments with MAPQ greater than that, default is 20.  
- `--fragment_length`: for single-end data an average expected fragment length to extend reads to fragments for bigwig creation and FRiP calculation, default is 250. That is only used if `--atacseq false` as for ATAC-seq data everything is based on the transposome cutsites (that is the 5' ends of the alignments).  
- `keep_merge`: logical, whether to keep the merged fastq files, else they're not published to the output directory  
- `keep_trim`: logical, whether to keep the trimmed fastq files, else they're not published to the output directory  

**Process options**  

- `--do_not_trim`: logical, whether to skip adapter and quality trimming  
- `--trim_additional`: additional arguments for the `fastp` trimming process beyond what is coded in the module definition, default `--dont_eval_duplication -z 6` to skip duplicate level assessment and to compress outputs  
- `--align_additional`: additional arguments for the `bowtie2` alignment process beyond what is coded in the module definition, default is `-X2000 --very-sensitive`, see `bowtie2` [manual](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
- `--sort_additional`: additional arguments for the `samtools` sorting process beyond what is coded in the module definition, default is `-l 6` to compress the resulting BAM file to that level  
- `--filter_additional`: additional arguments for the `samtools view` filtering process beyond what is described above and given with the `-q` and `-F` flags
- `--macs_additional`: additional arguments for the `macs2 callpeak`, default is in any case `--keep-dup=all` since we provide already deduplicated data to that process and if ATAC-seq data are processed (default) then `--nomodel --extsize 100 --shift -50 --min-length 250` to provide some smoothing when using the cutsites for peak calling.

## Resources

The [nextflow.config](nextflow.config) files contains hardcoded defaults towards resources for the individual processes, suitable for use on HPC environments. The most demanding process
is the alignment steps, requiring 16 threads and 16GB of RAM per sample.

## Schedulers

The [schedulers.config](configs/schedulers.config) file currently contains a single scheduler profile for SLURM as used on or HPC,
submitting jobs (if using `-profile slurm`) to a quere called `normal` with a maximum 8h of walltime. Custom profiles should be added to this config.

