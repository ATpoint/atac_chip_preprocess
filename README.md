# atac_chip_preprocess

This Nextflow pipeline is intended for preprocessing of DNA-seq
experiments such as ATAC-seq, ChIP-seq and CUT&RUN. It currently
performs trimming and alignment of fastq files, filtering of the
resulting BAM files, basic peak calling and a QC assessment by
calculating Fractions Of Reads Per Peaks (FRiPs). 

<br>

![CI](https://github.com/ATpoint/atac_chip_preprocess/actions/workflows/CI.yml/badge.svg)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A4%2021.04.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

<br>

## Usage

Creating an index for the mapping on our HPC, with 16 cores using SLURM and the `normal` queue.
Output will be in `atac_chip_preprocess_results/bowtie2_idx/` created in the location from which the run is launched.

```bash

NXF_VER=21.04.3 \
    nextflow run atpoint/atac_chip_preprocess -r 2.0 -profile singularity,slurm -with-trace -with-report report.html \
        --ref_genome /path/to/genome.fa.gz --idx_threads 16 --idx_mem '16.GB' --only_idx --queue 'normal' \
        -bg > report.log

```

Mapping a dataset. Note, for single-end data use `--mode single`, for paired-end data use `--mode paired` and for ATAC-seq
add additionally the `--atacseq` flag.

```bash

NXF_VER=21.04.3 \
    nextflow run atpoint/atac_chip_preprocess -r 2.0 -profile singularity,slurm -with-trace -with-report report.html \
        --idx path/to/idx/idxbasename \
        --fastq path/to/fastqfolder/*_{1,2}.fastq.gz' \
        --align_threads 12 --sort_threads 2 --sort_mem '4G' --queue 'normal' \
        -bg > report.log

```

## Pipeline

- trim reads with `cutadapt`, either the TruSeq adapter (default) or 
the Nextera adapter (default if using `--atacseq`)
- map trimmed reads with `bowtie2` and mark duplicates with `samblaster`
- remove unmapped reads, alignments with MAPQ < 20, mitochondrial alignments (for ATAC-seq), alignments to non-primary chromosomes,
non-primary and supplementary alignments, PCR duplicates
- call peaks per sample with `macs2`
- calculate FRiPs per sample based on the per-sample peaks
- check insert sizes (for paired-end data)

## Containerization

Use either `-profile docker,singularity` to run via the hardcoded container. 
Using `-profile conda` is possible but untested, it will build the environment based on the `environment.yml` file which the container is based on.

