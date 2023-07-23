# atac_chip_preprocess

  
![CI](https://github.com/ATpoint/atac_chip_preprocess/actions/workflows/CI.yml/badge.svg)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with apptainer/singularity](https://img.shields.io/badge/run%20with-apptainer/singularity-1d355c.svg?labelColor=000000&logo=data%3Aimage%2Fjpeg%3Bbase64%2C%2F9j%2F4AAQSkZJRgABAQAAkACQAAD%2F2wBDABwcHBwcHDAcHDBEMDAwRFxEREREXHRcXFxcXHSMdHR0dHR0jIyMjIyMjIyoqKioqKjExMTExNzc3Nzc3Nzc3Nz%2F2wBDASIkJDg0OGA0NGDmnICc5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ub%2FwAARCACAAHgDASIAAhEBAxEB%2F8QAGgAAAgMBAQAAAAAAAAAAAAAABAUAAgMBBv%2FEADMQAAIBAwIDBgYBAwUAAAAAAAECAAMEESExEpHRBRVBUVJxEyIyYYGxQhQzwUNTgpLh%2F8QAGAEAAwEBAAAAAAAAAAAAAAAAAAIDAQT%2FxAAfEQEBAQACAgMBAQAAAAAAAAAAAQIDESExEkFRIjL%2F2gAMAwEAAhEDEQA%2FADe5rX1PzHSTua19T8x0jYbTsAUdzWvqfmOknc1r6n5jpG8yq1kpDLn8eMAW9zWvqfmOko3ZNmmrOw9yOk0qXdR9F%2BUfbeCHU5OspOP9SvJ%2BIbHs4f6jn2x0lk7OsnVmDVMLvt0lYZb%2FANqr7CG8yS0Z3bQw7OsD%2FNx746TVeybNvpdz%2BR0kk%2B4nLOSn7W7mtfU%2FMdJO5rX1PzHSbJcOujfMIalRagyplM6lb2WdzWvqfmOknc1r6n5jpG8kZpR3Na%2Bp%2BY6SRvJAODadnBtMq9UUU4vHwEILelLi4FEYGrGKWZnPExyTOMxYlmOSZtRoNWbTQDcy8kzHPbdVkqM54UGTDqdid6h%2FA6w6nTSmvCgxLxLu%2FSk459h1taA%2Fjn31mop01BCqADvLznEvmIndp%2BoyNCkf449pg9r4ofwYZkHYzsS4lHRQyspwwwZFYqeJTgxo6K4wwzF1Wk1I66jwMjrHXllg2jWFQYOjTeJwSpBG4jOlUFRc%2BPjKY334rZWskkko1wbRPc1fi1DjYaCMq7%2FDolhvsPzEsrxz7S5L9NKdM1XCL4x0iLTUIuwgtlT4U%2BId2%2FUNi7vd6bjPU7Ud0poXc4UakmILjtaq5K244F8zqeXhO9r1yai242UcR9ztEsMw9rZ69aocu7H8n%2FGJjgeIEk0p06lVxTpqWY%2BAjMZj5dV09tP1N0ubin9FRh%2Bc%2FvM2fs%2B8QZamSPsQYHAHFn2jd1KyUWw%2FEcZIwceJ0noGUOpVtjEHY9HiqPXP8Rwj3OpnoZPTYUuhRipl6L8DjyOhhNymU4xuP1AJzWfGsOZJlRbjpgyS8vZgt63yIvmcxbvp5w69yXT2MFpqTUUHzE6M%2BnPvzo7VQqhR4CWkkkHQ8l2ln%2BtqZ%2B2PbEBjzti3IIul20Vv8RJKy%2BC1yF2d0bSr8Th4gRgjx%2FEEkgHsaF7bXGlNhnyOh5TC9sEuV40wKg2Pn9jPKw%2B37QuaGnFxr5N13i%2FH8b29BYUDb2qIww27e5hkGtbqndoXp5GDgg%2BBhMStcYcSlT4xRtpHEU1Bh2H3MlyRlF2p0ZfIySlr9Te0kbHoRtVdkxjxEzFZiQDiXrj5QfKCx2mckqp4lB85aAUcIylKgBDaYPjEdz2QQS1qdPS3%2BDL9sViPh0VODniOPttNLPtRKgFO5PC%2B3F4HoY079sIKlKpRbhqqUP36zOe3dadVOFwGU%2BeonjrhaaV6iUjlAcCNL2yxjJJJGAuzuGt7hXB0JCsPsek9hPCHY4nuV%2Bke0TTYtF7XNQMQMYz5Q52CqWPgIonNy6666V4537H29V6hIbGkk5aD5S3mZI%2BP8%2BS79%2BBDrxIRF8ZDaB1k4WyNjHK1oNkcB8IRFqsVIYeEPVhUXIgHkr6r8a6dxsDwj2H%2FALmCR1V7GqDWjUDfZtDzHSBP2fepvTJ9iDKSwoMEgcIJA8gTicmzW9wv1UnH%2FEyvwq3%2B2%2F8A1MYM5IStpdP9NJvyMfuH0eyKznNdgg8hqekzsA7G3NzcKP4qQzH22H5nr5jRoUrdPh0hgfv3l6jimvEZPWvs0ga6fAFMeOpgMszF2LNuZtb0%2BN8nZZyW%2FLTon8wdSTgphZJpJOmTrw564NpV1Drgyw2nZoLmUqcGdRyhyIZUphxrv5wJlKHDQA5KiuNOUvFoJGomy12H1awAySYCuh3yJf4qFS2dBvANJIMbqmNsmYPdO2i%2FL%2B4l5Mw8xaMqVVpjXfyi2pUao2W5ShJJydZZEZzwqJDW7rwrnMy4ql2CruY1poKahRK0qS0h5k7may3HjrzUt67SSSSUI4Np2J%2B%2Bbb0vyHWd75tvS%2FIdYA3lWUMMMMxV3zbel%2BQ6yd823pfkOsAMegRqmswKlfqGJl3zbel%2BQ6znfNqf4PyHWAazVf7VT2gR7VsjvTbkOsnetngrwPg77dZlnhsWnVVm%2BkZmY7TsRtTbkOs075tRsj8h1kZw%2FtVvJ%2BCUtWOtQ4%2BwhqoqDCjAirvm29L8h1k75tvS%2FIdZXOJPSd1abyRR3zbel%2BQ6yd823pfkOsYpvJE%2FfNt6X5DrJAP%2F2Q%3D%3D)](https://sylabs.io/docs/)  

--| [Overview](overview)  
--| [Usage](usage)  
--| [Options](options)  
--| [Output](output)  
--| [Resources](resources)  
--| [Schedulers](schedulers)  
--| [Software](software)  
## Overview

**atac_chip_preprocess** is a containerized Nextflow pipeline for preprocessing of ATAC-seq and ChIP-seq data.  

The pipeline consists of these steps, starting from a [samplesheet](https://github.com/ATpoint/atac_chip_preprocess/blob/dev/test/samplesheet.csv) to read the fastq files:  

- Validation of the samplesheet to ensure that fastq files are not duplicated and all file paths exist. If any validation fails then the process will return an informative error for debuggung.
- Initial fastq QC with `fastqc`. 
- Merging of lane/technical replicates per sample, see the samplesheet section below on how to indicate technical replicates (=multiple fastq files) per sample. 
- Adapter and quality trimming with `fastp`.
- Mapping of reads with `bowtie2`, by default with the `--very-sensitive` (and `-X 2000` for paired-end data) flags.
The pipeline does not include a dedicated indexing step for the reference genome since this is a one-liner via `bowtie2-build`. The user needs to build the index upfront and then provide the path to the folder with the `*.bt2` index files via `--index`. The files will then automatically be found in that folder. It is expected that only a single set of index file is in that folder.
- Duplicate marking with `samblaster`.
- Removal of alignments with MAPQ below 20, removal of non-primary or supplementary alignments, removal of reads not mapped to primary chromosomes (the regex to keep alignments is `chr[1-9,X,Y]`), removal of mitochondrial alignments and duplicate reads with `samtools` and combinations of GNU tools.
- For paired-end data fetching of insert size metrics with `picard CollectInsertSizeMetrics`.
- For ATAC-seq data extraction of transposome insertion events (cutsites, that is the 5' end of alignments) using GNU tool combinations,
output in gzipped BED format.
- Peak calling with `macs2`. For ATAC-seq the options `--keep-dup=all --nomodel --extsize 100 --shift -50 --min-length 250` are used by default, else it is only `--keep-dup=all` since the deduplicated BAM files are used as input. Then, peaks are filtering against NGS blacklists (ENCODE + mitochondrial homologs in the nuclear genome, the latter for ATAC-seq only) using `bedtools`. This is species-specific.
Currently, human and mouse is supported via the flag `--species` with either `hs` or `mm`.
- Calculation of Fractions Of Reads in Peaks (FRiPs) as a signal-to-noise QC metric with `featureCounts`.
- Creation of bigwig tracks (raw counts) for quick visual inspection of data quality with `bedtools genomecov`.
- Summary report collecting fastqc, trimming, alignment metrics and FRiP/featureCounts metrics with `MultiQC`.
- Output of all used software versions and the exact command lines per process and sample using custom scripts.

Execute this command to run the pipeline on a tiny test dataset with minimal resources to explore outputs:

```bash
NXF_VER=23.04.0 nextflow run atpoint/atac_chip_preprocess -r main -profile docker,test --keep_merge --keep_trim
```

An overview of current software versions and exact command lines when using default settings of the pipeline can be found [here](misc/).

## Usage

The minimal parameters the user has to provide are the following ones:

- `--samplesheet`: path to a [samplesheet csv file](test/samplesheet.csv) with three columns, being `sample` (the sample name), `r1` (path to R1) and `r2` (path to R2), where r2 can be empty. If empty, then the sample is considered single-end.  
- `--index`: path to a folder containing a `bowtie2` index with the typical `*.bt2` files. Note, it is the path to the folder, not the path to the index basename, as the pipeline will find the bt2 files automatically.  
- `--species`: either of `mm` or `hs` to let the peak caller know whether mouse or human data are provided, so it gets the effective genome length right.  

Note that the `bowtie2` index must be produced upfront, we did not include that into the pipeline as it is trivially just `bowtie2-build genome.fa idx`.

On our HPC we typically use this command below, with Apptainer (currently the profiles `docker`, `singularity` and `apptainer` are supported) as container engine and SLURM as executor. If any [other executor](https://www.nextflow.io/docs/latest/executor.html) shall be used then the user needs to add it to the [scheduler config file](./configs/schedulers.config) file.

```bash
# Example for mouse ATAC-seq data
NXF_VER=23.04.0 nextflow run atpoint/atac_chip_preprocess -r main -profile apptainer,slurm --samplesheet path/to/samplesheet.csv --index path/to/index_folder --species mm

# Example for mouse ChIP-seq data
NXF_VER=23.04.0 nextflow run atpoint/atac_chip_preprocess -r main -profile apptainer,slurm --samplesheet path/to/samplesheet.csv --index path/to/index_folder --species mm --atacseq false
```

Use either of `-profile docker/singularity/apptainer` to use any of these container engines.
## Options

The pipeline uses (in our opinion) reasonable defaults for all processing steps. Still, the following options exist for customization:  

**General options**  

- `--outdir`, path to desired output folder collecting all results. By default, it is `./atac_chip_preprocess_results/` in the directory from which the pipeline is launched. 
- `--atacseq`, a logical, set to `false` if processing something like ChIP-seq data, by default `true` for ATAC-seq data.

**Filtering options**  

- `--blacklist`: path to a BED file to filter peaks against. By default when `--species` is `mm` then the provided mm10 blacklist is used, for `hs` the hg38 one is used.  
- `--filter_blacklist`: logical, set to `false` to turn off any blacklist filtering, default `true`. If any species other than human and mouse is used this can be set to false so no filtering takes place and both `--species` and `--blacklist` have no effect, hence can be left at defaults.  
- `--flag_remove`: a numeric flag to be used with `samtools view -F`, so indicating which alignments to remove. Default is 3332, so discard unmapped, not primary, supplementary alignments and duplicates . See [here](https://broadinstitute.github.io/picard/explain-flags.html) for details.  
- `--chr_regex`: a groovy-compatible regex to indicate which chromosomes to keep in the BAM alignments. Default is `chr[1-9,X,Y]` which means keep everything starting with `chr` and then a number or X/Y. That in turn removes typical decoys (`chrEBV`) and unplaced/random contigs such as `chrU...`. As a result it keeps only the primary autosomes and sex chromosomes.  
- `--min_mapq`: an integer, keep only alignments with MAPQ greater than that, default is 20.  
- `--fragment_length`: for single-end data an average expected fragment length to extend reads to fragments for bigwig creation and FRiP calculation, default is 250. That is only used if `--atacseq false` as for ATAC-seq data everything is based on the transposome cutsites (that is the 5' ends of the alignments).  
- `keep_merge`: logical, whether to keep the merged fastq files, else they're not published to the output directory.  
- `keep_trim`: logical, whether to keep the trimmed fastq files, else they're not published to the output directory.  

**Process options**  

- `--do_not_trim`: logical, whether to skip adapter and quality trimming.  
- `--trim_additional`: additional arguments for the `fastp` trimming process beyond what is coded in the module definition, default `--dont_eval_duplication -z 6` to skip duplicate level assessment and to compress outputs  
- `--align_additional`: additional arguments for the `bowtie2` alignment process beyond what is coded in the module definition, default is `-X 2000 --very-sensitive`, see `bowtie2` [manual](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml).  
- `--sort_additional`: additional arguments for the `samtools sort` process beyond what is coded in the module definition, default is `-l 6` to compress the resulting BAM file to that level. Do not add `-m` or `-@` here, as resources are hardcoded in the `nextflow.config` file.   
- `--filter_additional`: additional arguments for the `samtools view` filtering process beyond what is described above and given with the `-q` and `-F` flags
- `--macs_additional`: additional arguments for the `macs2 callpeak`, default is in any case `--keep-dup=all` since we provide already deduplicated data to that process and if ATAC-seq data are processed (default) then `--nomodel --extsize 100 --shift -50 --min-length 250` to provide some smoothing to the pileup when using the cutsites for peak calling.  

## Output

By default, all outputs will be collected in `./atac_chip_preprocess_results/` relative to the directory from which the pipeline is launched. Use `--outdir` to change this. Outputs are:

- **alignments_filtered:** Sorted bam, bam index and flagstats for the filtered alignments.
- **alignments_unfiltered:** Sorted bam, bam index and flagstats for the unfiltered alignments.
- **bed_files:** In ATAC-seq mode gzipped BED files with the cutsites (5' end of filtered alignments) per sample. These are used for peak calling in the pipeline.
- **bigwig:** Bigwig files of filtered alignments, without any normalization, for a quick visual inspection of data quality.
- **fastq_merged:** The merged fastq files in case there were technical replicates per sample and `--keep_merge` was used.
- **fastq_trimmed:** The trimmed fastq files and trimming stats from `fastp` in case trimming was activated (default yes) and `--keep-trim` was used.
- **fastqc:** The fastqc per-sample outputs.
- **frips:** A file `frips_all.txt` with the FRiP score per sample. If a sample had a FRiP of zero then currently it is not listed but ignored.
- **misc:** Contains the chromsizes extracted from the BAM files and the insert sizes per sample in case of paired-end data.
- **multiqc:** The multiQC summary report summarizing all fastqc, trim, alignment and FRiP stats.
- **peaks:** The narrowPeak and summit BED files per sample from `macs2`.
- **pipeline_info:** A file `command_lines.txt` summarizing all used command lines per process and sample, and a file `sortware_versions.txt` summarizing all software versions that were used in the pipeline.

## Resources

The [nextflow.config](nextflow.config) files contains hardcoded defaults towards resources for the individual processes, suitable for use on HPC or workstation environments. The most demanding process is the alignment steps, requiring 16 threads and 16GB of RAM per sample to finish in a reasonable amount of time.  

## Schedulers

The [schedulers.config](configs/schedulers.config) file currently contains a single scheduler profile for SLURM as used on or HPC,
submitting jobs (if using `-profile slurm`) to a quere called `normal` with a maximum 8h of walltime. Custom profiles should be added to this config. Users can add custom configurations here.

## Software

For reproducibility we recommend to use the container options (Docker, Singularity, Apptainer) for `-profile docker,singularity,apptainer` to take care of all software, using a provided [Docker-based image](https://hub.docker.com/r/atpoint/atac_chip_preprocess/tags). By default, the pipeline does not support conda/mamba. However, the user can create such a software environment with the `environment.yml` file, and then simply run the pipeline without specifying any of the above container engines, so the software available in the current environment will be used.