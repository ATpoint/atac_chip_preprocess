---
editor_options: 
  markdown: 
    wrap: 72
---

# atac_chip_preprocess

This Nextflow pipeline is intended for preprocessing of DNA-seq
experiments such as ATAC-seq, ChIP-seq and CUT&RUN. It currently
performs trimming and alignment of fastq files, filtering of the
resulting BAM files, basic peak calling and a QC assessment by
calculating Fractions Of Reads Per Peaks (FRiPs). Specific CUT&RUN
features might be added in later versions, the defaults should work well
for most input data though.

## Usage

Up front: Do not use `.` as delimiters in file- or folder names. We
expect the dot to separate file suffixes from its basenames, not as
delimiter anywhere else!

It is recommended to first create an index and then in a second run
perform all other steps. All parameters are explained below. We provide
an `environment.yml` file which can be used to create a conda or mamba
environment that should solve on both Linux and macOS (Mac M1 chips were
not tested). There is also a Docker image based on this environment
available at the [Docker
Hub](https://hub.docker.com/r/atpoint/atac_chip_preprocess). One can use
the `-profile` argument of Nextflow to take care of the software
dependencies:  
There is `-profile conda/docker/singularity` which automatically creates
the conda environment or pulls the respective image. Submission via the
SLURM scheduler is possible with `-profile slurm`. If no profile is set
then the software is expected locally in `$PATH`. For testing there are
two profiles `-profile test_single/test_paired` which use some minimal
example data in the `./test` folder, both can be used with or without
the `--atacseq` option, see below for more details.

<br>

### Indexing

This process indexes the genome with `bowtie2-build`. The minimal
indexing command is:

``` bash
nextflow run main.nf --skip_align --ref_genome '/path/to/genome.fa.gz' --idx_name '/path/to/index_directory'
```

This will run `bowtie2-build` on the reference genome (expected
`*.fa.gz`) and then move the index files (named `idx.*`) to the folder
specified in the `--idx_name` argument, e.g.   
`--idx_name /scratch/username/index_dir/bowtie2_index_mm10/`.  
Do not use `.` as delimiter of the index name as we expect that the dot
separates the basename from the `.bt2` suffix.

**Options with defaults:**

-   `--idx_threads 1`  
    =\> number of threads
-   `--idx_mem '8.GB'`  
    =\> memory allocation for this process
-   `--idx_additional ''`  
    =\> additional arguments for `bowtie2-build`, escape first one with
    `\` e.g.  
    `--idx_additional '\--large-index'`
-   `--idx_name 'bowtie2_idx'`  
    =\> name (full path) of the output folder to contain the bt2 index
    files, default is `bowtie2_index` in the `$launchDir` of the
    Nextflow run
-   `--idx_pubmode 'relink'`  
    =\> publishDir mode for Nextflow

Alternatively, if a bowtie2 index has already been made with this (or
any other) command or downloaded from a repository such as the [bowtie2
SourceForge
website](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml), then
it can be specified via the `--idx` argument, e.g.   
`--idx /path/to/index/prefix` so basically the `-x` argument of bowtie2,
e.g. `--idx /scratch/user/index/genome` given that there are index files
witn `genome.*bt2`.

Note that `--idx_pubmode 'move'` is not compatible with using the
created index during the same Nextflow run by the align process. If
using 'move' then one has to specify the path to the index and invoke a
second run. It is therefore preferred to first build an index with
'move' (so it can be long-term stored somewhere) and then run the
alignment in a second run. Alternatively, if indexing and alignment
should happen in a single run then use the pubmode 'rellink' or any
other of the link types
[compatible](https://www.nextflow.io/docs/latest/process.html#publishdir)
with Nextflow caching.

The minimal commands for the processes downstream of the indexing are
shown below:

``` bash
#/ minimal command for paired-end data:
nextflow run main.nf --idx '/full/path/to/index_basename' --fastq '/path/to/fastq_folder/*{1,2}_.fastq.gz' --mode paired

#/ minimal command for single-end data:
nextflow run main.nf --idx '/full/path/to/index_basename' --fastq '/path/to/fastq_folder/*.fastq.gz' --mode single
```

If this is ATAC-seq data then add `--atacseq` to the command. This will
then trigger some defaults specific to ATAC-seq experiments such as the
using the Nextera adapter for fastq trimming, extraction of transposase
cutting sites and slightly different peak calling parameters for the QC

Below we list all available params that can be used to customize the
run.

<br>

### Trimming & Alignment

This process trims the input fastq files, expected `*_{1,2}.fastq.gz`
for paired-end and `*.fastq.gz` for single-end data, with `cutadapt`.
This is currently not optional as trimming and alignment is hardcoded as
a Unix pipe. The output is streamed into `bowtie2` for alignment and
`salblaster` for duplicate marking (not removal, this optionally comes
in the next process). The resulting file is then sorted, indexed and
flagstated with `samtools`. The output name is `*_raw.bam`. The
following params can be used to customize the run:

**Options with defaults:**

-   `--skip_align 'false'`  
    =\> whether to skip all processes and only run the indexing
-   `--mode paired`  
    =\> "paired" or "single" for paired- or single-end data
-   `--trim_adapter ''`  
    =\> adapter sequence to trim via `cutadapt`, the default (when
    argument is not set) is the TruSeq adapter (`AGATCGGAAGAGC)`. When
    `--atacseq` is set then defaults to the Nextera adapter
    (`CTGTCTCTTATACACATCT`). Enter a sequence to use a custom adapter.
-   `--trim_additional ''`  
    =\> additional arguments for `cutadapt` beyond `--quiet -j -a -A -o`
    which are already set
-   `--align_threads 1`  
    =\> threads for alignment
-   `--align_mem  '8.GB'`  
    =\> allocated memory for alignment, the defaults should be enough
    for human and mouse, expected in GB.
-   `--align_additional '\--very-sensitive'`=\> additional parameters
    for bowtie2 alignment beyond `-q --threads --rg-id -x`
-   `--align_dir $(realpath ./bam_raw/)`  
    =\> output directory for the resulting bam file, index and flagstat
-   `--align_pubmode 'rellink'`=\> publish mode for
    [publishDir](https://www.nextflow.io/docs/latest/process.html#publishdir)
-   `--sort_threads 1`  
    =\> threads for `samtools sort`
-   `--sort_mem '1G'`  
    =\> memory per thread for `samtools sort`. This must be in the
    format recognized by the `-@` argument of `samtools`, so <value><G>.
    This is expected to be in Gigabytes.
-   `--sort_additional ''` =\> additional arguments for samtools sort
    beyond `-@ -m -o --write-index` which are already set

Note that the total memory required for this process is:

`sort_mem * sort_threads + align_mem`

<br>

### Filtering

This process performs filtering on the sorted `*_raw.bam` file. This
typically involves removal of non-primary chromosomes (unplaced
scaffolds/contigs) and chrM, MAPQ filtering and duplicate reloval. There
are four options that can be used to define the filtering which wrap
around the `-q`, `-f` and `-F` options of `samtools view` to filter MAPQ
and by presence/absence of the SAM flags, and an option to specify which
chromosomes (=alignments to that chromosome) to remove from the BAM. For
appropriate flags see
<https://broadinstitute.github.io/picard/explain-flags.html>

**Options with defaults:**

-   `--keep_chr 'chr[1-9,X,Y,EBV]'`  
    =\> a regex or search string for chromosomes to keep. Technically
    this string is parsed to a `grep` command on the list of
    chromosomes. The matched chromosomes will be retained, all others
    will be removed. It can be anything that Unix `grep` can accept. The
    default regex keeps all chromosomes prefixed with `chr` followed by
    any number as well as X and Y. It also keeps (if present) the EBV
    decoy. This is sufficient for GENCODE-formatted reference genomes
    where unplaced contigs etc are called either chrU or something like
    GL(...) or JH(...). Change this regex in case dhromosome identifiers
    are different in your reference. Use `grep '^>'` on your reference
    to get all present identifiers. One can set `--keep_chr ''` to keep
    all chromosomes without any filtering.
-   `--bamfilter_mapq 20`  
    =\> integer, keep only alignments with MAPQ above this value. Set to
    0 to keep all everything regardless of MAPQ.
-   `--bamfilter_flag_keep ''`  
    =\> a SAM flag for `samtools view -f`, so keep alignments with this
    flag set. Default is `''` and in this case means keep all mapped
    reads for single-end data and all mapped and paired reads for
    paired-end data. Set to 'nofilter' to deactivate this filter. That
    would be 0 for single-end data and 1 for paired-end data.
-   `--bamfilter_flag_remove ''`  
    =\> a SAM flag for `samtools view -F`, so remove alignments with
    this flag set. Default is `''` and in this case means a flag of 3332
    soremove all non-primary and supplementary alignments, optical/PCR
    duplicates and unmapped reads. Set to `nofilter` to deactivate this
    filter.
-   `--skip_bamfilter`  
    =\> boolean, when set skips the filtering step. In this case the run
    stops after the alignment. If no filtering is desired on the BAM
    file while downstream processes still should run then use:
    `--bamfilter_mapq 0 --bamfilter_flag_keep 'nofilter' --bamfilter_flag_remove 'nofilter' --bamfilter_keepchr "''"`,
    but this is quite unusual when processing ATAC/ChIP-seq etc data so
    there is no shortcut for it.
-   `--bamfilter_additional`  
    =\> additional arguments to give to `samtools view` which runs the
    actual filtering beyond options `-@ -m -f -F --write-index -q -o.`
    As usual for additional options the first one must be escaped, e.g.
    `--bam_additional '\--verbosity 2`
-   `--bamfilter_dir $(realpath ./bam_filtered/)`  
    =\> output directory for the files
-   `--bamfilter_pubmode 'rellink'`=\> publish mode for
    [publishDir](https://www.nextflow.io/docs/latest/process.html#publishdir)

<br>

### InsertSizes

If in paired-end mode this process will use `CollectInsertSizeMetrics`
from Picard to collect the paired-end insert sizes (TLENs) which can be
used as a QC e.g. in ATAC-seq which should show the characteristic
banding pattern.

**Options with defaults:**

-   `--skip_isizes`  
    =\> boolean, turn this process of
-   `--isizes_mem '4.GB'`  
    =\> probably no need to ever change this, process should have little
    memory footprint
-   `--isizes_dir $(realpath ./insert_sizes/)`  
    =\> output folder for files
-   `--isizes_pubmode 'rellink'`  
    =\> publish mode for
    [publishDir](https://www.nextflow.io/docs/latest/process.html#publishdir)

<br>

### Cutsites extraction

This process runs if the `--atacseq` param is set. It takes the filtered
BAM file from above and extracts the transposase integration or
"cutting" events which are then 5'-ends of the reads, shifted by +4/-5bp
to account for the 9bp duplication event that the Tn5 creates when
integrating into the target site in the genome. This file can be useful
when plotting insertion frequencies around transcription factor motifs
but in this case requires additional normalization for sequencing depth
and composition. Here it is only used as input for for peak calling/QC.

**Options with defaults:**

-   `--atacseq`  
    =\> boolean, whether to turn on ATAC-seq mode and run this process
-   `--cutsites_threads 3`  
    =\> threads for this process. Defaults to 3 as the process
    internally is a pipe of three tools. If providing more than 3 then
    the additional threads will be used for the sorting step of the
    resulting file. In the testing profiles
-   `--cutsites_mem '1G'`  
    =\> memory for the GNU sort step, in a format that the `-S` option
    of GNU sort accepts, e.g. `1G`.
-   `--cutsites_dir $(realpath ./bed_cutsites/)`  
    =\> output directory for the `cutsites.bed.gz` file
-   `--cutsites_pubmode 'rellink'`  
    =\> publish mode for
    [publishDir](https://www.nextflow.io/docs/latest/process.html#publishdir)

<br>

### FRiPs

This process calls peaks with `macs2` and then uses `featureCounts` from
the [subread package](http://subread.sourceforge.net/) on these peaks to
calculate the Fractions Of Reads Per Peak (FRiPs) as a proxy for data
quality. It is basically a measure of signal to noise ratio, as reads
overlapping peaks are signal and all other reads are noise. For ChIP-seq
this can have a wide range depending on protein abundance, antibody
quality, moon phase and the current mood of the ChIP gods. It can well
range from 0.01 to something like 0.2 - 0.3. It is encouraged to always
check data on a genome browser and see by eye whether there is a good
separation between peaks and noise. For ATAC-seq this should well be \>
0.1. On fresh ex vivo and cell line material from mice and human (cells
similar to hematopoietic progenitors from bone marrow or PBMCs) we
usually get FRiPs up to 0.5 using the
[OmniATAC](https://www.nature.com/articles/nmeth.4396) protocol.

This is celltype-dependent and might be notably different in other
celltypes/organs/tissues.

**Options with defaults:**

-   `--macs_additional ''`  
    =\> additional parameters given to `macs2 callpeak`. For `--atacseq`
    there is a default when this option is left empty which is:  
    `--nomodel --extsize 100 --shift -50 --keep-dup=all --min-length 150 -q 0.005`
    using the extracted cutsite BED file for peak calling. If
    `--atacseq` is not set then it defaults to:  
    `--keep-dup=all --min-length 150`.  
    Allowed parameters for this argument is everything beyond
    `-t -c -n -f -g`. For a broad-peak ChIP-seq dataset one could use
    e.g.:  
    `--macs_additional '\--broad --keep-dup=all --min-length`
-   `--macs_gflag 'mm'`  
    =\> parameter `-g` in `macs2 callpeak` so either an in-built genome
    size flag (hs, mm, ce, dm) or the effective genome size as a plain
    number or scientific, e.g. human (if not using hs preset) would be
    2.7e9
-   `--macs_format ''`  
    =\> the `-f` option of `macs2 callpeak`, is automatically set and
    should not require manual change
-   `--macs_suffix ''`  
    =\> an optional suffix to append to the output files. Default would
    be e.g. `<basename>_peaks.narrowPeak`, and with
    `--suffix '_coolsuffix'` it would be
    `<basename>_coolsuffix_peaks.narrowPeak`. Mind that the delimiter
    must be provided so here the underscore in `_coolsuffix`.
-   `--macs_control ''`  
    =\> optional control file, currently experimental, don't use this.
    The peak calling done here is mainly for per-sample QC (FRiPs) and
    the final peak list should be created outside. For ATAC-seq we
    generally find the Genrich peak caller useful. For ChIP-seq one
    often uses an input/IgG as control.
-   `--macs_dir $(realpath ./macs2/)`  
    =\> output directory for peak files
-   `--macs_mem '4.GB'`  
    =\> allocated memory, should not be necessary to change that
-   `--macs_pubmode 'rellink'`  
    =\> publish mode for
    [publishDir](https://www.nextflow.io/docs/latest/process.html#publishdir)
-   `--frips_threads 1`  
    =\> threads for `featureCounts` to calculate FRiPs based on the
    per-sample count matrix it builds using the peaks from `macs2`
-   `--frips_additional ''`  
    =\> additional arguments for `featureCounts`, leave it empty, will
    be set automatically depending on `--mode` and `--atacseq`
-   `--frips_mem '4.GB'`=\> memory allocation, probably no need to
    change that
-   `--frips_dir $(realpath ./frips/)`  
    =\> output folder for files
-   `--frips_pubmode 'rellink'`   
    =\> publish mode for
    [publishDir](https://www.nextflow.io/docs/latest/process.html#publishdir)  

<br>

## Citations

-   [nf-core project](https://nf-co.re/)

-   [Ewels et al (2020) The nf-core paper. Nature Biotechnology volume
    38, pages
    276–278](https://www.nature.com/articles/s41587-020-0439-x)

-   [Nextflow Docs](https://www.nextflow.io/docs/latest/index.html#)

-   [Seqera Training](https://seqera.io/training/)

-   [https://github.com/nextflow-io/rnaseq-nf -- The Seqera Labs DSL2
    proof-of-concept workflow](https://github.com/nextflow-io/rnaseq-nf)

-   [Merkel, D (2014). Docker: lightweight linux containers for
    consistent development and deployment. Linux
    Journal](https://dl.acm.org/doi/10.5555/2600239.2600241)

-   [Kurtzer et al (2017) Singularity: Scientific containers for
    mobility of compute. PLoS
    ONE](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0177459)

-   [Grüning et al (2018) Bioconda: sustainable and comprehensive
    software distribution for the life sciences. Nat Methods
    15:475-476](https://www.nature.com/articles/s41592-018-0046-7)
