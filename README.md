# atac_chip_preprocess

This Nextflow pipeline is intended for preprocessing of DNA-seq
experiments such as ATAC-seq, ChIP-seq and CUT&RUN. It currently
performs trimming and alignment of fastq files, filtering of the
resulting BAM files, basic peak calling and a QC assessment by
calculating Fractions Of Reads Per Peaks (FRiPs). Specific CUT&RUN
features might be added in later versions, the defaults should work well
for most input data though.

<br>

![CI](https://github.com/ATpoint/atac_chip_preprocess/actions/workflows/basic_test.yml/badge.svg)

<br>

## Usage

...documentation will follow...

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
