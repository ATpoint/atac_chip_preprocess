// produce a BED file with cutsites:

process InsertSizes {

    tag "$sample_id"

    cpus 1
    memory params.isizes_mem

    publishDir params.isizes_dir, mode: params.isizes_pubmode

    input:
    tuple val(sample_id), path(bam), path(bai)                  

    output:
    path("${sample_id}_InsertSizes.txt")

    script:

    """
    picard CollectInsertSizeMetrics \
        --INPUT ${bam} \
        --OUTPUT ${sample_id}_InsertSizes.txt \
        --Histogram_FILE /dev/null \
        --QUIET true --VERBOSITY ERROR --VALIDATION_STRINGENCY LENIENT 2> /dev/null
    """

}