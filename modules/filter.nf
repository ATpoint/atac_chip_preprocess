// Filter BAM files post alignment

process FilterBam {

    tag "$sample_id"

    cpus params.bamfilter_threads 
    memory '1.GB' // this process should not have any notable footprint

    publishDir params.bamfilter_dir, mode: params.bamfilter_pubmode

    input:
    tuple val(sample_id), path(bam), path(bai)                  

    output:
    tuple val("${sample_id}"), path("*.bam"), path("*.bai"), emit: bam    

    script:

    // catch edge case when empty params.bamfilter_keepchr was provided (--params.bamfilter_keepchr '')
    // to turn off chr filtering but nf interpreted it as boolean:
    if(params.bamfilter_keepchr == true) { regex = "''" } else { regex = params.bamfilter_keepchr }

    """

    samtools idxstats $bam \
    | cut -f1 | grep $regex | grep -v '*' \
    | xargs samtools view --write-index $params.flag_keep $params.flag_remove \
        $params.bamfilter_additional -@ $task.cpus -o ${sample_id}_filtered.bam##idx##${sample_id}_filtered.bam.bai $bam

    samtools flagstat -@ $params.align_threads ${sample_id}_filtered.bam > ${sample_id}_filtered.flagstat        
    
    """

}
