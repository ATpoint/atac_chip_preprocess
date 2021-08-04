// produce a BED file with cutsites:

process Cutsites {

    tag "$sample_id"

    cpus params.threads 
    memory params.cutsites_mem 

    publishDir params.cutsites_dir, mode: params.cutsites_pubmode

    input:
    tuple val(sample_id), path(bam), path(bai)                  

    output:
    tuple val("${sample_id}"), path("${sample_id}_cutsites.bed.gz"), emit: bed

    script:

    // depending on total allocated CPUs give the sort some extra
    if(task.cpus == 2) { thready = 1 } else { thready = (task.cpus - 2) }

    """

    #/ this is the same for paired-end and single-end as we count reads and not pairs/fragments:

    bedtools bamtobed -i $bam \
    | $baseDir/bin/shift_reads.sh /dev/stdin \
    | sort -k1,1 -k2,2n -k3,3n -k6,6 -S ${params.cutsites_mem} --parallel=${thready} \
    | bgzip -@ 1 > ${sample_id}_cutsites.bed.gz
    
    """

}