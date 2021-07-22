// Trim reads with cutadapt

process Trim {

    tag "$sample_id"

    cpus   params.threads
    memory params.mem

    publishDir params.outdir, mode: params.pubmode

    input:
    tuple val(sample_id), path(reads)
            
    output:
    tuple val(sample_id), path("*.fastq.gz"), emit: trimmed
        
    script: 

    def readfiles_in = (params.mode=='single') ? "$reads" : "${reads[0]} ${reads[1]}"

    def readfiles_out = (params.mode=='single') ? "-o ${sample_id}_trimmed.fastq.gz" : "-o ${sample_id}_trimmed_1.fastq.gz -p ${sample_id}_trimmed_2.fastq.gz"

    def adapter = (params.mode=='single') ? "-a $params.adapter_seq" : "-a $params.adapter_seq -A $params.adapter_seq"

    """    
    cutadapt --quiet \
        --max-n 0.1 -m 25 $params.additional -j $task.cpus $adapter $readfiles_out $readfiles_in
    """                

}