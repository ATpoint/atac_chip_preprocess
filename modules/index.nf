// Index a genome with bowtie2

process Bowtie2Idx {

    cpus   params.idx_threads
    memory params.idx_mem

    errorStrategy 'finish'

    publishDir params.idx_name, mode: params.publishmode

    input:
    path(genome)
        
    output:
    path 'idx.*', emit: idx
        
    script: 
    // inspired from https://github.com/nf-core/modules/blob/master/modules/bowtie2/build/main.nf
    """
    bowtie2-build -q --seed 1234 --threads $task.cpus $params.idx_additional $genome idx
    """                

}