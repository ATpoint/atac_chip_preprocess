// Create genome-decoyed index with salmon and a tx2gene mapping table

process Tximport {

    cpus   1
    memory params.mem

    publishDir params.outdir, mode: params.pubmode

    input:
    tuple val(sample_id), path(quants)
    path(tx2gene)
        
    output:
    path("*_counts.txt")
    path("*_lengths.txt")
    path("*_infreps.txt.gz"), optional: true
            
    script: 
    """

    Rscript --vanilla $baseDir/bin/tximport.R $quants $sample_id $tx2gene

    """      

}