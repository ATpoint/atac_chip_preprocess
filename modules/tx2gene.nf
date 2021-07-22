// tx2gene map for tximport

process TX2GENE {

    cpus   1
    memory '1G'

    publishDir params.outdir, mode: params.pubmode

    input:
    path(gtf)
            
    output:
    path("tx2gene.txt")
        
    script:    
    """
    Rscript --vanilla ${baseDir}/bin/tx2gene.R \
        annot.gtf.gz tx2gene.txt params.transcript_id params.transcript_name params.gene_id params.gene_name params.gene_type
    """                

}