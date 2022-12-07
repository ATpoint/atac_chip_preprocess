// This process extracts chromosome sizes from BAM file headers

process Chromsizes {

    cpus 1 
    memory 500.MB

    errorStrategy 'finish'
    
    publishDir = [
        path: params.outdir,
        mode: params.publishmode,
        saveAs: { filename -> filename.equals("versions.txt") || filename.equals("command_lines.txt") ? null : filename } 
    ]

    container params.container

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    path("chromsizes.txt"), emit: chromsizes
    tuple path("versions.txt"), path("command_lines.txt"), emit: versions
    
    script:

    """
    samtools view -H $bam | grep "^@SQ" | cut -f2,3 | awk 'OFS="\\t" {gsub("SN:|LN:","");print}' | sort -k1,1 > chromsizes.txt

    echo ${task.process}:${meta.id} > command_lines.txt
    cat .command.sh | grep -vE '^#!/bin|versions.txt\$|command_lines.txt\$|cat \\.command.sh' | sed 's/  */ /g' | awk NF >> command_lines.txt

    echo 'awk:' \$(awk -W version 2>&1 | head -n 1 | cut -d " " -f2,3 | tr " " "-") > versions.txt
    echo 'grep:' \$(grep --version 2>&1 | head -n1 | cut -d " " -f4) >> versions.txt
    echo 'samtools:' \$(samtools --version 2>&1 | head -n 1 | cut -d " " -f2) >> versions.txt
    echo 'sort:' \$(sort --version 2>&1 | head -n1 | cut -d " " -f4) >> versions.txt
    """

}