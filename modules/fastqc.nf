// This process runs fastqc on fastq files

process Fastqc {

    tag "$meta.id"

    cpus   1
    memory 1.GB

    errorStrategy 'finish'

    publishDir = [
        path: params.outdir,
        mode: params.publishmode,
        saveAs: { filename -> filename.equals("versions.txt") || filename.equals("command_lines.txt") ? null : filename } 
    ]

    container params.container

    input:
    tuple val(meta), path(reads, stageAs: "?/*")
            
    output:
    path("*.html"), emit: html
    path("*.zip") , emit: zip
    tuple path("versions.txt"), path("command_lines.txt"), emit: versions
    
    script: 
    if(!meta.single_end){

        """
        fastqc --threads 1 -o ./ -q ${reads[0]}
        fastqc --threads 1 -o ./ -q ${reads[1]}

        echo ${task.process}:${meta.id} > command_lines.txt
        cat .command.sh | grep -vE '^#!/bin|versions.txt\$|command_lines.txt\$|cat \\.command.sh' | sed 's/  */ /g' | awk NF >> command_lines.txt

        echo 'FastQC:' \$(fastqc --version | cut -d " " -f2) > versions.txt
        """

    } else {

        """
        fastqc --threads 1 -o ./ -q $reads

        echo ${task.process}:${meta.id} > command_lines.txt
        cat .command.sh | grep -vE '^#!/bin|versions.txt\$|command_lines.txt\$|cat \\.command.sh' | sed 's/  */ /g' | awk NF >> command_lines.txt

        echo 'FastQC:' \$(fastqc --version | cut -d " " -f2) > versions.txt
        """

    }

}
