// This process runs fastp for adapter and quality trimming

process Trim {

    tag "${meta.id}"   

    label 'process_trim'

    errorStrategy 'finish'
    
    publishDir = [
        path: params.outdir,
        mode: params.publishmode,
        saveAs: { filename -> 
            (filename.endsWith("fq.gz") & params.keep) || 
            (filename.endsWith(".html")) ||
            (filename.endsWith(".json")) ? filename : null
        }
    ]

    container params.container

    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("${meta.id}_trimmed*.fq.gz"), emit: tuple_trimmed
    path("${meta.id}_fastp.html"), emit: html
    path("${meta.id}_fastp.json"), emit: json
    tuple path("versions.txt"), path("command_lines.txt"), emit: versions
    
    script:

    def indata  = !meta.single_end ? "--in1 ${reads[0]} --in2 ${reads[1]}" : "--in1 ${reads}"
    def outdata = !meta.single_end ? "--out1 ${meta.id}_trimmed_R1.fq.gz --out2 ${meta.id}_trimmed_R2.fq.gz" : "--out1 ${meta.id}_trimmed.fq.gz"

    """
    fastp $params.args $indata $outdata --thread $task.cpus --json "${meta.id}_fastp.json" --html "${meta.id}_fastp.html"
        
    echo ${task.process}:${meta.id} > command_lines.txt
    cat .command.sh | grep -vE '^#!/bin|versions.txt\$|command_lines.txt\$|cat \\.command.sh' | sed 's/  */ /g' | awk NF >> command_lines.txt

    echo 'fastp:' \$(fastp --version 2>&1 | cut -d " " -f2) > versions.txt
    """

}
