// This process collects insert size metrics of paired-end alignments

process Isizes {

    tag "${meta.id}"

    cpus 1
    memory 1.GB

    errorStrategy 'finish'
    
    publishDir = [
        path: params.outdir,
        mode: params.publishmode,
        saveAs: { filename -> 
                    (filename.equals("versions.txt") || 
                     filename.equals("command_lines.txt")) ? null : filename } 
    ]

    container params.container

    input:
    tuple val(meta), path(bam), path(bai) 

    output:
    path("${meta.id}_isizes.txt"), emit: isizes
    tuple path("versions.txt"), path("command_lines.txt"), emit: versions
    
    script:
    """
    picard CollectInsertSizeMetrics --INPUT ${bam} --OUTPUT ${meta.id}_isizes.txt --Histogram_FILE /dev/null --QUIET true --VERBOSITY ERROR --VALIDATION_STRINGENCY LENIENT 2> /dev/null

    echo ${task.process}:${meta.id} > command_lines.txt
    cat .command.sh | grep -vE '^#!/bin|versions.txt\$|command_lines.txt\$|cat \\.command.sh' | sed 's/  */ /g' | awk NF >> command_lines.txt

    echo 'picard:' \$(picard CollectInsertSizeMetrics --version 2>&1 | grep -v setlocale | cut -d ":" -f2) > versions.txt
    """

}