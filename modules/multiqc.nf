// This process runs MultiQC

process Multiqc {

    label 'process_multiqc'

    errorStrategy 'finish'

    publishDir = [
        path: params.outdir,
        mode: params.publishmode,
        saveAs: { filename -> filename.equals("versions.txt") || filename.equals("command_lines.txt") ? null : filename } 
    ]

    container params.container

    input:
    path(everything)
            
    output:
    path("multiqc*")
    tuple path("versions.txt"), path("command_lines.txt"), emit: versions
    
    script: 
    """
    multiqc .

    echo ${task.process}: > command_lines.txt
    cat .command.sh | grep -vE '^#!/bin|versions.txt\$|command_lines.txt\$|cat \\.command.sh' | sed 's/  */ /g' | awk NF >> command_lines.txt

    echo 'MultiQC:' \$(multiqc --version 2>&1 | cut -d " " -f3)  > versions.txt
    """     

}
