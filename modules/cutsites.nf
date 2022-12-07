// This process extracts ATAC-seq insertion events from alignments

process Cutsites {

    tag "$meta.id"

    label 'process_cutsites'

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
    tuple val(meta), path("${meta.id}_cutsites.bed.gz"), emit: tuple_cutsites
    tuple path("versions.txt"), path("command_lines.txt"), emit: versions
    
    script:
    def threads_sort = (task.cpus==1) ? 1 : (task.cpus/2).round(0)
    
    """
    bedtools bamtobed -i $bam \
    | mawk 'OFS=\"\t\" {if (\$6 == \"+\") print \$1, \$2+4, \$2+5, \".\", \".\", \$6} {if (\$6 == \"-\") print \$1, \$3-5, \$3-4, \".\", \".\", \$6}' \
    | sort -k1,1 -k2,2n -k3,3n -k6,6 -S 2G --parallel=$threads_sort > ${meta.id}_cutsites.bed
    
    bgzip -l 6 -@ $task.cpus ${meta.id}_cutsites.bed

    echo ${task.process}:${meta.id} > command_lines.txt
    cat .command.sh | grep -vE '^#!/bin|versions.txt\$|command_lines.txt\$|cat \\.command.sh' | sed 's/  */ /g' | awk NF >> command_lines.txt

    echo 'bedtools:' \$(bedtools --version | head -n 1 | cut -d " " -f2) > versions.txt
    echo 'bgzip:' \$(bgzip --version | head -n 1 | cut -d " " -f3) >> versions.txt
    echo 'mawk:' \$(mawk -W version 2>&1 | head -n 1 | cut -d " " -f2,3 | tr " " "-") >> versions.txt
    echo 'sort:' \$(sort --version 2>&1 | head -n1 | cut -d " " -f4) >> versions.txt
    """

}