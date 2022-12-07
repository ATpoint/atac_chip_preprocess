// This process runs filtering on sorted BAM files

process Filter {

    tag "$meta.id"

    label 'process_filter'

    errorStrategy 'finish'
    
    publishDir = [
        path: params.outdir,
        mode: params.publishmode,
        saveAs: { filename -> filename.equals("versions.txt") || filename.equals("command_lines.txt") ? null : filename } 
    ]

    container params.container

    input:
    tuple val(meta), path(bam), path(bai)
    val(flag_remove)
    val(chr_regex)

    output:
    tuple val(meta), path("${meta.id}_filtered.bam"), path("${meta.id}_filtered.bam.bai"), emit: tuple_bam
    tuple val(meta), path("${meta.id}_filtered.flagstat"), emit: tuple_flagstat
    tuple path("versions.txt"), path("command_lines.txt"), emit: versions
    
    script:
    
    """
    samtools idxstats $bam \
    | cut -f1 | grep $chr_regex | grep -v '*' \
    | xargs samtools view $params.args --write-index -q $params.min_mapq -F $flag_remove -@ $task.cpus -o ${meta.id}_filtered.bam##idx##${meta.id}_filtered.bam.bai $bam

    samtools flagstat ${meta.id}_filtered.bam > ${meta.id}_filtered.flagstat

    echo ${task.process}:${meta.id} > command_lines.txt
    cat .command.sh | grep -vE '^#!/bin|versions.txt\$|command_lines.txt\$|cat \\.command.sh' | sed 's/  */ /g' | awk NF >> command_lines.txt

    echo 'cut:' \$(cut --version 2>&1 | head -n1 | cut -d " " -f4) > versions.txt
    echo 'grep:' \$(grep --version 2>&1 | head -n1 | cut -d " " -f4) >> versions.txt
    echo 'samtools:' \$(samtools --version 2>&1 | head -n 1 | cut -d " " -f2) >> versions.txt
    echo 'xargs: \$(xargs --version 2>&1 | head -n1 | cut -d " " -f4) >> versions.txt
    """
    
}
    
