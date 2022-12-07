// This process reads a tuple with meta map and fastq files, runs alignment, duplicate marking, sorts the BAM file and creates a flagstat.

process Align {

    tag "${meta.id}"

    label 'process_align'

    errorStrategy 'finish'
    
    publishDir = [
        path: params.outdir,
        mode: params.publishmode,
        saveAs: { filename -> filename.equals("versions.txt") || filename.equals("command_lines.txt") ? null : filename } 
    ]

    container params.container

    input:
    tuple val(meta), path(reads)
    path(idx)     

    output:
    tuple val(meta), path("${meta.id}_raw.bam"), path("${meta.id}_raw.bam.bai"), emit: tuple_bam
    tuple val(meta), path("${meta.id}_raw.flagstat"), emit: tuple_flagstat
    tuple path("versions.txt"), path("command_lines.txt"), emit: versions
    
    script:
    def sambamba_memory = task.memory.toString().replaceAll(".GB", "G").replaceAll(".MB", "M")
    def readin = !meta.single_end ? "-1 ${reads[0]} -2 ${reads[1]}" : "-U $reads"
    """
    use_idx=\$(find ${idx}/ -name '*.bt2' | awk -F \".\" '{print \$1 | \"sort -u\"}')

    bowtie2 -x \${use_idx} $params.args --threads $task.cpus --rg-id ${meta.id} $readin \
    | samblaster --quiet --ignoreUnmated \
    | samtools view -O bam,level=1 -@ 1 -o ${meta.id}_raw_unsorted.bam
        
    sambamba sort $params.args2 --tmpdir ./ -t $task.cpus -m $sambamba_memory -o ${meta.id}_raw.bam ${meta.id}_raw_unsorted.bam
       
    samtools flagstat -@ $task.cpus ${meta.id}_raw.bam > ${meta.id}_raw.flagstat
        
    echo ${task.process}:${meta.id} > command_lines.txt
    cat .command.sh | grep -vE '^#!/bin|versions.txt\$|command_lines.txt\$|cat \\.command.sh' | sed 's/  */ /g' | awk NF >> command_lines.txt

    echo 'bowtie2:' \$(bowtie2 --version | head -n1 | cut -d " " -f3) > versions.txt
    echo 'sambamba:' \$(sambamba 2>&1 | head -n2 | awk NF | cut -d " " -f2) >> versions.txt
    echo 'samblaster:' \$(samblaster -h 2>&1 | head -n1 | cut -d " " -f3) >> versions.txt
    echo 'samtools:' \$(samtools --version 2>&1 | head -n 1 | cut -d " " -f2) >> versions.txt
    """

}
    
