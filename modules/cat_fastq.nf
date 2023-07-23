// This process concatenates lane replicates / technical replicates of fastq files into a single file

process CatFastq {

    tag "$meta.id"

    cpus   1
    memory 1.GB

    errorStrategy 'finish'

    publishDir = [
        path: params.outdir,
        mode: params.publishmode,
        saveAs: { filename -> 
            (filename.endsWith("fq.gz") & params.keep) ? filename : null
        } 
    ]

    container params.container

    input:
    tuple val(meta), path(reads, stageAs: "?/*")
            
    output:
    tuple val(meta), path("*merged*.fq.gz"), emit: fastq_tuple
    tuple path("versions.txt"), path("command_lines.txt"), emit: versions
    
    script: 

    // reads can either be a BlankSeparatedList or TaskPath
    def readList = reads instanceof List ? reads.collect{ it.toString() } : [reads.toString()]
    def rs = readList.size()

    // 
    if(meta.single_end){

        r1x = readList
        r2x = ''

    } else {

        if(rs>2){

            r1x = readList[0..(rs/2-1)]
            r2x = readList[(rs/2)..(rs-1)]

        } else {

            r1x = readList[0]
            r2x = readList[1]

        }

    }

    r1 = r1x.toString().replaceAll(",|\\[|\\]", "")
    r2 = r2x.toString().replaceAll(",|\\[|\\]", "")

    r1_cat = "${meta.id}_merged_R1.fq.gz"
    r2_cat = "${meta.id}_merged_R2.fq.gz"

    if(meta.single_end){
        
        """
        cat ${r1} > ${r1_cat}

        echo ${task.process}:${meta.id} > command_lines.txt
        cat .command.sh | grep -vE '^#!/bin|versions.txt\$|command_lines.txt\$|cat \\.command.sh' | sed 's/  */ /g' | awk NF >> command_lines.txt

        echo 'cat:' \$(cat --version | head -n 1 | cut -d ")" -f2) > versions.txt
        """     

    } else {

        """
        cat ${r1} > ${r1_cat}
        cat ${r2} > ${r2_cat}
        
        echo ${task.process}:${meta.id} > command_lines.txt
        cat .command.sh | grep -vE '^#!/bin|versions.txt\$|command_lines.txt\$|cat \\.command.sh' | sed 's/  */ /g' | awk NF >> command_lines.txt

        echo 'cat:' \$(cat --version | head -n 1 | cut -d ")" -f2) > versions.txt
        """     

    } 

}
