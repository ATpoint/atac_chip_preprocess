// This process produces unscaled (that is non-normalized) bigwig tracks

process Bigwigs {

    tag "$meta.id"

    label 'process_bigwigs'

    errorStrategy 'finish'
    
    publishDir = [
        path: params.outdir,
        mode: params.publishmode,
        saveAs: { filename -> filename.equals("versions.txt") || filename.equals("command_lines.txt") ? null : filename } 
    ]

    container params.container

    input:
    tuple val(meta), path(data) // meta map and then either BAM or cutsites BED
    path(chromsizes)

    output:
    path("${meta.id}_extended_unscaled.bigwig"), emit: bigwig
    tuple path("versions.txt"), path("command_lines.txt"), emit: versions
    
    script:
    def pc = !params.atacseq & !meta.single_end ? '-pc' : ''
    def fs = !params.atacseq & meta.single_end ? "-fs ${params.fragment_length}" : ''

    if(params.atacseq){

        // For ATAC-seq first extend cutting sites to 100bp (smoothing), then make the pipelup
        """
        bedtools slop -b 50 -g $chromsizes -i $data | bedtools genomecov -bga -i - -g $chromsizes | sort -k1,1 -k2,2n --parallel $task.cpus > "${meta.id}_extended_unscaled.bedGraph"

        bedGraphToBigWig "${meta.id}_extended_unscaled.bedGraph" $chromsizes "${meta.id}_extended_unscaled.bigwig"

        echo ${task.process}:${meta.id} > command_lines.txt
        cat .command.sh | grep -vE '^#!/bin|versions.txt\$|command_lines.txt\$|cat \\.command.sh' | sed 's/  */ /g' | awk NF >> command_lines.txt

        echo 'bedGraphToBigWig' \$(bedGraphToBigWig 2>&1 | head -n1 | cut -d " " -f3) > versions.txt
        echo 'bedtools:' \$(bedtools --version | head -n 1 | cut -d " " -f2) >> versions.txt
        """

    } else {
        
        // For non ATAC-seq extend reads to fragments or use paired-end fragment size directly
        """
        bedtools genomecov -bga -ibam $data $pc $fs | sort -k1,1 -k2,2n --parallel $task.cpus > "${meta.id}_extended_unscaled.bedGraph"

        bedGraphToBigWig "${meta.id}_extended_unscaled.bedGraph" $chromsizes "${meta.id}_extended_unscaled.bigwig"

        echo ${task.process}:${meta.id} > command_lines.txt
        cat .command.sh | grep -vE '^#!/bin|versions.txt\$|command_lines.txt\$|cat \\.command.sh' | sed 's/  */ /g' | awk NF >> command_lines.txt

        echo 'bedGraphToBigWig:' \$(bedGraphToBigWig 2>&1 | head -n1 | cut -d " " -f3) > versions.txt
        echo 'bedtools:' \$(bedtools --version | head -n 1 | cut -d " " -f2) >> versions.txt
        """

    }

}