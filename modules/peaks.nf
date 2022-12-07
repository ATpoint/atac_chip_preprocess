// This process calls peaks with macs2 and filters resulting peaks against provided blacklists

process Peaks {

    tag "$meta.id"

    label 'process_peaks'

    errorStrategy 'finish'
    
    publishDir = [
        path: params.outdir,
        mode: params.publishmode,
        saveAs: { filename -> filename.equals("versions.txt") || filename.equals("command_lines.txt") ? null : filename } 
    ]

    container params.container

    input:
    tuple val(meta), path(data)
    path(blacklist)

    output:
    tuple val(meta), path("${meta.id}_peaks.narrowPeak"), path("${meta.id}_peaks.bed"), emit: tuple_peaks
    tuple path("versions.txt"), path("command_lines.txt"), emit: versions
    
    script:

    if(params.atacseq){ format = 'BED' } else { format = meta.single_end ? 'BAM' : 'BAMPE' }
    def species = params.species
    def args2 = params.atacseq ? '--nomodel --extsize 100 --shift -50 --min-length 250' : ' '
    
    if(params.filter_blacklist){
        
        """
        export OPENBLAS_NUM_THREADS=1 # should keep blas library from unwanted multithreading

        macs2 callpeak $params.args $args2 --tempdir ./ -f $format -g $species -t $data -n ${meta.id}_all

        bedtools intersect -v -a ${meta.id}_all_peaks.narrowPeak -b $blacklist > ${meta.id}_peaks.narrowPeak
        bedtools intersect -v -a ${meta.id}_all_summits.bed -b $blacklist > ${meta.id}_peaks.bed

        echo ${task.process}:${meta.id} > command_lines.txt
        cat .command.sh | grep -vE '^#!/bin|versions.txt\$|command_lines.txt\$|cat \\.command.sh' | sed 's/  */ /g' | awk NF >> command_lines.txt

        echo 'bedtools:' \$(bedtools --version | head -n 1 | cut -d " " -f2) > versions.txt
        echo 'macs2:' \$(macs2 --version | cut -d " " -f2) >> versions.txt
        """

    } else {

        """
        export OPENBLAS_NUM_THREADS=1 # should keep blas library from unwanted multithreading

        macs2 callpeak $params.args $args2 --tempdir ./ -f $format -g $species -t $data -n ${meta.id}

        echo ${task.process}:${meta.id} > command_lines.txt
        cat .command.sh | grep -vE '^#!/bin|versions.txt\$|command_lines.txt\$|cat \\.command.sh' | sed 's/  */ /g' | awk NF >> command_lines.txt

        echo 'bedtools:' \$(bedtools --version | head -n 1 | cut -d " " -f2) > versions.txt
        echo 'macs2:' \$(macs2 --version | cut -d " " -f2) >> versions.txt
        """

    }
    

}

