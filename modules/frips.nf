// This process calculates Fraction Of Reads in Peaks

process Frips {

    tag "${meta.id}"

    label 'process_frips'

    errorStrategy 'finish'
    
    publishDir = [
        path: params.outdir,
        mode: params.publishmode,
        saveAs: { filename -> 
                    (filename.equals("versions.txt") || 
                     filename.equals("command_lines.txt") ||
                     filename.contains("_frips.txt")) ? null : filename } 
    ]

    container params.container

    input:
    tuple val(meta), path(data) // meta map is first element, second one is a map with bam first and peaks second

    output:
    tuple val(meta), path("${meta.id}_frips.txt"), emit: frips
    tuple path("versions.txt"), path("command_lines.txt"), emit: versions
    
    script:
    def bam = data[0]
    def peaks = data[1]
    
    // If ATAC-seq then count cutsites
    args_atac_specific = params.atacseq ? '--read2pos 5 --readExtension3 4' : '' 
    
    // If non ATAC-seq and single-end extend reads to provided fragment length, default 250
    args_single_specific = (!params.atacseq & meta.single_end & params.fragment_length != '') ? "--readExtension3 $params.fragment_length" : ""

    // If paired-end then count fragments instead of reads
    args_paired_specific = (!params.atacseq & !meta.single_end) ? "-p" : ""

    """
    awk 'OFS="\\t" {print \$1":"\$2"-"\$3, \$1, \$2, \$3, "."}' $peaks > saf

    featureCounts $args_atac_specific $args_single_specific $args_paired_specific -a saf -F SAF -T $task.cpus -o fc_out $bam

    signal=\$(grep -w 'Assigned' fc_out.summary | cut -f2)
    noise=\$(grep -w 'Unassigned_NoFeatures' fc_out.summary | cut -f2)
    paste <(echo ${meta.id}) <(bc <<< "scale=6;\${signal}/(\${signal}+\${noise})") > ${meta.id}_frips.txt

    echo ${task.process}:${meta.id} > command_lines.txt
    cat .command.sh | grep -vE '^#!/bin|versions.txt\$|command_lines.txt\$|cat \\.command.sh' | sed 's/  */ /g' | awk NF >> command_lines.txt

    echo 'bc:' \$(bc --version | head -n1 | cut -d " " -f2) > versions.txt
    echo 'featureCounts:' \$(featureCounts 2>&1 | head -n2 | tail -n1 | cut -d " " -f2) >> versions.txt
    """
    
}

// This process collects FRiPs per sample and merges them into a single file, sorted by decreasing FRiP

process FripsCollect {

    cpus 1
    memory 100.MB
    
    publishDir = [
        path: params.outdir,
        mode: params.publishmode,
        saveAs: { filename -> 
                    (filename.equals("versions.txt") || 
                     filename.equals("command_lines.txt")) ? null : filename } 
                     
    ]

    container params.container

    input:
    path(frips)

    output:
    path("frips_all.txt")
    
    script:
    """
    cat $frips | sort -k2,2 -r > frips_all.txt
    """

}  
