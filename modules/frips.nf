// Use featureCounts to calculate FRiPs:

process FRiPs {

    tag "$sample_id"

    cpus    params.frips_threads
    memory  params.frips_mem

    publishDir params.frips_dir, mode: params.frips_pubmode

    input:
    tuple val(sample_id), path(bam)
    path(saf)

    output:
    path("*.counts"),       emit: counts
    path("*.summary"),      emit: summary
    path("*_frips.txt"),    emit: frips

    script:

    """

    featureCounts $params.frips_additional -a $saf -F SAF -T $task.cpus -o ${saf.simpleName}.counts $bam

    $baseDir/bin/calc_frips.sh ${saf.simpleName}.counts.summary | paste <(echo ${saf.simpleName}) <(cat /dev/stdin) > ${saf.simpleName}_frips.txt

    """

}