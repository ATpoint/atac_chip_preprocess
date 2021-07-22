// Single-sample peak calling to calculate FRiPs:

process MacsPeaks {

    tag "$sample_id"

    cpus 1
    memory params.macs_mem

    publishDir params.macs_dir, mode: params.macs_pubmode

    input:
    tuple val(sample_id), path(infile)
    path(optional_control)

    output:
    path("*.{narrowPeak,broadPeak}"), emit: peaks
    path("*.bed"),  optional: true
    path("*.r"),    optional: true
    path("*.xls"),  optional: true
    path("*.bdg"),  optional: true
    path("*.saf"),  emit: saf

    script:

    // tweak to make optional inputs work: https://github.com/nextflow-io/nextflow/issues/1694#issuecomment-789334369
    if(optional_control instanceof List && optional_control.isEmpty()){
        ctrl = optional_control
    } else {
        if(optional_control == infile) { ctrl = '' } else ctrl = "-c $optional_control"
    }

    """

    #/ call peaks:
    macs2 callpeak $params.macs_additional $params.macs_format $params.macs_gflag \
        -t $infile $ctrl -n ${sample_id}${params.macs_suffix}

    #/ write as SAF:
    $baseDir/bin/write_saf.sh ${sample_id}${params.macs_suffix}_peaks.{narrowPeak,broadPeak} > ${sample_id}${params.macs_suffix}_peaks.saf

    """

}

