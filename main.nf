#! /usr/bin/env nextflow

nextflow.enable.dsl=2

ANSI_RESET   = "\u001B[0m"
ANSI_RED     = "\u001B[31m"
ANSI_GREEN   = "\u001B[32m"
ANSI_YELLOW  = "\u001B[33m"
DASHEDDOUBLE = "=".multiply(121)

// Check that NF version is 21.04.1 since we still have include statements in the workflow definitions
// which is an error in later Nextflow versions.
if( nextflow.version.matches(">21.04.3") ) {
    println "$ANSI_RED" + "$DASHEDDOUBLE"
    println "[VERSION ERROR] Please make sure you use a Nextflow version not greater 21.04.3"
    println "=> You are running version $nextflow.version."
    println "=> We recommend using NXF_VER=21.04.3 nextflow run (...)"
    println "$DASHEDDOUBLE ${ANSI_RESET}"
    System.exit(1)
}

//------------------------------------------------------------------------
// Intro message
//------------------------------------------------------------------------

Date date = new Date()
String datePart = date.format("yyyy-dd-MM -- ")
String timePart = date.format("HH:mm:ss")
def start_date = datePart + timePart

println ""
println "\u001B[33m$DASHEDDOUBLE"
println "Pipeline:      atac_chip_preprocess"
println "GitHub:        https://github.com/ATpoint/atac_chip_preprocess/"
println "Documentation: https://github.com/ATpoint/atac_chip_preprocess/README.md"
println "Author:        Alexander Toenges (@ATpoint)"
println "Runname:       $workflow.runName"
println "Profile:       $workflow.profile"
println "Start:         $start_date"
println "$DASHEDDOUBLE\u001B[0m"

// Print summary of all params:
def max_char = params.keySet().collect { it.length() }.max()
println ''
println "\u001B[33m$DASHEDDOUBLE"
params.each { nm, en -> 
        
        def use_length = max_char - nm.length()
        def spacer = ' '.multiply(use_length)
        println "${nm} ${spacer}:: ${en}" 

    }
println "$DASHEDDOUBLE\u001B[0m"
println ''

// if only_idx then imply skipping everything:
if(params.only_idx){
    skip_align       = true
    skip_bamfilter   = true
    skip_isizes      = true
    skip_qc          = true
    fastq            = null
} else {
    skip_align       = params.skip_align
    skip_bamfilter   = params.skip_bamfilter
    skip_isizes      = params.skip_isizes
    skip_qc          = params.skip_qc
    fastq            = params.fastq
}

// make sure the fastq channel is handled properly:
// Check for presence of the fastq files, throw error if empty or set to null if only indexing shall be executed:
if(!skip_align) {

    if(fastq.isEmpty()){
        println("[Error] --fastq is empty!")
        System.exit(1)
    }

    if(params.mode == "paired"){

        ch_fastq = Channel
                    .fromFilePairs(fastq, checkIfExists: true)

    } else if(params.mode == "single"){
        ch_fastq = Channel
                    .fromPath(fastq, checkIfExists: true)
                    .map { file -> tuple(file.simpleName, file) }
    }

} else ch_fastq = null

//-----------------------------------------------------------------------------------------------------------------------------------//

// Define the final workflow:
workflow ATAC_CHIP {

    //-------------------------------------------------------------------------------------------------------------------------------//
    // Indexing

    include{ Bowtie2Idx } from './modules/index'
    
    // Either make a new index from scratch or use provided one if exists:
    if(params.idx == ''){
        
        Bowtie2Idx(params.ref_genome)
        use_index = Bowtie2Idx.out.idx

    } else use_index = params.idx     

    //-------------------------------------------------------------------------------------------------------------------------------//
    // Trim/Alignment

    if(!skip_align) {

    // calculate total threads required (bowtie2+samtools-sort+cutadapt + cutadapt). Ignore samblaster as it runs on
    // just a few %CPU (<<< 100% per core) so it should not merit an extra core:
        align_threads_overall = (params.align_threads + params.sort_threads + 1)

        // calculate total memory required for align+sort:
        m1 = params.align_mem.replaceAll(".GB|GB|G", "").toInteger()
        m2 = params.sort_mem.replaceAll(".GB|GB|G", "").toInteger() * params.sort_threads
        total_mem = "${(m1+m2)}.GB"

        if(params.trim_adapter == ''){
            if(params.atacseq){
                use_adapter = 'CTGTCTCTTATACACATCT' // Nextera
            } else use_adapter = 'AGATCGGAAGAGC'    // TruSeq
        } else use_adapter = params.trim_adapter    // custom user-provided

        include{ Bowtie2Align } from './modules/align' addParams(   threads:        align_threads_overall,
                                                                    memory:         total_mem,
                                                                    trim_adapter:   use_adapter             )  

        Bowtie2Align(ch_fastq, use_index)
        use_bam_raw = Bowtie2Align.out.bam

    } else use_bam_raw = ''

    //-------------------------------------------------------------------------------------------------------------------------------//
    // Filtering of the aligned BAM file.

    if(!skip_bamfilter){

        // set defaults for -f and -F option of samtools view in case the params from the config file are empty
        f = params.bamfilter_flag_keep
        if(params.bamfilter_flag_keep == ''){
            if(params.mode == "paired") f = '-f 1'
            if(params.mode == "single") f = '-f 0'
        } else if(params.bamfilter_flag_keep == "nofilter") f = ''
            
        F = params.bamfilter_flag_remove     
        if(params.bamfilter_flag_remove == ''){
            F = '-F 3332'
        } else if(params.bamfilter_flag_remove == "nofilter") F = ''

        include{ FilterBam } from './modules/filter'    addParams(  flag_keep:    f,
                                                                    flag_remove:  F)

        FilterBam(use_bam_raw)
        use_filtered_bam = FilterBam.out.bam

    } else use_filtered_bam = ''

    //-------------------------------------------------------------------------------------------------------------------------------//
    // Insertsize metrics:
    if(params.mode == 'paired' && !skip_align && !skip_isizes){

        include { InsertSizes } from './modules/isizes'

        InsertSizes(use_filtered_bam)

    }

    //-------------------------------------------------------------------------------------------------------------------------------//
    // In ATAC-seq mode produce a bed.gz with the cutting sites:
    if(params.atacseq & !skip_bamfilter){

        // need at least 3 threads, (bedtools+awk+sort, and then later for bgzip which starts when sort finishes)
        if(params.cutsites_threads < 3){
            if(!params.testing) { cutsites_threads = 3 } else { cutsites_threads = 2 } // there are only two CPUs on the Linux VMs for CI testing
        } else { cutsites_threads = params.cutsites_threads }

        // the requested memory should be the sort memory (--cutsites_mem) with 10% extra to avoid
        // requested memory being near limit and oom kill on SLURM
        cutsites_mem = "${params.cutsites_mem.replaceAll(".GB|GB|G", "").toInteger() * 1.1}GB"

        include{ Cutsites } from './modules/cutsites' addParams(threads: cutsites_threads,
                                                                memory:  cutsites_mem)

        Cutsites(use_filtered_bam)

    }

    //-------------------------------------------------------------------------------------------------------------------------------//
    // Basic QC, for this call peaks and calculate FRiPs:
    if(!skip_bamfilter && !skip_qc){

        // experimental feature, README recommends not to use
        if(params.macs_control == '') { macs_ctrl = [] } else macs_ctrl = params.macs_control
        
        // preset for ATAC-seq:
        if(params.atacseq){

            macs_input       = Cutsites.out.bed
            macs_format      = '-f BED'

            if(params.macs_additional == ''){
                m_additional = '--nomodel --extsize 100 --shift -50 --keep-dup=all --min-length 150 -q 0.005'
            } else m_additional = params.macs_additional

        // and non ATAC-seq preset:
        } else {

            macs_input = use_filtered_bam.map { file -> tuple(file[0], file[1]) }
            if(params.mode == 'single') { macs_format = '-f BAM' } else { macs_format = '-f BAMPE' }

            if(params.macs_additional == ''){
                m_additional = '--keep-dup=all --min-length 150'
            } else m_additional = params.macs_additional

        }
        
        include { MacsPeaks } from './modules/macs_peaks'   addParams(  macs_additional:    m_additional,
                                                                        macs_format:        macs_format)

        MacsPeaks(macs_input, macs_ctrl)

        // FRiP calculation with featureCounts:
        if(params.atacseq) frips_additional = '--read2pos 5' 
        if(!params.atacseq && params.mode == 'paired') frips_additional = '-p'
        if(!params.atacseq && params.mode == 'single') frips_additional = ''

        include { FRiPs } from './modules/frips' addParams(frips_additional: frips_additional)

        FRiPs(use_filtered_bam.map { file -> tuple(file[0], file[1]) }, MacsPeaks.out.saf)

        }
    
}             

workflow { 
    
    ATAC_CHIP() 
    
    def od = params.outdir
    workflow.onComplete {
        Date date2 = new Date()
        String datePart2 = date2.format("yyyy-dd-MM -- ")
        String timePart2 = date2.format("HH:mm:ss")
        def end_date = datePart2 + timePart2
        println ""
        println "Pipeline completed!"
        println "End: $end_date"
        println "Results are in:"
        println od
        println ""
    }
    
}
