#! /usr/bin/env nextflow

nextflow.enable.dsl=2

println ''
println '|------------------------------------------------------------------'
println "[Info] This is atac_chip_preprocess version ::: $params.version"
println '|------------------------------------------------------------------'
println ''

// if no alignment then turn off everything else (that is only the case if only-indexing shall run)
if(params.skip_align) params.skip_bamfilter = true

// make sure the fastq channel is handled properly:
// Check for presence of the fastq files, throw error if empty or set to null if only indexing shall be executed:
if(!params.skip_align) {

    if(params.fastq.isEmpty()){
        println("[Error] --fastq is empty!")
        System.exit(1)
    }

    if(params.mode == "paired"){

        ch_fastq = Channel
                    .fromFilePairs(params.fastq, checkIfExists: true)

    } else if(params.mode == "single"){
        ch_fastq = Channel
                    .fromPath(params.fastq, checkIfExists: true)
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

    if(!params.skip_align) {

    // calculate total threads required (bowtie2+samtools-sort+cutadapt + cutadapt). Ignore samblaster as it runs on
    // just a few %CPU (<<< 100% per core) so it should not merit an extra core:
        align_threads_overall = (params.align_threads + params.sort_threads + 1)

        include{ Bowtie2Align } from './modules/align' addParams(   threads:    align_threads_overall,
                                                                    memory:     params.align_mem_total    )                                                                        

    
        Bowtie2Align(ch_fastq, use_index)
        use_bam_raw = Bowtie2Align.out.bam
    } else use_bam_raw = ''

    //-------------------------------------------------------------------------------------------------------------------------------//
    // Filtering of the aligned BAM file.

    if(!params.skip_bamfilter){

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
    if(params.mode == 'paired' && !params.skip_align && !params.skip_isizes){

        include { InsertSizes } from './modules/isizes'

        InsertSizes(use_filtered_bam)

    }

    //-------------------------------------------------------------------------------------------------------------------------------//
    // In ATAC-seq mode produce a bed.gz with the cutting sites:
    if(params.atacseq & !params.skip_align){

        // process needs three threads at least (one for awk, one for bedtools and one for sort)
        cutsite_threads = (2 + params.cutsites_threads)
        include{ Cutsites } from './modules/cutsites' addParams(threads: cutsite_threads)

        Cutsites(use_filtered_bam)

    }

    //-------------------------------------------------------------------------------------------------------------------------------//
    // Basic QC, for this call peaks and calculate FRiPs:
    if(!params.skip_align){

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

workflow { ATAC_CHIP() }