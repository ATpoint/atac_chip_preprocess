#! /usr/bin/env nextflow

nextflow.enable.dsl=2

//------------------------------------------------------------------------
// Intro message
//------------------------------------------------------------------------

def longline="========================================================================================================================="

ANSI_RESET   = "\u001B[0m"
ANSI_RED     = "\u001B[31m"
ANSI_GREEN   = "\u001B[32m"
ANSI_YELLOW  = "\u001B[33m"
DASHEDDOUBLE = "=".multiply(121)
DASHEDSINGLE = "-".multiply(121)
    
// Function for printing consistent error messages:
include{ ErrorMessenger } from './functions/validate_schema_params.nf'

// Date for intro message
Date date = new Date()
String datePart = date.format("yyyy-dd-MM -- ")
String timePart = date.format("HH:mm:ss")
def start_date = datePart + timePart

println ""
println "\u001B[33m$longline"
println ""
println "Pipeline:          rnaseq_preprocess"
println "GitHub:            https://github.com/ATpoint/rnaseq_preprocess"
println "Revision:          $workflow.revision"
println "Commit:            $workflow.commitId"
println "Documentation:     https://github.com/ATpoint/rnaseq_preprocess/README.md"
println "Author:            Alexander Toenges (@ATpoint)"
println "Runname:           $workflow.runName"
println "Profile:           $workflow.profile"
if(workflow.containerEngine!=null) println "Container Engine:  $workflow.containerEngine"
println "Command line:      $workflow.commandLine"
println "Start:             $start_date"
println ""
println "$longline\u001B[0m"

//------------------------------------------------------------------------
// Validate input params via schema.nf
//------------------------------------------------------------------------

evaluate(new File("${baseDir}/functions/validate_schema_params.nf"))

//------------------------------------------------------------------------
// Load the modules and pass params
//------------------------------------------------------------------------

include{ ValidateSamplesheet } from './modules/validatesamplesheet'

//------------------------------------------------------------------------

include{ CatFastq } from './modules/cat_fastq' addParams(outdir: params.merge_dir, keep: params.keep_merge)

//------------------------------------------------------------------------

include{ Fastqc } from './modules/fastqc' addParams(outdir: params.fastqc_dir)

//------------------------------------------------------------------------

include{ Trim } from './modules/trim' addParams(outdir: params.trim_dir, args: params.trim_additional, keep: params.keep_merge)

//------------------------------------------------------------------------

include{ Align } from './modules/align' addParams(outdir: params.align_dir, args: params.align_additional, args2: params.sort_additional)

//------------------------------------------------------------------------

include{ Chromsizes } from './modules/chromsizes' addParams(outdir: params.misc_dir)

//------------------------------------------------------------------------

include{ Filter } from './modules/filter' addParams(outdir: params.filter_dir, args: params.filter_additional)

//------------------------------------------------------------------------

include{ Isizes } from './modules/isizes' addParams(outdir: params.misc_dir)

//------------------------------------------------------------------------

include{ Cutsites } from './modules/cutsites' addParams(outdir: params.bed_dir)

//------------------------------------------------------------------------

// If ATAC-seq => set default, if not ATAC-seq => set different defaults, or overwrite completely if the user uses --macs_additional
args_peaks = params.macs_additional=='' ? (params.atacseq ? '--keep-dup=all --nomodel --extsize 100 --shift -50 --min-length 250' : '--keep-dup=all') : params.macs_additional

include{ Peaks } from './modules/peaks' addParams(outdir: params.peaks_dir, args: args_peaks)

//------------------------------------------------------------------------

include{ Frips } from './modules/frips' addParams(outdir: params.frip_dir)

//------------------------------------------------------------------------

include{ FripsCollect } from './modules/frips' addParams(outdir: params.frip_dir)

//------------------------------------------------------------------------

include{ Bigwigs } from './modules/bigwigs' addParams(outdir: params.bigwig_dir)

//------------------------------------------------------------------------

include{ Multiqc } from './modules/multiqc' addParams(outdir: params.multiqc_dir)

//------------------------------------------------------------------------

include{ CommandLines } from './modules/commandline' addParams(outdir: params.pipeline_dir)                                                              

//------------------------------------------------------------------------
// Define workflows
//------------------------------------------------------------------------

workflow VALIDATESSAMPLESHEET {

    take: 
        samplesheet_unvalidated

    main:
        ValidateSamplesheet(samplesheet_unvalidated)

    samplesheet = ValidateSamplesheet.out.samplesheet
           .splitCsv(header:true)
           .map {
               
                // Samplesheet allows Nextflow variables to be used, replace by absolute path
                r1 = it['r1']
                        .replaceAll('\\$baseDir|\\$\\{baseDir\\}', new String("${baseDir}/"))
                        .replaceAll('\\$launchDir|\\$\\{launchDir\\}', new String("${launchDir}/"))
                        .replaceAll('\\$projectDir|\\$\\{projectDir\\}', new String("${projectDir}/"))

                rx = it['r2']
                        .replaceAll('\\$baseDir|\\$\\{baseDir\\}', new String("${baseDir}/"))
                        .replaceAll('\\$launchDir|\\$\\{launchDir\\}', new String("${launchDir}/"))
                        .replaceAll('\\$projectDir|\\$\\{projectDir\\}', new String("${projectDir}/"))

                r2 = rx.toString()=='' ? '.' : rx   

                // Make sure fastq file paths exist
                is_error = false

                r1_file = file(r1).exists()           
                if(!r1_file){
                    ErrorMessenger("$r1 does not exist")
                    is_error = true
                }

                if(r2!=".") { 
                    r2_file = file(r2).exists()
                    if(!r1_file){
                        ErrorMessenger("$r2 does not exist") 
                        is_error = true
                    }
                }

                se = rx.toString()=='' ? true : false

                // meta map inspired by nf-core
                meta = [id:it['sample'], group:it['group'], single_end: se]     
                reads = [r1: r1, r2: r2]      
                counter = [1]

                if(!is_error){
                    return tuple(meta, reads, counter)
                } else {
                    return null
                }
                
            }
            .groupTuple(by:0)
            .map{ meta, grouped_reads, counter -> [ meta, grouped_reads.flatten(), counter ] }
            .map{ meta, grouped_reads, counter -> 

                    if(meta.single_end){
                        reads = grouped_reads['r1'].flatten()
                    } else {
                        reads = [grouped_reads['r1'].flatten(), grouped_reads['r2'].flatten()].flatten()
                    }

                    // counter is > 1 if a sample has more than one fastq file per read, requiring a merge
                    // before actual fastq processing
                    [meta, reads, counter.size()]

            }

    emit:
        samplesheet = samplesheet
        versions = ValidateSamplesheet.out.versions

}

// THE MAIN WORKFLOW

workflow {

   // ----------------------------------------------------------------------------------------
   // Validate the provided samplesheet and merge fastq if necessary
   // ----------------------------------------------------------------------------------------

    sx = file(params.samplesheet, checkIfExists: true)
        
    VALIDATESSAMPLESHEET(params.samplesheet)
    validatesamplesheet_versions = VALIDATESSAMPLESHEET.out.versions

    // Samples with > 1 fastq per read
    VALIDATESSAMPLESHEET.out.samplesheet
    .map {meta, reads, counter -> 
        
        if(counter>1) [meta, reads]
 
    }.set { ch_needMerge }

    CatFastq(ch_needMerge)
    ch_merged = CatFastq.out.fastq_tuple
    cat_versions = CatFastq.out.versions  

    // Samples with 1 fastq per read
    VALIDATESSAMPLESHEET.out.samplesheet
       .map {meta, reads, counter -> 
            
            if(counter==1) [meta, reads]

        }.set { ch_noMerge }

    // This channel is now [meta, reads] and can go in all downstream processes that require fastq
    ch_fastq = ch_noMerge.concat(ch_merged)

    // ----------------------------------------------------------------------------------------
    // FastQC
    // ----------------------------------------------------------------------------------------

    Fastqc(ch_fastq)
    fastqc_versions = Fastqc.out.versions
    fastqc_for_multiqc = Fastqc.out.zip

    // ----------------------------------------------------------------------------------------
    // Trim (or not)
    // ----------------------------------------------------------------------------------------

    if(!params.do_not_trim){

        Trim(ch_fastq)
        trim_versions = Trim.out.versions
        ch_for_align = Trim.out.tuple_trimmed
        trim_for_multiqc = Trim.out.json
    
    } else {

        trim_versions = Channel.empty()
        ch_for_align = ch_fastq
        trim_for_multiqc = Channel.empty()

    }
    
    // ----------------------------------------------------------------------------------------
    // Align, mark dulicates, sort
    // ----------------------------------------------------------------------------------------

    Align(ch_for_align, params.index)
    align_versions = Align.out.versions
    align_for_multiqc = Align.out.tuple_flagstat.map {it[1]}

    // ----------------------------------------------------------------------------------------
    // Chromsizes
    // ----------------------------------------------------------------------------------------

    // This sort makes the channel order deterministic as channel order itself is unsorted/random
    // and hence a simple first() would result in different files to fetch chromsizes from so
    // a -resume would still cause rerun of that process and as such all downstream processes
    // that depend on the chromsizes output
    ch_for_chromsizes = Align.out.tuple_bam
                        .map { [it[1].getName(), it[1]] }
                        .toSortedList( { a, b -> b[1] <=> a[1] } )
                        .map{it[0]}
                        .map{it[1]}
                        
    Chromsizes(ch_for_chromsizes)
    chromsizes_versions = Chromsizes.out.versions
    
    // ----------------------------------------------------------------------------------------
    // Filtering of alignments
    // ----------------------------------------------------------------------------------------

    // Get the default filtering flags if the user does not provide explicit ones
    // remove unmapped, not primary, duplicate or supplementary, default is 3332
    Filter(Align.out.tuple_bam, params.flag_remove, params.chr_regex)
    filter_versions = Filter.out.versions
    filter_for_multiqc = Filter.out.tuple_flagstat.map {it[1]}

    // ----------------------------------------------------------------------------------------
    // Extraction of cutsites for ATAC-seq data
    // ----------------------------------------------------------------------------------------

    // Get 5' end of alignments, shifted by -4/+5 to get transposome insertion position
    if(params.atacseq){

        Cutsites(Filter.out.tuple_bam)
        cutsites_versions = Cutsites.out.versions

        input_peak_calling = Cutsites.out.tuple_cutsites

    } else {

        cutsites_versions = Channel.empty()
        input_peak_calling = Filter.out.tuple_bam.map { [meta: it[0], reads: it[1]] }

    }

    // ----------------------------------------------------------------------------------------
    // Insert sizes for paired-end data
    // ----------------------------------------------------------------------------------------

    ch_for_isizes = Filter.out.tuple_bam.map { if(!it[0].single_end) return it }
    Isizes(ch_for_isizes)
    isizes_version = Isizes.out.versions
    isizes_for_multiqc = Isizes.out.isizes

    // ----------------------------------------------------------------------------------------
    // Peak calling and filtering against blacklist
    // ----------------------------------------------------------------------------------------
    
    // Decide for blacklist based on species string
    if(params.filter_blacklist){

        if(params.atacseq){
            if(params.species=="hs") {
                blacklist = Channel.fromPath("$baseDir/assets/hg38/hg38_combined_blacklist.bed")
            }
            if(params.species=="mm") {
                blacklist = Channel.fromPath("$baseDir/assets/mm10/mm10_combined_blacklist.bed")
            }
        }

        if(!params.atacseq){
            if(params.species=="hs") {
                blacklist = Channel.fromPath("$baseDir/assets/hg38/hg38_encode_blacklist_v2.bed")
            }
            if(params.species=="mm") {
                blacklist = Channel.fromPath("$baseDir/assets/mm10/mm10_encode_blacklist_v2.bed")
            }
        }

    } else {
        blacklist = "."
    }

    Peaks(input_peak_calling, blacklist.collect())
    peaks_versions = Peaks.out.versions
    //peaks_for_multiqc = Peaks.out.xls

    // ----------------------------------------------------------------------------------------
    // FRiP calculation
    // ----------------------------------------------------------------------------------------

    // Make sure that only samples that produced any peaks are submitted to this process
    tmp_bam   = Filter.out.tuple_bam.map { [it[0], it[1]] }
    
    tmp_peaks = Peaks.out.tuple_peaks.map { meta, peaks, summits -> 
                    
                    lns = peaks.countLines()
                    if(lns>0){
                        return [meta, peaks]
                    } else return null
                }

    ch_for_frip = tmp_bam.concat(tmp_peaks)
                  .groupTuple(by:0, remainder: false)
                  .map {if(it[1][1]!=null) return it}

    // This writes sample names with no peaks to a file in params.peaks_dir
    no_peaks_file = file("${params.peaks_dir}/no_peaks.txt")
    if(no_peaks_file.exists()) no_peaks_file.delete()
    no_peaks = tmp_bam.concat(tmp_peaks)
               .groupTuple(by:0, remainder: false)
               .map { if(it[1][1]==null) no_peaks_file.append("${it[0].id}\n") }                    
    
    // Samples that produced peaks now go into the FRiP calculation process
    Frips(ch_for_frip)

    FripsCollect(Frips.out.frips.map {it[1]}.collect())
    
    frips_versions = Frips.out.versions

    // ----------------------------------------------------------------------------------------
    // Bigwig creation
    // ----------------------------------------------------------------------------------------

    // Create browser tracks for visual inspection of the data -- without normalization
    // as we later (manually as it is hard to automate that) scale the tracks using the
    // scale factors from the differential analysis, not per-million,
    // here is why: https://www.biostars.org/p/413626/#414440

    ch_for_bigwigs = params.atacseq ? Cutsites.out.tuple_cutsites : Filter.out.tuple_bam.map { 
                                                                        [it[0], it[1]]
                                                                    }
    Bigwigs(ch_for_bigwigs, Chromsizes.out.chromsizes)
    bigwigs_versions = Bigwigs.out.versions

    // ----------------------------------------------------------------------------------------
    // MultiQC
    // ----------------------------------------------------------------------------------------

    Multiqc(fastqc_for_multiqc.concat(trim_for_multiqc, align_for_multiqc, filter_for_multiqc, 
            isizes_for_multiqc).collect()) 
    multiqc_versions = Multiqc.out.versions

    // ----------------------------------------------------------------------------------------
    // Command lines and software versions
    // ----------------------------------------------------------------------------------------
    
    x_commands = validatesamplesheet_versions.concat(cat_versions, fastqc_versions, trim_versions, align_versions, 
                                                     chromsizes_versions, bigwigs_versions,
                                                     filter_versions, cutsites_versions,
                                                     peaks_versions, frips_versions,
                                                     isizes_version, multiqc_versions)
                 .map {it [1]}.flatten().collect()

    x_versions = validatesamplesheet_versions
                 .concat(cat_versions, align_versions, fastqc_versions, trim_versions, chromsizes_versions, 
                         filter_versions, bigwigs_versions, cutsites_versions,
                         peaks_versions, frips_versions, isizes_version, multiqc_versions)
                 .map {it [0]}
                 .flatten()
                 .collect()

    CommandLines(x_commands, x_versions)

}

def od = params.outdir
workflow.onComplete {
    Date date2 = new Date()
    String datePart2 = date2.format("yyyy-dd-MM -- ")
    String timePart2 = date2.format("HH:mm:ss")
    def end_date = datePart2 + timePart2
    println ""
    println "\u001B[33m$longline"
    println "Pipeline completed!"
    println "End: $end_date"
    println "Results are in:"
    println od
    println "$longline\u001B[0m"
    println ""
}
