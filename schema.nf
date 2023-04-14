#! /usr/bin/env nextflow

// --------------------------------------------------------------------------------------------------------------------------------------------------------------------
//
// SCHEMA DEFINITION FOR PARAMS VALIDATION
//
// --------------------------------------------------------------------------------------------------------------------------------------------------------------------

def Map schema = [:] // don't change this line

// --------------------------------------------------------------------------------------------------------------------------------------------------------------------

// General options:
schema.title1         = [title: 'GENERAL OPTIONS']
schema.min_nf_version = [value: '23.04.0', type: 'string', mandatory: true, allowed: '']
schema.publishmode    = [value: 'copy', type: 'string', mandatory: true, allowed: ['symlink', 'rellink', 'link', 'copy', 'copyNoFollow', 'move']]
schema.container      = [value: "atpoint/atac_chip_preprocess:v1.2.0", type: 'string', mandatory: false]

// Output directories
schema.title1         = [title: 'OUTPUT DIRECTORIES']
schema.outdir         = [value: params.outdir, type: 'string', mandatory: true]
schema.pipeline_dir   = [value: "${params.outdir}/pipeline_info/", type: 'string']
schema.merge_dir      = [value: "${params.outdir}/fastq_merged/", type: 'string']
schema.fastqc_dir     = [value: "${params.outdir}/fastqc/", type: 'string']
schema.multiqc_dir    = [value: "${params.outdir}/multiqc/", type: 'string']
schema.trim_dir       = [value: "${params.outdir}/fastq_trimmed/", type: 'string']
schema.align_dir      = [value: "${params.outdir}/alignments_unfiltered/", type: 'string']
schema.filter_dir     = [value: "${params.outdir}/alignments_filtered/", type: 'string']
schema.bed_dir        = [value: "${params.outdir}/bed_files/", type: 'string']
schema.peaks_dir      = [value: "${params.outdir}/peaks/", type: 'string']
schema.frip_dir       = [value: "${params.outdir}/frips/", type: 'string']
schema.bigwig_dir     = [value: "${params.outdir}/bigwig/", type: 'string']
schema.misc_dir       = [value: "${params.outdir}/misc/", type: 'string']

// Input options
schema.title3         = [title: 'INPUT OPTIONS']
schema.samplesheet    = [value: '', type: 'string', mandatory: true]
schema.index          = [value: '', type: 'string', mandatory: true]
schema.atacseq        = [value: true, type: 'logical', mandatory: true, pattern: /.*\.csv$/]
schema.species        = [value: 'mm', type: 'string', mandatory: true, allowed: ['hs', 'mm']]
schema.blacklist      = [value: '', type: 'string']

// Merging fastq, trim align, duplicate marking, sort
schema.title4           = [title: 'MERGE/TRIM/ALIGN/MARK DUPLICATES/SORT OPTIONS']
schema.keep_merge       = [value: false, type: 'logical']
schema.keep_trim        = [value: false, type: 'logical']
schema.do_not_trim      = [value: false, type: 'logical', mandatory: true]
schema.trim_additional  = [value: '--dont_eval_duplication -z 6', type: 'string']
schema.align_additional = [value: '-X 2000 --very-sensitive', type: 'string']
schema.sort_additional  = [value: '-l 6', type: 'string']

// Filtering options
schema.title5            = [title: 'FILTERING OPTIONS']
schema.flag_remove       = [value: 3332, type: 'numeric']
schema.chr_regex         = [value: "chr[1-9,X,Y]", type: 'string']
schema.min_mapq          = [value: 20, type: 'numeric', mandatory: true]
schema.filter_additional = [value: '', type: 'string']

// Peak calling, FRiPs and Bigwigs
schema.title6           = [title: 'PEAK CALLING / QC OPTIONS']
schema.macs_additional  = [value: '', type: 'string']
schema.filter_blacklist = [value: true, type: 'logical']
schema.fragment_length  = [value: 250, type: 'numeric']

// --------------------------------------------------------------------------------------------------------------------------------------------------------------------

return schema // don't change this line