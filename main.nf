#!/usr/bin/env nextflow

// https://github.com/bioinfo-pf-curie/scRNA-SmartSeq3/tree/master
// https://www.protocols.io/view/smart-seq3xpress-yxmvmk1yng3p/v3?step=18


nextflow.enable.dsl=2


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { EXTRACT_UMI } from './modules/umi_tools/extract'
include { TRIM_POLYTAIL } from './modules/cutadapt/trimpolytail'
include { ALIGN_READS } from './modules/star/alignreads'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
  Channel
    .fromPath(params.input)
    .splitCsv(header: true)
    .map {
      row ->
        def fastq = row.findAll { it.key != 'sample' }
        return [[id:row.sample], [fastq.values()]]
    }
    .groupTuple(by: [0])
    .map {meta, fastq -> [meta, fastq.flatten()] }
    .set { ch_fastq }
    
  EXTRACT_UMI(ch_fastq)
  TRIM_POLYTAIL(EXTRACT_UMI.out.reads)
  ALIGN_READS(TRIM_POLYTAIL.out.reads)
  
}

