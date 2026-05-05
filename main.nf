#!/usr/bin/env nextflow

// https://github.com/bioinfo-pf-curie/scRNA-SmartSeq3/tree/master
// https://www.protocols.io/view/smart-seq3xpress-yxmvmk1yng3p/v3?step=18


nextflow.enable.dsl=2


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CAT_FASTQ } from './modules/cat/fastq'
include { FASTQC } from './modules/fastqc/fastqc'
include { MULTIQC } from './modules/fastqc/fastqc'
include { EXTRACT_UMI } from './modules/umi_tools/extract'
include { RECONSTRUCT_ISOFORMS } from './modules/umi_tools/isoform'
include { TRIM_POLYTAIL } from './modules/cutadapt/trimpolytail'
include { ALIGN_READS } from './modules/star/alignreads'
include { ASSIGN_READS } from './modules/subread/assignreads'
include { CREATE_COUNTMTX } from './modules/custom/countmtx'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
  Channel
    .fromPath(params.input)
    .splitCsv(header: true)
    .map { row  -> [ [id: row.sample], [row.fastq1, row.fastq2 ] ] }
    .groupTuple(by: [0])
    .map {meta, fastq -> [meta, fastq.flatten()] }
    .branch { meta, fastq ->
      single: fastq.size() == 2
      multiple: fastq.size() > 2
    }
    .set { ch_fastq }
  
  CAT_FASTQ(ch_fastq.multiple)
    .mix(ch_fastq.single)
    .set { ch_catfastq }
    
  FASTQC(ch_catfastq) | collect |  MULTIQC

  ch_reads = params.umi ? EXTRACT_UMI(ch_catfastq).out.reads : ch_catfastq
  TRIM_POLYTAIL(ch_reads)
  ALIGN_READS(TRIM_POLYTAIL.out.reads)
  ASSIGN_READS(ALIGN_READS.out.bam).set { assign_reads_ch }
  
  if (params.umi) {
    RECONSTRUCT_ISOFORMS(assign_reads_ch.bam)
  }
  
  assign_reads_ch.counts
    .map { meta, counts -> counts }
    .collect()
    .set { ftcount_ch }

  CREATE_COUNTMTX(ftcount_ch)
}

