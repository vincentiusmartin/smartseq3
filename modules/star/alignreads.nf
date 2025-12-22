

process ALIGN_READS {
  tag "$meta.id"
  label 'process_med'
  publishDir "${params.outdir}/star", mode: 'copy'
  
  module 'star/2.7.11'
  
  input:
  tuple val(meta), path(reads)
  
  output:
  tuple val(meta), path("*.bam"), emit: reads
  tuple val(meta), path("*.tab"), emit: spjunc
  tuple val(meta), path("*.{out,yml}"), emit: log
  
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def cores = 1
  if (task.cpus) {
      cores = (task.cpus as int) - 1
      if (cores < 1) cores = 1
      if (cores > 8) cores = 8
  }
  
  """
  for type in umi nonumi; do
    STAR \
      --genomeDir $params.star_index_path \
      --sjdbGTFfile $params.gtf_path \
      --readFilesIn ${prefix}_\${type}_trimmed.1.fastq.gz ${prefix}_\${type}_trimmed.2.fastq.gz \
      --readFilesCommand zcat \
      --outFilterMultimapNmax 10 \
      --limitSjdbInsertNsj 2000000 \
      --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
      --clip3pAdapterSeq CTGTCTCTTATACACATCT CTGTCTCTTATACACATCT \
      --clip3pAdapterMMp 0.1 0.1 \
      --outSAMtype BAM SortedByCoordinate \
      --outFileNamePrefix ${prefix}_\${type} \
      --runThreadN $cores
  done 
  
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      star: \$(echo \$(umi_tools --version 2>&1) | sed 's/^.*version //; s/Last.*\$//')
  END_VERSIONS
  """
}