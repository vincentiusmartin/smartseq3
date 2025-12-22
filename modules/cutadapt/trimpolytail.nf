

process TRIM_POLYTAIL {
  tag "$meta.id"
  label 'process_med'
  publishDir "${params.outdir}/cutadapt", mode: 'copy'
  
  module 'cutadapt/5.2'
  
  input:
  tuple val(meta), path(reads)
  
  output:
  tuple val(meta), path("*_trimmed*.fastq.gz"), emit: reads
  tuple val(meta), path("versions.yml"), emit: log
  
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  for type in umi nonumi; do
      cutadapt -g T{20} -G T{20} -a A{20} -A A{20} -a G{20} -A G{20} \
          --minimum-length 20 -j 8 \
          -o ${prefix}_\${type}_trimmed.1.fastq.gz \
          -p ${prefix}_\${type}_trimmed.2.fastq.gz \
          ${prefix}_\${type}.1.fastq.gz ${prefix}_\${type}.2.fastq.gz
  done
  
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      cutadapt: \$(echo \$(cutadapt --version 2>&1) | sed 's/^.*version //; s/Last.*\$//')
  END_VERSIONS
  """
}
