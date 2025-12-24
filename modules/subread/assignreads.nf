
/*
 * ASSIGN_READS
 *
 * Assign aligned reads to genes using featureCounts. Bam is sorted immediately.
 *
 */
process ASSIGN_READS {
    tag "$meta.id"
    label 'process_med'
    publishDir "${params.outdir}/subread", mode: 'copy'
    
    module 'subread/2.1.1'
    module 'samtools/1.22.1'

    input:
    tuple val(meta), path(bamfiles)

    output:
    tuple val(meta), path("*bam"), emit: bam, optional: true
    tuple val(meta), path("*_counts.txt"), emit: counts
    tuple val(meta), path("*.summary"), emit: log
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cores = 1
    if (task.cpus) {
        cores = (task.cpus as int) - 1
        if (cores < 1) cores = 1
        if (cores > 8) cores = 8
    }
    def bamFiles = bamfiles.join(' ')
    
    """
    valid_bams=()
    for bam in $bamFiles; do
        if [[ -s "\$bam" ]] && [[ \$(samtools view -c "\$bam") -gt 0 ]]; then
            valid_bams+=("\$bam")
        fi
    done
    
    if [[ \${#valid_bams[@]} -eq 0 ]]; then
        echo "No valid BAMs for ${prefix}, generating empty counts"
        awk 'BEGIN{OFS="\t"; print "Geneid","Length","nonumi_count","umi_count"}' > ${prefix}_counts.txt
        touch ${prefix}_counts.summary
    else
      featureCounts -p \
        -a ${params.gtf_path} \
        -o ${prefix}_counts.txt \
        -T $cores \
        -g gene_name \
        -R BAM \
        "\${valid_bams[@]}"
  
      awk 'BEGIN{OFS="\t"}
      NR==2{
          genecol=1
          lengthcol=6
          nonumi=0
          umi=0
          for(i=7;i<=NF;i++){
              if(\$i ~ /nonumiAligned/){ nonumi=i }
              if(\$i ~ /umiAligned/){ umi=i }
          }
          print "Geneid","Length","nonumi_count","umi_count"
          next
      }
      NR>2{
          nonumi_val = (nonumi>0 ? \$nonumi : 0)
          umi_val    = (umi>0 ? \$umi : 0)
          print \$genecol,\$lengthcol,nonumi_val,umi_val
      }' ${prefix}_counts.txt > ${prefix}_counts.tmp \
      && mv ${prefix}_counts.tmp ${prefix}_counts.txt
        
      mv ${prefix}_counts.txt.summary ${prefix}_counts.summary
    
      for bam in *.featureCounts.bam; do
        if [[ "\$bam" == *"_umiAligned"* ]]; then
          out="${prefix}_umi.featureCounts.bam"
        elif [[ "\$bam" == *"_nonumiAligned"* ]]; then
          out="${prefix}_nonumi.featureCounts.bam"
        else
          continue
        fi
        mv "\$bam" "\$out"
      done
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        featureCounts: \$(echo \$(featureCounts -v 2>&1) | sed 's/^.*version //; s/Last.*\$//') 
        samtools: \$(samtools --version | head -n 1 | awk '{print \$2}')
    END_VERSIONS
    """
    
}


