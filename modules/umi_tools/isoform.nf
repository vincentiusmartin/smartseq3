
process RECONSTRUCT_ISOFORMS {
    tag "$meta.id"
    label 'process_med'
    publishDir "${params.outdir}/umi_tools/isoform", mode: 'copy'

    conda '/research/groups/northcgrp/home/common/Vincentius/envs/umitools'
    module 'samtools/1.22.1'

    input:
    tuple val(meta), path(bamfiles)

    output:
    tuple val(meta), path("*_isoforms.bam"), emit: isoforms_bam
    tuple val(meta), path("*_isoforms.tsv"), emit: isoforms_table

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bamFiles = bamfiles.join(' ')

    """
    out="${prefix}_isoforms.bam"
    samtools merge -@ ${task.cpus} merged.bam ${bamFiles}
    samtools sort -@ ${task.cpus} -o sorted.bam merged.bam
    samtools index sorted.bam
    
    # Group by UMI
    umi_tools group \
        --stdin=sorted.bam \
        --stdout=\$out \
        --output-bam \
        --extract-umi-method=tag \
        --umi-tag=RX
    
    # Sort and index grouped BAM
    samtools sort -@ ${task.cpus} -o sorted.bam \$out
    samtools index sorted.bam
    
    # Generate per-isoform counts
    umi_tools count \
        --per-gene \
        --gene-tag=GN \
        --stdin=sorted.bam \
        --stdout=${prefix}_isoforms.tsv
    """
}

