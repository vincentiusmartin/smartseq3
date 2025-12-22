
process EXTRACT_UMI {
    tag "$meta.id"
    label 'process_med'
    publishDir "${params.outdir}/umi_tools", mode: 'copy'
    
    module 'seqkit/2.10.1'
    conda '/research/groups/northcgrp/home/common/Vincentius/envs/umitools'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*umi*.fastq.gz"), emit: reads
    tuple val(meta), path("*_umistats.csv"), emit: umistats
    tuple val(meta), path("*.{log,yml}"), emit: log

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    def args1_list = [
        "--extract-method=regex", 
        "--bc-pattern='(?P<discard_1>.*$params.sequence_tag)(?P<umi_1>.{$params.umi_size}(?P<discard_2>.{3}))'",
        "--stdin=${reads[0]}",
        "--stdout=${prefix}_umi.1.fastq.gz",
        "--read2-in=${reads[1]}",
        "--read2-out=${prefix}_umi.2.fastq.gz",
        "--filtered-out=${prefix}_nonumi.1.fastq.gz",
        "--filtered-out2=${prefix}_nonumi.2.fastq.gz",
        "--log=${prefix}_umi.log"
    ]

    """
    [ ! -f ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
    [ ! -f ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz

    umi_tools extract ${args1_list.join(' ')}
    #seqkit seq -j ${task.cpus} -n -i ${prefix}_nonumi.1.fastq.gz -o ${prefix}_nonUmiIds.txt
    
    n_reads=\$(( \$(zcat ${prefix}_1.fastq.gz | wc -l) / 4 ))
    n_umis=\$((\$(zcat ${prefix}_umi.1.fastq.gz| wc -l) / 4 ))
    percent_umi=\$(echo "scale=3; \$n_umis*100/\$n_reads" | bc)
    echo "sampleid,n_reads,n_umis,percent_umi" > ${prefix}_umistats.csv
    echo "${prefix},\${n_reads},\${n_umis},\${percent_umi}" >> ${prefix}_umistats.csv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        umi_tools: \$(echo \$(umi_tools --version 2>&1) | sed 's/^.*version //; s/Last.*\$//')
    END_VERSIONS
    """
}
