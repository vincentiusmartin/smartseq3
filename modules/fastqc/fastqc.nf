process FASTQC {
    tag "$meta.id"
    label 'process_med'
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    module "fastqc/0.11.8"

    input:
    tuple val(meta), path(reads)

    output:
    path "*_fastqc.{zip,html}"

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Make list of old name and new name pairs to use for renaming in the bash while loop
    def old_new_pairs = reads instanceof Path || reads.size() == 1 ? [[ reads, "${prefix}.${reads.extension}" ]] : reads.withIndex().collect { entry, index -> [ entry, "${prefix}_${index + 1}.${entry.extension}" ] }
    def rename_to = old_new_pairs*.join(' ').join(' ')
    def renamed_files = old_new_pairs.collect{ old_name, new_name -> new_name }.join(' ')
    """
    printf "%s %s\\n" $rename_to | while read old_name new_name; do
        [ -f "\${new_name}" ] || ln -s \$old_name \$new_name
    done

    fastqc $renamed_files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
    END_VERSIONS
    """
}

process MULTIQC {
    label 'process_med'
    publishDir "${params.outdir}/multiqc", mode: "copy"

    module "multiQC/1.15"

    input:
    path 'fastqc/*'

    output:
    path "multiqc*"
    path "versions.yml"

    script:
    """
    multiqc .
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version )
    END_VERSIONS
    """
}