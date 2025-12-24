

// NOTE: INPUT/OUTPUT PATH SHOULD NEVER CONTAIN PATTERN, e.g. *counts.txt
// Otherwise, nextflow will rename files for 1.counts.txt, ...
process CREATE_COUNTMTX {
    label 'process_med'
    publishDir "${params.outdir}/counts", mode: 'copy'

    module 'R/4.2.2'

    input:
    path 'counts/*'

    output:
    tuple path("matrix.mtx"), path("barcodes.tsv"), path("features.tsv")

    script:
    """
    countmtx.R counts/*counts.txt
    """

}
