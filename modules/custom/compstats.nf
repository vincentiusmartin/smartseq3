

process SUMMARIZE_RUNSTATS{
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
