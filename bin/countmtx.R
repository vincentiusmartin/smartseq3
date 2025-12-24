#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

args <- commandArgs(trailingOnly=TRUE)
allpaths <- Sys.glob(args)

countmat <- lapply(allpaths, function(path){
  df <- fread(path) %>%
    mutate(cellid = sub("_counts.*$", "", basename(path)),
           count = nonumi_count + umi_count)
}) %>% bind_rows() %>%
  dplyr::select(Geneid, cellid, count) %>%   
  tidyr::pivot_wider(
    names_from = cellid,
    values_from = count,
    values_fill = 0 
  ) %>%
  tibble::column_to_rownames("Geneid") %>%
  as.matrix()

countmat <- Matrix::Matrix(countmat, sparse = TRUE)
tmp <- Matrix::writeMM(countmat, "matrix.mtx")
write.table(colnames(countmat), "barcodes.tsv", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(rownames(countmat), "features.tsv", quote = FALSE, row.names = FALSE, col.names = FALSE)
