#!/usr/bin/env Rscript

library(dplyr)
library(purrr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: get_counts.r <feature_counts_files>", call. = FALSE)
}

files <- unlist(strsplit(args, split = " "))

map(files, ~ {

  sample <- gsub("\\.ReadsPerGene.out.tab", "", .x)
  feature_counts <- read.table(.x, header = TRUE, skip = 4,
                               col.names = c("gene", sample, "forward", "reverse")) %>%
    select(1, 2)
}) %>%
  reduce(inner_join, by = "gene") -> count_table

write.table(count_table, file = "count_table_star.txt", row.names = FALSE, quote = FALSE)

colnames(count_table) <- sapply(strsplit(colnames(count_table), split = "_"), "[[", 1)
rownames(count_table) <- count_table$gene
count_table$gene <- NULL

count_table <- as.matrix(count_table)

saveRDS(count_table, file = "count_table_star.rds")
