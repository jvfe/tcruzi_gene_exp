#!/usr/bin/env Rscript

library(dplyr)
library(purrr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: get_counts_hisat2.R <feature_counts_files>", call. = FALSE)
}

files <- unlist(strsplit(args, split = " "))

map(files, ~ {
  feature_counts <- read.table(.x, header = TRUE) %>%
    select(1, 7)
}) %>%
  reduce(inner_join, by = "Geneid") -> count_table

write.table(count_table, file = "count_table_hisat2.txt", row.names = FALSE, quote = FALSE)

colnames(count_table) <- sapply(strsplit(colnames(count_table), split = "_"), "[[", 1)
rownames(count_table) <- count_table$Geneid
count_table$Geneid <- NULL

count_table <- as.matrix(count_table)

saveRDS(count_table, file = "count_table_hisat2.rds")
