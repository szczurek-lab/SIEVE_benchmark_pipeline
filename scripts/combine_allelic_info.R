library(dplyr)

results <- lapply(snakemake@input, readRDS) %>% bind_rows()

saveRDS(results, snakemake@output[[1]])