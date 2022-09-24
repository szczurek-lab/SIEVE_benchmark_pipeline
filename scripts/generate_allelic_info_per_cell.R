library(dplyr)

cell.names <- read.table(snakemake@input$cellNames)[[1]]

size.factors <- matrix(read.table(snakemake@input$trueSizeFactors)[[1]], nrow = 1)
colnames(size.factors) <- cell.names

allelic.seq.info.ori <- read.table(snakemake@input$trueAllelicSeqInfo)
colnames(allelic.seq.info.ori) <- c("t", "v")
allelic.seq.info <- lapply(allelic.seq.info.ori, function(x) matrix(x, ncol = 1))

allelic.info <- list(
    t = allelic.seq.info$t %*% size.factors,
    v = allelic.seq.info$v %*% (size.factors * size.factors)
)

file.list <- read.table(snakemake@input$sieveVarCallFiles)
colnames(file.list) <- c("file_name")
estimated.allelic.info.file <- filter(file.list, grepl(".+allelic_info", file_name))$file_name[1]
sieve.cell.names.file <- filter(file.list, grepl(".+cell_names", file_name))$file_name[1]

estimated.allelic.info <- read.table(paste(dirname(snakemake@input$sieveVarCallFiles), estimated.allelic.info.file, sep = "/"))
sieve.cell.names <- read.table(paste(dirname(snakemake@input$sieveVarCallFiles), sieve.cell.names.file, sep = "/"))
estimated <- tibble(
    cell_names = sieve.cell.names[[1]],
    t_estimated = estimated.allelic.info[[1]],
    v_estimated = estimated.allelic.info[[2]]
)

results <- tibble(
    cell_num = as.integer(snakemake@params$cellNum),
    coverage_mean = snakemake@params$covMean,
    coverage_variance = snakemake@params$covVariance,
    eff_seq_err_rate = snakemake@params$effSeqErrRate,
    ado_rate = snakemake@params$adoRate,
    deletion_rate = snakemake@params$deletionRate,
    insertion_rate = snakemake@params$insertionRate,
    tool = snakemake@params$tool,
    tool_setup = snakemake@params$toolSetup,
    fine_tune_type = snakemake@params$fineTuneType,
    data_type = snakemake@params$dataType,
    dataset = snakemake@params$datasetLabel,
    cell_names = cell.names,
    t_mean = apply(allelic.info$t, 2, mean),
    t_median = apply(allelic.info$t, 2, median),
    t_0.025 = apply(allelic.info$t, 2, function(x) quantile(x, 0.025)),
    t_0.975 = apply(allelic.info$t, 2, function(x) quantile(x, 0.975)),
    v_mean = apply(allelic.info$v, 2, mean),
    v_median = apply(allelic.info$v, 2, median),
    v_0.025 = apply(allelic.info$v, 2, function(x) quantile(x, 0.025)),
    v_0.975 = apply(allelic.info$v, 2, function(x) quantile(x, 0.975))
)

results <- results %>%
    left_join(estimated, by = c("cell_names")) %>%
    mutate(t_diff_wrt_median = t_estimated - t_median) %>%
    mutate(t_rela_diff_wrt_median = t_diff_wrt_median / t_median) %>%
    mutate(t_in_0.95 = (t_0.025 <= t_estimated & t_estimated <= t_0.975)) %>%
    mutate(v_diff_wrt_median = v_estimated - v_median) %>%
    mutate(v_rela_diff_wrt_median = v_diff_wrt_median / v_median) %>%
    mutate(v_in_0.95 = (v_0.025 <= v_estimated & v_estimated <= v_0.975))

saveRDS(results, file = snakemake@output[[1]])