suppressPackageStartupMessages(library(dplyr))

files <- data.frame(
    dataset = sort(as.character(snakemake@params$datasets)),
    trueSNVSites = sort(snakemake@input$trueSNVSites),
    parallelMut = sort(snakemake@input$parallelMuts)
) %>%
split(., .$dataset)

info <- bind_rows(
    sapply(
        names(files),
        function(x) {
            .trueSNVSites <- read.table(as.character(files[[x]][[2]]))

            .parallelMut <- tryCatch(
                read.table(as.character(files[[x]][[3]])),
                error = function(cond) return(NULL)
            )

            data.frame(
                cell_num = snakemake@params$cellNum,
                coverage_mean = snakemake@params$covMean,
                coverage_variance = snakemake@params$covVariance,
                eff_seq_err_rate = snakemake@params$effSeqErrRate,
                ado_rate = snakemake@params$adoRate,
                deletion_rate = snakemake@params$deletionRate,
                insertion_rate = snakemake@params$insertionRate,
                dataset = as.integer(x),
                total_snv_sites = dim(.trueSNVSites)[1L],
                total_snv_entries = dim(.trueSNVSites)[1L] * snakemake@params$cellNum,
                parallel_mu_sites = ifelse(
                    is.null(.parallelMut),
                    0L,
                    dim(.parallelMut)[1L]
                ),
                parallel_mu_entries = ifelse(
                    is.null(.parallelMut),
                    0L,
                    sum(rowSums(select(.parallelMut, -V1)))
                )
            )
        },
        simplify = FALSE,
        USE.NAMES = FALSE
    )
) %>%
mutate(
    prop_parallel_mu_sites = parallel_mu_sites / total_snv_sites,
    prop_parallel_mu_entries = parallel_mu_entries / total_snv_entries
)

write.table(
    x = info,
    file = snakemake@output[[1]],
    row.names = FALSE
)