repository <- "https://stat.ethz.ch/CRAN/"

if(!"optparse" %in% installed.packages())
  install.packages("optparse", repos = repository)
library(optparse)

if(!"ape" %in% installed.packages())
  install.packages("ape", repos = repository)
library(ape)

if(!"phangorn" %in% installed.packages())
  install.packages("phangorn", repos = repository)
library(phangorn)

if(!"tools" %in% installed.packages())
  install.packages("tools", repos = repository)
library(tools)

# function: compute tree distances
compute.tree.distance <- function(row.data) {
  real.tree <- read.tree(file = row.data["real_tree"])

  if (grepl("sieve", row.data["tool"], ignore.case = TRUE, fixed = TRUE) || grepl("sciphi", row.data["tool"], ignore.case = TRUE, fixed = TRUE)) {
    inferred.tree <- read.nexus(file = row.data["inferred_tree"], force.multi = FALSE)
  } else {
    inferred.tree <- read.tree(file = row.data["inferred_tree"])
  }

  results <- numeric()

  # trees from Cellphy are not rooted
  is.rooted <- !grepl("cellphy", row.data["tool"], ignore.case = TRUE, fixed = TRUE)

  #print(paste("Working on:", row.data["real_tree"], "and", row.data["inferred_tree"], sep = " "))
  rf <- RF.dist(real.tree, inferred.tree, rooted = is.rooted, check.labels = TRUE)
  nrf <- RF.dist(real.tree, inferred.tree, rooted = is.rooted, normalize = TRUE, check.labels = TRUE)
  results <- c(results, rf / nrf, rf, nrf)

  if (grepl("sciphi", row.data["tool"], ignore.case = TRUE, fixed = TRUE)) {
    results <- c(results, NA, NA, NA)
  } else {
    wrf <- wRF.dist(real.tree, inferred.tree, rooted = is.rooted, check.labels = TRUE)
    nwrf <- wRF.dist(real.tree, inferred.tree, rooted = is.rooted, normalize = TRUE, check.labels = TRUE)
    rbc <- KF.dist(real.tree, inferred.tree, rooted = is.rooted, check.labels = TRUE)
    results <- c(results, wrf, nwrf, rbc)
  }

  names(results) <- c("max_clades", "RF_distance", "normalized_RF_distance", "weighted_RF_distance", "normalized_weighted_RF_distance", "rooted_branch_score_difference")

  return(results)
}

# function: update tree table with tree distances
get.tree.info <- function(input.file) {
  tree.info <- read.table(input.file, header = TRUE, sep = "\t")
  tree.info <- cbind(tree.info, t(apply(tree.info, 1, compute.tree.distance)))
  return(tree.info)
}

# global variables
tree.info <- data.frame()

# command line arguments
option.list <- list(
    make_option(c("-i", "--input"), type = "character", default = NULL, metavar = "INPUT TREE INFORMATION",
                help = "input tree information file in tsv format"),
    make_option(c("-o", "--output"), type = "character", default = NULL, metavar = "OUTPUT TREE INFORMATION",
                help = "output tree distance file in tsv format")
)
opt.parser <- OptionParser(option_list = option.list)
opt <- parse_args(opt.parser)

# command line arguments sanity check and variable initialization
if (is.null(opt$input) || is.null(opt$output)) {
  print_help(opt.parser)
  stop("Both input and output file should be specified.")
}

tree.info <- get.tree.info(opt$input)
write.table(tree.info, file = opt$output, sep = "\t", row.names = FALSE)
