# Title     : Utils
# Objective : Provide utensil functions for R
# Created by: senbaikang
# Created on: 07.03.21

if(!"scales" %in% installed.packages()){
  install.packages("scales", repos = repository)
}
library(scales)

source(file = "scripts/magic_color.R")

load.data <- function(
  file,
  sep = "\t",
  data.label.name = NA,
  data.label = NA,
  mu.rate.name = NA,
  mu.rate = NA
) {
  data <- read.table(file = file, header = TRUE, sep = sep)

  if (!is.na(data.label.name) && !is.na(data.label))
    data[data.label.name] <- data.label

  if (!is.na(mu.rate.name) && !is.na(mu.rate))
    data[mu.rate.name] <- mu.rate

  return(data)
}

load.rds <- function(
  file,
  data.label.name = NA,
  data.label = NA,
  mu.rate.name = NA,
  mu.rate = NA
) {
  data <- readRDS(file = file)

  if (!is.na(data.label.name) && !is.na(data.label))
    data[[data.label.name]] <- data.label

  if (!is.na(mu.rate.name) && !is.na(mu.rate))
    data[[mu.rate.name]] <- mu.rate

  return(data)
}

scientific <- function(vals){
  sapply(
    vals,
    function(x) {
      if (x == 0)
        return("0")

      x <- ifelse(
        grepl("[eE]", x),
        x,
        label_scientific()(x)
      )

      if (grepl("^1[eE]", x)) {
        x <- gsub("^1[eE]", "10^", x)
      } else {
        x <- gsub("[eE]", "%*%10^", x)
      }

      parse(text = gsub(
        "[+]",
        "",
        x
      ))
    }
  )
}

define.legend <- function(
  use.common.legend,
  common.legend,
  original.legend
) {
  if (use.common.legend) {
    common.legend
  } else {
    original.legend
  }
}
