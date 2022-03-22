repository <- "https://stat.ethz.ch/CRAN/"

if (!"stringr" %in% installed.packages())
  install.packages("stringr", repos = repository)
library(stringr)


prettify.colors <- function(tool.names) {
  sorted.tool.names <- character()
  sorted.tool.fills <- character()
  sorted.tool.colors <- character()

  tool.cellphy.names <- character()
  tool.cellphy.fills <- character()
  tool.cellphy.colors <- character()

  tool.sifit.names <- character()
  tool.sifit.fills <- character()
  tool.sifit.colors <- character()

  tool.monovar.names <- character()
  tool.monovar.fills <- character()
  tool.monovar.colors <- character()

  tool.sciphi.names <- character()
  tool.sciphi.fills <- character()
  tool.sciphi.colors <- character()

  tool.sieve.names <- character()

  predefined.colors.cellphy <- c(
    "#a0bc1d",
    "#d0e766",
    "#CBE453",
    "#d52eff",
    "#a12eff"
    )
  predefined.colors.sifit <- c(
    "#149af2",
    "#72c2f7",
    "#7FC8F8",
    "#27d81b",
    "#096000"
    )
  predefined.colors.monovar <- c(
    "#3c5488",
    "#7c93c5",
    "#6681bb",
    "#2897ff",
    "#0000ef"
    )
  predefined.colors.sciphi <- c(
    "#c06c00",
    "#ffab40",
    "#ff9d20",
    "#FF9914",
    "#ebc90b"
    )
  predefined.colors.sieve <- c(
    "#e64b35",
    "#f09285",
    "#ed8071",
    "#F21B3F",
    "#ff3d3d",
    "#c70000")

  # cellphy, sifit, monovar, sciphi, sieve
  color.index <- c(1, 1, 1, 1, 1)
  for (tool in tool.names) {
    if (grepl("^cellphy", tool, ignore.case = TRUE)) {
      tool.cellphy.names <- c(tool.cellphy.names, tool)
      tool.cellphy.colors <- c(tool.cellphy.colors, predefined.colors.cellphy[color.index[1]])
      tool.cellphy.fills <- c(tool.cellphy.fills, predefined.colors.cellphy[color.index[1] + 1])
      color.index[1] <- color.index[1] + 2
    } else if (grepl("^sifit", tool, ignore.case = TRUE)) {
      tool.sifit.names <- c(tool.sifit.names, tool)
      tool.sifit.colors <- c(tool.sifit.colors, predefined.colors.sifit[color.index[2]])
      tool.sifit.fills <- c(tool.sifit.fills, predefined.colors.sifit[color.index[2] + 1])
      color.index[2] <- color.index[2] + 2
    } else if (grepl("^monovar", tool, ignore.case = TRUE)) {
      tool.monovar.names <- c(tool.monovar.names, tool)
      tool.monovar.colors <- c(tool.monovar.colors, predefined.colors.monovar[color.index[3]])
      tool.monovar.fills <- c(tool.monovar.fills, predefined.colors.monovar[color.index[3] + 1])
      color.index[3] <- color.index[3] + 2
    } else if (grepl("^sciphi", tool, ignore.case = TRUE)) {
      tool.sciphi.names <- c(tool.sciphi.names, tool)
      tool.sciphi.colors <- c(tool.sciphi.colors, predefined.colors.sciphi[color.index[4]])
      tool.sciphi.fills <- c(tool.sciphi.fills, predefined.colors.sciphi[color.index[4] + 1])
      color.index[4] <- color.index[4] + 2
    } else if (grepl("^sieve", tool, ignore.case = TRUE)) {
      tool.sieve.names <- c(tool.sieve.names, tool)
    }
  }

  sorted.tool.names <- c(sorted.tool.names, tool.cellphy.names, tool.sifit.names, tool.monovar.names, tool.sciphi.names)
  sorted.tool.colors <- c(sorted.tool.colors, tool.cellphy.colors, tool.sifit.colors, tool.monovar.colors, tool.sciphi.colors)
  sorted.tool.fills <- c(sorted.tool.fills, tool.cellphy.fills, tool.sifit.fills, tool.monovar.fills, tool.sciphi.fills)

  if (length(tool.sieve.names) > 0) {
    tool.sieve.names <- sort.sieve.names(tool.sieve.names)
    
    if (length(tool.sieve.names) > length(predefined.colors.sieve)) {
      stop("Short of predefined colors. Please report this bug.")
    }
    
    step <- 2 * floor(length(predefined.colors.sieve) / length(tool.sieve.names))
    
    for (index in seq_along(tool.sieve.names)) {
      sorted.tool.names <- c(sorted.tool.names, tool.sieve.names[index])
      sorted.tool.colors <- c(sorted.tool.colors, predefined.colors.sieve[(index - 1) * step + 1])
      sorted.tool.fills <- c(sorted.tool.fills, predefined.colors.sieve[(index - 1) * step + 2])
    }
  }
  
  names(sorted.tool.colors) <- sorted.tool.names
  names(sorted.tool.fills) <- sorted.tool.names
  
  return(list(
    tool = sorted.tool.names,
    color = sorted.tool.colors,
    fill = sorted.tool.fills
  ))
}


sort.sieve.names <- function(sieve.names) {
  # sorted.sieve.names.1 <- character()
  # sorted.sieve.names.2 <- character()
  # sorted.sieve.names.3 <- character()
  # sorted.sieve.names.4 <- character()
  
  # for (index in 1:length(splited.sieve.names)) {
  #   if (splited.sieve.names[[index]][2] == "n") {
  #     if (grepl("-", sieve.names[index], fixed = TRUE)) {
  #       sorted.sieve.names.4 <- c(sorted.sieve.names.4, sieve.names[index])
  #     }
  #     else {
  #       sorted.sieve.names.1 <- c(sorted.sieve.names.1, sieve.names[index])
  #     }
  #   } 
  #   else if(splited.sieve.names[[index]][2] == "f") {
  #     if (grepl("[uF]", splited.sieve.names[[index]][3], fixed = FALSE)) {
  #       sorted.sieve.names.2 <- c(sorted.sieve.names.2, sieve.names[index])
  #     } else {
  #       sorted.sieve.names.3 <- c(sorted.sieve.names.3, sieve.names[index])
  #     }
  #   }
  # }
  
  # return(c(sorted.sieve.names.1, sorted.sieve.names.2, sorted.sieve.names.3,
  # sorted.sieve.names.4))
  
  sorted.sieve.names.1 <- character()
  sorted.sieve.names.2 <- character()
  sorted.sieve.names.3 <- character()
  sorted.sieve.names.4 <- character()

  for (name in str_sort(sieve.names)) {
    if (grepl("none", name, ignore.case = TRUE, fixed = TRUE)) {
      sorted.sieve.names.1 <- c(sorted.sieve.names.1, name)
    } else if (grepl("fels", name, ignore.case = TRUE, fixed = TRUE)) {
      sorted.sieve.names.2 <- c(sorted.sieve.names.2, name)
    } else if (grepl("2stage", name, ignore.case = TRUE, fixed = TRUE)) {
      sorted.sieve.names.3 <- c(sorted.sieve.names.3, name)
    } else {
      sorted.sieve.names.4 <- c(sorted.sieve.names.4, name)
    }
  }

  return(c(sorted.sieve.names.1, sorted.sieve.names.2, sorted.sieve.names.3, sorted.sieve.names.4))
}


rename.sieve <- function(tool.name, tool.setup) {
  if (!is.na(tool.setup)) {
    if (grepl("-", tool.setup, fixed = TRUE)) {
      return(paste(tool.name, "2stage", sep = "_"))
    } else if(grepl("^none_univ", tool.setup, ignore.case = TRUE, fixed = FALSE) || grepl("^none_univ", tool.setup, ignore.case = TRUE, fixed = FALSE)) {
      return(paste(tool.name, "none_univ", sep = "_"))
    } else if (grepl("^fels_fixed", tool.setup, ignore.case = TRUE, fixed = FALSE) || grepl("^fels_fixed", tool.setup, ignore.case = TRUE, fixed = FALSE)) {
      return(paste(tool.name, "fels_fixed", sep = "_"))
    } else if (grepl("^fels_univ", tool.setup, ignore.case = TRUE, fixed = FALSE)) {
      return(paste(tool.name, "fels_univ", sep = "_"))
    } else {
      stop(paste("Unhandled condition encountered:", tool.name, "and", tool.setup, sep = " "))
    }
  }
}


rename.tool <- function(tool.name, snv.type, tool.setup) {
  renamed.tool.setup <- ""
  if (!is.null(tool.setup) && !is.na(tool.setup)) {
    if (grepl("true_parameters$", tool.setup, ignore.case = TRUE)) {
      renamed.tool.setup <- "tp"
    } else if (grepl("true_parameters_zero$", tool.setup, ignore.case = TRUE)) {
      renamed.tool.setup <- "tp0"
    } else if (grepl("^ep$", tool.setup, ignore.case = TRUE)) {
      renamed.tool.setup <- "EP"
    } else if (grepl("^pl$", tool.setup, ignore.case = TRUE)) {
      renamed.tool.setup <- "PL"
    } else {
      stop(paste("Unhandled condition encountered:", tool.name, ",", snv.type, ",", tool.setup, sep = " "))
    }
  }

  renamed.snv.type <- ""
  if (!is.null(snv.type) && !is.na(snv.type)) {
    if (grepl("true_monovar_snvs$", snv.type, ignore.case = TRUE)) {
      renamed.snv.type <- "mnv"
    } else if (grepl("sieve_snvs$", snv.type, ignore.case = TRUE)) {
      renamed.snv.type <- "sv"
    } else {
      stop(paste("Unhandled condition encountered:", tool.name, ",", snv.type, ",", tool.setup, sep = " "))
    }
  }

  if (renamed.tool.setup != "" && renamed.snv.type != "") {
    return(paste(tool.name, renamed.snv.type, renamed.tool.setup, sep = "_"))
  } else if (renamed.tool.setup != "") {
    return(paste(tool.name, renamed.tool.setup, sep = "_"))
  } else if (renamed.snv.type != "") {
    return(paste(tool.name, renamed.snv.type, sep = "_"))
  } else {
    return(tool.name)
  }
}


rename.tools.internal <- function(row.data, rename.sieve, rename.cellphy, rename.sifit.monovar) {
  if (rename.sieve > 0 && grepl("sieve", row.data["tool"], ignore.case = TRUE)) {
    return(
      rename.sieve(
        row.data["tool"], 
        row.data["tool_setup"]
        )
      )
  } else if (rename.cellphy > 0 && grepl("cellphy", row.data["tool"], ignore.case = TRUE)) {
    return(
      rename.tool.helper(
        rename.cellphy, 
        row.data
        )
      )
  } else if (rename.sifit.monovar > 0 && grepl("sifit", row.data["tool"], ignore.case = TRUE)) {
    return(
      rename.tool.helper(
        rename.sifit.monovar, 
        row.data
        )
      )
  } else if (rename.sifit.monovar > 0 && grepl("monovar", row.data["tool"], ignore.case = TRUE)) {
    return(
      rename.tool.helper(
        rename.sifit.monovar, 
        row.data
        )
      )
  } else {
    return(row.data["tool"])
  }
}


rename.tool.helper <- function(rename.level, row.data) {
  if (rename.level == 0) {
    return(row.data["tool"])
  } else if (rename.level == 1) {
    return(
      rename.tool(
        row.data["tool"], 
        row.data["snv_type"], 
        NULL
        )
      )
  } else if (rename.level == 2) {
    return(
      rename.tool(
        row.data["tool"], 
        NULL, 
        row.data["tool_setup"]
        )
      )
  } else if (rename.level == 3) {
    return(
      rename.tool(
        row.data["tool"], 
        row.data["snv_type"], 
        row.data["tool_setup"]
        )
      )
  }
}


rename.tools <- function(data) {
  rename.sieve <- determine.rename.level(data[data$tool == "SIEVE",])
  rename.cellphy <- determine.rename.level(data[data$tool == "CellPhy",])
  rename.sifit.monovar <- determine.rename.level(data[data$tool == "SiFit",]) + determine.rename.level(data[data$tool == "Monovar",])

  return(
    apply(
      data, 1, rename.tools.internal, 
      rename.sieve = rename.sieve, 
      rename.cellphy = rename.cellphy, 
      rename.sifit.monovar = rename.sifit.monovar
      )
  )
}


# 0: do not rename
# 1: rename with "snv_type"
# 2: rename with "tool_setup"
# 3: rename with both "snv_type" and "tool_setup"
determine.rename.level <- function(data) {
  if (nrow(data) == 0) {
    return(0)
  }

  length.snv.type <- length(levels(factor(data$snv_type)))
  length.tool.setup <- length(levels(factor(data$tool_setup)))

  if (length.snv.type <= 1 && length.tool.setup <= 1) {
    return(0)
  } else if (length.snv.type > 1 && length.tool.setup > 1) {
    return(3)
  } else if (length.snv.type <= 1 && length.tool.setup > 1) {
    return(2)
  } else if (length.snv.type > 1 && length.tool.setup <= 1) {
    return(1)
  }
}
