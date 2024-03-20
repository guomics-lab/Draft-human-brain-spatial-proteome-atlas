# ---
# title: "PulseDIA combine"
# author: "Weigang Ge"
# author's unit: "Westlake Omics Ine"
# date: "2020/12/17"
# ---

## library package
#------------------------------------------------------
req.pcg.install <- function(pcg) {
  new <- pcg[!(pcg %in% installed.packages()[, "Package"])]
  if (length(new))
    
    install.packages(new, dependencies = T, repos
                     = "http://cran.rstudio.com")
  sapply(pcg, require, ch = T)
  
}
req.pcg.install("getopt")


spec <- matrix(
  c(
    'input',
    'i',
    1,
    'character',
    'required, File path',
    'out',
    'o',
    1,
    'character',
    'required, The output path and result file name',
    'combine_method',
    "m",
    2,
    "character",
    "optional,the combined value method ,mean or max, default mean",
    "DIA_method",
    "s",
    2,
    "character",
    "optional, the search software ,DIANN or Spectronaut or OpenSWATH, default DIANN",
    "outmatrix",
    "x",
    2,
    "character",
    "optional, the type of exported matrix ,peptide or protein, default protein",
    "fdr",
    "f",
    2,
    "double",
    "optional, if DIA_method is DIANN,the cutoff of FDR in the result extraction , default 0.01",
    "spn.pep.col",
    "p",
    2,
    "integer",
    "optional, if DIA_method is Spectronaut, specify which column the peptide is in, default 2",
    "spn.prot.col",
    "t",
    2,
    "integer",
    "optional, if DIA_method is Spectronaut, specify which column the protein is in, default 1",
    "spn.sample.col.start",
    "a",
    2,
    "integer",
    "optional, if DIA_method is Spectronaut, need to specify which column the quantitative result of the sample starts from, default 3",
    "pulse.bool",
    "b",
    2,
    "logical",
    "optional, Whether do PulseDIA merge, default T",
    
    "help",
    "h",
    0,
    "logical",
    "This is Help"
  ),
  byrow = TRUE,
  ncol = 5
)

opt <- getopt(spec = spec, debug = T)

if (!is.null(opt$help)) {
  cat(paste(getopt(spec = spec, usage = T), "\n"))
  quit()
}

install.pcg <-
  c("reshape2",
    "stringr")

req.pcg.install(install.pcg)


## Set default parameters
#----------------------------------------------------------------
# TMT original matrix
input_path <- opt$input

# Information correspondence table
out_name <- opt$out

# The combined value method
combine_method <-
  ifelse(!is.null(opt$combine_method), opt$combine_method, "mean")# mean or max

# Whether do PulseDIA merge
pulse.bool <-
  ifelse(!is.null(opt$pulse.bool), opt$pulse.bool, T)     # T or F

# the search software
DIA_method <-
  ifelse(!is.null(opt$DIA_method), opt$DIA_method , "DIANN")    # 'DIANN' or 'Spectronaut' or 'OpenSWATH'

# The type of exported matrix
outmatrix <-
  ifelse(!is.null(opt$outmatrix), opt$outmatrix , "protein")  # 'peptide' or 'protein'

#the cutoff of FDR in the result extraction
fdr <- ifelse(!is.null(opt$fdr), opt$fdr , 0.01)

spn.pep.col <-
  ifelse(!is.null(opt$spn.pep.col), opt$spn.pep.col , 2)  # specify which column the peptide is in

spn.prot.col <-
  ifelse(!is.null(opt$spn.prot.col), opt$spn.prot.col , 1) # specify which column the protein is in

spn.sample.col.start <-
  ifelse(!is.null(opt$spn.sample.col.start),
         opt$spn.sample.col.start ,
         3) # specify which column the quantitative result of the sample starts from


out_name <- gsub("//", "\\/", out_name)
input_path <- gsub("//", "\\/", input_path)
# print(paste0("test:", input_path))

pulse.combine <- function(data, method, filename) {

  #nm <- gsub("^\\d","",nm)
  
  if (outmatrix == "peptide") {
    nm <- as.character(sapply(colnames(data)[-(1:2)], function(v) {
      pa = str_split(v, "_part")[[1]]
      if (length(pa) > 2) {
        print("WARNING: Sample naming convention: The keyword '_part' can only appear once")
      }
      if (length(pa) == 1) {
        print("WARNING: Sample naming convention: The keyword '_part' must be included")
      } else{
        return(pa[length(pa) - 1])
      }
    }))
    names(data) <- c("peptide_group_label", "prot", nm)
    
    df0 <-
      data.frame(t(as.data.frame(lapply(
        data[, -c(1, 2)], as.numeric
      ))))
    df0$label <- names(data)[-c(1:2)]
    
    print("Data combining...")
    if (method == "mean") {
      aggdata <- aggregate(df0,
                           by = list(df0$label),
                           FUN = mean,
                           na.rm = TRUE)
    } else if (method == "max") {
      aggdata <- aggregate(df0,
                           by = list(df0$label),
                           FUN = max,
                           na.rm = TRUE)
    } else{
      print("WARNING: Other merging methods are not supported")
    }
    result <- data.frame(t(aggdata[, -c(1, ncol(aggdata))]))
    
    colnames(result) <- aggdata[, 1]
    result$peptide_group_label <- data$peptide_group_label
    result$prot <- data$prot
    result <- result[, c('peptide_group_label', 'prot', aggdata[, 1])]
    rownames(result) <- 1:dim(result)[1]
  } else if (outmatrix == "protein") {
    nm <- as.character(sapply(colnames(data)[-(1)], function(v) {
      pa = str_split(v, "_part")[[1]]
      if (length(pa) > 2) {
        print("WARNING: Sample naming convention: The keyword '_part' can only appear once")
      }
      if (length(pa) == 1) {
        print("WARNING: Sample naming convention: The keyword '_part' must be included")
      } else{
        return(pa[length(pa) - 1])
      }
    }))
    names(data) <- c("prot", nm)
    df0 <-
      data.frame(t(as.data.frame(lapply(
        data[, -c(1)], as.numeric
      ))))
    df0 <- data.frame(t(data[, -c(1)]))
    df0$label <- names(data)[-c(1)]
    
    print("Data combining...")
    if (method == "mean") {
      aggdata <- aggregate(df0,
                           by = list(df0$label),
                           FUN = mean,
                           na.rm = TRUE)
    } else if (method == "max") {
      aggdata <- aggregate(df0,
                           by = list(df0$label),
                           FUN = max,
                           na.rm = TRUE)
    } else{
      print("WARNING: Other merging methods are not supported")
    }
    result <- data.frame(t(aggdata[, -c(1, ncol(aggdata))]))
    
    colnames(result) <- aggdata[, 1]
    result$prot <- data$prot
    result <- result[, c('prot', aggdata[, 1])]
    rownames(result) <- 1:dim(result)[1]
  }
  
  
  print("Writing result matrix")
  write.table(
    result,
    file = filename,
    sep = "\t",
    col.names = T,
    row.names = F,
    quote = F
  )
  print("Completed ! ")
}


###############################################
#diann

read.diann.file <- function(path, outmatrix, fdr = 0.01) {
  sum.data.list <- list()
  sum.data <- c()
  file <- list.files(path, ".tsv", full.names = T)
  n = 1
  for (fl in file) {
    print(paste0("reading file:  ", fl))
    df <- read.delim(fl)
    bb <- which(df$Precursor.Quantity == 0)
    if (length(bb)) {
      df <- df[-bb, ]
    }
    df <- df[df$Q.Value < fdr & df$Protein.Q.Value < fdr, ]
    if (outmatrix == "peptide") {
      df1 <-
        tryCatch(
          data.frame(
            sample = df$Run,
            pep = df$Modified,
            prot = df$Protein.Ids,
            intesity = df$Precursor.Quantity
          ),
          error = function(e) {
            print("The specified column name does not exist")
            stop()
          }
        )
      sum.data.list[[n]] <-
        dcast(df1, pep + prot ~ sample, value.var = 'intesity', mean)
      
      if (n == 1) {
        sum.data <- sum.data.list[[1]]
      } else if (n > 1) {
        sum.data <-
          merge(sum.data,
                sum.data.list[[n]],
                by.y = c("pep", "prot"),
                all = T)
      }
    } else if (outmatrix == "protein") {
      if (length(which(df$Protein.Ids == ""))) {
        df <- df[-which(df$Protein.Ids == ""), ]
      }
      df1 <-
        tryCatch(
          data.frame(
            sample = df$Run,
            prot = df$Protein.Ids,
            intesity = df$PG.Quantity
          ),
          error = function(e) {
            print("The specified column name does not exist")
            stop()
          }
        )
      sum.data.list[[n]] <-
        dcast(df1, prot ~ sample, value.var = 'intesity', mean)
      
      if (n == 1) {
        sum.data <- sum.data.list[[1]]
      } else if (n > 1) {
        sum.data <-
          merge(sum.data,
                sum.data.list[[n]],
                by.y = c("prot"),
                all = T)
      }
    }
    
    n = n + 1
  }
  return(sum.data)
}




####################################################
#OpenSWATH

read.OSW.file <- function(path, outmatrix) {
  sum.data <- c()
  file <- list.files(path, ".txt", full.names = T)
  n = 1
  for (fl in file) {
    print(paste0("reading file:  ", fl))
    df <- read.delim(fl)
    if (outmatrix == "peptide") {
      if (n == 1) {
        sum.data <- df
      } else if (n > 1) {
        sum.data <-
          merge(sum.data,
                df,
                by.y = c(names(df)[1], names(df)[2]),
                all = T)
      }
    } else if (outmatrix == "protein") {
      print("Only OSW peptides combining is supported!")
      stop()
    }
    n = n + 1
  }
  return(sum.data)
}

###############################################
#spn

read.spn.file <- function(path, outmatrix) {
  sum.data <- c()
  file <- list.files(path, ".xls", full.names = T)
  n = 1
  for (fl in file) {
    print(paste0("reading file:  ", fl))
    df <- read.delim(fl)
    df[df == "Filtered"] <- NA
    df[df == "NaN"] <- NA
    if (outmatrix == "peptide") {
      df1 <-
        data.frame(pep = df[, spn.pep.col], prot = df[, spn.prot.col], df[, spn.sample.col.start:ncol(df)])
      if (n == 1) {
        sum.data <- df1
      } else if (n > 1) {
        sum.data <- merge(sum.data,
                          df1,
                          by.y = c("pep", "prot"),
                          all = T)
      }
      if (length(strsplit(names(sum.data)[3], '\\.')[[1]]) > 3) {
        names(sum.data)[-c(1:2)] <-
          sapply(names(sum.data)[-c(1:2)], function(v) {
            strsplit(v, '\\.')[[1]][4]
          })
      }
      
    } else if (outmatrix == "protein") {
      df1 <-
        data.frame(prot = df[, spn.prot.col], df[, spn.sample.col.start:ncol(df)])
      if (n == 1) {
        sum.data <- df1
      } else if (n > 1) {
        sum.data <- merge(sum.data,
                          df1,
                          by.y = c("prot"),
                          all = T)
      }
      if (length(strsplit(names(sum.data)[3], '\\.')[[1]]) > 3) {
        names(sum.data)[-c(1)] <-
          sapply(names(sum.data)[-c(1)], function(v) {
            strsplit(v, '\\.')[[1]][4]
          })
      }
    }
    
    n = n + 1
  }
  return(sum.data)
}



if (DIA_method == "DIANN") {
  sum.data <- read.diann.file(input_path, outmatrix, fdr)
} else if (DIA_method == "OpenSWATH" |
           DIA_method == "openswath" | DIA_method == "Openswath") {
  sum.data <- read.OSW.file(input_path, outmatrix)
  if (outmatrix == "protein") {
    print(
      "It is not recommended to directly use the protein matrix merger, please combine the peptide matrix, use proteomeExpert's pep2prot function merger"
    )
  }
} else if (DIA_method == "Spectronaut" |
           DIA_method == "spectronaut" | DIA_method == "spn") {
  sum.data <- read.spn.file(input_path, outmatrix)
  if (outmatrix == "protein") {
    print(
      "It is not recommended to directly use the protein matrix merger, please combine the peptide matrix, use proteomeExpert's peptoprot function merger"
    )
  }
}

if (pulse.bool) {
  pulse.combine(sum.data, combine_method, out_name)
} else{
  print("Writing result matrix")
  write.table(
    sum.data,
    file = out_name,
    sep = "\t",
    col.names = T,
    row.names = F,
    quote = F
  )
  print("Completed ! ")
}


