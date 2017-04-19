#' Read SDA output into R list
#'
#' \code{load_results} simulate digital gene expression matrix containing count data from given parameters
#'
#'
#' @param results_folder string; path to folder containing estimates, same as the string you passed to '--out' when running SDA
#' 
#' @param iteration integer; iteration number of save to load
#' 
#' @return A list of matrices
#'
#' @examples
#' 
#' 
#' @export
#' @import data.table

load_results <- function(results_folder, iteration) {
  
  folder <- paste0(results_folder, "/it", iteration)
  stopifnot(dir.exists(folder))
  
  out <- list()
  
  # read each file in folder as item in list
  files <- list.files(folder)
  for (file in files) {
    out[[file]] <- as.matrix(fread(paste0(folder, "/", file)))
  }
  
  out1 <- reformat_data(out) # collect X matrices together, same for S, B

  out1$n$iterations <- iteration
  
  out1$free_energy <- out$free_energy
  out1$miss <- out$miss
  
  # if log saved, extract PIP fraction < 0.5
  if (file.exists(paste0(results_folder, "/log.txt"))){
    out1$pip_fraction <- fread(paste0("grep -o '[0]) : [0-9].*[0-9]*' ", results_folder, "/log.txt", sep=" "))[1:iteration]$V3
  }

  # if command.txt saved, extract arguments used
  if (file.exists(paste0(results_folder, "/command.txt"))){
    command <- readLines(paste0(results_folder, "/command.txt"), warn=F)
    command <- strsplit(command, "--")[[1]][-1] # split into argument value string and remove first string (location of sda)
    command <- regmatches(command, regexpr(" ", command), invert = TRUE) # split argument name from value
    
    values <- unlist(command)[c(FALSE,TRUE)]
    values <- gsub(" $","", values) # remove trailing space
    names(values) <- unlist(command)[c(TRUE,FALSE)]
    
    out1$command_arguments <- as.data.table(t(values))
  }
  
  return(out1)
}

reformat_data <- function(out) {
  matrix_names <- names(out)
  
  est <- list()
  est$scores <- out$A
  
  est$n <- list()
  est$n$individuals <- nrow(est$scores) # number individuals
  est$n$components <- ncol(est$scores) # number components
  est$n$omics <- length(matrix_names[grep("X[0-9]?", matrix_names)]) # num_X_mats 
  
  stopifnot(length(matrix_names[grep("S[0-9]?", matrix_names)]) == est$n$omics) # check all pips are loaded
  stopifnot(est$n$omics != 0) # at least one loading matrix
  est$n$context_matrices <- length(matrix_names[grep("B[0-9]?", matrix_names)]) # number of B matrices
  
  est$loadings <- list()
  out$pips <- list()
  out$context_scores <- list()
  
  for (d in seq_len(est$n$omics)) {
    est$loadings[[d]] <- out[[paste0("X", d)]]
    est$pips[[d]] <- out[[paste0("S", d)]]
    if(est$n$context_matrices > 0){
      est$context_scores[[d]] <- out[[paste0("B", d)]]
    }
  }

  return(est)
}
