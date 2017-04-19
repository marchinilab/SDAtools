#' Simulate 2D Data
#'
#' \code{simulate_2D_data} Simulate 2D Data
#'
#' @param n_individuals integer; the number of individuals to simualte
#' 
#' @param n_variables integer vector; vector of integers specifying the number of variables (e.g. genes) in each data type
#' 
#' @param n_components integer; number of true components to simulate from
#' 
#' @param sparsity numeric; value between 0 and 1 for the sparsity
#' 
#' @return list of Y, A, X and noise matrices
#'
#' @examples
#' 
#' 
#' @export
simulate_2D_data <- function(n_individuals=100, n_variables=500, n_components=10, sparsity=0.9) {
  A <- matrix(rnorm(n_individuals * n_components), n_individuals, n_components) # individual scores
  
  X <- matrix(rnorm(n_components * n_variables), n_components, n_variables) # loadings
  X[sample(1:length(X), sparsity * length(X))] <- 0  # make loadings sparse

  noise <- matrix(rnorm(n_individuals * n_variables, 0, 1), n_individuals, n_variables)

  Y <- A %*% X + noise

  data <- list()
  data$Y <- Y
  data$A <- A
  data$X <- X
  data$noise <- noise

  return(data)
}


#' check_simulation_scores
#'
#' \code{check_simulation_scores} check_simulation_scores
#' NB use nmf.options(grid.patch=TRUE) to stop blank pages appearing!
#' 
#' @param data character; location of data saved by simulate_and_save e.g. simulated_true/yourname_data.Rdata"
#' 
#' @param results list; output of load_results()
#' 
#' @return Panel of Heatmaps and scatter plots comparing the truth (data) to SDA inference (results)
#'
#' @examples
#' 
#' 
#' @export
#' @import data.table
#' @importFrom NMF aheatmap

check_simulation_scores <- function(data, results) {

  layout(matrix(c(1,1, 2,2, 3,5,4,6), ncol=4))
  
	colnames(data$A) = paste("True", seq_len(ncol(data$A)))
	colnames(results$scores) = paste("Est.", toupper(letters[26:1])[seq_len(ncol(results$scores))])
  
  aheatmap(data$A,
           Rowv=NA,
           Colv=NA,
           sub="Components",
           main="Scores (truth)")
  
  aheatmap(results$scores,
           Rowv=NA,
           Colv=NA,
           sub="Components",
           main="Scores (estimated)")
  
  aheatmap(cor(data$A, results$scores),
           Rowv= FALSE,
           Colv= FALSE,
           sub="Components",
           breaks=0, # always use beige as 0 correlation
           main="Corr. between true\n& estimated scores",
		   cexRow=0.5)
  
  # plot scatter plot for most correlated true for each estimated component
  # for each 3 components, choose estimated which has highest correlation
  for (y in sample(1:ncol(data$A),3)){
    for (i in 1:ncol(results$scores)){
      if (sd(results$scores[,i])!=0){ # skip this component estimate if all 0
        if (max(abs(cor(data$A[,y],results$scores)),na.rm=TRUE)==abs(cor(data$A[,y],results$scores[,i]))){
          plot(data$A[,y],results$scores[,i],
               xlab=paste(dimnames(data$A)[2][[1]][y]),
               ylab=paste(dimnames(results$scores)[2][[1]][i]))
          if(cor(data$A[,y],results$scores[,i])>0){
            abline(0, 1)
          }else{
            abline(0, -1)
          }
        }
      }else{
      }
    }
	}

}


#' Plot loadings matrices
#'
#' \code{plot_loadings} plot_loadings
#' NB use nmf.options(grid.patch=TRUE) to stop blank pages appearing!
#' 
#' @param data character; location of data saved by simulate_and_save e.g. simulated_true/yourname_data.Rdata"
#' 
#' @param data_type numeric; index of component to plot loadings for
#' @return Panel of Heatmaps of loadings
#'
#' @examples
#' 
#' 
#' @export
#' @importFrom NMF aheatmap
plot_loadings <- function(data, data_type) {
  layout(matrix(c(1,2,3), 3, 1, byrow = TRUE))
  aheatmap(data$X[[data_type]],
           Rowv=NA,
           Colv=NA,
           breaks=0, # always use beige as 0 correlation
           sub="L",
           main="Loadings")
}
