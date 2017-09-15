#' run_SDA
#'
#' \code{run_SDA} Run SDA
#'
#' @param sda_location string; location of the sda executable e.g. "../SDA/build/sda"
#' the default value 'sda' assumes the location of sda is in your $PATH variable
#' 
#' @param All the other parameters are identical to those for the command line bash arguments. Please see the SDA manual for the descriptions.
#' NB Not all parameters are implemented here.
#' 
#' @return Nothing, the results are saved on disk in the directory specified in the parameter 'out'.
#'
#' @examples
#' 
#' 
#' @export

run_SDA <- function(sda_location = "sda",
	out = NULL,
	data = NULL,
	N = NULL,
	num_comps = 10,
	max_iter = 5000,
	save_freq = 5001,
	free_freq = 50,
	set_seed = NULL,
	ignore_missing = FALSE,
	remove_zero_comps = TRUE,
	num_blocks = 1,
	num_openmp_threads = 1) {
  
	if(is.null(sda_location) | is.null(out) | is.null(data) | is.null(N)) stop("Must Specify all required parameters")
  
	if(dir.exists(out)) stop("That results directory already exists")
	
	dir.create(out)

	if(!is.null(set_seed)){
		seed <- paste0(" --set_seed ",set_seed)
	}else{
		seed <- ""
	}
	
	command <- paste0(sda_location,
		' --data ',data,
		' --N ',N,
		' --out ',out,
		' --num_comps ',num_comps,
		' --max_iter ',max_iter,
		' --save_freq ',save_freq,
		' --free_freq ',free_freq,
		seed,
		' --ignore_missing ',ignore_missing,
		' --remove_zero_comps ',remove_zero_comps,
		' --num_blocks ',num_blocks,
		' --num_openmp_threads ',num_openmp_threads,
		' > ',out,'/log.txt')
	
	system(command)
  
}

#' Export data and labels
#'
#' \code{export_data} Export data and labels
#'
#' @param matrix numeric matrix; the matrix you wish to save in a format readable by SDA
#' 
#' @param path character; the directory in which you want the data to be saved
#' 
#' @param name character; the name of the data, will be included in the file name
#' 
#' @return Nothing, the data is saved on disk in the directory specified in the parameter 'path'.
#'
#' @examples
#' 
#' 
#' @export

export_data <- function(matrix, path=NULL, name){

	# Create directory if it doesn't exist
	if(!is.null(path)){
		if(!dir.exists(path)){
			dir.create(path)
		}
	}

	# Save a copy of the dimension labels
	saveRDS(dimnames(matrix), paste0(path,name,"_dimnames.rds"))

	# Create file for input to SDA
	fwrite(data.table(signif(matrix, 3)),
	       sep = " ",
	       paste0(path,name,".data"),
	       row.names = FALSE,
	       col.names = FALSE,
	       quote = FALSE)

}