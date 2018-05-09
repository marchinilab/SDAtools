#' run_SDA
#'
#' \code{run_SDA} Run SDA
#'
#' @param sda_location string; location of the sda executable e.g. "../SDA/build/sda"
#' the default value 'sda' assumes the location of sda is in your $PATH variable
#'
#' @param out string; name of folder to be created to put results in
#'
#' @param data string; name of exported file with data in, should be equal to the value you set as the 'name' parameter for export_data.
#'
#' @param N int; number of individuals, will be inferred automatically by numer of lines in data file if not provided
#' but this will only work if the data is 2D / a single tissue type / data / omic type!
#'
#' @param num_comps int; how many components to infer (maximum).
#' @param max_iter int; maximum number of iterations to run for.
#' @param save_freq int; results will be saved every save_freq iterations. Defaults to max_iter/2
#' @param free_freq int; free energy will be calculated every free_freq ierations. Defaults to max_iter/20
#' @param set_seed string; random numbers to initiate starting position. e.g. "79151 17351". Can be left blank and a seed will be generated.
#' @param ignore_missing logical; Should sda ignore samples with completely missing data.
#' @param remove_zero_comps logical; should components with all PIPs < 0.05 be removed
#' @param eigen_parallel logical; use multi-threading options in eigen.
#' @param num_blocks int; number of blocks to split loadings matrix into for parallel updates.
#' @param num_openmp_threads int; number of threads to use in parallel mode
#' @param debug logical; if true will calculate free energy at every iteration (slower)
#' @param save_everything logical; if true will save all parameter estimated including for example noise precision.
#'
#' @description The parameters are identical to those for the command line bash arguments. Please see the SDA manual for the descriptions.
#' NB Not all parameters are implemented here.
#'
#' @return Nothing, the results are saved on disk in the directory specified in the parameter 'out'.
#'
#'
#' @export

run_SDA <- function(sda_location = "sda",
	out = NULL,
	data = NULL,
	N = NULL,
	num_comps = 10,
	max_iter = 500,
	save_freq = NULL,
	free_freq = NULL,
	set_seed = NULL,
	ignore_missing = FALSE,
	remove_zero_comps = TRUE,
	eigen_parallel = FALSE,
	num_blocks = 1,
	num_openmp_threads = 1,
	save_everything = FALSE,
	debug = FALSE) {

	if(is.null(sda_location) | is.null(out) | is.null(data)) stop("Must Specify all required parameters")

  if(!file.exists(data)) stop("data file not found")

  if(is.null(N)){
    warning("Inferring number of individuals automatically - for 2D data only!")
    N <- as.numeric(system(paste0("wc -l < ",data), intern = T))
  }

	if(dir.exists(out)) stop("That results directory already exists")

	dir.create(out)

	if(is.null(save_freq)){
	  save_freq <- max_iter/2
	}

	if(is.null(free_freq)){
	  free_freq <- max_iter/20
	}

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
		' --eigen_parallel ',eigen_parallel,
		' --save_everything ',save_everything,
		' --debug ',debug,
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
#' @param name character; the name of the data file
#'
#' @return Nothing, the data is saved on disk in the directory specified in the parameter 'path'.
#'
#' @examples
#' data <- simulate_2D_data()
#' export_data(data$Y, name = "simulated.data")
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
	saveRDS(dimnames(matrix), paste0(path,tools::file_path_sans_ext(name),"_dimnames.rds"))

	# Create file for input to SDA
	fwrite(data.table(signif(matrix, 3)),
	       sep = " ",
	       paste0(path,name),
	       row.names = FALSE,
	       col.names = FALSE,
	       quote = FALSE)

}
