#' Plot free energy over time
#'
#' \code{plot_free_energy} Plot free energy over time
#'
#'
#' @param results list; object containing result  (output of read_output)
#' 
#' @return A ggplot2 object
#'
#' @examples
#' 
#' 
#' @export
#' @import ggplot2
plot_free_energy <- function(results) {
qplot(seq_len(ncol(results$free_energy)) * as.numeric(results$command$free_freq),
      results$free_energy[1, ],
      ylab="Free Energy",
      xlab="Iteration")
}


#' Plot change in free energy over time
#'
#' \code{plot_free_energy_change} Plot change in free energy over time
#'
#'
#' @param results list; object containing result  (output of read_output)
#' 
#' @return A ggplot2 object
#'
#' @examples
#' 
#' 
#' @export
#' @import ggplot2
plot_free_energy_change <- function(results) {
qplot(seq_along(diff(results$free_energy[1,]))  * as.numeric(results$command$free_freq),
      diff(results$free_energy[1,]),
      ylab="Change in Free Energy",
      xlab="Iteration") +
    scale_y_log10()
}


#' Plot overall distribution of gene loadings
#'
#' \code{loading_distribution} Plot overall distribution of gene loadings
#'
#'
#' @param results list; object containing result  (output of read_output)
#' 
#' @param omic numeric or character; index of the omic type to plot, or name of omic type if you assigned a name
#' 
#' @return A ggplot2 object
#'
#' @examples
#' 
#' 
#' @export
#' @import ggplot2
loading_distribution <- function(results, omic=1){
	qplot(as.numeric(results$loadings[[omic]]),
	      binwidth = 0.01,
	      main="Overall distribution of gene loadings",
	      xlab="Gene Loading")
}


#' Plot overall distribution of individual scores
#'
#' \code{scores_distribution} Plot overall distribution of individual scores
#'
#'
#' @param results list; object containing result  (output of read_output)
#' 
#' @return A ggplot2 object
#'
#' @examples
#' 
#' 
#' @export
#' @import ggplot2
scores_distribution <- function(results){
	qplot(as.numeric(results$scores),
	      binwidth = 0.01,
	      main = "Overall distribution of individual scores",
	      xlab="Score")
}


#' Plot maximum score and loading by component
#'
#' \code{plot_maximums} Plot maximum score and loading by component
#'
#'
#' @param results list; object containing result  (output of read_output)
#' 
#' @param omic numeric or character; index of the omic type to plot loading for, or name of omic type if you assigned a name
#' 
#' @return A ggplot2 object
#'
#' @examples
#' 
#' 
#' @export
#' @import data.table ggplot2 ggrepel
plot_maximums <- function(results, omic=1){
  
  if(is.null(results$component_statistics$max_score)){ # skip if previously calculated
    
  component_stats_temp <- data.table(
    Component = 1:results$n$components,
    Component_name = dimnames(results$scores)[[2]],
    max_score = apply(abs(results$scores), 2, max),
    max_loading = apply(abs(results$loadings[[omic]]), 1, max),
    mean_score = apply(abs(results$scores), 2, mean),
    mean_loading = apply(abs(results$loadings[[omic]]), 1, mean)
    )[order(-Component)]
  
  if(is.null(results$component_statistics)){
  results$component_statistics <<- component_stats_temp
  results$component_statistics <- component_stats_temp
  }else{
    results$component_statistics <<- cbind(results$component_statistics,
                                           component_stats_temp[, -c("Component","Component_name"), with=FALSE])
    results$component_statistics <- cbind(results$component_statistics,
                                           component_stats_temp[, -c("Component","Component_name"), with=FALSE])
  }
  
  }
  
ggplot(results$component_statistics, aes(max_score, max_loading, label = Component_name)) +
  geom_point() +
  geom_label_repel()
}


#' Scree Plot
#'
#' \code{plot_scree} Plot scree
#'
#'
#' @param results list; object containing result  (output of read_output)
#'
#' @param omic numeric or character; index of the omic type to calculate contributions for, or name of omic type if you assigned a name
#'  
#' @return A ggplot2 object
#'
#' @examples
#' 
#' 
#' @export
#' @import ggplot2
plot_scree <- function(results, omic=1){
  
  if(is.null(results$component_statistics$product_sum)){ # skip if previously calculated
    
  # Calculate product of score and loadings vectors for each component
  # get outer product statistics for each component
  summarise_component_matrix <- function(x){
    temp <- results$scores[ , x, drop=F] %*% results$loadings[[omic]][x, , drop=F] # i.e. outer product
    c(sum(abs(temp)), sd(temp))
  }
  
  x <- matrix(nrow = nrow(results$loadings[[omic]]), ncol=2)
  for (i in 1:nrow(results$loadings[[omic]])){
    x[i, ] <- summarise_component_matrix(i)
  }
  
  component_stats_temp <- data.table(
    Product_sd = x[,2],
    Product_sum = x[,1],
    Component = seq_len(results$n$components),
    Component_name = dimnames(results$scores)[[2]]
  )[order(-Component)]
  
  if(is.null(results$component_statistics)){
    results$component_statistics <<- component_stats_temp
    results$component_statistics <- component_stats_temp
  }else{
    results$component_statistics <<- cbind(results$component_statistics,
                                           component_stats_temp[, -c("Component","Component_name"), with=FALSE])
    results$component_statistics <- cbind(results$component_statistics,
                                           component_stats_temp[, -c("Component","Component_name"), with=FALSE])
  }
  
  }
  
ggplot(results$component_statistics[order(-Product_sd)], aes(factor(Component, levels = Component), Product_sd)) +
  geom_bar(stat="identity") +
  labs(x="Component")
}


#' Plot PIP fraction <0.5 by iteration
#'
#' \code{plot_PIP} Plot PIP fraction <0.5 by iteration
#'
#'
#' @param results list; object containing result  (output of read_output)
#' 
#' @param burn_in integer; number of iterations to skip plotting at the start
#' 
#' @return A ggplot2 object
#'
#' @examples
#' 
#' 
#' @export
#' @import ggplot2
plot_PIP <- function(results, burn_in = 50) {
	# excluding first 10 as burn in
	qplot(seq_along(results$pip_fraction[-seq_len(burn_in)]),
	      results$pip_fraction[-seq_len(burn_in)],
	      geom="line",
	      ylab="% PIP < 0.5",
	      xlab="Iteration")
}

#' Plot PIP distribution
#'
#' \code{PIP_distribution} Plot final PIP distribution over all components
#'
#'
#' @param results list; object containing result  (output of read_output)
#' 
#' @param omic numeric or character; index of the omic type to calculate contributions for, or name of omic type if you assigned a name
#' 
#' @return A ggplot2 object
#'
#' @examples
#' 
#' 
#' @export
#' @import ggplot2
PIP_distribution <- function(results, omic=1){
	qplot(as.numeric(results$pips[[omic]]), geom="histogram", binwidth=0.005) +
    xlab("PIP")
}

#' Plot PIP distribution for a single component
#'
#' \code{PIP_component_distribution} Plot final PIP distribution for a single component
#'
#'
#' @param results list; object containing result  (output of read_output)
#' 
#' @param component integer or character; index or name of the component you want to plot
#' 
#' @param omic integer or character; index of the omic type to calculate contributions for, or name of omic type if you assigned a name
#' 
#' @return A ggplot2 object
#'
#' @examples
#' 
#' 
#' @export
#' @import ggplot2
PIP_component_distribution <- function(results, component, omic=1){
  qplot(results$pips[[omic]][component,],
        geom="histogram",
        binwidth=0.005) +
    xlab("pip")
}

#' Plot faction PIP < 0.5 distribution
#'
#' \code{PIP_threshold_distribution} Plot faction PIP < 0.5 distribution
#'
#'
#' @param results list; object containing result  (output of read_output)
#' 
#' @param omic integer or character; index of the omic type to calculate contributions for, or name of omic type if you assigned a name
#' 
#' @return A ggplot2 object
#'
#' @examples
#' 
#' 
#' @export
#' @import ggplot2
PIP_threshold_distribution <- function(results, omic=1){
  qplot(apply(results$pips[[omic]], 1, function(x) sum(x<0.5)) / ncol(results$pips[[omic]]),
	      main="Distribution of fraction of PIPs < 0.5",
	      xlim=c(0,1),
	      binwidth=0.01,
	      xlab="fraction of PIPs < 0.5 for a component")
}

##### Functions for looking at Loadings matrix

#' Which component has the highest loading for a given gene
#'
#' \code{which_component} plot component loadings for a given gene
#'
#'
#' @param results list; object containing result  (output of read_output)
#' 
#' @param omic integer or character; index of the omic type to calculate contributions for, or name of omic type if you assigned a name
#' 
#' @param variable_name character string; name or index of variable e.g. gene
#' 
#' @param top_n integer; the number of components to label
#' 
#' @return A ggplot2 object
#'
#' @examples
#' 
#' 
#' @export
#' @import ggplot2 scales
# 
highest_components <- function(results, variable_name, top_n=5, omic=1){

  reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, 
              log_breaks(base = base), 
              domain = c(1e-100, Inf))
    }
  
  temp <- data.table(component = seq_len(results$n$components),
                   loading = results$loadings[[omic]][,variable_name])
  
  for (c in seq_len(results$n$components)){
    temp[component==c, rank := which(names(sort(-abs(results$loadings[[omic]][c,])))==variable_name)]
  }

  ggplot(temp, aes(loading, rank)) +
    geom_point() +
    scale_y_continuous(trans=reverselog_trans(10)) +
    annotation_logticks(side="l") + 
    geom_label_repel(data = temp[order(rank)][seq_len(top_n)], 
                     aes(label = component), 
                     size = 3, 
                     box.padding = unit(0.5, "lines"), 
                     point.padding = unit(0.1, "lines"), 
                     force=1,
                     segment.size=0.2,
                     segment.color="blue") +
    ggtitle(paste("Components with highest loading for variable:", variable_name))
}

#' Which gene has the highest loading for a given component
#'
#' \code{which_gene} plot gene loadings for a given component
#'
#'
#' @param results list; object containing result  (output of read_output)
#' 
#' @param omic integer or character; index of the omic type to calculate contributions for, or name of omic type if you assigned a name
#' 
#' @param component character string; name or index of component
#' 
#' @param top_n integer; the number of components to label
#' 
#' @return A ggplot2 object
#'
#' @examples
#' 
#' 
#' @export
#' @import ggplot2
# 
# 
highest_genes <- function(results, component, omic=1, max.items = 20, label.size = 3, label.repel = 1){

temp <- data.table(value = results$loadings[[omic]][component,],
                   name = colnames(results$loadings[[omic]]))

temp[, ranking := rank(value, ties.method = "first")]
temp[, gene_index := seq_along(temp$value)]


ggplot(temp, aes(gene_index, value)) +
  geom_point() +
  geom_label_repel(data = temp[abs(value) > 0.1][order(-abs(value))][1:max.items], 
                   aes(label = gsub("^[XY0-9:-]+:", "", name)),
                   size = label.size,
                   box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.1, "lines"),
                   force = label.repel,
                   segment.size=0.2,
                   segment.color="blue") +
  ggtitle(paste("Gene Loading for Component:",component)) +
  ylab("Loading") +
  xlab("Gene Index") +
  expand_limits(y = c(-0.1,0.1))
}


#' Load genomic locations from biomart
#'
#' \code{genome_loadings} Load genome locations for a list of genes
#'
#' @param gene.symbols character vector; Vector of gene names to look up genomic coordinates for
#' 
#' @param data.set character; which dataset from Ensembl should be used, e.g. "mmusculus_gene_ensembl" or "hsapiens_gene_ensembl"
#'
#' @param gene_name character; which biomart value to use for matching e.g. "external_gene_name"
#' 
#' @return stuff
#'
#' @examples
#' 
#' 
#' @import biomaRt data.table
get.location <- function(gene.symbols, data_set, gene_name){
  
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org", dataset = data_set)

mapTab <- getBM(attributes = c(gene_name,'chromosome_name','start_position'),
				filter = gene_name, values = gene.symbols, mart = ensembl, uniqueRows=TRUE)

mapTab <- as.data.table(mapTab)

setnames(mapTab, gene_name,"gene_symbol")

# Remove duplicate genes!!
# first which genes are duplicated
duplicates <- mapTab[duplicated(mapTab, by="gene_symbol")]$gene_symbol
# order by chr so patch versions to go bottom, then choose first unique by name
unduplicated <- unique(mapTab[gene_symbol %in% duplicates][order(chromosome_name)], by="gene_symbol")

# remove duplicates and replace with unique version
mapTab <- mapTab[!gene_symbol %in% duplicates]
mapTab <- rbind(mapTab, unduplicated)

# change all patch chr names to NA
mapTab[!chromosome_name %in% c(1:22,"X","Y","MT")]$chromosome_name <- NA # mice actually 19
mapTab[is.na(chromosome_name)]$start_position <- NA
return(mapTab)
}


#' Load gene coordinates from Ensembl biomart or Cached data
#'
#' \code{load_gene_locations} Load gene coordinates from Ensembl biomart or Cached data
#'
#' @param genes character vector; Vector of gene names to look up genomic coordinates for
#' 
#' @param cache logical; If TRUE a cache of the downloaded locations will be saved for quick loading next time
#' 
#' @param path character; string of path where the cached location file should be saved
#' 
#' @param organism character; which dataset from Ensembl should be used, e.g. "mmusculus_gene_ensembl" or "hsapiens_gene_ensembl"
#' 
#' @param name character; descriptive name of dataset, this will be appended to the standard cache file name when saving / loading the cache
#' 
#' @return A data table of genes chromosome and coordinate. An object named chromosome.lengths will also be created in the global environment.
#'
#' @examples
#' 
#' 
#' @export
#' @import data.table
load_gene_locations <- function(genes=NULL, cache=TRUE, path="", organism="mmusculus_gene_ensembl", name="mouse"){
if (file.exists(paste0(path, "SDAtools_gene_locations_", name, ".rds"))==FALSE){
	gene_locations <- get.location(genes, data_set = organism, gene_name = "external_gene_name")
	saveRDS(gene_locations, paste0(path,"SDAtools_gene_locations_",name,".rds"))
}
# assign("rna.location", readRDS(paste0(path,"SDAtools_gene_locations_",name,".rds")), envir=globalenv())

# set chromosome lengths for calculating plotting position
if (organism=="mmusculus_gene_ensembl"){
# http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/mouse/data/
chromosome.lengths <- data.table(chromosome=factor(c(1:19,"X","Y","MT","Un")),
                                 length=c(196283350,182113224,160039680,157219549,153573022,149736546,145617427,129401213,124595110,130694993,122082543,
                                          120129022,120421639,124902244,104043685,98207768,94987271,90702639,61431566,171368232,92500857,18000,803895))

}else if(organism=="hsapiens_gene_ensembl"){
# http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/data/ GRCh38.p7
chromosome.lengths <- data.table(chromosome=factor(c(1:22,"X","Y","MT","Un")),
                                 length=c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,
                                          114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468,156040895,57227415,22222222,4485509))
}else{
  print("Error: Organism chromosme lengths not found / not currently supported")
}
chromosome.lengths[,length_padded:=length+5000000]
chromosome.lengths[length<5e7, length_padded:=length_padded+5e7]
chromosome.lengths[,genomic_offset:=cumsum(as.numeric(length_padded))-(length_padded)]
chromosome.lengths[,center := genomic_offset + length/2]
setkey(chromosome.lengths,chromosome)

assign("chromosome.lengths", chromosome.lengths, envir=globalenv())

return(readRDS(paste0(path,"SDAtools_gene_locations_",name,".rds")))
}


#' Plot Component Loadings by Genomic Location
#'
#' \code{genome_loadings} Plot Component Loadings by Genomic Location
#'
#' @param component numeric; A named numeric vector containing the gene loadings and their names as the attribute "names"
#' 
#' @param gene_locations data.table; The output of load_gene_locations() containing the gene chromosome location coordinates
#' 
#' @param label_both logical; If TRUE the top N genes will be labeled for positive and negative loadings independently. 
#' If FALSE only genes with either positive OR negative loadings will be labeled, depending on which side has the highest loading.
#' 
#' @param max.items integer; Number of genes to label each side of 0.
#' I.e. if label_both is TRUE then the total number of labels will be twice max.items.
#' 
#' @return A data table of statistics and their values
#'
#' @examples
#' genome_loadings(results$loadings[[1]][48,])
#' 
#' @export
#' @import ggplot2 ggrepel
genome_loadings <- function(component = NULL,
                            max.items = 20,
                            label.size = 3,
                            label.repel = 1,
                            label_both = FALSE,
                            label_X = FALSE,
                            min_loading = 0.01,
                            gene_locations=rna_locations){

temp <- data.table(loading = component, gene_symbol = names(component))

setkey(temp, gene_symbol)
setkey(gene_locations, gene_symbol)
temp <- merge(temp, gene_locations, all.x = TRUE)

temp$chromosome <- factor(temp$chromosome_name, levels = c(1:22,"X","Y","MT","Un"))

temp[is.na(chromosome)]$chromosome <- "Un"
temp[chromosome=="Un"]$start_position <- sample(1:chromosome.lengths[chromosome=="Un"]$length, nrow(temp[chromosome=="Un"]))

setkey(temp,chromosome)
temp <- chromosome.lengths[temp]

temp[, genomic_position := genomic_offset + start_position]


if (label_both == TRUE){
	label_data <- temp[loading > min_loading][order(-loading)][1:max.items] # positives
	label_data <- rbind(label_data,
	temp[loading < (-min_loading)][order(loading)][1:max.items]) # add negatives
} else {
	label_data <- temp[abs(loading) > min_loading][order(-abs(loading))][1:max.items]
}

if (label_X == TRUE){
	label_data <- rbind(label_data,
	temp[abs(loading) > min_loading][chromosome_name == "X"])
}


P <- ggplot(temp, aes(genomic_position, loading)) +
	geom_point(size = 0.5, aes(color = chromosome)) +
	scale_colour_manual(values = c(rep_len(c("black", "cornflowerblue"), length(levels(temp$chromosome))), "grey")) +
	xlab("Genomic Coordinate") +
	ylab("Loading") +
	theme(legend.position = "none") +
	scale_x_continuous(breaks = chromosome.lengths$center, labels = chromosome.lengths$chromosome, minor_breaks=NULL) +
	geom_label_repel(data = label_data, 
		aes(label = gene_symbol), 
		size = label.size, 
		box.padding = unit(0.5, "lines"), 
		point.padding = unit(0.1, "lines"), 
		force= label.repel,
		segment.size=0.2,
		segment.color="blue")
	#	expand_limits(y = c(-0.2,0.2)) +, 
	return(P)
}

#' Plot Individual Scores
#'
#' \code{plot_scores} Plot Individual Scores
#'
#' @param results list; object containing result  (output of read_output)
#' 
#' @param component integer; The component index to plot scores for
#' 
#' @return A ggplot2 object
#'
#' @examples
#' 
#' 
#' @export
#' @import ggplot2
plot_scores <- function(results, component){
  tmp <- data.table(cell_index = seq_len(nrow(results$scores)), score = results$scores[,component])
  ggplot(tmp, aes(cell_index, score)) +
    geom_point(size=0.5) +
    xlab("Cell Index") +
    ylab("Score")
}

#' Check Convergence
#'
#' \code{check_convergence} Check Convergence
#'
#' @param results list; object containing result  (output of read_output)
#' 
#' @return A ggplot2 object
#'
#' @examples
#' 
#' 
#' @export
#' @import ggplot2
#' @importFrom cowplot ggdraw draw_plot
check_convergence <- function(results){
  ggdraw() +
    draw_plot(plot_free_energy(results), 0, 0.5, .5, .5) +
    draw_plot(plot_free_energy_change(results), .5, 0.5, .5, .5) +
    draw_plot(plot_PIP(results), 0, 0, 1, .5)
  # plot_grid(plot_free_energy(results),
  #           plot_free_energy_change(results),
  #           plot_PIP(results), nrow = 1)
}