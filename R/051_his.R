

#' Hierarchical Interaction Score -- read results from matlab and write into R object
#'
#' @param project
#' @return nothing
#' @export
step_051_his <- function( project, nem_method, perturbed_genes, prepared_dir, nems_dir) {
    matching <- dir(prepared_dir, pattern=project)
    matching <- Filter(function(x) grepl(paste("\\.", "csv", sep=""), x), matching)
    matching <- Filter(function(x) grepl(nem_method, x), matching)
    for (input_file_name in matching) {
        d <- read.csv(file.path(prepared_dir, input_file_name), header=FALSE)
        rownames(d) <- perturbed_genes
        colnames(d) <- perturbed_genes
        l <- list()
        l$graph <- d
        output_file_name <- gsub("csv", "Rds", input_file_name)
        saveRDS(l, file.path(nems_dir, output_file_name))
    }

}

