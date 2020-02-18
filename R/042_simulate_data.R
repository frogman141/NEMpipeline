
step_042_simulate_data <- function(project, alpha, beta, perturbed_genes, prepared_dir, lfc=FALSE, plots=FALSE) {
    m <- matrix( c( 0, 1, 1, 0, 0, 0, 0,
                    0, 0, 0, 1, 1, 0, 0,
                    0, 0, 0, 0, 0, 1, 0,
                    0, 0, 0, 0, 1, 0, 1,
                    0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0  ), byrow=TRUE, nrow=7, ncol=7)
    rownames(m) <- perturbed_genes
    colnames(m) <- perturbed_genes

    generated_data <- make_data(m)
    if (plots) {
        heatmap(as.matrix(generated_data), scale="none", Rowv=NA, Colv=NA)
    }
    generated_data <- set_fp_rate(generated_data, alpha)
    generated_data <- set_fn_rate(generated_data, beta)
    if (plots) {
        heatmap(as.matrix(generated_data), scale="none", Rowv=NA, Colv=NA)
    }
    
    prep_method <- paste(alpha, beta, sep="_")
    if (lfc) {
        lfc_data <- make_lfc(generated_data)
        output_file_name <- paste("lfc", prep_method, project, "Rds", sep=".")
        saveRDS(lfc_data, file.path(prepared_dir, output_file_name))
        output_file_name <- paste("lfc", prep_method, project, "csv", sep=".")
        write.csv(lfc_data, file.path(prepared_dir, output_file_name))
    }
    output_file_name <- paste("binary", prep_method, project, "Rds", sep=".")
    saveRDS(generated_data, file.path(prepared_dir, output_file_name))
}

make_lfc <- function(generated_data) {
    # convert binary data to lfc 
    
    len <- 3 * nrow(generated_data) * ncol(generated_data)

    d <- rnorm(len, mean=0, sd=1)
    effects <- d[d > 0.5 | d< -0.5]
    noneffects <- d[-0.5 < d & d < 0.5]

    v <- as.matrix(generated_data)

    idx <- which(generated_data == 1)
    v[idx] <- effects[1:length(idx)]

    idx <- which(generated_data == 0)
    v[idx] <- noneffects[1:length(idx)] 

    d <- as.data.frame(v)
    return(d)
}


set_fp_rate <- function(generated_data, alpha) {
    idx <- which(generated_data == 0)
    fps <- sample(idx, ceiling(alpha*length(idx)), replace=FALSE)

    v <- as.matrix(generated_data)
    v[fps] <- 1
    d <- as.data.frame(v)
    return(d)
}

set_fn_rate <- function(generated_data, beta) {
    idx <- which(generated_data == 1)
    fns <- sample(idx, ceiling(beta*length(idx)), replace=FALSE)

    v <- as.matrix(generated_data)
    v[fns] <- 0
    d <- as.data.frame(v)
    return(d)
}

make_data <- function(m) {
    data_cols <- nrow(m)
    data_rows <- 11*data_cols
    reporter_names <- make.unique(c(colnames(m), rep(colnames(m), each=10)))
    generated_data <- as.data.frame(matrix(0, nrow=data_rows, ncol=data_cols))
    colnames(generated_data) <- colnames(m)
    rownames(generated_data) <- reporter_names

    complexes <- list()
    for (node_idx in 1:length(rownames(m))) {
        complexes[[node_idx]] <- walk_tree(node_idx, m)
    }

    # assign values to the reporter genes
    for (path in unlist(complexes)) {
        affected <- c()
        path_array <- strsplit(path, "_")[[1]]
        for (node in as.numeric(path_array)) {
            node_name <- colnames(m)[[node]]
            affected <- c(affected, node_name)
            # paths climb up tree -- genes at top of path are affected by those underneath
            for (a in affected) {
                rows <- grep(paste(node_name, ".", sep=""), reporter_names)  # XXX only works with 26 nodes or fewer and no underscores in names
                generated_data[rows,a] <- 1
            }
        }
    }
    
    # assign values to the s-genes as reporters
    sgenes <- colnames(generated_data)
    for (sgene in sgenes) {
        generated_data[sgene, sgene] <- 1
    }

    # format of this object is perturbations in columns, reporters in rows
    # rows with same name as columns are looking at whether the knockdown is successful
    # n.# rownames are the egenes attached to the respective sgene

    return(generated_data)
}


walk_tree <- function(node_idx, m) {
    complexes <- wt_h(as.character(node_idx), m)
    return(complexes)
}

wt_h <- function(path, m) {
    x <- strsplit(path, "_")[[1]]
    last <- x[[length(x)]]
    first <- x[[1]]
    if (last == first && length(path) > 1) {
        print(paste0("cycle detected at ", paste(path, collapse=" ")))
        return(complexes)
    }
    newnodes <- as.integer(which(m[,as.numeric(first)]==1)) # find parents of this node
    complexes <- list() # track any complexes that branch from this node
    if (length(newnodes) > 0) {
        for (i in 1:length(newnodes)) {
            newpath <- paste(newnodes[[i]], path, sep="_")
            complexes[[i]] <- wt_h(newpath, m)
        }
    } else {
        complexes <- path
    }
    return(complexes)
}

