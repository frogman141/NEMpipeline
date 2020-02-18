# pipeline
# do nem calculations



attach_egenes <- function(nem_model, egenes) {
    # we are going to get $graph, an adjacency matrix and $mappos, a list of lists of egenes
    # and want to turn this into a single matrix

    m <- as(nem_model$graph, "matrix")
    sgenes <- rownames(m)
    egenes <- unique(unlist(nem_model$mappos))
    genes <- append(sgenes, egenes)
    
    bigger <- matrix(0L, nrow = length(genes), ncol= length(genes))
    colnames(bigger) <- genes
    rownames(bigger) <- genes
    bigger[sgenes, sgenes] <- m[sgenes, sgenes]
    
    for (sg in sgenes) {
        at <- nem_model$mappos[[sg]]
        bigger[sg, at] <- 1
        bigger[at, sg] <- 1
    }
    nem_model$graph <- bigger
    return(nem_model)
}


# functions
run_nems <- function(nem_method, expr_data, prepared_dir, nems_dir, egenes_dir, selected.genes, nem_method_compat) {
    print("starting running nems\n")
    start_time <- Sys.time()
    if (nem_method == "bnem") {
        if (requireNamespace("bnem")) {
            # requires data prepared by nem_method boolean
            NEMlist <- list()
            NEMlist$exprs <- NULL
            NEMlist$fc <- expr_data

            if (project == "fly") {
                stimuli <- c("LPS")
            } else {
                stimuli <- c("Ctrl")
            }
            inhibitors <- selected.genes
            CNOlist <- CellNOptR::dummyCNOlist(stimuli = stimuli, inhibitors = inhibitors, maxStim = 2, maxInhibit = 1)

            Sgenes <- c(stimuli, inhibitors)

            sifMatrix <- numeric()
            for (i in Sgenes) {
              for (j in Sgenes) {
                if (i %in% j) { next() }
                sifMatrix <- rbind(sifMatrix, c(i, 1, j))
                sifMatrix <- rbind(sifMatrix, c(i, -1, j)) # if you want negative edges
              }
            }
            write.table(sifMatrix, file = file.path(nems_dir, "temp.sif"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
            PKN <- readSIF(file.path(nems_dir, "temp.sif"))

            model <- bnem::preprocessing(CNOlist, PKN, maxInputsPerGate=1)

            initBstring <- rep(0, length(model$reacID)) # start with the empty network
            parallel <- 8
            output_file_name <- paste(nem_method, input_file_name, "pdf", sep=".")
            print(output_file_name)
            pdf(output_file_name)
            locRun <- bnem::localSearch(
                CNOlist=CNOlist,
                NEMlist=NEMlist,
                model=model,
                parallel=parallel,
                initSeed=initBstring,
                draw = TRUE # FALSE does not draw the network evolution and can be faster
            )
            dev.off()
            resString <- locRun$bStrings[1, ]
            end_time <- Sys.time()
            print("Warning: output not written to file")
            return(locRun)

        } else {
            print("To install bnem:")
            print("devtools::install_github('MartinFXP/B-NEM')")
        }
    } else if (nem_method == "lem") {
        if (requireNamespace("lem")) {
            print("nem_method: lem")
            if (project == "ps") {
                bin_exp_mat <- readRDS(file.path(prepared_dir, paste("bin_exp_mat", project, "Rds", sep=".")))
            } else {
                bin_exp_mat <- diag(length(selected.genes))
            }

            nocontrols <- FALSE
            if (nocontrols) {
                bin_exp_mat <- bin_exp_mat[,5:14]
            }

            # ensure that the matricies can be multiplied
            common <- intersect(rownames(bin_exp_mat), colnames(expr_data))
            bin_exp_mat <- bin_exp_mat[common,]
            expr_data <- expr_data[,common]
            lem_result <- lem::lem(t(as.matrix(expr_data)), as.matrix(bin_exp_mat), inference="greedy", parameter.estimation="linear.reg", verbose=TRUE)

            # extract which sgenes control which egenes (many:many)
            beta <- t(lem_string$beta)
            rownames(beta) <- rownames(expr_data)
            #table(apply(beta, 1, function(x) {order(x)})[9,])
            beta[beta < 0] <- 0
            n <- list()
            for (r in 1:nrow(beta)) {
                row <- beta[r,]
                ord <- row[rev(order(row))]
                egene <- rownames(beta)[r]
                n[[egene]] <- ord[ord > 0]
            }
            # write the egene attachments before we lose data in changing the format
            # TODO: HACK this violates the one script one output dir paradigm
            output_file_name <- paste("egenes", nem_method, project, "Rds", sep=".")
            sink(file.path(egenes_dir, output_file_name))
            print(n)
            sink()

            # format these results the same as the results from the nem package
            r <- list()
            r$graph <- l$graph
            # note: an E-gene can be attached to multiple S-genes (linear effects model...)
            mappos <- list()
            sgene_idx <- 1
            for (sg in colnames(expr_data)) {
                attached_genes <- rownames(expr_data)[which(l$beta[sgene_idx,]>0)]
                mappos[[sg]] <- attached_genes
                sgene_idx <- sgene_idx + 1
            }
            r$mappos <- mappos
            return(r)
        } else {
            print("To install lem:")
            print("install.packages('https://www.mimuw.edu.pl/~szczurek/lem/lem_1.0.tar.gz',repos=NULL, type='source')")
        }
    } else if (nem_method == "mnem") {
        if (requireNamespace("mnem")) {
            k <- 3
            result <- mnem::mnem(expr_data, k = k, starts = 10, search="greedy") # could do this with nem_method="disc" for discrete data
            ret_list <- list()
            for (k_idx in 1:k) {
                # format these results the same as the results from the nem package
                r <- list()
                sgenes <- colnames(result$data)
                egenes <- rownames(result$data)
                r$graph <- result$comp[[1]]$phi
                colnames(r$graph) <- unique(sgenes)
                rownames(r$graph) <- unique(sgenes)
                attachments <- result$comp[[1]]$theta
                mappos <- list()
                sgene_idx <- 1
                for (sg in sgenes) {
                    attached_egenes <- egenes[which(attachments == sgene_idx)]
                    mappos[[sg]] <- attached_egenes
                    sgene_idx <- sgene_idx + 1
                }
                r$mappos <- mappos
                ret_list[[k_idx]] <- r
                # TODO? return weighted probabilities to be incorporated into the file name
                # TODO: guarantee that all egenes are in $mappos
            }
            return(ret_list)
        }
    } else if (nem_method == "pcnem") {
        if (requireNamespace("pcnem")) {
            nem_method <- "AdaSimAnneal"
            type <- "mLL"
            control <- pcnem::set.default.parameters(selected.genes, type="mLL", pcombi=TRUE, trans.close=FALSE)
            control$map <- as.matrix(filtered) # works with binary data??
            b <- pcnem::nem(as.matrix(filtered), inference=nem_method, control=control, verbose=TRUE)
        } else {
            print("To install pcnem:")
            print("devtools::install_github('cbg-ethz/pcNEM')")
        }
    } else if (nem_method == "depn") {
        #not run
        data(SahinRNAi2008)
        control = set.default.parameters(setdiff(colnames(dat.normalized),"time"), map=map.int2node, type="depn",debug=FALSE) # set mapping of interventions to perturbed nodes
        net = nem(dat.normalized, control=control) # greedy hillclimber to find most probable network structure
    } else if (nem_method == "fgnem") {
        if (requireNamespace("fgnem")) {
            # now comes from prepared
            if (project == "fly") {
                input_file_name <- file.path(prepared_dir, paste(prep_method, project, "tsv", sep="."))
            } else {
                input_file_name <- file.path(prepared_dir, paste(prep_method, diffexp_method, aligner, project, "Rds", sep="."))
            }
            eg <- read.egene.tab(input_file_name)
            paramMean <- RCommandArgDouble("MEAN", default=1.5, gt=0, errorMsg=usage)
            paramSD <- RCommandArgDouble("SD", default=1, gt=0)
            params <- paramGen(paramMean, paramSD)
            results <- scoreBestModelEstimate(eg, params=params, doTransitivity=FALSE, summarization=marginop)
        } else {
            print("To install fgnem:")
            print("Follow instructions at https://sysbio.soe.ucsc.edu/projects/fgnem/")
        }
    } else {
        contr <- c(0.15,0.05)
        if (nem_method %in% unique(c(nem_method_compat[['binary']], "search", "nem.greedy", "triples", "pairwise", "ModuleNetwork", "ModuleNetwork.orig"))) {
            type <- "mLL" 
        } else {
            type <- "CONTmLLBayes"
        }
        hyper <- nem::set.default.parameters(selected.genes, 
                                   para=contr, 
                                   type = type)
        b <- nem::nem.bootstrap(as.matrix(expr_data),
            inference=nem_method,
            control=hyper, 
            verbose=F,
            nboot=1 # 
        )
        return(b)
    }
}


step_050_nems <- function( project, aligner, diffexp_method, prep_method, nem_method, nem_method_compat, prepared_dir, nems_dir, egenes_dir, benchmark_file, report_attached_egenes, selected.genes) {
    # main
    # read all expression data that has been written into data dir, and calculate nems for that, by each method
    timing <- data.frame(input=character(0), nem_method=character(0), seconds=numeric(0), date=character(0), stringsAsFactors=FALSE)
    matching <- dir(prepared_dir, pattern=project)
    matching <- Filter(function(x) grepl(paste("\\.", "Rds", sep=""), x), matching)
    for (input_file_name in matching) {
        #TODO: loading RDS needs to not happen for fgnem
        expr_data <- readRDS(file.path(prepared_dir, input_file_name))
        if (nem_method %in% nem_method_compat[[prep_method]]) {
            print(paste0("nem method: ", nem_method))
            tryCatch ({
                start_time <- Sys.time()
                nem_model <- run_nems(nem_method, expr_data, prepared_dir, nems_dir, egenes_dir, selected.genes, nem_method_compat)
                end_time <- Sys.time()
                # save models
                if (! ("graph" %in% names(nem_model))) {
                    # for models that return multiple networks, save them all
                    idx <- 1
                    for (nm in nem_model) {
                        if (report_attached_egenes) {
                            nem_model <- attach_egenes(nem_model)
                            output_file_name <- paste(nem_method, "with_egenes", input_file_name, idx, sep=".")
                        }
                        output_file_name <- paste(nem_method, input_file_name, idx, sep=".")
                        saveRDS(nm, file.path(nems_dir, output_file_name))
                        idx <- idx + 1
                    }
                } else {
                    if (report_attached_egenes) {
                        nem_model <- attach_egenes(nem_model)
                        output_file_name <- paste(nem_method, "with_egenes", input_file_name, sep=".")
                    }
                    output_file_name <- paste(nem_method, input_file_name, sep=".")
                    saveRDS(nem_model, file.path(nems_dir, output_file_name))
                }
                # save timing info
                timing[nrow(timing)+1,] <- c(input_file_name, nem_method, end_time - start_time, as.character(as.POSIXct(end_time, format='%Y-%m-%d-%h:%m:%s')))
            }, warning = function(w) {
              message(w)
            }, error = function(e) {
              message(e)
            })
            # TODO: make sure all errors are reported
        }
    }
    output_file_name <- file.path(nems_dir, benchmark_file)
    if (file.exists(output_file_name)) {
        write.table(timing, output_file_name, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE, append=TRUE)
    } else {
        write.table(timing, output_file_name, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
    }
}

