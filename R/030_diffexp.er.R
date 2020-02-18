#!/usr/bin/env Rscript

# pipeline
# select differentially expressed genes from the parsed bam lfc data

# functions

preprocess <- function(fc, experiment_definitions, expr.cutoff, samples) {
    dedup <- grep("bam\\.1", colnames(fc$counts)) # some filenames have been included twice in the list of files, remove them here
    x <- edgeR::DGEList(counts=fc$counts[,-dedup], genes=fc$annotation[,c("GeneID","Length")])

    id <- as.character(experiment_definitions$Bam.File)  
    id <- vapply(strsplit(id,"/"),"[",4, FUN.VALUE=character(1))
    id <- vapply(strsplit(id,"-"),"[",2, FUN.VALUE=character(1))

    # do we have an id for each column?
    names <- grep(paste(id,collapse = "|"),colnames(x$counts),value = T)
    stopifnot(all.equal(names,colnames(x$counts)[1:471])) # here something could potentially go wrong 
    x$counts <- x$counts[,which(colnames(x$counts) %in% names)]
    x$samples <- x$samples[which(rownames(x$samples) %in% names),]


    # Look up the gene names in the UCSC hg38 annotation
    #print("Are the rownames of the counts matrix in the same order as the gene IDs?")
    stopifnot(all.equal(as.integer(rownames(x$counts)),x$genes$GeneID))
    rownames(x$counts) <- x$genes$GeneID
    #print("Look Up gene names in UCSC hg38 annotation")
    require(org.Hs.eg.db, quietly = TRUE)  # sorry for this hack
    rownames(x$counts) <- annotate::lookUp(rownames(x$counts), 'org.Hs.eg', 'SYMBOL')       
    
    # rename colnames(x$counts)
    samples_all <- paste(experiment_definitions$Pool, experiment_definitions$Barcode, experiment_definitions$Cell.Line, experiment_definitions$Gene, experiment_definitions$shRNA, experiment_definitions$Biological.Rep)
    # don't overwrite samples from config
    colnames(x$counts) <- samples_all
  
    # print dimensions of the raw counts matrix
    #print("Dimensions of the raw counts matrix:")
    #dim(x)
    #print("Range of absolute unnormalized expression values [cpm]: ")
    #range(edgeR::cpm(x$counts))
    
    # filter out genes that don't vary by more than expr.cutoff in more than 2 experiments
    keep <- rowSums(edgeR::cpm(x) > expr.cutoff) >= 2
    x <- x[keep, , keep.lib.sizes=FALSE]

    # filtering out strange one "RARA 8C B" 
    #get.rid.off <- c("RARA 8C B")
    #x <- x[,-which(colnames(x$counts) %in% get.rid.off), keep.lib.sizes=TRUE]

    # don't try to reorder stuff
#    ctrlcols <- grep("CTRL", colnames(x$counts))
#    treatcols <- grep("CTRL", colnames(x$counts), invert=TRUE)
#    x$counts[,c(treatcols,ctrlcols)]
#
    controls <- grep("CTRL",colnames(x$counts),value = T)
#    index <- colnames(x$counts[,order(colnames(x$counts))])[-which(colnames(x$counts) %in% controls)]
#
#    index <- c(controls,samples)
#    sel <- x
#    sel <- sel[,which(colnames(sel$counts) %in% index)]

    sel <- x$counts
    return(list(x$counts, controls))
}
    
normalize_counts <- function(pos.norm) {
    # Normalize by dividing all samples INCLUDING pos.control by the pos.control
      controls <- grep("CTRL",colnames(pos.norm$counts),value = T)
      non_nas <- sel$counts[ -which(0 == rowMeans(sel$counts[,controls])), ]
      pos.norm$counts <- non_nas
      pos.norm$counts <- pos.norm$counts/rowMeans(pos.norm$counts[,controls], na.rm=TRUE)
      return(pos.norm)
      #return(pos.norm$counts)
}

process_with_edger <- function(sel, controls, selected.genes) {
    # based on: https://www.bioconductor.org/help/workflows/RnaSeqGeneEdgeRQL/
    data <- sel
      
    # adjusting design matrix (because we removed some samples above!)
    group <- vapply(strsplit(colnames(data$counts)," "),"[",1, FUN.VALUE=character(1))
    group <- factor(group)
    data$samples$group <- group
  
    # new design matrix for the conditions we want to contrast
    design.succ <- model.matrix(~0+group)
    colnames(design.succ) <- levels(group)
    # This design matrix simply links each group to the samples that belong to it. Each row of the design matrix        corresponds to a sample whereas each column represents a coefficient corresponding to one of the six groups

    # calculate normalization factors
    norm.factors <- edgeR::calcNormFactors(data, method = "TMM")
    
    # obtain normalized counts for downstream purposes
    edgeR.cpm <- edgeR::cpm(norm.factors, normalized.lib.sizes=F)
    
    # Dispersion estimation
    dispersion <- edgeR::estimateDisp(norm.factors, 
                               design.succ, 
                               robust=TRUE
                               )

    # NB model extended with quasi-likelihood estimation
    fit <- edgeR::glmQLFit(dispersion, 
                    design.succ, 
                    robust=TRUE
                    )

    # generate the vector of contrasted genes
    targets <- colnames(design.succ)[-which(colnames(design.succ)=="POS.CONTROL")]
    contrasted.genes <- list()
    for (i in 1:length(targets)) {
        contrasted.genes[[i]] <- paste0(targets[[i]],"vsPOS.CONTROL = ",targets[[i]],"-POS.CONTROL")
    }
    contrasted.genes <- as.vector(unlist(contrasted.genes))
    
    # selecting DGEs for downstream modeling
    # cutoff will be drawn by padj=0.05 
    
    # looping over every single condition to perform the tests (padj=0.05)
    contrast.list <- list()
    glmQLFTest.list <- list()
    for (i in 1:length(contrasted.genes)){
        contrast.list[[i]] <- edgeR::makeContrasts(contrasts=contrasted.genes[[i]], levels = design.succ)
    }
            
    for (i in 1:length(contrast.list)){    
        glmQLFTest.list[[i]] <- edgeR::glmQLFTest(fit, contrast=contrast.list[[i]]) # An object of class "DGELRT"
    }
    return.val <- list()
    for (i in 1:length(glmQLFTest.list)) {
        table <- glmQLFTest.list[[i]]$table
        return.val[[i]] <- table[,c(1,4)]
        colnames(return.val[[i]]) <- c("log2FoldChange", "padj")
    }
    names(return.val) <- selected.genes
    return(return.val)
}

reformat <- function(sel, controls) {
    counts <- sel$counts
    condition <- factor(vapply(strsplit(colnames(counts)," "),"[",1, FUN.VALUE=character(1)))
    experiments <- grep(paste0(condition[-which(condition %in% grep("CTRL",condition,value = T))], 
                               collapse = "|"),
                               colnames(counts),
                               value = T
                               )
    sort.idx <- c(controls,sort(experiments))
    counts <- counts[,sort.idx]
    
    # adjusting condition before generating dds since counts are sorted differently now
    condition <- factor(vapply(strsplit(colnames(counts)," "),"[",1, FUN.VALUE=character(1)))
    return(list(counts, condition))
}

process_with_deseq <- function(counts, design, selected.genes) {
    #coldata <- data.frame(row.names=colnames(counts), design)
    rownames(design) <- colnames(counts)
    design$shRNA <- factor(design$shRNA, levels=c("Renilla", levels(design$shRNA)[-which(levels(design$shRNA) == "Renilla")]))
        
    # actual DESeq2 analysis
    dds <- DESeq2::DESeqDataSetFromMatrix(countData=counts, 
                                  colData=design, 
                                  design=~ Cell.line + shRNA + Cell.line:shRNA
                                  )
    
    # just pulling out the normalizes counts (+devided by sizeFactor)
    dds <- DESeq2::estimateSizeFactors(dds)
    DESeq.norm.counts <- DESeq2::counts(dds, normalized=TRUE)
    
    # DGE analysis
    dds <- DESeq2::DESeq(dds)

    guides <- unique(design$shRNA)
    guides <- as.character(guides[ -which(guides == "Renilla") ])
    # to make pairwise comparisons
    # has to be done by looping 
    # over the contrasts (pairwise)
    res.list <- list()
    for(i in 1:length(guides)){
        res.list[[i]] <- DESeq2::results(dds, 
                        alpha = 0.05, # alpha refers to FDR cutoff
                        contrast = c("shRNA","Renilla",paste0(guides[[i]]))) 
    }
    names(res.list) <- guides

    guides <- gsub("-", ".", guides)
    res.list.t47d <- list()
    for (i in 7:length(guides)) {
        res.list.t47d[[i]] <- DESeq2::results(dds, list( c(paste("shRNA_", guides[[i]], "_vs_Renilla", sep=""), paste("Cell.lineT47D.shRNA", guides[[i]], sep=""))), alpha=0.05)
    }
    names(res.list.t47d) <- guides

    return(list(res.list, res.list.t47d))

    
    # not implemented for designs with interactions
    shrink.list <- list()
    for(i in 1:length(guides)){  
        shrink.list[[i]] <- DESeq2::lfcShrink(dds, 
                                  contrast = c("shRNA","Renilla",paste0(guides[i])),
                                  res = res.list[[i]] # this is important, because otherwise the settings above will be neglected - e.g. the p.adjust <0.05 setting
                                  )
                               }
    names(shrink.list) <- guides
    return(shrink.list)
}

sgene_heatmap <- function(diffexp, diffexp_method, input_file_name, selected.genes, diffexp_dir) {
    stat <- "log2FoldChange"
    d <- data.frame()
    for (sg in selected.genes) {
        d <- rbind(d, diffexp[[sg]][selected.genes,stat])
    }
    rownames(d) <- paste("ko", selected.genes)
    colnames(d) <- selected.genes
    ifn <- gsub("Rds", "pdf", input_file_name)
    output_file_name <- paste("sgene_heatmap", diffexp_method, ifn, sep=".")
    pdf(file.path(diffexp_dir, output_file_name))
    pheatmap::pheatmap(d, cluster_rows=FALSE, cluster_cols=FALSE, color=color, breaks=breaks)
    # rows are knockdowns
    # cols are lfc
}

guide_heatmap <- function(diffexp, diffexp_method, input_file_name, selected.genes, guides, cellline, diffexp_dir) {
    stat <- "log2FoldChange"
    selected.genes <- toupper(selected.genes)
    names(diffexp) <- gsub("-", ".", names(diffexp))
    d <- data.frame()
    for (g in guides) {
        d <- rbind(d, diffexp[[g]][selected.genes,stat])
    }
    colnames(d) <- selected.genes
    rownames(d) <- guides

    e <- experiment_definitions %>% filter(Gene != "CTRL") %>% group_by(Gene,shRNA) %>% summarize() %>% print(n=200)
    e$shRNA <- gsub("-", ".", e$shRNA)

    #group genes and guides
    f <- d[as.character(e$shRNA), sort(colnames(d))]

    # log2foldchange
    breaksdn <- seq(-1.5, -0.5, 0.5)
    breaksup <- seq(0.5, 3, 0.5)
    breaksdn <- seq(-3, -0.5, 0.5)
    breaksup <- seq(0.5, 2, 0.5)
    colordn <- rev(colorRampPalette(RColorBrewer::brewer.pal(n=7, name="Reds"))(length(breaksdn)))
    colorup <- colorRampPalette(RColorBrewer::brewer.pal(n=7, name="Blues"))(length(breaksup))

    ifn <- gsub("Rds", "pdf", input_file_name)
    output_file_name <- paste("guide_heatmap", cellline, diffexp_method, ifn, sep=".")
    pdf(file.path(diffexp_dir, output_file_name))
    pheatmap::pheatmap(f, cluster_rows=FALSE, cluster_cols=FALSE, color=c(colordn, "#FFFFFF", colorup), breaks=c(breaksdn, 0, breaksup))
    dev.off()
}

pval_heatmap <- function(diffexp, diffexp_method, input_file_name, selected.genes, guides, cellline, diffexp_dir) {
    # padj
    stat <- "padj"

    selected.genes <- toupper(selected.genes)
    names(diffexp) <- gsub("-", ".", names(diffexp))
    d <- data.frame()
    for (g in guides) {
        d <- rbind(d, diffexp[[g]][selected.genes,stat])
    }
    colnames(d) <- selected.genes
    rownames(d) <- guides

    e <- experiment_definitions %>% filter(Gene != "CTRL") %>% group_by(Gene,shRNA) %>% summarize() %>% print(n=200)
    e$shRNA <- gsub("-", ".", e$shRNA)

    #group genes and guides
    f <- d[as.character(e$shRNA), sort(colnames(d))]

    breaksup <- seq(2.5, 5, 0.5)
    breaksdn <- seq(0, 2, 0.5)
    colordn <- rep("#FFFFFF", length(breaksdn))
    colorup <- colorRampPalette(RColorBrewer::brewer.pal(n=7, name="Blues"))(length(breaksup))

    ifn <- gsub("Rds", "pdf", input_file_name)
    output_file_name <- paste("pval_heatmap", cellline, diffexp_method, ifn, sep=".")
    pdf(file.path(diffexp_dir, output_file_name))
    pheatmap::pheatmap(-log(f,10), cluster_rows=FALSE, cluster_cols=FALSE, color=c(colordn, colorup), breaks=c(breaksdn, breaksup))
    dev.off()
    # rows are knockdowns
    # cols are lfc
}





# main

step_030_diffexp <- function(project, aligner, diffexp_method, lfc_dir, diffexp_dir, experiment_definitions, expr.cutoff, samples, selected.genes) {
    print("030_diffexp")
    input_file_name <- paste(aligner, project, "Rds", sep=".")
    fc <- readRDS(file.path(lfc_dir, input_file_name))
    foo <- preprocess(fc, experiment_definitions, expr.cutoff, samples)
    sel <- foo[[1]] #x$counts
    controls <- foo[[2]]
    #nsel <- normalize_counts(sel)
    nsel <- sel

    # for multifactor design, pass dataframe as condition

    design <- experiment_definitions[,c(1,6,7,8)]


    #l <- reformat(nsel, controls)
    #counts <- l[[1]]
    #condition <- l[[2]]
    counts <- nsel
    if (diffexp_method == "DESeq") {
        # diffexp is now a list
        diffexp <- process_with_deseq(counts, design, selected.genes)
    } else if (diffexp_method == "edgeR") {
        #diffexp <- process_with_edger(sel, controls, selected.genes)
    }
    # generate sanity check heatmap
    #sgene_heatmap(diffexp, diffexp_method, input_file_name, selected.genes, diffexp_dir)
    cellline <- c("MCF7", "T47D")
    for (idx in 1:length(celllines)) {
        output_file_name <- paste(diffexp_method, cellline[[idx]], input_file_name, sep=".")
        saveRDS(diffexp[[idx]], file.path(diffexp_dir, output_file_name))
        guide_heatmap(diffexp[[idx]], diffexp_method, input_file_name, selected.genes, guides, cellline[[idx]], diffexp_dir)
        pval_heatmap(diffexp[[idx]], diffexp_method, input_file_name, selected.genes, guides, cellline[[idx]], diffexp_dir)
    }

    # save control data for bnem -- abs values, not lfc
    #ctrl <- as.data.frame(rowMeans(sel[!grepl("^NA[0-9]*$",rownames(sel$counts)),controls]$counts))
    #rownames(ctrl) <- rownames(sel[!grepl("^NA[0-9]*$",rownames(sel$counts)),]$counts)
    #colnames(ctrl) <- c("ctrl")
    #output_file_name <- paste("controls", input_file_name, sep=".")
    #saveRDS(ctrl, file.path(diffexp_dir, output_file_name))
}

