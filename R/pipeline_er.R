#!/usr/bin/env Rscript

# execute the analysis pipeline for the ER data
# each step will write intermediate output files to disk

library(NEMpipeline)

run_ER_pipeline <- function() {

    ################################################################################
    #
    # Definitions
    #
    ################################################################################

    # file locations
    base_input_dir <- "~/NEMpipelineBAM"

    # output locations
    base_output_dir <- "~/projects/NEMpipelineoutput"
    source_dir <- file.path(base_output_dir, "000_source")
    lfc_dir <- file.path(base_output_dir, "010_lfc")
    diffexp_dir <- file.path(base_output_dir, "030_diffexp")
    #e2_regulon_dir <- file.path(base_output_dir, "039_e2_regulon")
    prepared_dir <- file.path(base_output_dir, "040_prepared")
    heatmap_dir <- file.path(base_output_dir, "045_heatmap")
    nems_dir <- file.path(base_output_dir, "050_nems")
    consensus_dir <- file.path(base_output_dir, "060_consensus")
    plots_dir <- file.path(base_output_dir, "070_plots")
    benchmark_dir <- file.path(base_output_dir, "075_benchmark")
    egenes_dir <- file.path(base_output_dir, "080_egenes")
    peaks_dir <- file.path(base_output_dir, "090_chip")

    # and ensure they exist
    for (output_dir in c(base_output_dir, lfc_dir, diffexp_dir, prepared_dir, heatmap_dir, nems_dir, consensus_dir, plots_dir, benchmark_dir, egenes_dir, peaks_dir)) {
        if( ! dir.exists(output_dir)) {
            dir.create(output_dir)
        }
    }

    benchmark_file <- "timing.log"

    # project name and definitions from included files from included files
    project <- "er"
    data(ER_experiment_definitions, package="NEMpipeline")

    # regulon to limit differentially expressed genes to
    # gives 1591 genes
    data(ER_regulon, package="NEMpipeline")

    # read in the sample.table excel sheet specifying the experiment details
    #csv.file <- file.path(base_input_dir, project, projects_definition[[project]]$experiment)
    #if (file.exists(csv.file)) {
        #experiment_definitions <- read.csv(csv.file)
        #sampleTable <- experiment_definitions
    #}



    # aligners: hisat2, bowtie
    aligner <- "hisat2" 

    # diffexp_methods: DESeq, edgeR
    diffexp_method <- "DESeq" 

    # prep_methods: binary, lfc, pvalue, boolean, fgnem, decompose, progressive, bootstrap, random, geneset, cluster_bin
    # binary - discretize expression data to 0 or 1 (there is an effect or not)
    # lfc - use log fold change of differential expression
    # pvalue - use pvalue of differential expression
    # boolean - model inhibition and activation with bnem
    # fgnem - use factor graph nem
    # decompose - subset sgenes into a set of smaller nems, then evaluate with search to see if the models are compositional
    # progressive - take N most significant egenes, then 2N most significant, then 3N etc
    # bootstrap - randomly sample from the egenes
    # random - randomly generate data, to compare the bootstraps against
    # geneset - use gsea to group the egenes into genesets to run the nems on
    # cluster_bin - cluster binary expression data to generate ensemble egenes to run nems on
    prep_method <- "binary"

    # Which NEM methods should be applied?
    # The way the data is prepared determines which NEM methods are compatible with which prepared data
    #
    # nem_methods: search, triples, pairwise, BN.exhaustive, BN.greedy, ModuleNetwork.orig, ModuleNetwork, mc.eminem, depn, lem
        # search -- Markowetz 2005 "Non-transcriptional pathway features reconstructed from secondary effects of RNA interference"
        # triples -- Markowetz 2007 "Nested effects models for high-dimensional phenotyping screens"
        # pairwise -- Markowetz 2007 "Nested effects models for high-dimensional phenotyping screens"
        # BN.exhaustive -- Zeller 2008 "A Bayesian Network View on Nested Effects Models"
        # BN.greedy -- Zeller 2008 "A Bayesian Network View on Nested Effects Models"
        # ModuleNetwork.orig -- Frohlich 2007 "Large scale statistical inference of signaling pathways from RNAi and microarray data"
        # ModuleNetwork -- Frohlich 2008, "Estimating large-scale signaling networks through nested effect models with intervention effects from microarray data"
        # nem.greedyMAP -  Tresch 2008 "Structure Learning in Nested Effects Models"
        # nem.greedy -- H. Fröhlich, A. Tresch, T. Beißbarth, Nested Effects Models for Learning Signaling Networks from Perturbation Data, Biometrical Journal, 2(51):304 - 323, 2009
        # mc.eminem -- Niederberger 2012 "MC EMiNEM Maps the Interaction Landscape of the Mediator"
        # depn -- Frohlich 2009 "Deterministic Effects Propagation Networks for reconstructing protein signaling networks from multiple interventions"
        # lem -- Szczurek 2016 "Linear effects models of signaling pathways from combinatorial perturbation data"
    nem_method_compat <- list("binary" = c("greedy", "search", "nem.greedy", "triples", "pairwise", "ModuleNetwork", "ModuleNetwork.orig"), 
                         "cluster_bin" = c("greedy", "search", "nem.greedy", "triples", "pairwise", "ModuleNetwork", "ModuleNetwork.orig"), 
                         "lfc" = c("mnem", "bnem", "lem", "mnem", "mc.eminem"),
                         "pvalue" = c("mc.eminem", "bnem", "lem", "mnem", "mc.eminem"),
                         "boolean" = c("bnem"), 
                         "geneset" = c("mc.eminem"),
                         "bootstrap" = c("nem.greedy"),
                         "random" = c("nem.greedy"),
                         "decompose" = c("search"),
                         "progressive" = c("nem.greedy"),
                         "sc" = c("lem")
                         )

    nem_method <- "triples"
    nem_methods <- c("greedy", "nem.greedy", "triples", "pairwise", "ModuleNetwork")


    # successful screens
    samples <- c("EP300 12BR3 A", "EP300 12BR3 B", "EP300 12BR3 C", "ESRRA 2C A", "ESRRA 2C B", "ESRRA 2C C", "NCOA3 1DR3 A", "NCOA3 1DR3 B", "NCOA3 1DR3 C", "NCOA3 2BR A", "NCOA3 2BR B", "NCOA3 2BR C", "NR2F2 17AR3 A", "NR2F2 17AR3 B", "NR2F2 17AR3 C", "NRIP1 5C A", "NRIP1 5C B", "NRIP1 6A C", "RARA 7B A", "RARA 8C A", "RARA 8C C", "SUMO1 M1 A", "SUMO1 M1 B",     "SUMO1 M1 C",    "SUMO1 M35 A",    "SUMO1 M35 B",    "SUMO1 M35 C", "SUMO3 24AR3 A", "SUMO3 24AR3 B", "SUMO3 24AR3 C", "SUMO3 M16 A", "SUMO3 M16 C", "TRIM33 13CR3 A", "TRIM33 13CR3 B", "TRIM33 13CR3 C",     "TRIM33 M2 B", "TRIM33 M2 C", "ZMIZ1 19AR3 A", "ZMIZ1 19AR3 B", "ZMIZ1 19AR3 C", "ZMIZ1 22m35b A", "ZMIZ1 22m35b B", "ZMIZ1 22m35b C")

    selected.genes <- unique(vapply(strsplit(samples," "),"[",1, FUN.VALUE=character(1)))

    # 2 experiments must be above this lfc
    expr.cutoff <- 0.5


    # threshold for determining whether lfc is significant and there is an effect or not
    adjusted_pvalue_cutoff <- 0.05

    
    report_attached_egenes <- FALSE

    ################################################################################
    #
    # Pipeline
    #
    ################################################################################

    print("Running ER pipeline")
    print(paste("project: ", project))
    print(paste("aligner: ", aligner))

    # turn bam files into log fold change
    #step_010_lfc(project, aligner, ER_experiment_definitions, base_input_dir, lfc_dir)

    # calculate differential expression
    #step_030_diffexp(project, aligner, diffexp_method, lfc_dir, diffexp_dir, ER_experiment_definitions, expr.cutoff, samples, selected.genes)

    # prepare data in correct format to use with NEMs
    step_040_prepare_data(project, aligner, diffexp_method, prep_method, diffexp_dir, prepared_dir, adjusted_pvalue_cutoff, ER_regulon)

    # nems
    for (nem_method in nem_methods) {
        step_050_nems(project, aligner, diffexp_method, prep_method, nem_method, nem_method_compat, prepared_dir, nems_dir, egenes_dir, benchmark_file, report_attached_egenes, selected.genes)
    }
}

run_ER_pipeline()

