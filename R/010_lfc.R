#!/usr/bin/env Rscript

# pipeline
# go from aligned bam files to log fold change counts


#' Compute log fold change from bam files
#'
#' Output of this function is also stored in data/hisat2.project.rda
#'
#' @param project
#' @param aligner
#' @param experiment_definitions
#' @param base_input_dir
#' @return nothing. writes output to file
#' @export
step_010_lfc <- function(project, aligner, experiment_definitions, base_input_dir, lfc_dir) {
    print("010_lfc")
    # set up the input/output locations

    # loading bam.files
    filenames <- as.character(experiment_definitions$Bam.File)
    filenames <- paste(base_input_dir, filenames, sep="/")
    if (all(file.exists(filenames))) {
        bamfiles <- Rsamtools::BamFileList(filenames)
        seqinfo(bamfiles[1])
              
        # counting reads
        ER_lfc <- Rsubread::featureCounts(filenames,annot.inbuilt="hg38",isPairedEnd=FALSE)

        # store data if we have regenerated it
        #usethis::use_data(ER_lfc) # store in package
        saveRDS(ER_lfc, file.path(lfc_dir, paste(aligner, project, "Rds", sep=".")))
    }
    else {
        print("using cached file. To run these steps, Download BAM files from ... and put in the dir specified by base_input_dir")
        # else get data from pregenerated data
        data(ER_lfc, package="NEMpipeline")
        # and save to local data
        saveRDS(ER_lfc, file.path(lfc_dir, paste(aligner, project, "Rds", sep=".")))
    }
}
