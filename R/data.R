
#' @name ER_experiment_definitions
#' @title defines ER experiments
#' @format A data frame with variables:
#' Pool,Barcode,Sequence,Sample name,Bam File,Gene,shRNA,Biological Rep
"ER_experiment_definitions"


#' @name ER_regulon
#' @title filter for genes that are under the regulation of estrogen
#' @format A data frame with variables:
#' Affy Probe Id,Gene,Ensembl gene ID,gene.index,RP/Rsum,FC:(class1/class2),pfp,P.value
"ER_regulon"


#' @name ER_lfc
#' @title lfc data for ER experiments
#' @format Output from featureCounts()
#' counts, annotation, targets, stat
"ER_lfc"

