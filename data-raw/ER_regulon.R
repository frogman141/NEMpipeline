
infile <- "Jagannathan&Robinson-Rechavi2011.csv"
if (file.exists(infile)) {
    meta <- read.csv(infile)
    ER_regulon <- unique(meta$Gene)
}
usethis::use_data(ER_regulon)
