library(rtracklayer)
library(tidyverse)
#https://www.biostars.org/p/464924/

fix_names <- function(nms) {
  w0 <- gsub("[\t\n\r\v\f\a\b ]","",x=gsub("(?<=[\t\n\r\v\f\a\b ]).+","",x = nms, perl = T))
  w <- gsub(".+\\|","", x=w0)
  w1 <- gsub("#.+","",x = w)
  w1
}

te_fa <- import(snakemake@input$fasta, format = "fasta")

names(te_fa) <- names(te_fa) %>% fix_names()

Biostrings::writeXStringSet(te_fa, snakemake@output$fasta, format = "fasta")