library(rtracklayer)

#gtffile <- "results/flybase-anno/transcriptome.gtf.gz"
gtffile <- snakemake@input[["host_gtf"]]

tx <- import(gtffile)

#te_gtf_file <- "results/tidal-anno/Tidalbase_transposon_sequence.dmel.fixed.sensechecked.consensus.gtf.gz"
te_gtf_file <- snakemake@input$te_gtf
te_gtf <- import(te_gtf_file)

te_gtf$type <- "mRNA"
te_gtf$gene_symbol <- te_gtf$gene_id
te_gtf$transcript_symbol <- te_gtf$transcript_id

# make final combined gr
combined <- c(tx, te_gtf)

combined$`#` <- NULL

# write to disk
export(combined,snakemake@output[["gtf"]])