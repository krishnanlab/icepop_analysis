suppressPackageStartupMessages({
    library(zellkonverter)
    library(SingleCellExperiment)
})

args <- commandArgs(trailingOnly = TRUE)

infile  <- args[1]
outdir  <- args[2]

cat("Processing:", infile, "\n")

sce <- suppressWarnings(readH5AD(infile))

# set gene names
rownames(sce) <- rowData(sce)$entrez

# set logcounts
assay(sce, "logcounts") <- assay(sce, "X")

# output file
fname <- basename(infile)
fname <- sub("\\.h5ad$", ".rds", fname)
outfile <- file.path(outdir, fname)

saveRDS(sce, file = outfile)

cat("Saved:", outfile, "\n")