library(seismicGWAS)

args <- commandArgs(trailingOnly = TRUE)
f <- args[1]
indir <- args[2]
outdir <- args[3]

# mkdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# get infile and outfile name
base <- sub("\\.genes.out$", "", basename(f))
outfile <- paste0(outdir, "/", base, ".tsv")

# load
sce <- readRDS(file.path(indir, "expr.rds"))

# calc spec score with minimal filtering
mmscore <- calc_specificity(
  sce,
  ct_label_col='cell_type',
  min_uniq_ct = 2,
  min_ct_size = 1,
  min_cells_gene_exp = 1,
  min_avg_exp_ct = 0.000001
)

# convert across species
gene_mapping_table <- read.table(
    file.path(indir, "mm2hs.tsv"),
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
)
sscore <- translate_gene_ids(mmscore, from='mm', to='hs', gene_mapping_table=gene_mapping_table)

# get association
magma_z <- read.table(f, header = TRUE, sep = "", stringsAsFactors = FALSE)
ra <- get_ct_trait_associations(sscore, magma_z)

# save
write.table(
    ra,
    file = outfile,
    sep = "\t",              # tab-separated
    quote = TRUE,           # no quotes around strings
    row.names = FALSE        # do not save row names
)

