suppressPackageStartupMessages(library(seismicGWAS))

argv <- commandArgs(trailingOnly = TRUE)
indir <- argv[1]
sc_file <- argv[2]
outdir <- argv[3]
setting <- basename(indir)

# load
sce <- readRDS(sc_file)

# mkdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# run seismic on all traits
files <- list.files(indir, pattern = "\\.tsv$", full.names = TRUE)

t_start <- Sys.time()
n_total <- length(files)

pb <- txtProgressBar(min = 0, max = length(files), style = 3)
for (i in seq_along(files)) {
    f <- files[i]
    base <- sub("\\.tsv$", "", basename(f))
    base <- sub(".*__", "", base)
    setting_name <- paste0(setting, "__run-", base)
    outfile = paste0(outdir, "/", setting_name, ".tsv")

    if (file.exists(outfile)) {
        next
    }

    # calc spec score
    sscore <- calc_specificity(
        sce,
        ct_label_col=paste0('run_', base),
        min_uniq_ct = 2,
        min_ct_size = 2,
        min_cells_gene_exp = 1,
        min_avg_exp_ct = 0.000001
    )

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

    # update progress bar
    setTxtProgressBar(pb, i)

    # print timing every 10 files (and at the end)
    if (i %% 10 == 0 || i == n_total) {
        elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
        avg_per_file <- elapsed / i
        remaining <- avg_per_file * (n_total - i)

        message(sprintf(
            "[%d/%d] elapsed: %.1f min | ETA: %.1f min",
            i, n_total, elapsed / 60, remaining / 60
        ))
    }
}
close(pb)
