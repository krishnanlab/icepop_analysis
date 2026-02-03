library(zellkonverter)
library(SingleCellExperiment)

# sce <- readH5AD("../data/TM_FACS/TM_FACS_normed.h5ad")
# rownames(sce) <- rowData(sce)$entrez
# assay(sce, "logcounts") <- assay(sce, "X")
# saveRDS(sce, file = "../data/TM_FACS/seismic/expr.rds")

sce <- readH5AD("../data/mouse_colon/mouse_colon_normed.h5ad")
rownames(sce) <- rowData(sce)$entrez
assay(sce, "logcounts") <- assay(sce, "X")
saveRDS(sce, file = "../data/mouse_colon/seismic/expr.rds")
