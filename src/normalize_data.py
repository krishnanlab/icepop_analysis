import scanpy as sc

# =============================
# Tabula Muris (TM_FACS)
# =============================

# load data
adata = sc.read('../data/TM_FACS/TM_FACS_cnt.h5ad')

# add gene identifier column
adata.var['entrez'] = adata.var_names

# standard preprocessing
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.tl.pca(adata, n_comps=50, use_highly_variable=True)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# save processed data
adata.write('../data/TM_FACS/TM_FACS_normed.h5ad')


# =============================
# Mouse colon dataset
# =============================

# load data
adata = sc.read('../data/mouse_colon/mouse_colon_cnt.h5ad')

# add gene identifier column
adata.var['entrez'] = adata.var_names

# standard preprocessing
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.tl.pca(adata, n_comps=50, use_highly_variable=True)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# save processed data
adata.write('../data/mouse_colon/mouse_colon_normed.h5ad')