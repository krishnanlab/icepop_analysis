import scanpy as sc
adata = sc.read('../data/TM_FACS/TM_FACS_cnt.h5ad')
adata.var['entrez'] = adata.var_names
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.tl.pca(adata, n_comps=50, use_highly_variable=True)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
adata.write('../data/TM_FACS/TM_FACS_normed.h5ad')

adata = sc.read('../data/mouse_colon/mouse_colon_cnt.h5ad')
adata.var['entrez'] = adata.var_names
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.tl.pca(adata, n_comps=50, use_highly_variable=True)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
adata.write('../data/mouse_colon/mouse_colon_normed.h5ad')

adata = sc.read('../data/mouse_ileum/mouse_ileum_cnt.h5ad')
adata.var['entrez'] = adata.var_names
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.tl.pca(adata, n_comps=50, use_highly_variable=True)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
adata.write('../data/mouse_ileum/mouse_ileum_normed.h5ad')
