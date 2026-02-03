import scanpy as sc
from pathlib import Path
import pandas as pd
import anndata as ad
import numpy as np
from scipy.sparse import csr_matrix
from icepop.convert_score import CrossSpeciesScoreConverter


if __name__ == '__main__':
    min_cells = 10
    min_cell_in_gene = 10

    outdir = '../data/simulation_input_subset'
    Path(outdir).mkdir(exist_ok=True)

    # load data
    adata = sc.read('../data/sc_data/TM_FACS_MetaQ.h5ad')

    # randomizer
    rng = np.random.default_rng(1)

    # sample 1/10 of cells
    n = adata.n_obs
    k = int(n * 0.1)
    idx = rng.choice(n, size=k, replace=False)
    adata_sub = adata[idx, :].copy()

    # remove cell type with few cells
    counts = adata_sub.obs['cell_type'].value_counts()
    valid_types = counts[counts >= min_cells].index
    adata_sub = adata_sub[adata_sub.obs['cell_type'].isin(valid_types)].copy()

    # remove genes that don't have expression
    expr_counts = (adata_sub.X > 0.0).sum(axis=0)
    expr_counts = np.array(expr_counts).ravel()
    genes_keep = expr_counts >= min_cell_in_gene
    adata_sub = adata_sub[:, genes_keep].copy()

    # get count and convert across species
    score_converter = CrossSpeciesScoreConverter(adata_sub, sp='10090')
    score_converter.generate_cross_sp_matrix()
    exp = score_converter.convert_score_across_species(
        pd.DataFrame(
            adata_sub.layers['raw'].toarray(),
            index=adata_sub.obs_names,
            columns=adata_sub.var['entrez']),
        normed=False,
    )
    # make new adata
    X = csr_matrix(np.asarray(exp, dtype=np.float32))
    varnames = pd.DataFrame(zip(
        exp.columns,
        exp.columns),
        columns=['name', 'entrez']
    )
    varnames = varnames.set_index('name')
    new_adata = ad.AnnData(
        X=X,              # matrix (cells × genes)
        obs=adata_sub.obs,      # keep same metadata
        var=varnames
    )
    new_adata.obs = new_adata.obs.drop(columns=['SEACell'])
    new_adata.layers['raw'] = X.copy()

    # preprocess
    sc.pp.normalize_total(new_adata, target_sum=1e4)
    sc.pp.log1p(new_adata)
    sc.pp.highly_variable_genes(new_adata, n_top_genes=2000)
    sc.tl.pca(new_adata, n_comps=50, mask_var="highly_variable")
    new_adata.write(f'{outdir}/TM_subset.h5ad')
