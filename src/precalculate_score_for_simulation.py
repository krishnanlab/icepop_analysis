import scanpy as sc
import numpy as np
import pandas as pd
from icepop.specificity_score import specificity_score


if __name__ == '__main__':
    mcsize = 30
    infile = '../data/simulation/TM_subset_cnt.h5ad'
    mc_assign_infile = f'../data/simulation/mc-{mcsize}/mc_assign.csv'
    outdir = f'../data/simulation/mc-{mcsize}'

    # load data and metacell assignments
    adata = sc.read(infile)
    adata.obs['metacell'] = pd.read_csv(mc_assign_infile, header=None)[0].values

    # store raw counts
    adata.layers['raw'] = adata.X.copy()

    # standard preprocessing
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.tl.pca(adata, n_comps=50, use_highly_variable=True)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    # compute metacell specificity scores
    spec = specificity_score(adata, n_jobs=20)
    if "spec_score" not in adata.uns:
        spec.get_metacell_spec_score()
        np.savez_compressed(
            f'{outdir}/mc_spec_score.npz',
            score=np.asarray(adata.uns['spec_score'], np.float32),
            mc=adata.uns['spec_score'].index.values,
            genes=adata.uns['spec_score'].columns.values,
        )
    del adata.uns['spec_score']

    # add required gene identifier column
    adata.var['entrez'] = adata.var_names

    # save processed object
    adata.write(f'{outdir}/TM_subset__score_calc.h5ad')
