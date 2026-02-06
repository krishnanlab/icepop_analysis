import scanpy as sc
import numpy as np
import pandas as pd
from icepop.specificity_score import specificity_score


def get_centroids(adata):
    """
    get metacell centroids
    """
    # metacell-level PCA
    X = adata.obsm["X_pca"]  # shape N × D
    mcs = adata.obs['metacell'].values  # length N

    # unique metacells
    unique_mcs = np.array(sorted(np.unique(mcs)))
    M = len(unique_mcs)

    # centroid PCA for each metacell
    centroids = np.zeros((M, X.shape[1]))
    for i, mc in enumerate(unique_mcs):
        idx = np.where(mcs == mc)[0]
        centroids[i] = X[idx].mean(axis=0)
    adata.uns['centroids'] = pd.DataFrame(
        centroids,
        index=unique_mcs,
        columns=[f'PC{i}' for i in np.arange(centroids.shape[1])]
    )


if __name__ == '__main__':
    indir = '../data/simulation'
    outdir = '../data/simulation/mc-75'

    # load data
    adata = sc.read(f'{indir}/TM_subset_cnt.h5ad')
    adata.obs['metacell'] = pd.read_csv(f'{outdir}/mc_assign.csv', header=None, index_col=None)[0].values

    # save raw
    adata.layers['raw'] = adata.X.copy()

    # normalize
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.tl.pca(adata, n_comps=50, use_highly_variable=True)

    # get spec score
    spec = specificity_score(adata, n_jobs=20)
    if "spec_score" not in adata.uns.keys():
        spec.get_metacell_spec_score()
        np.savez_compressed(
            f'{outdir}/mc_spec_score.npz',
            score=np.asarray(adata.uns['spec_score'], np.float32),
            mc=adata.uns['spec_score'].index.values,
            genes=adata.uns['spec_score'].columns.values,
        )
        del adata.uns['spec_score']

    # get centroids of metacells
    get_centroids(adata)

    # save
    adata.write(f'{outdir}/TM_subset__score_calc.h5ad')
