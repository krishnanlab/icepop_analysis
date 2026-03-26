import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import trange

# load gene list from dataset
adata = sc.read('../data/simulation/TM_subset_cnt.h5ad')
genes = adata.var_names.values

# collect MAGMA result files
files = Path('../data/TM_FACS/magmaz').glob('*.genes.out')
files = [str(i) for i in files]

# output directory for null simulations
outdir = '../data/simulation/null_simulation'
Path(outdir).mkdir(exist_ok=True)

# random generator for reproducibility
rng = np.random.default_rng(42)

# generate null datasets by permuting gene scores
for idx in trange(10000):
    file = rng.choice(files, 1)[0]

    # load MAGMA gene-level z-scores
    score_df = pd.read_csv(file, sep=r'\s+')
    score_df['GENE'] = score_df['GENE'].astype(str)
    score_df = score_df.set_index('GENE')

    # restrict to genes present in dataset
    shared_genes = score_df.index.intersection(genes)
    score_df = score_df.loc[shared_genes, :].copy()

    # permute z-scores to create null
    score_df['ZSTAT'] = rng.permutation(score_df['ZSTAT'].values)

    # save null dataset
    score_df['ZSTAT'].to_csv(f'{outdir}/{idx}.tsv', sep='\t')