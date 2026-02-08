import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import trange


adata = sc.read('../data/simulation/TM_subset_cnt.h5ad')
genes = adata.var_names.values
files = Path('../data/TM_FACS/magmaz').glob('*.genes.out')
files = [str(i) for i in files]
outdir = '../data/simulation/null_simulation'
Path(outdir).mkdir(exist_ok=True)

# randomizer
rng = np.random.default_rng(42)

for idx in trange(10000):
    file = rng.choice(files, 1)[0]
    score_df = pd.read_csv(file, header=0, index_col=None, sep=r'\s+')
    score_df['GENE'] = score_df['GENE'].astype(str)
    score_df = score_df.set_index('GENE')
    shared_genes = score_df.index.intersection(genes)
    score_df = score_df.loc[shared_genes, :].copy()
    score_df['ZSTAT'] = rng.permutation(score_df['ZSTAT'].values)
    score_df.loc[:, 'ZSTAT'].to_csv(f'{outdir}/{idx}.tsv', header=True, index=True, sep='\t')
