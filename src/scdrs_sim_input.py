import pandas as pd
import numpy as np
import warnings
from pathlib import Path
import scanpy as sc
from tqdm import tqdm
warnings.filterwarnings("ignore")


if __name__ == '__main__':
    top_n = 1000

    outdir = '../data/simulation/scdrs'
    causal_outdir = f'{outdir}/causal'
    null_outdir = f'{outdir}/null'
    processed_adata = f"{outdir}/expr.h5ad"
    Path(causal_outdir).mkdir(exist_ok=True, parents=True)
    Path(null_outdir).mkdir(exist_ok=True, parents=True)

    # load adata
    print('preprocess sc data')
    adata = sc.read_h5ad("../data/simulation/TM_subset_cnt.h5ad")

    # # save
    # adata.write(processed_adata)

    ####################################
    # Casusal
    ####################################
    # subset gene sets
    geneset_files = Path('../data/simulation/causal_simulation').rglob('gwasz__*.tsv')
    geneset_files = [str(geneset_file) for geneset_file in geneset_files]
    geneset_dict = {}
    for file in tqdm(geneset_files, total=len(geneset_files)):
        setting = Path(file).parent.name
        runid = Path(file).stem.split('__')[1]
        df = pd.read_csv(file, header=0, index_col=None, sep='\t')
        df = df.sort_values(by='ZSTAT', ascending=False, ignore_index=True)
        df['GENE'] = df['GENE'].astype(str)
        df = df[df['GENE'].isin(adata.var_names)].copy()
        df = df.iloc[0:top_n, :]
        geneset_dict[f'{setting}__run-{runid}'] = ','.join([f'{gene}:{score}' for gene, score in zip(df['GENE'], df['ZSTAT'])])
    df_gs = pd.DataFrame.from_dict(geneset_dict, orient='index')
    df_gs = df_gs.reset_index()
    df_gs.columns = ['TRAIT', 'GENESET']

    # output
    for idx, df in enumerate(np.array_split(df_gs, 600)):
        df.to_csv(f'{causal_outdir}/geneset_{idx}.gs', header=True, index=False, sep='\t')

    ####################################
    # Null
    ####################################
    geneset_files = Path('../data/simulation/null_simulation').glob('*.tsv')
    geneset_files = [str(i) for i in geneset_files]
    geneset_dict = {}
    for file in tqdm(geneset_files, total=len(geneset_files)):
        runid = Path(file).stem
        df = pd.read_csv(file, header=0, index_col=None, sep='\t')
        df = df.sort_values(by='ZSTAT', ascending=False, ignore_index=True)
        df['GENE'] = df['GENE'].astype(str)
        df = df[df['GENE'].isin(adata.var_names)].copy()
        df = df.iloc[0:top_n, :]
        geneset_dict[f'run-{runid}'] = ','.join([f'{gene}:{score}' for gene, score in zip(df['GENE'], df['ZSTAT'])])
    df_gs = pd.DataFrame.from_dict(geneset_dict, orient='index')
    df_gs = df_gs.reset_index()
    df_gs.columns = ['TRAIT', 'GENESET']

    # output
    for idx, df in enumerate(np.array_split(df_gs, 200)):
        df.to_csv(f'{null_outdir}/geneset_{idx}.gs', header=True, index=False, sep='\t')
