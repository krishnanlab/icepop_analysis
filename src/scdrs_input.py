import pandas as pd
import numpy as np
import warnings
from glob import glob
from pathlib import Path
import scanpy as sc
from icepop.data import HomologyData

warnings.filterwarnings("ignore")


if __name__ == '__main__':
    top_n = 1000

    outdir = '../data/TM_FACS/scdrs'
    processed_adata = f"{outdir}/expr.h5ad"
    Path(outdir).mkdir(exist_ok=True)

    # load adata
    print('preprocess sc data')
    adata = sc.read_h5ad("../data/TM_FACS/TM_FACS_cnt.h5ad")

    # translate mouse to human genes
    # get human to model sp genes map
    ortho_map = HomologyData(sp='mmusculus').load()
    m2g = {}
    for hgene, mgene in ortho_map.items():
        m2g[mgene[0]] = hgene

    # keep mouse genes with human orthologs only
    adata.var['hgene'] = [m2g[mgene] if mgene in m2g else None for mgene in adata.var_names]
    adata = adata[:, ~pd.isna(adata.var['hgene'])].copy()
    adata.var_names = np.asarray(adata.var['hgene'])

    # save
    adata.write(processed_adata)

    # subset gene sets
    print('preprocess gene set')
    geneset_files = glob('../data/TM_FACS/magmaz/*.genes.out')
    geneset_dict = {}
    for file in geneset_files:
        name = Path(file).stem.replace('.genes.out', '')
        df = pd.read_csv(file, header=0, index_col=None, sep=r'\s+')
        df['GENE'] = df['GENE'].astype(str)
        df = df.sort_values(by='ZSTAT', ascending=False, ignore_index=True)
        df = df[df['GENE'].isin(adata.var_names)].copy()
        df = df.iloc[0:top_n, :]
        geneset_dict[name] = ','.join([f'{gene}:{score}' for gene, score in zip(df['GENE'], df['ZSTAT'])])
    df_gs = pd.DataFrame.from_dict(geneset_dict, orient='index')
    df_gs = df_gs.reset_index()
    df_gs.columns = ['TRAIT', 'GENESET']

    for trait in df_gs['TRAIT']:
        df_gs[df_gs['TRAIT'] == trait].to_csv(
            f'{outdir}/processed_{trait}_genesets.gs',
            header=True,
            index=False,
            sep="\t",
        )
