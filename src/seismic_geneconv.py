from icepop.data import HomologyData
import pandas as pd
import scanpy as sc

ortho_map = HomologyData().load()
m2g = {}
for hgene, mgene in ortho_map.items():
    m2g[mgene[0]] = hgene

adata = sc.read_h5ad("../data/mouse_colon/mouse_colon_cnt.h5ad")
hgene = [m2g[mgene] if mgene in m2g else 'NA' for mgene in adata.var_names]
df = pd.DataFrame(zip(adata.var_names, hgene), columns=['mm', 'hs'])
df = df[df['hs'] != 'NA']

df.to_csv('../data/mouse_colon/seismic/mm2hs.tsv', header=True, index=False, sep='\t')
