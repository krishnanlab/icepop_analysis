from icepop.data import HomologyData
import pandas as pd
import scanpy as sc

# load human–mouse orthology mapping
ortho_map = HomologyData().load()

# build mouse → human gene map (take first match)
m2g = {}
for hgene, mgene in ortho_map.items():
    m2g[mgene[0]] = hgene

# load single-cell dataset
adata = sc.read_h5ad("../data/mouse_colon/mouse_colon_cnt.h5ad")

# map mouse genes to human genes (NA if not found)
hgene = [m2g[mgene] if mgene in m2g else 'NA' for mgene in adata.var_names]

# create mapping table and remove unmapped genes
df = pd.DataFrame(zip(adata.var_names, hgene), columns=['mm', 'hs'])
df = df[df['hs'] != 'NA']

# save mapping for downstream analysis
df.to_csv('../data/mouse_colon/seismic/mm2hs.tsv', header=True, index=False, sep='\t')