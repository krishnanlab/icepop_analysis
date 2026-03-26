import pandas as pd
from tqdm import tqdm
import pickle

# download gene annotation file (if not already downloaded)
# wget -O ../data/gene_info.gz https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz

sp = '10090'  # mouse taxonomy ID

# load and filter gene info for mouse
gene_info_df = pd.read_csv('../data/gene_info.gz', sep='\t', dtype=str)
gene_info_df = gene_info_df[gene_info_df['#tax_id'] == sp]

# build mapping: gene name / synonym → GeneID
name2id = {}

# add synonyms
for syn, id in tqdm(zip(gene_info_df['Synonyms'], gene_info_df['GeneID']), total=len(gene_info_df)):
    for i in syn.split('|'):
        name2id.setdefault(i, []).append(id)

# add official symbols (if not already included)
for sym, id in tqdm(zip(gene_info_df['Symbol'], gene_info_df['GeneID']), total=len(gene_info_df)):
    if sym not in name2id:
        name2id.setdefault(sym, []).append(id)

# save mapping for downstream use
with open('../data/mouse_name2id.pkl', 'wb') as f:
    pickle.dump(name2id, f)