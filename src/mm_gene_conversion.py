import pandas as pd
from tqdm import tqdm
import pickle

# download 
# wget -O ../data/gene_info.gz https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz

sp = '10090'
gene_info_df = pd.read_csv('../data/gene_info.gz', header=0, index_col=None, sep='\t', dtype=str)
gene_info_df = gene_info_df[gene_info_df['#tax_id'] == sp]

name2id = {}
for syn, id in tqdm(zip(gene_info_df['Synonyms'], gene_info_df['GeneID']), total=len(gene_info_df)):
    for i in syn.split('|'):
        name2id.setdefault(i, [])
        name2id[i].append(id)

for sym, id in tqdm(zip(gene_info_df['Symbol'], gene_info_df['GeneID']), total=len(gene_info_df)):
    if sym not in name2id:
        name2id.setdefault(sym, [])
        name2id[sym].append(id)

with open('../data/mouse_name2id.pkl', 'wb') as f:
    pickle.dump(name2id, f)
