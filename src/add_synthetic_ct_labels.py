import scanpy as sc
import numpy as np
import pickle
from pathlib import Path
from tqdm import tqdm

adata = sc.read('../data/simulation/TM_subset_cnt.h5ad')
sim_dir = '../data/simulation/causal_simulation_synthetic'
adata_outdir = '../data/simulation/TM_synthetic_adata'
adata_normed_outdir = '../data/simulation/TM_synthetic_adata_normed'
Path(adata_outdir).mkdir(exist_ok=True)
Path(adata_normed_outdir).mkdir(exist_ok=True)

# synthetic cell type name
label = 'synthetic_cell_type'

cts = np.asarray(list(adata.obs['cell_type']))

# get norm data
adata_norm = adata.copy()
sc.pp.normalize_total(adata_norm, target_sum=1e4)
sc.pp.log1p(adata_norm)
adata_norm.var['entrez'] = adata_norm.var_names

# load labels into adata
tc_files = Path(sim_dir).rglob('tc.pkl')
tc_files = [str(tc_file) for tc_file in tc_files]
for tc_file in tqdm(tc_files, total=len(tc_files)):
    adata_tmp = adata.copy()
    adata_norm_tmp = adata_norm.copy()

    setting = Path(tc_file).parent.name
    with open(tc_file, 'rb') as f:
        tc = pickle.load(f)

    for i in range(len(tc)):
        colname = f'run_{i}'
        active_cells = tc[i]['active_cells']
        inactive_cells = tc[i]['inactive_cells']
        all_cells = np.concatenate([active_cells, inactive_cells])
        tmps = cts.copy()
        tmps[adata_tmp.obs_names.get_indexer(all_cells)] = label
        adata_tmp.obs[colname] = tmps
        adata_norm_tmp.obs[colname] = tmps

    # save
    outfile = f'{adata_outdir}/{setting}.h5ad'
    adata_tmp.write(outfile)

    # normalize for seismic input
    normed_outfile = f'{adata_normed_outdir}/{setting}.h5ad'
    adata_norm_tmp.write(normed_outfile)
