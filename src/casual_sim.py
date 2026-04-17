import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
from time import time
from multiprocessing import Pool
from multiprocessing.shared_memory import SharedMemory
from scipy.sparse import csr_matrix
import pickle
import argparse

import warnings
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)


# ------------------------------------------
# simulate gwas signal from spec score
# ------------------------------------------
def choose_causal_genes(spec, signal_frac=0.2, rng=None):
    if rng is None:
        rng = np.random.default_rng()
    G = len(spec)
    spec = np.asarray(spec)

    # --- choose signal genes (π1 = signal_frac) ---
    n_signal = int(np.floor(signal_frac * G))
    if n_signal < 1:
        raise ValueError("signal_frac too small, no signal genes selected.")

    weights = np.log1p(np.exp(spec))
    # weights = np.where(spec <= 0, 0, spec)
    weights = weights / weights.sum()

    signal_idx = rng.choice(G, size=n_signal, replace=False, p=weights)

    is_causal = np.zeros(G, dtype=bool)
    is_causal[signal_idx] = True
    return is_causal


def simulate_latent_from_specificity(
    spec: np.ndarray,
    is_causal: np.ndarray,
    beta: float = 1.0,
    noise_sd: float = 1.0,
    rng: np.random.Generator | None = None,
):
    """
    spec: 1D array of specificity scores (higher = more cell-type specific)
    is_causal: boolean array for which genes are 'signal'
    beta: strength of coupling between specificity and latent score
    noise_sd: how noisy the mapping is (higher = weaker separation)
    """
    if rng is None:
        rng = np.random.default_rng()

    noise = rng.normal(loc=0.0, scale=noise_sd, size=len(spec))

    y = np.zeros(len(spec))
    # causal genes: spec + noise
    y[is_causal] = beta * spec[is_causal] + noise[is_causal]
    # noncausal genes: pure noise
    y[~is_causal] = noise[~is_causal]

    return y


def assign_z_from_latent(y: np.ndarray, z_pool: np.ndarray) -> np.ndarray:
    """
    Rank-match GWAS-like z-scores to latent 'importance' scores.

    - Sort y ascending
    - Sort z_pool ascending
    - Assign largest z to largest y, etc.
    """
    assert y.shape[0] == z_pool.shape[0]
    # get ranks of y (0 .. G-1), ties broken arbitrarily
    y_rank = np.argsort(np.argsort(y))

    z_sorted = np.sort(z_pool)  # real z distribution
    z_sim = z_sorted[y_rank]  # map by rank

    return z_sim


def simulate_gwas_from_specificity_celltype(
    spec: pd.DataFrame,
    z_pool: pd.DataFrame,
    signal_frac: float = 0.1,
    beta: float = 1.0,
    noise_sd: float = 1.0,
    rng: np.random.Generator | None = None,
):
    """
    Simulate GWAS-like z-scores that:
      - follow the real z_pool distribution
      - are positively coupled to specificity for a subset of genes
    """
    if rng is None:
        rng = np.random.default_rng()

    # --- align indices ---
    genes = z_pool.index.intersection(spec.index)
    z_pool = z_pool.loc[genes].copy()
    spec = spec.loc[genes].copy()
    G = np.asarray(genes)

    # pick causal genes at cell type level (major trend of specificity of a cell type)
    is_causal = choose_causal_genes(spec, signal_frac=signal_frac, rng=rng)

    # turn to np array
    spec = np.asarray(spec)
    z_pool = np.asarray(z_pool.iloc[:, 0])

    # perturb score
    y = simulate_latent_from_specificity(
        spec,
        is_causal=is_causal,
        beta=beta,
        noise_sd=noise_sd,
        rng=rng,
    )

    z_sim = assign_z_from_latent(y, z_pool)

    return z_sim, G, is_causal


def compute_score_serial(X, idx):
    Xg = X[idx]
    Xr = X[~idx]

    mean_g = Xg.mean(axis=0).A1 if hasattr(Xg, "A1") else np.asarray(Xg.mean(axis=0)).ravel()
    mean_r = Xr.mean(axis=0).A1 if hasattr(Xr, "A1") else np.asarray(Xr.mean(axis=0)).ravel()

    logfc = np.log2((mean_g + 1e-9) / (mean_r + 1e-9))
    return logfc


def compute_score(args):
    '''calculate spec score for every cell type'''
    # Unpack task parameters
    (idx,
     shm_data_name, shm_indices_name, shm_indptr_name,
     data_dtype, indices_dtype, indptr_dtype, data_shape,
     indices_shape, indptr_shape, csr_shape) = args

    # Attach to shared memory
    shm_data = SharedMemory(name=shm_data_name)
    shm_indices = SharedMemory(name=shm_indices_name)
    shm_indptr = SharedMemory(name=shm_indptr_name)

    try:
        # Reconstruct arrays from shared memory
        data = np.ndarray(data_shape, dtype=data_dtype, buffer=shm_data.buf)
        indices = np.ndarray(indices_shape, dtype=indices_dtype, buffer=shm_indices.buf)
        indptr = np.ndarray(indptr_shape, dtype=indptr_dtype, buffer=shm_indptr.buf)

        # Reconstruct CSR matrix
        X = csr_matrix((data, indices, indptr), shape=csr_shape)

        Xg = X[idx]
        Xr = X[~idx]

        mean_g = Xg.mean(axis=0).A1 if hasattr(Xg, "A1") else np.asarray(Xg.mean(axis=0)).ravel()
        mean_r = Xr.mean(axis=0).A1 if hasattr(Xr, "A1") else np.asarray(Xr.mean(axis=0)).ravel()

        logfc = np.log2((mean_g + 1e-9) / (mean_r + 1e-9))
        return logfc

    finally:
        # Close shared memory handles in worker
        shm_data.close()
        shm_indices.close()
        shm_indptr.close()


def get_spec_score(adata, group_key='cell_type'):
    """
    For every gene in every cell type, calculate prob of this gene is higher than the rest of cell types
    Use sparse matrix to calculate var and mean, make it faster
    """
    print('[INFO] Calculate specificity score')

    t0 = time()
    # transpose of input matrix and convert to sparse matrix
    data = adata.layers['raw']
    # Extract CSR components
    data_array = data.data
    indices_array = data.indices
    indptr_array = data.indptr
    csr_shape = data.shape

    # Create shared memory blocks
    shm_data = SharedMemory(create=True, size=data_array.nbytes)
    shm_indices = SharedMemory(create=True, size=indices_array.nbytes)
    shm_indptr = SharedMemory(create=True, size=indptr_array.nbytes)

    # Copy data to shared memory
    np_data_shm = np.ndarray(data_array.shape, dtype=data_array.dtype, buffer=shm_data.buf)
    np_data_shm[:] = data_array[:]

    np_indices_shm = np.ndarray(indices_array.shape, dtype=indices_array.dtype, buffer=shm_indices.buf)
    np_indices_shm[:] = indices_array[:]

    np_indptr_shm = np.ndarray(indptr_array.shape, dtype=indptr_array.dtype, buffer=shm_indptr.buf)
    np_indptr_shm[:] = indptr_array[:]

    # get idx
    cat_idx = [np.asarray(adata.obs[group_key] == cat) for cat in adata.obs[group_key].unique()]

    # Prepare task parameters
    tasks = [
        (
            idx,
            shm_data.name,
            shm_indices.name,
            shm_indptr.name,
            data_array.dtype,
            indices_array.dtype,
            indptr_array.dtype,
            data_array.shape,
            indices_array.shape,
            indptr_array.shape,
            csr_shape
        )
        for idx in cat_idx
    ]
    print('[INFO] Took %.2f min to parepare input' % ((time() - t0) / 60))

    t0 = time()
    try:
        with Pool(20) as pool:
            res = pool.map(compute_score, tasks)
    finally:
        # Cleanup shared memory in main process
        shm_data.close()
        shm_data.unlink()
        shm_indices.close()
        shm_indices.unlink()
        shm_indptr.close()
        shm_indptr.unlink()
    print('[INFO] Took %.2f min to calculate specificity score' % ((time() - t0) / 60))

    return pd.DataFrame(np.vstack(res), index=adata.obs[group_key].unique(), columns=adata.var_names)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=str, required=True)
    parser.add_argument("--logfc", type=str, required=True)
    parser.add_argument("--signal_frac", type=float, default=0.05)
    parser.add_argument("--beta", type=float, default=1.0)
    parser.add_argument("--noise_sd", type=float, default=1.0)
    parser.add_argument("--sample_rate", type=float, default=1.0)
    parser.add_argument("--min_ncell", type=int, default=30)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--n_run", type=int, default=100)
    parser.add_argument("--outdir", type=str, default="results")
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    # load parameters
    input_data = args.data
    celltype_spec_file = args.logfc
    outdir = args.outdir
    signal_frac = args.signal_frac
    beta = args.beta
    noise_sd = args.noise_sd
    sample_rate = args.sample_rate
    min_ncell = args.min_ncell
    seed = args.seed
    n_run = args.n_run

    # load sc
    adata = sc.read(input_data)

    # mk outdir
    Path(outdir).mkdir(exist_ok=True)

    # cell type count
    ct_cnts = adata.obs['cell_type'].value_counts()

    # load disease
    zstat_indir = '../data/TM_FACS/magmaz'
    files = Path(zstat_indir).glob('*.genes.out')
    files = [str(i) for i in files]

    # get spec score
    if not Path(celltype_spec_file).exists():
        celltype_spec = get_spec_score(adata, group_key='cell_type')
        with open(celltype_spec_file, 'wb') as f:
            pickle.dump(celltype_spec, f)
    else:
        # load spec score
        with open(celltype_spec_file, 'rb') as f:
            celltype_spec = pickle.load(f)

    # randomizer
    if not seed:
        rng = np.random.default_rng()
        seed = rng.integers(1e9)  # save actually used seed for reproduction
    rng = np.random.default_rng(seed)

    # mk outdir for current run
    setting = f'sf-{signal_frac}__b-{beta}__ns-{noise_sd}__sr-{sample_rate}__seed-{seed}'
    suboutdir = f'{outdir}/{setting}'
    Path(suboutdir).mkdir(exist_ok=True)

    # get cell types with enough purity
    cell_types = adata.obs['cell_type'].unique()

    n_finished = 0
    while 1:
        if n_finished == n_run:
            break

        # get spec score for a random chosen cell type
        target_ct = rng.choice(cell_types, 1)[0]

        # skip if number of selected cells are too small
        n_target = int(ct_cnts.loc[target_ct] * sample_rate)
        if n_target < min_ncell:
            continue

        # select a random disease
        score_df = pd.read_csv(rng.choice(files, 1)[0], header=0, index_col=None, sep=r'\s+')
        score_df['GENE'] = score_df['GENE'].astype(str)
        score_df = score_df.loc[:, ['GENE', 'ZSTAT']].copy()
        score_df = score_df.set_index('GENE')

        # select causal genes at cell type level
        z_sim, shared_genes, is_causal = simulate_gwas_from_specificity_celltype(
            celltype_spec.loc[target_ct, :], score_df,
            signal_frac=signal_frac, beta=beta, noise_sd=noise_sd,
            rng=rng
        )
        sim_df = pd.DataFrame({'GENE': shared_genes, 'ZSTAT': z_sim, 'is_causal': is_causal})

        # save
        gwasz_outfile = f"{suboutdir}/gwasz__{n_finished}.tsv"
        sim_df.to_csv(gwasz_outfile, header=True, index=False, sep='\t')
        tc_outfile = f"{suboutdir}/tc__{n_finished}.csv"
        with open(tc_outfile, 'w') as f:
            f.write(f'{target_ct}\n')

        n_finished += 1

        if n_finished % 10 == 0:
            print(f'[INFO] {n_finished} finished')
