import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
from time import time
from scipy.sparse import csr_matrix
import argparse


# ------------------------------------------
# ACTIVE / INACTIVE SPLIT + LOGFC
# ------------------------------------------
def compute_active_inactive_logfc(
    adata,
    target_ct,
    n_target,
    rng,
    celltype_key="cell_type",
    umap_key="X_umap",
):
    ct_mask = adata.obs[celltype_key] == target_ct
    ct_cells = np.array(adata.obs.index[ct_mask])

    if len(ct_cells) < n_target * 2:
        raise ValueError("Not enough cells for split")

    umap = adata.obsm[umap_key][ct_mask]
    x, y = umap[:, 0], umap[:, 1]

    theta = rng.uniform(0, np.pi)
    a, b = np.cos(theta), np.sin(theta)

    proj = a * x + b * y
    order = np.argsort(proj)

    active_cells = ct_cells[order[-n_target:]]
    inactive_cells = ct_cells[order[:n_target]]

    active_idx = adata.obs_names.get_indexer(active_cells)
    inactive_idx = adata.obs_names.get_indexer(inactive_cells)

    logfc_active = compute_logfc(adata.layers["raw"], active_idx)
    logfc_inactive = compute_logfc(adata.layers["raw"], inactive_idx)

    return logfc_active, logfc_inactive, active_cells


# ------------------------------------------
# LOGFC
# ------------------------------------------
def compute_logfc(X, idx):
    Xg = X[idx]
    Xr = X[~idx]

    mean_g = Xg.mean(axis=0).A1 if hasattr(Xg, "A1") else np.asarray(Xg.mean(axis=0)).ravel()
    mean_r = Xr.mean(axis=0).A1 if hasattr(Xr, "A1") else np.asarray(Xr.mean(axis=0)).ravel()

    return np.log2((mean_g + 1e-9) / (mean_r + 1e-9))


# ------------------------------------------
# CAUSAL GENES (ACTIVE-BIASED)
# ------------------------------------------
def choose_causal_genes(spec_contrast, signal_frac=0.1, rng=None):
    if rng is None:
        rng = np.random.default_rng()

    G = len(spec_contrast)
    n_signal = int(np.floor(signal_frac * G))
    if n_signal < 1:
        raise ValueError("signal_frac too small")

    weights = np.exp(spec_contrast)
    weights /= weights.sum()

    idx = rng.choice(G, size=n_signal, replace=False, p=weights)

    is_causal = np.zeros(G, dtype=bool)
    is_causal[idx] = True
    return is_causal


# ------------------------------------------
# LATENT SIGNAL
# ------------------------------------------
def simulate_latent(spec_active, spec_inactive, is_causal, beta, noise_sd, rng):
    noise = rng.normal(0, noise_sd, size=len(spec_active))

    spec_contrast = spec_active - spec_inactive
    spec_contrast = (spec_contrast - spec_contrast.mean()) / (spec_contrast.std() + 1e-6)

    y = noise.copy()
    y[is_causal] += beta * spec_contrast[is_causal]

    return y


# ------------------------------------------
# RANK MATCH
# ------------------------------------------
def assign_z(y, z_pool):
    y_rank = np.argsort(np.argsort(y))
    z_sorted = np.sort(z_pool)
    return z_sorted[y_rank]


# ------------------------------------------
# MAIN SIMULATION
# ------------------------------------------
def simulate_gwas(logfc_active, logfc_inactive, z_df, signal_frac, beta, noise_sd, rng):
    genes = z_df.index
    z_pool = np.asarray(z_df.iloc[:, 0])

    spec_active = np.asarray(logfc_active)
    spec_inactive = np.asarray(logfc_inactive)

    spec_contrast = spec_active - spec_inactive

    is_causal = choose_causal_genes(spec_contrast, signal_frac, rng)

    y = simulate_latent(spec_active, spec_inactive, is_causal, beta, noise_sd, rng)

    z_sim = assign_z(y, z_pool)

    return z_sim, genes, is_causal


# ------------------------------------------
# ARGS
# ------------------------------------------
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--data", required=True)
    p.add_argument("--signal_frac", type=float, default=0.05)
    p.add_argument("--beta", type=float, default=1.0)
    p.add_argument("--noise_sd", type=float, default=1.0)
    p.add_argument("--sample_rate", type=float, default=1.0)
    p.add_argument("--min_ncell", type=int, default=30)
    p.add_argument("--seed", type=int, default=None)
    p.add_argument("--n_run", type=int, default=100)
    p.add_argument("--outdir", default="results")
    p.add_argument("--save_active_cells", action="store_true")
    return p.parse_args()


# ------------------------------------------
# MAIN
# ------------------------------------------
if __name__ == "__main__":
    args = parse_args()

    adata = sc.read(args.data)
    Path(args.outdir).mkdir(exist_ok=True)

    ct_counts = adata.obs["cell_type"].value_counts()
    cell_types = adata.obs["cell_type"].unique()

    zstat_files = list(Path("../data/TM_FACS/magmaz").glob("*.genes.out"))

    rng = np.random.default_rng(args.seed)

    setting = f"sf-{args.signal_frac}__b-{args.beta}__ns-{args.noise_sd}__seed-{args.seed}"
    outdir = Path(args.outdir) / setting
    outdir.mkdir(exist_ok=True)

    n_finished = 0

    while n_finished < args.n_run:
        target_ct = rng.choice(cell_types)

        n_target = int(ct_counts.loc[target_ct] * args.sample_rate)
        if n_target < args.min_ncell:
            continue

        try:
            logfc_active, logfc_inactive, active_cells = compute_active_inactive_logfc(
                adata,
                target_ct,
                n_target,
                rng,
            )
        except:
            continue

        score_df = pd.read_csv(
            rng.choice(zstat_files),
            sep=r"\s+"
        ).set_index("GENE")[["ZSTAT"]]

        z_sim, genes, is_causal = simulate_gwas(
            logfc_active,
            logfc_inactive,
            score_df,
            args.signal_frac,
            args.beta,
            args.noise_sd,
            rng,
        )

        # --------------------------------
        # SAVE GWAS
        # --------------------------------
        sim_df = pd.DataFrame({
            "GENE": genes,
            "ZSTAT": z_sim,
            "is_causal": is_causal
        })

        sim_df.to_csv(outdir / f"gwasz__{n_finished}.tsv", sep="\t", index=False)

        # --------------------------------
        # SAVE CELL TYPE
        # --------------------------------
        with open(outdir / f"tc__{n_finished}.txt", "w") as f:
            f.write(f"{target_ct}\n")

        # --------------------------------
        # OPTIONAL: SAVE ACTIVE CELLS
        # --------------------------------
        if args.save_active_cells:
            pd.Series(active_cells).to_csv(
                outdir / f"active_cells__{n_finished}.txt",
                index=False,
                header=False
            )

        n_finished += 1

        if n_finished % 10 == 0:
            print(f"[INFO] {n_finished} finished")