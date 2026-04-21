import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
import argparse
import pickle


def pick_active_inactive_cells(
    adata,
    n_active,
    n_inactive,
    rng,
    leiden_key="leiden",
    pca_key="X_pca",
    low_q=0.3,
    high_q=0.5
):
    if rng is None:
        rng = np.random.default_rng()

    leiden = adata.obs[leiden_key].astype(str)
    pca = adata.obsm[pca_key]

    cluster_sizes = leiden.value_counts()

    # --------------------------------
    # ACTIVE: pick ONE cluster large enough
    # --------------------------------
    valid_clusters = cluster_sizes[cluster_sizes >= n_active].index.tolist()

    if len(valid_clusters) == 0:
        raise ValueError("No leiden cluster large enough for active")

    active_cluster = rng.choice(valid_clusters)

    # --------------------------------
    # ACTIVE CELLS: PCA-local patch inside cluster
    # --------------------------------
    cluster_mask = np.asarray(leiden == active_cluster)
    cluster_cells = np.array(adata.obs.index[cluster_mask])
    cluster_pca = pca[cluster_mask]

    # centroid
    centroid = np.median(cluster_pca, axis=0)

    # find seed = closest to centroid
    centroid_dists = np.linalg.norm(cluster_pca - centroid, axis=1)
    seed_idx = np.argmin(centroid_dists)
    seed_coord = cluster_pca[seed_idx]

    # local neighborhood around seed
    dists = np.linalg.norm(cluster_pca - seed_coord, axis=1)
    order = np.argsort(dists)

    active_cells = cluster_cells[order[:n_active]]

    inactive_cells = []
    if n_inactive > 0:
        # --------------------------------
        # INACTIVE: from moderate-distance clusters (MODIFIED)
        # --------------------------------
        cluster_idx_all = adata.obs_names.get_indexer(cluster_cells)
        active_center = pca[cluster_idx_all].mean(axis=0)

        other_clusters = [c for c in cluster_sizes.index if c != active_cluster]

        cluster_dist = {}
        for c in other_clusters:
            c_cells = np.array(adata.obs.index[leiden == c])
            c_idx = adata.obs_names.get_indexer(c_cells)
            center = pca[c_idx].mean(axis=0)
            cluster_dist[c] = np.linalg.norm(center - active_center)

        # -------- select clusters in moderate band --------
        cluster_names = list(cluster_dist.keys())
        dist_values = np.array([cluster_dist[c] for c in cluster_names])

        low = np.quantile(dist_values, low_q)
        high = np.quantile(dist_values, high_q)

        selected_clusters = [
            c for c in cluster_names
            if (cluster_dist[c] >= low) and (cluster_dist[c] <= high)
        ]

        # fallback if no cluster selected
        if len(selected_clusters) == 0:
            selected_clusters = cluster_names

        # shuffle to avoid ordering bias
        rng.shuffle(selected_clusters)

        # --------------------------------
        # sample cells from selected clusters (unchanged logic)
        # --------------------------------

        for c in selected_clusters:
            c_cells = np.array(adata.obs.index[leiden == c])

            c_idx = adata.obs_names.get_indexer(c_cells)
            c_pca = pca[c_idx]

            # sort inside cluster by distance to active center
            dists = np.linalg.norm(c_pca - active_center, axis=1)
            order = np.argsort(dists)
            c_cells = c_cells[order]

            need = n_inactive - len(inactive_cells)
            take = c_cells[:need]
            inactive_cells.extend(take.tolist())

            if len(inactive_cells) >= n_inactive:
                break

        if len(inactive_cells) < n_inactive:
            raise ValueError("Not enough cells for inactive")

    # --------------------------------
    # compute logfc
    # --------------------------------
    active_idx = adata.obs_names.get_indexer(active_cells)
    logfc_active = compute_logfc(adata.layers["raw"], active_idx)

    if len(inactive_cells) > 0:
        inactive_idx = adata.obs_names.get_indexer(inactive_cells)
        logfc_inactive = compute_logfc(adata.layers["raw"], inactive_idx)
    else:
        logfc_inactive = np.zeros(logfc_active.size)

    return logfc_active, logfc_inactive, np.array(active_cells), np.array(inactive_cells)


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

    # gate on positive active-enriched signal
    spec_contrast = np.clip(spec_contrast, -50, 50)
    # weights = np.log1p(np.exp(spec_contrast))
    weights = np.exp(spec_contrast)
    weights /= weights.sum()
    idx = rng.choice(G, size=n_signal, replace=False, p=weights)

    is_causal = np.zeros(G, dtype=bool)
    is_causal[idx] = True
    return is_causal


# ------------------------------------------
# LATENT SIGNAL
# ------------------------------------------
def simulate_latent(spec_contrast, is_causal, beta, noise_sd, rng):
    noise = rng.normal(0, noise_sd, size=len(spec_contrast))

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
def simulate_gwas(logfc_active, logfc_inactive, z_df, signal_frac, beta, noise_sd, sample_rate, rng):
    genes = z_df.index
    z_pool = np.asarray(z_df.iloc[:, 0])

    spec_active = np.asarray(logfc_active)
    spec_inactive = np.asarray(logfc_inactive)

    # keep your original structure, only stabilize the contrast a bit
    w = 1.0 - sample_rate

    # require active enrichment to be positive; only penalize positive inactive enrichment
    spec_contrast = spec_active - w * spec_inactive

    is_causal = choose_causal_genes(spec_contrast, signal_frac, rng)

    y = simulate_latent(spec_contrast, is_causal, beta, noise_sd, rng)

    z_sim = assign_z(y, z_pool)

    return z_sim, genes, is_causal


# ------------------------------------------
# ARGS
# ------------------------------------------
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--data", required=True)
    p.add_argument("--signal_frac", type=float, default=0.01)
    p.add_argument("--beta", type=float, default=1.0)
    p.add_argument("--noise_sd", type=float, default=1.0)
    p.add_argument("--sample_rate", type=float, default=0.5)
    p.add_argument("--min_ncell", type=int, default=100)
    p.add_argument("--low_q", type=float, default=0.3)
    p.add_argument("--high_q", type=float, default=0.5)
    p.add_argument("--seed", type=int, default=None)
    p.add_argument("--n_run", type=int, default=100)
    p.add_argument("--outdir", default="results")
    return p.parse_args()


# ------------------------------------------
# MAIN
# ------------------------------------------
if __name__ == "__main__":
    args = parse_args()

    adata = sc.read(args.data)
    Path(args.outdir).mkdir(exist_ok=True)

    ct_counts = adata.obs["cell_type"].value_counts()

    # get candiate cts
    cell_types = ct_counts[ct_counts >= args.min_ncell].index.values

    zstat_files = list(Path("../data/magmaz").glob("*.genes.out"))

    # leiden clustering
    leiden_key = 'leiden'
    if leiden_key not in adata.obs.columns:
        raise ValueError('Cannot find leiden clustering')

    # randomizer
    seed = args.seed
    if seed is None:
        rng = np.random.default_rng()
        seed = rng.integers(1e9)
    rng = np.random.default_rng(seed)

    setting = f"sf-{args.signal_frac}__b-{args.beta}__ns-{args.noise_sd}__sr-{args.sample_rate}__seed-{seed}"
    outdir = Path(args.outdir) / setting
    outdir.mkdir(exist_ok=True)

    n_finished = 0
    save_cells = []
    print('start running')
    while n_finished < args.n_run:
        target_ct = rng.choice(cell_types)
        n_active = int(ct_counts.loc[target_ct] * args.sample_rate)
        n_inactive = ct_counts.loc[target_ct] - n_active

        try:
            (
                logfc_active, logfc_inactive,
                active_cells, inactive_cells
            ) = pick_active_inactive_cells(
                adata,
                n_active, n_inactive,
                rng,
                leiden_key="leiden",
                low_q=args.low_q,
                high_q=args.high_q,
            )
        except Exception as e:
            print(f"[ERROR] target_ct={target_ct}: {e}")
            continue

        score_df = pd.read_csv(
            rng.choice(zstat_files),
            sep=r"\s+"
        ).set_index("GENE")[["ZSTAT"]]

        # get shared genes
        score_df.index = score_df.index.astype(str)
        shared_genes = score_df.index.intersection(adata.var_names)
        score_df = score_df.loc[shared_genes, :].copy()
        gene_idx = adata.var_names.get_indexer(shared_genes)
        logfc_active = logfc_active[gene_idx]
        logfc_inactive = logfc_inactive[gene_idx]

        z_sim, genes, is_causal = simulate_gwas(
            logfc_active,
            logfc_inactive,
            score_df,
            args.signal_frac,
            args.beta,
            args.noise_sd,
            args.sample_rate,
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
        # SAVE ACTIVE / INACTIVE CELL INDICES TO PKL
        # --------------------------------
        save_cells.append(
            {
                "active_cells": np.asarray(active_cells),
                "inactive_cells": np.asarray(inactive_cells),
            }
        )

        n_finished += 1

        if n_finished % 10 == 0:
            print(f"[INFO] {n_finished} finished")

        with open(outdir / "tc.pkl", "wb") as f:
            pickle.dump(save_cells, f)
