# src Code Overview

This directory contains scripts for:
1. preprocessing single-cell data,
2. generating null and causal GWAS-like simulation z-scores,
3. preparing benchmark-tool inputs (ICePop/scDRS/Seismic),
4. running Seismic association workflows.

Note: this folder does not provide end-to-end launch wrappers for ICePop/scDRS themselves; most scripts prepare inputs, while Seismic runners are included.

## Data preprocessing and utilities
- `normalize_data.py`: Normalizes/log-transforms TM_FACS and mouse_colon count data, computes PCA/neighbors/UMAP, and writes normalized `.h5ad`.
- `prepare_sim_data_subset.py`: Builds `TM_subset_cnt.h5ad` for simulation by subsampling cells, filtering too small cell types/low-expression genes, and converting expression scores to human space.
- `mm_gene_conversion.py`: Builds a mouse gene symbol/synonym -> GeneID mapping from NCBI `gene_info.gz`.
- `seismic_geneconv.py`: Writes mouse->human ortholog mapping (`mm2hs.tsv`) for Seismic cross-species translation.
- `h5ad2rds.R`: Converts normalized TM_FACS/mouse_colon `.h5ad` files to Seismic-compatible `SingleCellExperiment` RDS (`expr.rds`).
- `h5ad2rds_causal.R`: Generic `.h5ad` -> `.rds` converter used in synthetic cell type simulation Seismic workflows.

## Null simulation
- `null_sim.py`: Creates null MAGMA-like z-score tables by permuting `ZSTAT` values over genes shared with simulation data.

## Causal simulation
- `casual_sim.py`: Main causal simulator using cell-type specificity (logFC-style signal):
  - chooses target cell type,
  - injects signal with `signal_frac`, `beta`, and `noise_sd`,
  - rank-matches to empirical MAGMA z-score distribution,
  - outputs `gwasz__*.tsv` and `tc__*.csv`.
  - no heterogeneity is modeled in this configuration. The `--sample_rate` is fixed at 1.0. Note that `sample_rate` is the complement of disease signal heterogeneity level (i.e., `heterogeneity = 1 − sample_rate`).
- `submit_causal.py`: Grouped/chunked Slurm submission wrapper for `casual_sim.py` (multiple simulation commands per Slurm job).

## Synthetic cell-type causal simulation
- `casual_sim_synthetic.py`: Causal simulator simulating disease heterogeneity setting using synthetic assoicated/non-assoicated cells.
  - `--sample_rate` can be varied in this configuration. Note that `sample_rate` is the complement of disease signal heterogeneity level (i.e., `heterogeneity = 1 − sample_rate`).
- `submit_causal_synthetic.py`: Grouped/chunked Slurm submission wrapper for synthetic-cell-type causal simulation.
- `submit_causal_synthetic_sanitized.py`: Sanitized Slurm submission wrapper for synthetic-cell-type causal simulation (no hardcoded private path/account; configurable via args/environment).
- `add_synthetic_ct_labels.py`: Adds per-run synthetic labels (`run_i`) to simulation adata and writes normalized/unnormalized synthetic-labeled `.h5ad`.

## Benchmark input preparation
### ICePop
- `precalculate_score_for_simulation.py`: Prepares cached specificity/PCA-related inputs for faster repeated ICePop simulation benchmarking.

### scDRS
- `scdrs_sim_input.py`: Converts null + causal simulation z-score outputs into chunked scDRS `.gs` files.
- `scdrs_sim_input_synthetic_ct.py`: Converts synthetic-cell-type causal outputs into scDRS `.gs` files.
- `scdrs_input.py`: Prepares TM_FACS scDRS inputs by mouse->human ortholog conversion and MAGMA geneset conversion.
- `scdrs_input_asd.py`: Prepares mouse_colon ASD scDRS inputs using ortholog conversion and ASD MAGMA geneset processing.

### Seismic
- `run_seismic.R`: Runs Seismic for TM_FACS/mouse_colon trait analyses (including optional species-mapping workflow).
- `run_seismic_null.R`: Seismic runner for null simulation.
- `run_seismic_causal.R`: Seismic runner for causal simulation.
- `run_seismic_causal_synthetic_ct.R`: Seismic runner for synthetic cell-type causal simulation.
