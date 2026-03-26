# src Code Overview

This directory contains scripts used to:
1. generate simulation GWAS z-scores (null + causal),
2. preprocess single-cell inputs,
3. prepare gene-set inputs for external GWAS tools (e.g., Seismic, scDRS),
4. wrappers to run Seismic on simulation/Tabula Muris/mouse colon data.

Note: There are no wrappers provided for running ICePop or scDRS, as both tools offer command-line interfaces and can be executed directly.

## Prepare simulation data
- `prepare_sim_data_subset.py`: Creates a smaller simulation-ready expression dataset `TM_subset_cnt.h5ad` by:
  - sampling a fraction of cells,
  - removing cell types whose size is too small,
  - dropping genes not expressed in enough cells,
  - convert expression across species to human space.

## Null simulation
- `null_sim.py`: Generates *null* MAGMA z-score tables by permuting MAGMA-derived `ZSTAT` values across genes (restricted to genes present in `TM_subset_cnt.h5ad`). The outputs are written as per-run TSV files.

## Causal simulation
- `casual_sim.py`: Core *causal* simulation generator.
  - Computes cell-type specificity scores (logFC) from raw counts (implemented as cell-type logFC-style scores).
  - Generates simulated MAGMA z-scores by mixing:
    - causal gene fraction (`signal_frac`),
    - signal strength (`beta`),
    - cell sampling heterogeneity (`sample_rate`),
    - Gaussian noise (`noise_sd`).
  - Saves per-run outputs:
    - `gwasz__*.tsv`: gene `ZSTAT` plus a causal indicator,
    - `tc__*.csv`: the target cell type used for the run.
- `submit_gwas_simulation.py`: submission wrapper that submits `casual_sim.py` jobs over grids of `signal_frac`, `noise_sd`, `beta`, and `sample_rate`.

## Benchmark inputs / tool execution
### ICePop
- `precalculate_score_for_simulation.py`: Precomputes expression specificity scores to speed up repeated simulation runs, stored in an `.h5ad`.

### scDRS
- `scdrs_sim_input.py`: Converts simulated MAGMA z-score tables (causal and null) into scDRS-ready `.gs` gene-set files (chunked for many runs).
- `scdrs_input.py`: Prepares scDRS inputs for the TM FACS dataset by:
  - converting mouse genes to human orthologs,
  - converting MAGMA gene-set outputs into scDRS `.gs` files.
- `scdrs_input_asd.py`: Prepares scDRS mouse colon inputs by:
  - mapping mouse genes to human orthologs (using `HomologyData`),
  - convert expression from mouse to human
  - writing processed `.h5ad` plus trait `.gs` gene-set files.

### Seismic
- `run_seismic_null.R`: Runs Seismic for each generated null gene-set input TSV.
- `run_seismic_casual.R`: Runs Seismic for causal simulation gene-set inputs.
- `run_seismic.R` Runs Seismic for TM FACS and mouse colon data

## Others - preparation codes
- `normalize_data.py`: Normalizes and log-transforms raw count datasets; saves normalized `.h5ad` objects for both TM FACS and mouse colon workflows.
- `h5ad2rds.R`: Converts processed `h5ad` expression objects into Seismic-compatible `SingleCellExperiment` RDS objects (`expr.rds`) for mouse colon and TM_FACS.
- `mm_gene_conversion.py`: Downloads/reads NCBI `gene_info.gz` and builds a mouse gene-name → GeneID mapping (saved as a pickle).
- `seismic_geneconv.py`: Writes a mouse→human ortholog mapping file (`mm2hs.tsv`) for Seismic cross-species translation.
