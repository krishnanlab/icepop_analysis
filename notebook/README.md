# Notebook Overview

This directory contains Jupyter notebooks used for preprocessing, simulations, analyses, and downstream visualization.

## Process data
- `process_mouse_colon_scdata.ipynb`: Preprocesses mouse colon single-cell data (filters lowly expressed genes, renames cell types, converts gene symbols to Entrez IDs) and saves the cleaned object.
- `process_TM_FACS.ipynb`: Preprocesses Tabula Muris (TM) FACS input data (loads data, checks outlier cells, converts symbols to IDs) and saves the prepared dataset(s) for downstream ICePop analyses.

## Null simulation
- `null_sim.ipynb`: Evaluates ICePop under a null setting by aggregating metacell-level p-values from many runs and producing empirical calibration/QQ-style plots.
- `null_sim_ct.ipynb`: Tests null behavior stratified by cell-type size and compares ICePop against other methods (e.g., Seismic, scDRS).
- `null_sim_icepop_mc_size.ipynb`: Studies how ICePop metacell size (`mc-30/50/75`) changes null calibration; includes both metacell-level and cell-type-level null summaries.

## Casual simulation
- `causal_sim.ipynb`: Runs causal simulation evaluations using datasets with known ground-truth signals and compares inferred association behavior across methods. Evaluates statistical power under three settings: (A) varying the fraction of causal genes, (B) varying the noise variance, and (C) varying the signal scaling factor.
- `causal_sim_heterogeneity.ipynb`: Runs causal simulation evaluations with known ground-truth signals and compares inferred association behavior across methods. Evaluates statistical power under heterogeneous disease–cell associations within cell types.
- `causal_sim_compare_mc.ipynb`: Compares causal simulation performance across different ICePop metacell sizes (mc-30/50/75). Evaluates the statistical power of ICePop across metacell sizes under three settings: (A) varying the fraction of causal genes, (B) varying the noise variance, and (C) varying the signal scaling factor.
- `causal_sim_compare_mc_heterogeneity.ipynb`: Compares causal simulation performance across different ICePop metacell sizes (mc-30/50/75). Evaluates the statistical power of ICePop across metacell sizes under heterogeneous disease–cell associations within cell types.

## TM analysis
- `TM_analysis.ipynb`: Main TM FACS analysis workflow; compares differential associations, provides example subpopulations tied to traits/diseases, and saves trait/cell-type mapping outputs.
- `TM_analysis_compare_mc.ipynb`: Companion comparison notebook to assess how metacell choices affect results.

## Disease clusterting
- `disease_disease_sim.ipynb`: Performs disease-trait clustering/separation using distance metrics (z-score based and metacell-association based) and visualizes separations such as leukocyte-count trait vs autoimmune disease and structure among hair-related traits.

## ASD mouse colon enteric neuron analysis
- `asd__mouse_colon.ipynb`: ASD-focused analysis in mouse colon enteric neurons; checks overall disease-associated cell types, zooms into sensory neuron mechanisms, and examines if influential genes contributing to ASD-enteric sensory neuron asociation show up in SFARI high-confidence genes.
