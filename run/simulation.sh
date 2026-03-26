# subset data
python ../src/prepare_sim_data_subset.py

# get metacell
poetry run icepop metacell \
    --h5ad ../data/simulation/TM_subset_cnt.h5ad \
    --outdir ../data/simulation/mc-120 \
    --ncell_per_mc 120 \
    --save_name sim

poetry run icepop metacell \
    --h5ad ../data/simulation/TM_subset_cnt.h5ad \
    --outdir ../data/simulation/mc-75 \
    --save_name sim

poetry run icepop metacell \
    --h5ad ../data/simulation/TM_subset_cnt.h5ad \
    --outdir ../data/simulation/mc-50 \
    --ncell_per_mc 50 \
    --save_name sim

poetry run icepop metacell \
    --h5ad ../data/simulation/TM_subset_cnt.h5ad \
    --outdir ../data/simulation/mc-30 \
    --ncell_per_mc 30 \
    --save_name sim

# normalize and get metacell pca centroid
python ../src/precalculate_score_for_simulation.py

#############################################
# casual simulation
#############################################
# generate logfc score from a random run
python ../src/casual_sim.py \
    --data ../data/simulation/mc-75/TM_subset__score_calc.h5ad \
    --logfc ../data/simulation/celltype_logfc.pkl \
    --signal_frac 0.01 \
    --beta 1.0 \
    --noise_sd 0.5 \
    --sample_rate 0.2 \
    --n_run 1 \
    --outdir ../data/simulation/causal_simulation \
    --seed 123

# simulate data
python ../src/submit_gwas_simulation.py \
    ../data/simulation/mc-75/TM_subset__score_calc.h5ad \
    ../data/simulation/celltype_logfc.pkl \
    ../data/simulation/causal_simulation_test

python ../src/run_gwas_simulation.py \
    ../data/simulation/mc-75/TM_subset__score_calc.h5ad \
    ../data/simulation/celltype_logfc.pkl \
    ../data/simulation/causal_simulation

# convert h5ad to rds for seismic
Rscript ../src/h5ad2rds.R

# prepare input for scdrs
python ../src/scdrs_sim_input.py

# run in parallel
sbatch ../run/run_scDRS_casual_sim.sb
sbatch ../run/run_scDRS_casual_sim_ds.sb
sbatch ../run/run_seismic_causal_sim.sb
sbatch ../run/run_icepop_causal_sim__mc-30.sb
sbatch ../run/run_icepop_causal_sim__mc-50.sb
sbatch ../run/run_icepop_causal_sim__mc-75.sb

#############################################
# null simulation
#############################################
# permute y to simulate null data
python ../src/null_sim.py

# run in parallel
sbatch ../run/run_scDRS_null_sim.sb
sbatch ../run/run_seismic_null_sim.sb
sbatch ../run/run_icepop_causal_sim__mc-30.sb
sbatch ../run/run_icepop_causal_sim__mc-50.sb
sbatch ../run/run_icepop_causal_sim__mc-75.sb