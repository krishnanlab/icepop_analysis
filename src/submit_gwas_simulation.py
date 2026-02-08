import subprocess
from pathlib import Path
from sys import argv

# --------------------------------------------
# 1. Parameter grids
# --------------------------------------------

# fixed parameters
NOISE_SD_FIXED = 1.0
BETA_FIXED = 1.0
FRAC_SIGNAL_FIXED = 0.01
SAMPLE_RATE_FIXED = 1.0

# experiment A: vary signal strength
FRAC_SIGNAL_LIST = [0.001, 0.005, 0.01, 0.05, 0.1]

# experiment B: vary noise_sd
NOISE_SD_LIST = [0.1, 0.5, 1.0, 2.0, 5.0]

# experiment C: vary beta
BETA_LIST = [0.1, 0.5, 1.0, 2.0, 5.0]

# experiment D: vary cell sampling rate
SAMPLE_RATE_LIST = [0.2, 0.4, 0.6, 0.8, 1.0]

# signal source
signal_source = 'cell_type'

# job paramters
n_run = 100
min_ncell = 30
n_batch = 10

# input data
data = argv[1]
logfc = argv[2]

# create base outdir
outdir = argv[3]
Path(outdir).mkdir(exist_ok=True)


# --------------------------------------------
# 2. Where to run your simulation code
# --------------------------------------------
MAIN_SCRIPT = "../src/casual_sim.py"

SLURM_TEMPLATE = """#!/bin/bash -l
#SBATCH --time=03:50:00
#SBATCH --nodes=1
#SBATCH --mem=200GB
#SBATCH --cpus-per-task=1
#SBATCH --job-name={jobname}
#SBATCH --account=general
#SBATCH --output=../slurms/{jobname}_%j.out

cd ${{SLURM_SUBMIT_DIR}}

/mnt/ufs18/rs-032/FishEvoDevoGeno/Hao/icepop_analysis/.venv/bin/python {main_script} \\
    --data {data} \\
    --logfc {logfc} \\
    --min_ncell {min_ncell} \\
    --signal_frac {frac_signal} \\
    --noise_sd {noise_sd} \\
    --beta {beta} \\
    --sample_rate {sample_rate} \\
    --n_run {n_run} \\
    --outdir {outdir}
"""


# --------------------------------------------
# 3. Function to submit one job
# --------------------------------------------
def submit_job(frac_signal, noise_sd, beta, sample_rate, batch):
    jobname = f"fs{frac_signal}_ns{noise_sd}_b{beta}_sr{sample_rate}_{batch}".replace('.', 'p')
    slurm_script_name = f"../run/slurm_simulate_gwasz_{jobname}.sb"

    script_text = SLURM_TEMPLATE.format(
        data=data,
        logfc=logfc,
        min_ncell=min_ncell,
        jobname=jobname,
        main_script=MAIN_SCRIPT,
        frac_signal=frac_signal,
        noise_sd=noise_sd,
        beta=beta,
        sample_rate=sample_rate,
        n_run=n_run,
        outdir=outdir
    )

    # Write slurm script
    with open(slurm_script_name, "w") as f:
        f.write(script_text)

    # Submit to Slurm
    subprocess.run(["sbatch", slurm_script_name])
    print(f"[SUBMITTED] {jobname}")


# --------------------------------------------
# 4. Submit all experiments (one by one)
# --------------------------------------------
if __name__ == "__main__":
    combos = []

    # A: signal strength
    for fs in FRAC_SIGNAL_LIST:
        combos.append((fs, NOISE_SD_FIXED, BETA_FIXED, SAMPLE_RATE_FIXED))

    # B: noise_sd
    for ns in NOISE_SD_LIST:
        combos.append((FRAC_SIGNAL_FIXED, ns, BETA_FIXED, SAMPLE_RATE_FIXED))

    # C: beta
    for b in BETA_LIST:
        combos.append((FRAC_SIGNAL_FIXED, NOISE_SD_FIXED, b, SAMPLE_RATE_FIXED))

    # D: sample rate
    for sr in SAMPLE_RATE_LIST:
        combos.append((FRAC_SIGNAL_FIXED, NOISE_SD_FIXED, BETA_FIXED, sr))

    unique_combos = set(combos)

    print("submitting jobs...")
    for fs, ns, b, sr in unique_combos:
        for batch in range(n_batch):
            submit_job(
                frac_signal=fs,
                noise_sd=ns,
                beta=b,
                sample_rate=sr,
                batch=batch
            )

    print("\nAll jobs submitted.")
