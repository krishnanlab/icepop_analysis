import subprocess
from pathlib import Path
from sys import argv
from itertools import islice

# --------------------------------------------
# 1. Parameter grids
# --------------------------------------------

NOISE_SD_FIXED = 1.0
BETA_FIXED = 1.0
FRAC_SIGNAL_FIXED = 0.01
SAMPLE_RATE_FIXED = 1.0

FRAC_SIGNAL_LIST = [0.001, 0.005, 0.01, 0.05, 0.1]
NOISE_SD_LIST = [0.1, 0.5, 1.0, 2.0, 5.0]
BETA_LIST = [0.1, 0.5, 1.0, 2.0, 5.0]

n_run = 100
min_ncell = 40
n_batch = 10
chunk_size = 20

data = argv[1]
logfc = argv[2]
outdir = argv[3]
Path(outdir).mkdir(exist_ok=True)

MAIN_SCRIPT = "path/to/causal_sim.py"  # update this

# --------------------------------------------
# 2. SLURM TEMPLATE (multi-job)
# --------------------------------------------

SLURM_TEMPLATE = """#!/bin/bash -l
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --mem=128GB
#SBATCH --cpus-per-task=1
#SBATCH --job-name={jobname}
#SBATCH --output=logs/{jobname}_%j.out

cd ${{SLURM_SUBMIT_DIR}}

{commands}
"""

# --------------------------------------------
# 3. Helpers
# --------------------------------------------

def chunked(iterable, size):
    it = iter(iterable)
    while True:
        chunk = list(islice(it, size))
        if not chunk:
            break
        yield chunk


def build_command(frac_signal, noise_sd, beta, sample_rate, batch):
    return f"""
python {MAIN_SCRIPT} \
    --data {data} \
    --logfc {logfc} \
    --min_ncell {min_ncell} \
    --signal_frac {frac_signal} \
    --noise_sd {noise_sd} \
    --beta {beta} \
    --sample_rate {sample_rate} \
    --n_run {n_run} \
    --outdir {outdir}
"""


# --------------------------------------------
# 4. Build all tasks
# --------------------------------------------

combos = []

# A: vary signal fraction
for fs in FRAC_SIGNAL_LIST:
    combos.append((fs, NOISE_SD_FIXED, BETA_FIXED, SAMPLE_RATE_FIXED))

# B: vary noise
for ns in NOISE_SD_LIST:
    combos.append((FRAC_SIGNAL_FIXED, ns, BETA_FIXED, SAMPLE_RATE_FIXED))

# C: vary beta
for b in BETA_LIST:
    combos.append((FRAC_SIGNAL_FIXED, NOISE_SD_FIXED, b, SAMPLE_RATE_FIXED))

unique_combos = list(set(combos))

all_tasks = []
for fs, ns, b, sr in unique_combos:
    for batch in range(n_batch):
        all_tasks.append((fs, ns, b, sr, batch))


# --------------------------------------------
# 5. Submit grouped jobs
# --------------------------------------------

print(f"Submitting grouped jobs ({chunk_size} simulations per job)...")

for i, group in enumerate(chunked(all_tasks, chunk_size)):
    jobname = f"simgrp_{i}"

    cmd_list = []
    for fs, ns, b, sr, batch in group:
        cmd = build_command(fs, ns, b, sr, batch)
        cmd_list.append(cmd)

    script_text = SLURM_TEMPLATE.format(
        jobname=jobname,
        commands="\n".join(cmd_list)
    )

    script_path = f"slurm_group_{i}.sb"
    with open(script_path, "w") as f:
        f.write(script_text)

    subprocess.run(["sbatch", script_path])
    print(f"[SUBMITTED] {jobname}")

print("All jobs submitted.")
