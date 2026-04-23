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

SAMPLE_RATE_LIST = [0.2, 0.4, 0.6, 0.8, 1.0]

n_run = 100
min_ncell = 100
n_batch = 10

# control how many runs per Slurm job
CHUNK_SIZE = 20

data = argv[1]
outdir = argv[2]
Path(outdir).mkdir(exist_ok=True)

MAIN_SCRIPT = "path/to/causal_sim_synthetic.py"  # update as needed

# --------------------------------------------
# 2. SLURM TEMPLATE (multi-run)
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

for sr in SAMPLE_RATE_LIST:
    combos.append((FRAC_SIGNAL_FIXED, NOISE_SD_FIXED, BETA_FIXED, sr))

unique_combos = list(set(combos))

all_tasks = []
for fs, ns, b, sr in unique_combos:
    for batch in range(n_batch):
        all_tasks.append((fs, ns, b, sr, batch))


# --------------------------------------------
# 5. Submit grouped jobs
# --------------------------------------------

print(f"Submitting jobs with CHUNK_SIZE={CHUNK_SIZE} ...")

for i, group in enumerate(chunked(all_tasks, CHUNK_SIZE)):
    jobname = f"simgrp_synthetic_{i}"

    cmd_list = []
    for fs, ns, b, sr, batch in group:
        cmd = build_command(fs, ns, b, sr, batch)
        cmd_list.append(cmd)

    script_text = SLURM_TEMPLATE.format(
        jobname=jobname,
        commands="\n".join(cmd_list)
    )

    script_path = f"slurm_group_synthetic_{i}.sb"
    with open(script_path, "w") as f:
        f.write(script_text)

    subprocess.run(["sbatch", script_path])
    print(f"[SUBMITTED] {jobname}")

print("All jobs submitted.")
