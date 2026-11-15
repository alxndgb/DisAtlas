#!/bin/bash
# Full Figure S16 a K-selection bootstrap (100 replicates x 15 K) as a 150-task array, about 10 min per task on 4 CPUs.
#   sbatch scripts/slurm/k_selection_boot_array.sh OUTDIR
#   then: python -m disatlas.rebuild.k_selection_boot OUTDIR --merge
#         DISATLAS_KSEL_RUN=OUTDIR python -m pytest tests/test_rebuild_k_selection_boot.py
#SBATCH -p rd-bds
#SBATCH -c 4
#SBATCH --mem=14G
#SBATCH -t 06:00:00
#SBATCH --array=0-149%35
#SBATCH -J ksel
cd "${SLURM_SUBMIT_DIR}"
PYTHONPATH=src .venv/bin/python -m disatlas.rebuild.k_selection_boot "$1" --task "$SLURM_ARRAY_TASK_ID"
