#!/bin/bash
# Full SI_degree_distributions power-law bootstrap: 6 fit settings x 20 chunks (1,000 replicates per setting).
#   python -m disatlas.rebuild.powerlaw_csn observed OUTDIR
#   sbatch scripts/slurm/powerlaw_array.sh OUTDIR
#   python -m disatlas.rebuild.powerlaw_csn collect OUTDIR
#   DISATLAS_POWERLAW_BOOT=OUTDIR python -m pytest tests/test_rebuild_powerlaw_csn.py
#SBATCH -p rd-bds
#SBATCH -c 1
#SBATCH --mem=3G
#SBATCH -t 08:00:00
#SBATCH --array=0-119%25
#SBATCH -J plboot
cd "${SLURM_SUBMIT_DIR}"
PYTHONPATH=src .venv/bin/python -m disatlas.rebuild.powerlaw_csn boot_task "$1" "$SLURM_ARRAY_TASK_ID" 20
