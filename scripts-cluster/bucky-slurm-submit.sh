#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user@email
#SBATCH --array=1-100
export PATH="/workspace/software/bin:$PATH"
/workspace/claudia/software/TICR/scripts/bucky-slurm.pl mbsum -q $SLURM_ARRAY_TASK_ID
