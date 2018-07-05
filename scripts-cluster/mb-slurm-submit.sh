#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user@email
#SBATCH --array=1-94
#SBATCH -p long
export PATH="/workspace/software/bin:$PATH"
/s/mrbayes-3.2.6-1/bin/mb $SLURM_ARRAY_TASK_ID.nex
