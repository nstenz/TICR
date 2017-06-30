#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=solislemus@wisc.edu
#SBATCH --array=1-100
export PATH="/workspace/software/bin:$PATH"
/workspace/claudia/potato/bucky-slurm/TICR/scripts/bucky-slurm.pl mbsum -q $SLURM_ARRAY_TASK_ID -o only100
