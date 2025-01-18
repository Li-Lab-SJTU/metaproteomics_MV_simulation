#!/bin/bash
#SBATCH -J simulation
#SBATCH -p 64c512g
#SBATCH -o ./err_out/%A_%a.out
#SBATCH -e ./err_out/%A_%a.err
#SBATCH -n 1
#SBATCH --mail-type=fail,end
#SBATCH --mail-user=your@email.com
#SBATCH --array=1-2000%100

echo "SLURM_ARRAY_JOB_ID:" $SLURM_ARRAY_JOB_ID
echo "SLURM_ARRAY_TASK_ID:" $SLURM_ARRAY_TASK_ID
echo "SLURM_JOBID:" $SLURM_JOBID

id=$SLURM_ARRAY_TASK_ID
arg=$(sed -n "${id}p" simu_args.txt)
echo "arg:" $arg

Rscript simulation.R $arg
