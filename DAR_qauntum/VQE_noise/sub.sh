#!/bin/sh
#An example for serial job.
#SBATCH -J h6_test
#SBATCH -o job-%j.log
#SBATCH -e job-%j.err
#SBATCH -p normal
#SBATCH -N 1 -n 1
####SBATCH --gres=gpu:4
echo Running on hosts
echo Time is `date`
echo Directory is $PWD
echo This job runs on the following nodes:
echo $SLURM_JOB_NODELIST

python main_pennylane.py 

#mpirun -n 4 vasp_std > runlog
