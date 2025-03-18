#!/bin/sh
#An example for serial job.
#SBATCH -J Nxzzeng_ts 
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

export PYTHONPATH=/data/home/xzzeng/soft/ambit-master/build/src:/data/home/xzzeng/soft/psi4/objdir/stage/lib:/data/home/xzzeng/soft/FORTE/forte-master

export PATH=/data/home/xzzeng/soft/ambit-master/include:/data/home/xzzeng/soft/ambit-master/build:/data/home/xzzeng/soft/ambit-master/build/src/:/data/home/xzzeng/soft/psi4/objdir/stage/bin:/data/home/xzzeng/soft/psi4/objdir/stage/lib:/data/home/xzzeng/.conda/envs/dsrg/bin:/data/soft/anaconda3/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin

psi4 < input.dat

#python main.py

#mpirun -n 4 vasp_std > runlog
