#!/bin/bash -l
#SBATCH -o DIMACS.out
#SBATCH --qos=normal
#SBATCH -J DIMACS
#SBATCH --nodes=1

#SBATCH --cpus-per-task=16
#SBATCH -t 24:00:00


module purge
module add matlab


matlab -nodesktop -nodisplay -nosplash -r "run('../startup.m'); dataset = string('DIMACS'); instances = string('all'); test_socp; exit" 
