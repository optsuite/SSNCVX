#!/bin/bash -l
#SBATCH -o qap.out
#SBATCH --qos=normal
#SBATCH -J qap
#SBATCH --nodes=1

#SBATCH --cpus-per-task=16
#SBATCH -t 24:00:00


module purge
module add matlab


matlab -nodesktop -nodisplay -nosplash -r "run('../startup.m'); dataset = string('qap'); instances = string('all'); test_sdp; exit" 
