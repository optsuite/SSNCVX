#!/bin/bash -l
#SBATCH -o mittelmann.out
#SBATCH --qos=normal
#SBATCH -J mittelmann
#SBATCH --nodes=1
#SBATCH --nodelist=c04
#SBATCH --cpus-per-task=16
#SBATCH -t 24:00:00


module purge
module add matlab



matlab -nodesktop -nodisplay -nosplash -r "run('../startup.m'); dataset = string('mittelmann'); instances = string('all'); test_sdp; exit" 
