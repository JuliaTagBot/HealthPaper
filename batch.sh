#!/bin/bash -l

#SBATCH -A snic2018-8-254
#SBATCH  -p node -n 20
#SBATCH -t 1:00:00
#SBATCH -J Medicare
#SBATCH --mail-type=ALL
#SBATCH --mail-user=panagiotis.margaris@hhs.se
source $HOME/bin/julia-1.0.1/julia-environment
cd /home/pmarg/Projects/Medicare/src
julia -L HealthPaper.jl
echo Job Completed
