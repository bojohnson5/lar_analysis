#!/bin/bash
#SBATCH -N 1
#SBATCH -J lar_anlys
#SBATCH -p general
#SBATCH -A r00891
#SBATCH --mem-per-cpu=15G
#SBATCH --time=10:00:00
#SBATCH -o %j.o
#SBATCH -e %j.e

date

cd ~/lar_anlys/notebooks
source ~/cohar750_sim/setenv.sh
./testing

date
