#!/bin/bash
# Script to run multiple simulations
# Usage, on Peregrine:
#
#   sbatch ./run_multiple_janas.sh
#
# Peregrine directives:
#SBATCH --time=1000
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --job-name=run_multiple_janas
#SBATCH --output=multiple_janas.log

module load foss/2020a
g++ -Wall -o CompSource individual.cpp Source.cpp -std=c++17

  for i in $(seq 1 30)
do
 sbatch jana.sh $i
done


