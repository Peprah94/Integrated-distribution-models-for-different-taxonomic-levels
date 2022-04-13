#!/bin/sh
#SBATCH --partition=CPUQ
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=120000
#SBATCH --job-name="data_species"
#SBATCH --account=share-ie-imf
 
 
module purge
module load R/4.1.0-foss-2021a
 
/usr/bin/time -v Rscript nimble_species_specific.R
