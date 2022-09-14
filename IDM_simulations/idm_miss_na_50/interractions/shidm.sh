#!/bin/sh
#SBATCH --partition=CPUQ
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=18000
#SBATCH --job-name=“iidm”
#SBATCH --account=share-ie-imf
 
 
module purge
module load R/4.1.0-foss-2021a
 
/usr/bin/time -v Rscript nimble_interractions.R
