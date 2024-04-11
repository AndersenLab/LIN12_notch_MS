#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --job-name="mummer"

#activate environment
source activate nucmer

# align with nucmer (will spit out a .delta file)
nucmer --maxgap=500 --prefix=CB_comp2 --coords N2.genome.fa ../wormbase/WS280/c_briggsae.PRJNA10731.WS280.genomic.fa

# get coordinate file
show-coords -r -l -T CB_comp2.delta | awk '$5 > 1000' > N2_transformed2 
