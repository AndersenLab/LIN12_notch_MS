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
nucmer --maxgap=500 --mincluster=100 --prefix=JUtQX --coords /projects/b1059/projects/Nicolas/c.briggsae/gene_predictions/genomes/QX1410/decompressed/QX1410.genome.v2.1.fa caenorhabditis_nigoni.PRJNA384657.WBPS19.genomic.fa

# get coordinate file
show-coords -r -l -T QXtQX.delta | awk '$5 > 1000' > JUtQX_transformed.tsv 
