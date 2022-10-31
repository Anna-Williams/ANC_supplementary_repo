#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N run-velocyto-multiple
#$ -cwd
#$ -l h_rt=15:00:00
#$ -l h_vmem=16G
#$ -pe sharedmem 4
#$ -M s1359339@ed.ac.uk
#$ -m beas

# command line args
ID=$1

if [-z "$ID"]
then
 "No id argument supplied."
exit 0
fi

# Initialise the environment modules
. /etc/profile.d/modules.sh

# Load Python
module load roslin/python/3.6.8
module load roslin/samtools/1.9


# Run the program

velocyto run -b /exports/eddie/scratch/s1359339/$ID/outs/filtered_feature_bc_matrix/barcodes.tsv -o /exports/eddie/scratch/s1359339/veloOuts/$ID -m /exports/eddie/scratch/s1359339/Training/hs38_rmsk.gtf /exports/eddie/scratch/s1359339/$ID/outs/possorted_genome_bam.bam /exports/eddie/scratch/s1359339/RefGenome/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf


