#!/bin/sh

#$ -N run-cellranger-multiple
#$ -cwd
#$ -l h_rt=30:00:00
#$ -l h_vmem=24G
#$ -pe sharedmem 16
#$ -M s1359339@ed.ac.uk   # enter your email address here
#$ -m beas



# command line args
ID=$1

if [-z "$ID"] 
then
 "No id argument supplied."
 "USAGE: ./runCellRangerV2.sh <ID>"
exit 0
fi

# Initialise the environment modules
. /etc/profile.d/modules.sh

module load igmm/apps/cellranger/3.1.0


cellranger count --id=$ID \
                   --transcriptome=/exports/eddie/scratch/s1359339/RefGenome/refdata-cellranger-GRCh38-3.0.0  \
                   --fastqs=/exports/eddie/scratch/s1359339/fastqs/$ID \
		   &> /exports/eddie/scratch/s1359339/run_`date +"%Y-%m-%d-%H-%M-%S"`.log



