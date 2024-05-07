#!/usr/bin/bash -l
#SBATCH -p general -N 1 -n 1 -c 32 --mem 96gb --out logs/panedta%A.log


CPU=2
if [ $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi

OUTDIR=panEDTA_output

mkdir -p $OUTDIR

conda activate EDTA

CDS=$(/nas/longleaf/home/taniak/taniak/TE_Analysis/Histoplasma_capsulatum_G186AR-cleanER.cds.fna)

bash /nas/longleaf/home/taniak/taniak/TE_Analysis/pipeline/07_panEDTA.sh -g genome_list.txt -c $CDS -t $CPU

