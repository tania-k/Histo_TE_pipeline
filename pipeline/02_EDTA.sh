#!/bin/bash -l
#SBATCH -N 1 -n 48 --mem 96gb --time 32:00:00 --out logs/EDTA.%a.NEW.log

if [[ -z ${SLURM_CPUS_ON_NODE} ]]; then
    CPUS=1
else
    CPUS=${SLURM_CPUS_ON_NODE}
fi

N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "need to provide a number by --array or cmdline"
        exit
    fi
fi

INDIR=genomes
OUTDIR=EDTA_output
SAMPLES=strains.csv

IFS=,
tail -n +2 $SAMPLES | sed -n ${N}p | while read BASE	SPECIES	STRAIN	PHYLUM	LOCUS
do

mkdir -p $OUTDIR/$BASE
pushd $OUTDIR/$BASE
GENOME=$(realpath ../../$INDIR/$BASE.scaffolds.fa)
CDS=$(realpath ../../Histoplasma_capsulatum_G186A.cds.CLEAN.fna)
conda activate EDTA

EDTA.pl --genome $GENOME --cds $CDS --species others --step all --sensitive 1 --force 1 --anno 1 --evaluate 1
done
