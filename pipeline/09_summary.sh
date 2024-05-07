#!/usr/bin/bash -l
#SBATCH -p batch -N 1 -n 1  --mem 64G --out logs/Sum_EDTA.%a.%A.log

hostname
CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi
if [ -z $SLURM_JOB_ID ]; then
  SLURM_JOB_ID=$$
fi

module load repeatmasker/4.0.7
conda activate R-tools
shopt -s extglob

GENOMEDIR=$(realpath genomes)
INDIR=$(realpath RM_EDTA)
OUTDIR=$(realpath RM_Summarize)
SAMPFILE=strains.csv

mkdir -p $OUTDIR

N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "need to provide a number by --array or cmdline"
        exit
    fi
fi
MAX=$(wc -l $SAMPFILE | awk '{print $1}')
if [ $N -gt $(expr $MAX) ]; then
    MAXSMALL=$(expr $MAX)
    echo "$N is too big, only $MAXSMALL lines in $SAMPFILE"
    exit
fi

IFS=,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read BASE	SPECIES	STRAIN	PHYLUM	LOCUS
do
    name=$(echo -n ${BASE} | perl -p -e 's/\s+/_/g')
    genome=$GENOMEDIR/$name.scaffolds.fa
    twoBit=$GENOMEDIR/$name.scaffolds.2bit
    outfile=$INDIR/$name.RM/$name.scaffolds.fa.out
    if [[ ! -s $twoBit || $genome -nt $twoBit ]]; then
	   /nas/longleaf/home/taniak/.conda/envs/faToTwoBit $genome $twoBit
    fi
    # sum all the contigs to get genome size
    genomeSize=$(twoBitInfo -noNs $twoBit stdout | cut -f2 | paste -s -d+ - | bc)

    Rscript scripts/summarize_outputs.R $outfile $genomeSize
done
