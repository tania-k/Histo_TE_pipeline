#!/usr/bin/bash -l
#SBATCH --time 96:00:00 --ntasks 16 --nodes 1 --mem 64G --out logs/eggnog.%a.log

conda activate funannotate

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

INDIR=assembly
OUTDIR=ann/eggnog
mkdir -p $OUTDIR
SAMPFILE=strains_copy.csv
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "need to provide a number by --array or cmdline"
        exit
    fi
fi
MAX=`wc -l $SAMPFILE | awk '{print $1}'`

if [ $N -gt $MAX ]; then
    echo "$N is too big, only $MAX lines in $SAMPFILE"
    exit
fi

export FUNANNOTATE_DB=/nas/longleaf/home/taniak/.conda/envs/funannotate/

IFS=,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read SPECIES STRAIN PHYLUM LOCUSTAG
do
	emapper.py -m diamond --itype CDS -i assembly/$SPECIES.scaffolds.fa -o $OUTDIR/$SPECIES --resume
done
