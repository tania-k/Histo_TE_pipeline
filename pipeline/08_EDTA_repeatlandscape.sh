#!/usr/bin/bash -l
#SBATCH -p batch -N 1 -n 1  --mem 64G --out logs/RL_EDTA.%a.%A.SIL.log

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
OUTDIR=$(realpath RM_EDTA_plots)
SAMPFILE=strains.csv

mkdir -p $OUTDIR

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
    if [[ ! -s $twoBit || $genome -nt $twoBit ]]; then
	   /nas/longleaf/home/taniak/.conda/envs/faToTwoBit $genome $twoBit
    fi
    # sum all the contigs to get genome size
    genomeSize=$(twoBitInfo -noNs $twoBit stdout | cut -f2 | paste -s -d+ - | bc)

    RUNDIR=$INDIR/${name}.RM
    if [ ! -d $RUNDIR ]; then
	echo "No repeatmasker dir $RUNDIR"
    fi
    TMPDIV=$RUNDIR/${name}.scaffolds.fa.align
    DIVFILE=""
    # test to see if we have either an uncompressed version of the file or compressed ones
    for div in $TMPDIV $TMPDIV.gz $TMPDIV.bz2
    do
	if [ -f $div ]; then
	    DIVFILE=$div
	    break
	fi
    done
    if [ ! -z $DIVFILE ]; then
	echo "using divfile $DIVFILE"
    else 
	echo "Cannot find a divfile in $RUNDIR ($TMPDIV)"
	exit
    fi
    DIVSUM=$OUTDIR/${name}.divsum
    if [[ ! -s $DIVSUM || $DIVFILE -nt $DIVSUM ]]; then
	/nas/longleaf/home/taniak/taniak/TE_Analysis/scripts/calcDivergenceFromAlign.pl -s $DIVSUM $DIVFILE
    fi
    ./scripts/createReapeatLandscape_EDTA.pl -div $DIVSUM -twoBit $twoBit > $OUTDIR/${name}.landscape.html

    tail -n 72 $DIVSUM > $DIVSUM.tbl
    Rscript scripts/plot_kimuradist_TE.R $DIVSUM.tbl $genomeSize

done
