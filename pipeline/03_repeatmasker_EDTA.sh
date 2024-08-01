#!/usr/bin/bash -l
#SBATCH -p general -n 32 --mem 64gb --out logs/RM.%a.NEW.log -N 1

module load repeatmasker/
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

SAMPFILE=strains.csv

IFS=,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read SPECIES STRAIN PHYLUM LOCUS
do
 OUTDIR=RM_EDTA
 GENOME=genomes/$SPECIES.scaffolds.fa
 library=/nas/longleaf/home/taniak/taniak/TE_Analysis/EDTA_output/$SPECIES/$SPECIES.scaffolds.fa.mod.EDTA.TElib.fa

 mkdir -p $OUTDIR/$SPECIES.RM
if [ ! -f $OUTDIR/$SPECIES.RM/$SPECIES.sorted.fasta.tbl ]; then
 RepeatMasker -lib $library -a -s -e rmblast -dir $OUTDIR/$SPECIES.RM $GENOME > $OUTDIR/$SPECIES.RM/$SPECIES.RepeatMasker.out
else
    echo "Skipping ${name} as masked already"
fi
done
