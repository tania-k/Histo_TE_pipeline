#!/usr/bin/bash -l
#SBATCH -p general -n 32 --mem 64gb --out logs/OneCode.EDTA.RM.2.%a.log -N 1

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
 INDIR=RM_EDTA
 OUTDIR=OneCode_EDTA_RM
 RMFILE=$INDIR/$SPECIES.RM/$SPECIES.scaffolds.fa.out

 mkdir -p $OUTDIR/$SPECIES.OC
 Onecodetofindthemall/build_dictionary.pl --rm $RMFILE --unknown > $OUTDIR/$SPECIES.OC/$SPECIES.LTR-dictionary_output.txt
 Onecodetofindthemall/one_code_to_find_them_all.pl --rm $RMFILE --unknown --ltr $OUTDIR/$SPECIES.OC/$SPECIES.LTR-dictionary_output.txt
 mv $INDIR/$SPECIES.RM/*.csv $OUTDIR/$SPECIES.OC
 mv $INDIR/$SPECIES.RM/*.length $OUTDIR/$SPECIES.OC
 mv $INDIR/$SPECIES.RM/*.out.log.txt $OUTDIR/$SPECIES.OC
done
