#!/usr/bin/bash -l
#SBATCH --time 76:00:00 --ntasks 16 --nodes 1 --mem 64G --out logs/starfish_annotate.%A.log
mamba activate starfish

starfish annotate -T 2 -x Histo_tyr -a ome2assembly.txt -g ome2gff.txt -p ../starfish/db/YRsuperfams.p1-512.hmm -P ../starfish/db/YRsuperfamRefs.faa -i tyr -o geneFinder --namefield 'Parent='
