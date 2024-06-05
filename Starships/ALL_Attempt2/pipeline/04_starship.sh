#!/usr/bin/bash -l
#SBATCH --time 48:00:00 --ntasks 16 --nodes 1 --mem 64G --out logs/starfish_insert.%A.log
mamba activate starfish

starfish insert -T 2 -a ome2assembly.txt -d blastdb/histo_all.assemblies -b geneFinder/Histo.tyr.bed -i tyr -x Histo -o elementFinder/
