#!/usr/bin/bash -l
#SBATCH --time 48:00:00 --ntasks 4 --nodes 1 --mem 120G --out logs/starfish_visual.log
mamba activate starfish

starfish genome-viz -m all -a ome2assembly.txt -b elementFinder/Histo.elements.bed -o genomeViz/
