conda activate starfish
../aux/seq-gc.sh -Nbw 1000 blastdb/histo_all.assemblies.fna > histo_all.assemblies.gcContent_w1000.bed
rm blastdb/histo_all.assemblies.fna
