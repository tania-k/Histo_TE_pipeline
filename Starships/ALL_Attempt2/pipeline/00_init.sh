conda activate starfish
realpath assembly/* | perl -pe 's/^(.+?([^\/]+?).fna)$/\2\t\1/' > ome2assembly.txt
realpath gff3/* | perl -pe 's/^(.+?([^\/]+?).gff3)$/\2\t\1/' > ome2gff.txt
cat gff3/*.gff3 > histo_all.gff3
mkdir blastdb
cut -f2 ome2assembly.txt | xargs cat > blastdb/histo_all.assemblies.fna
makeblastdb -in blastdb/histo_all.assemblies.fna -out blastdb/histo_all.assemblies -parse_seqids -dbtype nucl
