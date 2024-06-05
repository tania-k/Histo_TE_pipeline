cut -f1,9  ann/*.emapper.annotations | grep -v  '#' | grep -v -P '\t-' | perl -pe 's/\t/\tEMAP\t/' | grep -vP '\tNA' > ann/Histo_all.gene2emap.txt
cut -f1,5 ann/*.emapper.annotations | grep -v '#' | perl -pe 's/^([^\s]+?)\t([^\|]+).+$/\1\t\2/' > ann/Histo_all.gene2og.txt
perl ../starfish/aux/geneOG2mclFormat.pl -i Histo_all.gene2og.txt -o ann/
