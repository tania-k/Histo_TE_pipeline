cut -f1,12  ann/*emapper.annotations | grep -v  '#' | grep -v -P '\t-' | perl -pe 's/\t/\tEMAP\t/' | grep -vP '\tNA' > ann/histo_all.gene2emap.txt
