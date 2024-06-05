#!/usr/bin/bash -l
#SBATCH --time 48:00:00 --ntasks 16 --nodes 1 --mem 64G --out logs/starfish_annotate.%A.log
mamba activate starfish

starfish annotate -T 2 -x Histo_tyr -a ome2assembly.txt -g ome2gff.txt -p starfish/db/YRsuperfams.p1-512.hmm -P starfish/db/YRsuperfamRefs.faa -i tyr -o geneFinder --namefield 'Parent='
starfish consolidate -o ./ -g histo_all.gff3 -G geneFinder/Histo_tyr.filt_intersect.gff --namefield 'Parent='

realpath Histo_tyr.filt_intersect.consolidated.gff | perl -pe 's/^/Histo\t/' > ome2consolidatedGFF.txt
starfish sketch -m 10000 -q geneFinder/Histo_tyr.filt_intersect.ids -g ome2consolidatedGFF.txt -i s -x Histo -o geneFinder/ --nameField 'Parent='
grep -P '\ttyr\t' geneFinder/Histo.bed > geneFinder/Histo.tyr.bed
starfish insert -T 2 -a ome2assembly.txt -d blastdb/histo_all.assemblies -b geneFinder/Histo.tyr.bed -i tyr -x Histo -o elementFinder/
starfish flank -a ome2assembly.txt -b elementFinder/Histo.insert.bed -x Histo -o elementFinder/
starfish summarize -a ome2assembly.txt -b elementFinder/Histo.flank.bed -x Histo -o elementFinder/ -S elementFinder/Histo.insert.stats -f elementFinder/Histo.flank.singleDR.stats -g ome2consolidatedGFF.txt -A ann/Histo_all.gene2emap.txt -t geneFinder/Histo_tyr.filt_intersect.ids
starfish pair-viz -m all -t empty -T 2 -A nucmer -a ome2assembly.txt -b elementFinder/Histo.elements.bed -f elementFinder/Histo.flank.singleDR.stats -S elementFinder/Histo.elements.named.stats -o pairViz/
mmseqs easy-cluster geneFinder/Histo_tyr.filt_intersect.fas regionFinder/Histo_tyr regionFinder/ --threads 2 --min-seq-id 0.5 -c 0.25 --alignment-mode 3 --cov-mode 0 --cluster-reassign
perl ../starfish/aux/mmseqs2mclFormat.pl -i regionFinder/Histo_tyr_cluster.tsv -g fam -o regionFinder/
starfish sim -m element -t nucl -b elementFinder/Histo.elements.bed -x Histo -o regionFinder/ -a ome2assembly.txt
starfish group -m mcl -s regionFinder/Histo.element.nucl.sim -i hap -o regionFinder/ -t 0.05
grep -P '\tcap\t' elementFinder/Histo.elements.bed | cut -f4,7 > regionFinder/Histo.cap2ship.txt
../starfish/aux/searchReplace.pl -i regionFinder/Histo_tyr_cluster.mcl -r regionFinder/Histo.cap2ship.txt > regionFinder/Histo.element_cluster.mcl
../starfish/aux/mergeGroupfiles.pl -t regionFinder/Histo.element_cluster.mcl -q regionFinder/Histo.element.nucl.I1.5.mcl > regionFinder/Histo.element.navis-hap.mcl
grep -f <(comm -23 <(cut -f1 geneFinder/Histo_tyr.filt_intersect.ids | sort) <(grep -P '\tcap\t|\ttyr\t' elementFinder/Histo.elements.bed | cut -f4| sort)) geneFinder/Histo.tyr.bed > regionFinder/unaffiliated_tyrs.bed
../starfish/aux/filterOG.pl -O ann/Histo_all.gene2og.mcl -a 1 -c 5 -o ann/
