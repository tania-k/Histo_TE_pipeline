# Transposable Element Pipeline
## How to run TE Analysis using this pipeline
1. Build the proper environment to run your analysis. You might need to create separate environments as not all are compatable with each other.
 You will need: 

```
* Annotated genomes
* RepeatMasker v.4.0+
* RepeatModeler v.2.0+
* bcbio - to get python bio in your environment
* MAFFT v.7.0+
* Clipkit v.2.0+
* R (I used v.4.3.2)
* EDTA environment v2.1.3
* strains.csv file (observe my copy to establish your own)
```

2. Use repeatmodeler pipeline/00_repeatmodeler.sh and build repeat libraries for each genome. Scripts are written to run iteratively with variables pulled from your samples.csv file.

3. Follow scripts in order set up in pipeline.
```
* 01_repeatmasker.sh
* 02_EDTA.sh
* 03_repeatmasker_EDTA.sh
* EDTA.pl
```

What makes my pipeline a little unusual when compared to other TE pipelines, is the use of Repeatmodeler/RepeatMasker + EDTA + Repeatmasker.
There's reason for this madness, mostly I noticed there were DNA TEs missing from RepeatMasker runs vs RNA TEs missing from EDTA runs. 
Using the gff3 TE library generated from one to then be used as the input in the other helped establish a more full analysis. 
Many choose one or the other program, many leaning to EDTA as it's easier.

This is up to your descretion.

4. Once analysis and summaries are complete, look through and run `scripts/analyze_TE_tree_phylo_signal.R` on your RStudio (or however you analyze your data)
This document is my working horse, where I hard coded a lot of work. I apologize for not making this more streamlined but the gist of my analysis is here.
I will also send you my tree building repository, but if you have your method of tree generation and tree output, please use your method instead as your hypotheses questions may differ.
Inside you will find my TE pulling/adding script which helps identify, add, and calculate the percentage of each TE against the genome length for each isolate.
```
* First calculate TEs, collect all varieties per isolate, combine into large table
* I worked the final table on my computer, had the script read that in to do final TE analysis.
* Trees are built here too, but had two CDS and AA trees, went with the AA tree. Add your own tree here
* Created trees using Aspergillus fumigatus, Emmonsia and Blastomyces as the root species. Final tree uses Blastomyces as root.
* Also added Phylogenetic signal to each rooted species. Make sure to run the phylogenetic signal analysis for each TE type.
* I had 7 groups done, RNA TEs, DNA TEs, Simple Repeats, Low Complexity, MITEs, Satellites, and Unknowns. I've seen many disregard Low complexity, simple repeats and unknowns. Up to your discretion. 
```
Good luck! Please message me at taniakurbessoian @ gmail for questions on results.

-T
