## Assign Phylostratigraphy 
- step 1: run `sh 01PS_blast.sh` to find orthologs in nr database in all lineage
- step 2: run `sh 02PS_assign.sh`
#### Requirement:
- software: [diamond](https://github.com/bbuchfink/diamond), R
- R packages: [taxonomizr](https://cran.r-project.org/web/packages/taxonomizr/index.html), optparse, tidyr, dplyr
- files: `protein_sequence.fa`, `*GenematchID`(if the protein ID differ from gene/transcript ID)

## Calculate dNdS
- step 1: run `sh 03filter_longest_sequence.sh` to get the representative sequence.
- step 2: run `sh 04NS_assign.sh`
#### Requirement:
- software: [cgat](https://github.com/cgat-developers/cgat-apps), [seqtk](https://github.com/lh3/seqtk), R
- R packages: [orthologr](https://github.com/drostlab/orthologr), optparse, tidyr, dplyr
- files: `cds.fa`, `*GenematchID`(if the cds ID differ from gene/transcript ID)

## Transcriptome profiling
- run `Transcriptome_*.sh`

## Phylotranscriptomic pattern in fungi
- run `hourglass_*.R`

## Functional annotation on protein sequences
- run `sh function_anno_eggnog.sh`
#### Requirement:
- software: [eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper)
- files: `protein_sequence.fa`

## Cross functional comparison on animal, plant, and fungi
- run `KOG_animal_plant_fungi.R`
- run `waist_animal_plant_fungi_enrich.R` (high-low-high or low-high-low gene expression)
