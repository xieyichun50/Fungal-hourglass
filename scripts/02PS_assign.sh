## This script will determine the phylostratum of all the protein sequences of specific species (e.g. Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.proteins.fa) you feed in.

## input: protein.fa
## edit line "i=Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.proteins.fa" to add your input protein.fa
## edit line "s=5346" to confirm your taxon ID (NCBI)
## file $i.diamond.taxon.tab generated by 01PS_blast.sh

## output1: protein.fa.diamond.taxon.tab
## output2: protein.fa.taxonlist.txt 
## output3: proteins.fa.taxonlist.txt.full.level.txt and proteins.fa.taxonlist.txt.simplePS.level.txt
## output4: proteins.fa.PSfinal.txt (for TAI calculation)

## For the first time you use this script, in R run "prepareDatabase("nameNode.sqlite", getAccessions=FALSE)"

#i=Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.proteins.fa
#s=5346
#i=Fusarium_graminearum.proteins.fa
#s=5518

## The following will try to cut down the memory usage by splitting the $i.diamond.taxon.tab into several smaller files
## "-l N" will split the blast output into smaller files every N lines. Change N to fit with your computer memory, -l 1000000 will use ~1Gb memory when processing each subfile in < 1min.

split -d --additional-suffix $i.diamond.taxon.tab -l 1000000 $i.diamond.taxon.tab
echo "$i.diamond.taxon.tab splitted into n files" 

ls x*$i.diamond.taxon.tab | wc -l 
echo "\n"

ls x*$i.diamond.taxon.tab | while read j;
do
echo "Sorting $j"
cat $j | sort | uniq > $j.sorted
rm $j
done

echo "Generating $i.taxonlist.txt"
cat *$i.diamond.taxon.tab.sorted | cut -f2 | sed 's/;/\n/g' | sort | uniq > $i.taxonlist.txt

echo "Splitting $i.taxonlist.txt into smaller files every 5000 lines."
split -d --additional-suffix $i.taxonlist.txt -l 5000 $i.taxonlist.txt

echo "Mapping taxonomy ID to taxonlist"
ls x*$i.taxonlist.txt | while read j;
do
echo "Rscript assign_taxon_level_to_taxonID.R -i $j -s $s &" >> 02mapping.sh
done
echo "wait" >> 02mapping.sh
sh 02mapping.sh
cat x*$i.taxonlist.txt.simplePS.level.txt | sort | uniq > $i.simplePS.level.txt
cat x*$i.taxonlist.txt.full.level.txt | sort | uniq > $i.full.level.txt
rm x*$i.taxonlist.txt x*$i.taxonlist.txt.simplePS.level.txt x*$i.taxonlist.txt.full.level.txt
rm 02mapping.sh

ls x*$i.diamond.taxon.tab.sorted | while read j;
do 
echo "Rscript assign_taxon_level_to_proteins.R -i $j -p $i.simplePS.level.txt &" >> 03protPS.sh
done
echo "wait" >> 03protPS.sh
sh 03protPS.sh

cat x*$i.diamond.taxon.tab.sorted.PS.txt | sort | uniq > $i.tab.sorted.PS.txt
Rscript assign_taxon_level_to_genes.R -i $i.tab.sorted.PS.txt -o $i.PSfinal.txt
rm x*$i.diamond.taxon.tab.sorted.PS.txt x*$i.diamond.taxon.tab.sorted
rm sh 03protPS.sh $i.taxonlist.txt
echo "DONE! The Phylostrata level of genes saved in $i.PSfinal.txt ."
