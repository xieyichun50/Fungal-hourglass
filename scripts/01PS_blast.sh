## This script will find ortholog against all available NCBI nr records of all the protein sequences of specific species (e.g. Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.proteins.fa) you feed in.

## input: protein.fa
## edit line "i=Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.proteins.fa" to add your input protein.fa
## edit line "s=5346" to confirm your taxon ID (NCBI)

## output1: protein.fa.diamond.taxon.tab

## To accelerate the search, Diamond blastp will be used instead of ncbi-blast*/blastp. Makesure you have diamond well installed. In addition, R package "taxonomizr" will be used for taxon mapping.
## For the first time you use this script, download nr from NCBI and makedb in Diamond format (for example, I 'makedb' under /tools/diamond/, unmask the following and run)

#wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
#wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/swissprot.gz
#wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
#wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
#rm -rf *.dmp
#unzip taxdmp.zip
#mv nr.dmnd nr.dmnd.$$
#diamond makedb --in nr.gz -d nr --taxonmap prot.accession2taxid.gz  --taxonnodes nodes.dmp --taxonnames names.dmp
#diamond makedb --in swissprot.gz -d swissprot --taxonmap prot.accession2taxid.gz  --taxonnodes nodes.dmp --taxonnames names.dmp

#i=Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.proteins.fa
#s=5346
#i=Fusarium_graminearum.proteins.fa
#s=5518

## In diamond, -k0 will blast against all available taxon

## If you want further filtering according to blast result, unmask the following 
#/tools/diamond/diamond blastp -d /tools/diamond/nr --evalue 1e-3 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids -k0 --threads 240 --query $i --out $i.diamond.tab
#awk -F$'\t' '($4 > 30){print}' $i.diamond.tab | cut -f1,13 > $i.diamond.taxon.tab

## No further filtering
/tools/diamond/diamond blastp -d /tools/diamond/nr --evalue 1e-3 --outfmt 6 qseqid staxids -k0 --threads 240 --query $i --out $i.diamond.taxon.tab

echo "01PS_blast.sh finished!"
echo "Next step, run 02PS_assign.sh to get phylostrata assign to all genes."
