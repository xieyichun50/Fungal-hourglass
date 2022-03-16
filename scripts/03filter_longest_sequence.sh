## This script will get the longest protein and cds sequence of each gene.

## As the downloaded genomic resources can be different in format, we better check one by one. To get species.CDS file in the following format.
#JAAGWA010000012.1       Genbank CDS     1702694 1702847 0       +       0       gene_id "KAG2003616.1";
#JAAGWA010000012.1       Genbank CDS     1702904 1705027 0       +       2       gene_id "KAG2003616.1";
#JAAGWA010000012.1       Genbank CDS     1705087 1706020 0       +       2       gene_id "KAG2003616.1";
#JAAGWA010000012.1       Genbank CDS     1706104 1707436 0       +       1       gene_id "KAG2003616.1";
#JAAGWA010000003.1       Genbank CDS     1573862 1574125 0       -       0       gene_id "KAG2017906.1";
#JAAGWA010000011.1       Genbank CDS     1789977 1789979 0       +       0       gene_id "KAG2004876.1";
#JAAGWA010000011.1       Genbank CDS     1790280 1790400 0       +       0       gene_id "KAG2004876.1";
#JAAGWA010000011.1       Genbank CDS     1790456 1791008 0       +       2       gene_id "KAG2004876.1";
#JAAGWA010000011.1       Genbank CDS     1791063 1791195 0       +       1       gene_id "KAG2004876.1";
#JAAGWA010000001.1       Genbank CDS     4125116 4125186 0       +       0       gene_id "KAG2024141.1";

##For NCBI genband direct download file
perl gff2gtf.pl $i.genomic.gff | grep -w CDS | sed 's/_id ".*"://' | sed 's/transcript_id/\t/' | cut -f1,2,3,4,5,6,7,8,9  > $i.genomic.gff.CDS
sed 's/lcl.*protein_id=//' $i.cds.fa | sed 's/]/\t/' | cut -f1 > $i.cds-simple.fa
##For JGI  Filtered Models ("best"). *_GeneCatalog_genes_*.gff
cat $i.genomic.gff | grep -w "CDS" | sed 's/name "/ gene_id "/;s/;/\t/' | cut -f1,2,3,4,5,6,7,8,9 > $i.genomic.gff.CDS
cat $i.cds.fa | sed 's/>.*|/>/' > $i.cds-simple.fa
mv $i.protein.fa $i.protein.fa.bk
cat $i.protein.fa.bk | sed 's/>.*|/>/' > $i.protein.fa

##Exception
perl gff2gtf.pl Beauveria_bassiana.genomic.gff | grep -w CDS | sed 's/_id ".*"://' | sed 's/transcript_id/\t/;s/XP/\tgene_id "XP/' | cut -f1,2,3,4,5,6,7,8,10  > Beauveria_bassiana.genomic.gff.CDS
grep -w "CDS" Fusarium_graminearum.genomic.gff | sed 's/Parent=/gene_id "/;s/\;/"\;/' > Fusarium_graminearum.genomic.gff.CDS
perl gff2gtf.pl Trichoderma_reesei.genomic.gff | grep -w CDS | sed 's/_id ".*"://' | sed 's/transcript_id/\t/' | cut -f1,2,3,4,5,6,7,8,9 | grep -v "gene_id \"\"" > Trichoderma_reesei.genomic.gff.CDS

##Get representative sequence
ls *.CDS | sed 's/.genomic.gff.CDS//' | while read i;
do
cgat gtf2gtf -I $i.genomic.gff.CDS  --method=filter --filter-method=longest-gene > $i.longest-gene.gtf
grep -w "gene_id" $i.longest-gene.gtf |sed 's/.*gene_id "//g;s/".*//' |sort -u > $i.longest-gene.gtf.id
seqtk subseq $i.protein.fa  $i.longest-gene.gtf.id > $i.longest-gene.proteins.fa
seqtk subseq $i.cds-simple.fa  $i.longest-gene.gtf.id > $i.longest-gene.cds.fa
wc -l $i*
grep ">" $i.cds-simple.fa | wc -l
grep ">" $i.longest-gene.cds.fa | wc -l
grep ">" $i.protein.fa | wc -l
grep ">" $i.longest-gene.proteins.fa | wc -l
done