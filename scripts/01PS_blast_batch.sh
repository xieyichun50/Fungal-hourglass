cut -f1 links.txt | while read i;
do
echo "diamond blastp -d /mnt/content_176/yichun/tools/diamond_nr/nr --evalue 1e-5 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids -k0 --threads 76 --query $i.protein.fa --out $i.protein.fa.diamond.tab"

done

cut -f1 links.txt | while read i;
do
echo "awk -F\$'\\\t' '(\$4 > 30){print}' $i.protein.fa.diamond.tab | cut -f1,13 > $i.protein.fa.diamond.taxon.tab"
done
