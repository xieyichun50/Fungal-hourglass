## Coprinopsis cinerea
ls | grep ".longest-gene.cds.fa" | grep -v "Coprinopsis_cinerea.longest-gene.cds.fa" | while read i;
do

echo "/store/jelly/yichun/ccin/dnds_calculation.R -q Coprinopsis_cinerea.longest-gene.cds.fa -r $i -b /home/yichun/miniconda3/envs/orthofinder/bin/blastp -m /home/yichun/miniconda3/envs/orthofinder/bin/mafft"

Rscript /store/jelly/yichun/ccin/dnds_calculation.R -q Coprinopsis_cinerea.longest-gene.cds.fa -r $i -b /home/yichun/miniconda3/envs/orthofinder/bin/blastp -m /home/yichun/miniconda3/envs/orthofinder/bin/mafft

done

## Fusarium_graminearum
ls | grep ".longest-gene.cds.fa" | grep -v "Fusarium_graminearum.longest-gene.cds.fa" | while read i;
do

echo "/store/jelly/yichun/ccin/dnds_calculation.R -q Fusarium_graminearum.longest-gene.cds.fa -r $i -b /home/yichun/miniconda3/envs/orthofinder/bin/blastp -m /home/yichun/miniconda3/envs/orthofinder/bin/mafft"

Rscript /store/jelly/yichun/ccin/dnds_calculation.R -q Fusarium_graminearum.longest-gene.cds.fa -r $i -b /home/yichun/miniconda3/envs/orthofinder/bin/blastp -m /home/yichun/miniconda3/envs/orthofinder/bin/mafft

done