## Coprinopsis cinerea
ls | grep ".longest-gene.cds.fa" | grep -v "Coprinopsis_cinerea.longest-gene.cds.fa" | while read i;
do

echo "dnds_calculation.R -q Coprinopsis_cinerea.longest-gene.cds.fa -r $i -b /home/yichun/miniconda3/envs/orthofinder/bin/blastp -m /home/yichun/miniconda3/envs/orthofinder/bin/mafft"

Rscript dnds_calculation.R -q Coprinopsis_cinerea.longest-gene.cds.fa -r $i -b /home/yichun/miniconda3/envs/orthofinder/bin/blastp -m /home/yichun/miniconda3/envs/orthofinder/bin/mafft

done

## Fusarium_graminearum
ls | grep ".longest-gene.cds.fa" | grep -v "Fusarium_graminearum.longest-gene.cds.fa" | while read i;
do

echo "dnds_calculation.R -q Fusarium_graminearum.longest-gene.cds.fa -r $i -b /home/yichun/miniconda3/envs/orthofinder/bin/blastp -m /home/yichun/miniconda3/envs/orthofinder/bin/mafft"

Rscript dnds_calculation.R -q Fusarium_graminearum.longest-gene.cds.fa -r $i -b /home/yichun/miniconda3/envs/orthofinder/bin/blastp -m /home/yichun/miniconda3/envs/orthofinder/bin/mafft

done

## Neurospora_crassa
ls | grep ".longest-gene.cds.fa" | grep -v "Neurospora_crassa.longest-gene.cds.fa" | while read i;
do

echo "dnds_calculation.R -q Neurospora_crassa.longest-gene.cds.fa -r $i -b /home/yichun/miniconda3/envs/orthofinder/bin/blastp -m /home/yichun/miniconda3/envs/orthofinder/bin/mafft"

Rscript dnds_calculation.R -q Neurospora_crassa.longest-gene.cds.fa -r $i -b /home/yichun/miniconda3/envs/orthofinder/bin/blastp -m /home/yichun/miniconda3/envs/orthofinder/bin/mafft

done

## Rhizopus_delemar
ls | grep ".longest-gene.cds.fa" | grep -v "Rhizopus_delemar.longest-gene.cds.fa" | while read i;
do

echo "dnds_calculation.R -q Rhizopus_delemar.longest-gene.cds.fa -r $i -b /home/yichun/miniconda3/envs/orthofinder/bin/blastp -m /home/yichun/miniconda3/envs/orthofinder/bin/mafft"

Rscript dnds_calculation.R -q Rhizopus_delemar.longest-gene.cds.fa -r $i -b /home/yichun/miniconda3/envs/orthofinder/bin/blastp -m /home/yichun/miniconda3/envs/orthofinder/bin/mafft

done
