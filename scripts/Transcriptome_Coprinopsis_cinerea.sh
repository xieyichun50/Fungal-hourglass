######read mapping
#!/bin/bash
source /PARA/app/scripts/cn-module.sh
module load Bowtie2/2.2.6
module load hisat2/2.0.5
module load samtools/1.5
module load java/jdk1.8.0_11
module load picard/2.9.0
export LD_LIBRARY_PATH=/PARA/pp585/Yichun/lib:$LD_LIBRARY_PATH

#hisat2-build -f Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.scaffolds.fa Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.scaffolds.fa -p 24

cat /PARA/pp585/Yichun/filelist_RNA | while read i;
do
/PARA/pp585/Yichun/tools/fastp -i /PARA/pp585/Yichun/fastq/$i.r1.fq -I /PARA/pp585/Yichun/fastq/$i.r2.fq -o /PARA/pp585/Yichun/fastp/$i.r1.fastp.fq -O /PARA/pp585/Yichun/fastp/$i.r2.fastp.fq -W 4 -M 20 -l 15 -j /PARA/pp585/Yichun/fastp/$i.fastp.json -h /PARA/pp585/Yichun/fastp/$i.fastp.html -R "$i.fastp_report" -w 16
echo "hisat2 -x /PARA/pp585/Yichun/genome/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.scaffolds.fa -1 /PARA/pp585/Yichun/fastp/$i.r1.fastp.fq -2 /PARA/pp585/Yichun/fastp/$i.r2.fastp.fq -S /PARA/pp585/Yichun/bam/$i.sam --phred33 --dta-cufflinks --novel-splicesite-outfile /PARA/pp585/Yichun/bam/splice/$i.splicesite.txt  -p 24"
hisat2 -x /PARA/pp585/Yichun/genome/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.scaffolds.fa -1 /PARA/pp585/Yichun/fastp/$i.r1.fastp.fq -2 /PARA/pp585/Yichun/fastp/$i.r2.fastp.fq -S /PARA/pp585/Yichun/bam/$i.sam --phred33 --dta-cufflinks --novel-splicesite-outfile /PARA/pp585/Yichun/bam/splice/$i.splicesite.txt  -p 24

echo "samtools view -b -S /PARA/pp585/Yichun/bam/$i.sam > /PARA/pp585/Yichun/bam/$i.bam"
samtools view --threads 24 -b -S /PARA/pp585/Yichun/bam/$i.sam > /PARA/pp585/Yichun/bam/$i.bam
echo "/PARA/pp585/Yichun/bam/$i.bam generated"

rm /PARA/pp585/Yichun/bam/$i.sam

echo "samtools sort /PARA/pp585/Yichun/bam/$i.bam -o /PARA/pp585/Yichun/bam/$i.sorted.bam"
samtools sort --threads 24 /PARA/pp585/Yichun/bam/$i.bam -o /PARA/pp585/Yichun/bam/$i.sorted.bam

####Remove duplicate reads (generate by PCR, etc.)

echo "java -Xmx40g -jar /WORK/app/picard/2.9.0/picard.jar MarkDuplicates INPUT=/PARA/pp585/Yichun/bam/$i.sorted.bam OUTPUT=/PARA/pp585/Yichun/redubam/$i.sorted.redu.bam METRICS_FILE=/PARA/pp585/Yichun/redubam/$i.sorted.redu.metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true"

java -Xmx40g -jar /WORK/app/picard/2.9.0/picard.jar MarkDuplicates INPUT=/PARA/pp585/Yichun/bam/$i.sorted.bam OUTPUT=/PARA/pp585/Yichun/redubam/$i.sorted.redu.bam METRICS_FILE=/PARA/pp585/Yichun/redubam/$i.sorted.redu.metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true

samtools index /PARA/pp585/Yichun/redubam/$i.sorted.redu.bam
echo "/PARA/pp585/Yichun/redubam/$i.sorted.redu.bam indexed"
done

####Expression level
gff=/home/yichun/RNAmodification/genome/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.gff3
dir=DEG

cat samples_n_reads_decribed.txt | while read id1 id2;
do  
	#hisat2 -p 8 --dta -x $ref -1 $fq1 -2 $fq2 --rf -S $id2.sam &>$id2.hisat2.log 
	#samtools sort -@ 8 -o $id2.bam $id2.sam 
	#stringtie -p 8 -G $ref.gtf -o $id2.stringtie.gtf -l $id2 $id2.bam
	stringtie -e -B -p 14 -G $gff -A ${id2}_gene_count.xls -o ballgown/${id2}_ballgown/$id2.gtf /home/yichun/RNAmodification/redubam/$id2.sorted.redu.bam
	#rm $id2.sam $id2.bam
done
perl /home/yichun/RNAmodification/stringte_gene_count_matrix.pl *_gene_count.xls

prepDE.py -i ballgown/ -l 150
sed -i 's/,/\t/g;s/_ballgown//g' gene_count_matrix.csv transcript_count_matrix.csv

##file ready
##
Rscript /home/yichun/RNAmodification/edgeR_TMM_norm.R

cd ${dir};
/home/yichun/tools/trinityrnaseq/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix gene_count_matrix.csv --method edgeR --samples_file samples_n_reads_decribed.txt --output edgeR_gene.min_reps2.min_cpm1 --min_reps_min_cpm 2,1 --contrasts sample_pair.txt
cd ${dir}/edgeR_gene.min_reps2.min_cpm1;
/home/yichun/tools/trinityrnaseq/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix ../gene_count_matrix.csv -P 0.05 -C 1 --samples ../samples_n_reads_decribed.txt

find  ballgown/ -name '*.gtf' | while read i;
do 
	cp $i .;
done

perl /home/yichun/RNAmodification/stringte_gtf_transcript_count_matrix.pl $( cut -f 2 samples_n_reads_decribed.txt | while read i; do echo -ne "$i.gtf ";done)
sed -i 's/.gtf//g' transcript_count_matrix.*
