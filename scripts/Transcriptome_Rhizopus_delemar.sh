cat samples_n_reads_decribed.txt | while read id1 id2 SRR ;
do
/tools/sratoolkit.3.0.0-ubuntu64/bin/fasterq-dump --split-3 -e 38 -m 40GB $SRR
/tools/fastp -i ${SRR}_1.fastq -I ${SRR}_2.fastq -o $id2.r1.fastp.fastq -O $id2.r2.fastp.fastq -W 4 -M 20 -l 15 -j $id2.fastp.json -h $id2.fastp.html -R "$id2.fastp_report" -w 38
hisat2 -x /store/jelly/yichun/ccin/rdel/genome/Rhizopus_delemar.genomic.fa -1 $id2.r1.fastp.fastq -2 $id2.r2.fastp.fastq -S $id2.sam --phred33 --dta-cufflinks --novel-splicesite-outfile $id2.splicesite.txt  -p 38
samtools view --threads 38 -b -S $id2.sam > $id2.bam
samtools sort --threads 38 $id2.bam -o $id2.sorted.bam
samtools index $id2.sorted.bam
java -Xmx40g -jar /tools/gatk/picard.jar MarkDuplicates INPUT=${id2}.sorted.bam OUTPUT=${id2}.sorted.redu.bam METRICS_FILE=${id2}.sorted.redu.metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true
samtools index ${id2}.sorted.redu.bam
done

gff=/store/jelly/yichun/ccin/rdel/genome/Rhizopus_delemar.genomic.gff
dir=DEG
cat samples_n_reads_decribed.txt | while read id1 id2 SRR;
do
stringtie -e -B -p 38 -G $gff -A ${id2}_gene_count.xls -o ballgown/${id2}_ballgown/$id2.gtf $id2.sorted.redu.bam
done
perl /data/transcriptome/scripts/stringte_gene_count_matrix.pl *_gene_count.xls

prepDE.py -i ballgown/ -l 100
sed -i 's/,/\t/g;s/_ballgown//g' gene_count_matrix.csv transcript_count_matrix.csv

Rscript edgeR_TMM_norm.R
