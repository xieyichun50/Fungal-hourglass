## Germination Wang 2019
cat samples_n_reads_decribed.txt | while read id1 id2 SRR ;
do
/tools/sratoolkit.3.0.0-ubuntu64/bin/fasterq-dump --split-3 -e 38 -m 40GB $SRR
echo "/tools/fastp -i ${SRR}_1.fastq -I ${SRR}_2.fastq -o $id.r1.fastp.fastq -O $id.r2.fastp.fastq -W 4 -M 20 -l 15 -j $id.fastp.json -h $id.fastp.html -R "$id.fastp_report" -w 38"
/tools/fastp -i ${SRR}_1.fastq -I ${SRR}_2.fastq -o $id.r1.fastp.fastq -O $id.r2.fastp.fastq -W 4 -M 20 -l 15 -j $id.fastp.json -h $id.fastp.html -R "$id.fastp_report" -w 38
echo "hisat2 -x /store/jelly/yichun/ccin/ncra/genome/Neurospora_crassa.genomic.fa -1 $id.r1.fastp.fastq -2 $id.r2.fastp.fastq -S $id.sam --fr --phred33 --dta-cufflinks --novel-splicesite-outfile $id.splicesite.txt -p 38"
hisat2 -x Neurospora_crassa.genomic.fa -1 $id.r1.fastp.fastq -2 $id.r2.fastp.fastq -S $id.sam --fr --phred33 --dta-cufflinks --novel-splicesite-outfile $id.splicesite.txt -p 38
samtools view --threads 38 -b -S $id.sam > $id.bam
echo "samtools sort --threads 38 $id.bam -o $id.sorted.bam"
samtools sort --threads 38 $id.bam -o $id.sorted.bam
samtools index $id.sorted.bam
done

gff=Neurospora_crassa.genomic.gff
dir=DEG
cat samples_n_reads_decribed.txt | while read id1 id2 SRR;
do
stringtie -e -B --fr -p 38 -G $gff -A ${id2}_gene_count.xls -o ballgown/${id2}_ballgown/$id2.gtf $id2.sorted.bam
done
prepDE.py -i germ/ -l 75

## Sex Liu 2017
cat samples_n_reads_decribed.txt | while read id1 id2 SRR ;
do
/tools/sratoolkit.3.0.0-ubuntu64/bin/fasterq-dump --split-3 -e 38 -m 40GB $SRR
echo "/tools/fastp -i ${SRR}_1.fastq -I ${SRR}_2.fastq -o $id.r1.fastp.fastq -O $id.r2.fastp.fastq -W 4 -M 20 -l 15 -j $id.fastp.json -h $id.fastp.html -R "$id.fastp_report" -w 38"
/tools/fastp -i ${SRR}_1.fastq -I ${SRR}_2.fastq -o $id.r1.fastp.fastq -O $id.r2.fastp.fastq -W 4 -M 20 -l 15 -j $id.fastp.json -h $id.fastp.html -R "$id.fastp_report" -w 38
echo "hisat2 -x Neurospora_crassa.genomic.fa -1 $id.r1.fastp.fastq -2 $id.r2.fastp.fastq -S $id.sam --fr --phred33 --dta-cufflinks --novel-splicesite-outfile $id.splicesite.txt -p 38"
hisat2 -x Neurospora_crassa.genomic.fa -1 $id.r1.fastp.fastq -2 $id.r2.fastp.fastq -S $id.sam --fr --phred33 --dta-cufflinks --novel-splicesite-outfile $id.splicesite.txt -p 38
samtools view --threads 38 -b -S $id.sam > $id.bam
echo "samtools sort --threads 38 $id.bam -o $id.sorted.bam"
samtools sort --threads 38 $id.bam -o $id.sorted.bam
samtools index $id.sorted.bam
done

gff=Neurospora_crassa.genomic.gff
dir=DEG
cat samples_n_reads_decribed.txt | while read id1 id2 SRR;
do
stringtie -e -B --fr -p 38 -G $gff -A ${id2}_gene_count.xls -o ballgown/${id2}_ballgown/$id2.gtf $id2.sorted.bam
done
prepDE.py -i sexL/ -l 150

## Sex Wang 2012
cat samples_n_reads_decribed.txt | while read id1 id2 SRR ;
do
/tools/fastp -i $SRR.fastq -o $id2.fastp.fastq -W 4 -M 20 -l 15 -j $id2.fastp.json -h $id2.fastp.html -R "$id2.fastp_report" -w 38
hisat2 -x Neurospora_crassa.genomic.fa -U $id2.fastp.fastq -S $id2.sam --phred33 --dta-cufflinks --novel-splicesite-outfile $id2.splicesite.txt  -p 38
samtools view --threads 38 -b -S $id2.sam > $id2.bam
samtools sort --threads 38 $id2.bam -o $id2.sorted.bam
samtools index $id2.sorted.bam
done

gff=Neurospora_crassa.genomic.gff
dir=DEG
cat samples_n_reads_decribed.txt | while read id1 id2 SRR;
do  
stringtie -e -B -p 38 -G $gff -A ${id2}_gene_count.xls -o ballgown/${id2}_ballgown/$id2.gtf $id2.sorted.bam
done
perl /home/yichun/RNAmodification/stringte_gene_count_matrix.pl *_gene_count.xls
prepDE.py -i sexW/ -l 36

## Vegetative JGI 2016
cat samples_n_reads_decribed.txt | while read id1 id2 SRR ;
do
/tools/fastp -i $SRR.fastq -o $id2.fastp.fastq -W 4 -M 20 -l 15 -j $id2.fastp.json -h $id2.fastp.html -R "$id2.fastp_report" -w 38
hisat2 -x Neurospora_crassa.genomic.fa -U $id2.fastp.fastq -S $id2.sam --phred33 --dta-cufflinks --novel-splicesite-outfile $id2.splicesite.txt  -p 38
samtools view --threads 38 -b -S $id2.sam > $id2.bam
samtools sort --threads 38 $id2.bam -o $id2.sorted.bam
samtools index $id2.sorted.bam
done

gff=Neurospora_crassa.genomic.gff
dir=DEG
cat samples_n_reads_decribed.txt | while read id1 id2 SRR;
do  
stringtie -e -B -p 38 -G $gff -A ${id2}_gene_count.xls -o ballgown/${id2}_ballgown/$id2.gtf $id2.sorted.bam
done
perl /home/yichun/RNAmodification/stringte_gene_count_matrix.pl *_gene_count.xls
prepDE.py -i veg/ -l 100

#In R
{
sexcount<-read.csv("expression/ballgown/sex.gene_count_matrix.csv", header = T)
germcount<-read.csv("expression/ballgown/germ.gene_count_matrix.csv", header = T)
rawcount<-merge(germcount,sexcount, by = "gene_id", all = T)
names(rawcount)<-gsub("_ballgown", "", names(rawcount))
write.table(rawcount,"gene_count_matrix.csv", sep = "\t", quote = F, row.names = F)
sexcount<-read.csv("expression/ballgown/sex.transcript_count_matrix.csv", header = T)
germcount<-read.csv("expression/ballgown/germ.transcript_count_matrix.csv", header = T)
rawcount<-merge(germcount,sexcount, by = "transcript_id", all = T)
names(rawcount)<-gsub("_ballgown", "", names(rawcount))
write.table(rawcount,"transcript_count_matrix.csv", sep = "\t", quote = F, row.names = F)
}

Rscript edgeR_TMM_norm.R
