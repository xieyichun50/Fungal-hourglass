wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/240/135/GCF_000240135.3_ASM24013v3/GCF_000240135.3_ASM24013v3_genomic.fna.gz -O Fusarium_graminearum.genomic.fa.gz

####Mapping
hisat2-build -p 39 -f Fusarium_graminearum.genomic.fa Fusarium_graminearum.genomic.fa
cat samples_n_reads_decribed.txt | while read id1 id2 SRR ;
do
/tools/sratoolkit.3.0.0-ubuntu64/bin/fasterq-dump --split-3 -e 38 -m 40GB $SRR
/tools/fastp -i $SRR.fastq -o $id2.fastp.fastq -W 4 -M 20 -l 15 -j $id2.fastp.json -h $id2.fastp.html -R "$id2.fastp_report" -w 38
hisat2 -x /store/jelly/yichun/ccin/fgra/genome/Fusarium_graminearum.genomic.fa -U $id2.fastp.fastq -S $id2.sam --phred33 --dta-cufflinks --novel-splicesite-outfile $id2.splicesite.txt  -p 38
samtools view --threads 38 -b -S $id2.sam > $id2.bam
samtools sort --threads 38 $id2.bam -o $id2.sorted.bam
samtools index $id2.sorted.bam
#java -Xmx40g -jar /tools/gatk/picard.jar MarkDuplicates INPUT=$id.sorted.bam OUTPUT=$id.sorted.redu.bam METRICS_FILE=$id.sorted.redu.metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true
#samtools index $id.sorted.redu.bam

####Expression level
gff=/store/jelly/yichun/ccin/fgra/genome/Fusarium_graminearum.genomic.gff
dir=DEG
cat samples_n_reads_decribed.txt | while read id1 id2 SRR;
do  
stringtie -e -B -p 38 -G $gff -A ${id2}_gene_count.xls -o ballgown/${id2}_ballgown/$id2.gtf $id2.sorted.bam
done
perl /home/yichun/RNAmodification/stringte_gene_count_matrix.pl *_gene_count.xls

prepDE.py -i ballgown/ -l 50

sed -i 's/,/\t/g;s/_ballgown//g' gene_count_matrix.csv transcript_count_matrix.csv

Rscript /home/yichun/RNAmodification/edgeR_TMM_norm.R
