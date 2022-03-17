#load package
library(edgeR)

####Read in TMM expression matrix
rawcount<-read.delim("gene_count_matrix.csv", header = TRUE)
row.names(rawcount)<-rawcount$gene_id
TMM<-rawcount[,2:ncol(rawcount)]

data=as.matrix(TMM)

####Read in sample group
trait<-read.delim("samples_n_reads_decribed.txt", header = FALSE)
group <- factor(trait[,2])
y<-DGEList(counts=data,group=group)
y <- calcNormFactors(y, method = "TMM")
TMM.norm<-cpm(y)

write.table(TMM.norm, file = "gene_count_matrix.TMM.xls", sep = "\t")

TMM.norm.log2<-log2(TMM.norm+1)
write.table(TMM.norm.log2, file = "gene_count_matrix.TMM_log2.xls", sep = "\t")

rawcount<-read.delim("transcript_count_matrix.csv", header = TRUE)
row.names(rawcount)<-rawcount$gene_id
TMM<-rawcount[,2:ncol(rawcount)]

data=as.matrix(TMM)

####Read in sample group
trait<-read.delim("samples_n_reads_decribed.txt", header = FALSE)
group <- factor(trait[,2])
y<-DGEList(counts=data,group=group)
y <- calcNormFactors(y, method = "TMM")
TMM.norm<-cpm(y)

write.table(TMM.norm, file = "transcript_count_matrix.TMM.xls", sep = "\t")

TMM.norm.log2<-log2(TMM.norm+1)
write.table(TMM.norm.log2, file = "transcript_count_matrix.TMM_log2.xls", sep = "\t")
