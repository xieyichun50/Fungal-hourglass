##usage
#Rscript dNdS_calculation -q query.cds-transcripts.fa -r reference.cds-transcripts.fa -b /root/miniconda3/bin/blastp -m /root/miniconda3/bin/mafft

#Output the following
##Species.taxidN.txt (unique performed to reduce record numbers)

#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(orthologr))

##Create parameters
option_list <- list(
  make_option(c("-q","--query"), type="character", default=NULL,
              help="cds sequence of your targeted species' [default %default]",
              dest="qry"),
  make_option(c("-r","--reference"), type="character", default=NULL,
              help="cds sequence of your reference species' [default %default]",
              dest="ref"),
  make_option(c("-b","--blast"), type="character", default=NULL,
              help="which blastp' [default %default]",
              dest="blast"),
  make_option(c("-m","--mafft"), type="character", default=NULL,
              help="which mafft' [default %default]",
              dest="mafft"),
  make_option(c("-k","--kaks"), type="character", default=NULL,
              help="which KaKs_Calculator' [default %default]",
              dest="kaks"),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print out all parameter settings [default]")
)

options(error=traceback)

parser <- OptionParser(usage = "%prog -i orthologs.txt -o out [options]",option_list=option_list)
opt = parse_args(parser)

##manual add files
{
  #opt$qry<-"Coprinopsis_cinerea.longest-gene.cds.fa"
  #opt$ref<-"Fusarium_graminearum.longest-gene.cds.fa"
  #opt$blast<-"/root/miniconda3/bin/blastp"
  #opt$mafft<-"/root/miniconda3/bin/mafft"
  #opt$kaks<-"/root/miniconda3/bin/KaKs_Calculator"
}

opt$output<-paste0(gsub(".longest-gene.cds.fa","",opt$qry),
                   "_VS_",
                   gsub(".longest-gene.cds.fa","",opt$ref),
                   ".NSfinal.txt")
cat(paste0("File will be written to ",opt$output,"\n"))

# using the `aa_aln_path` or `blast_path` arguments

dnds.result<-dNdS(query_file = opt$qry,
                  subject_file = opt$ref,
                  ortho_detection = "RBH",
                  blast_path = opt$blast,
                  aa_aln_type = "multiple", 
                  aa_aln_tool = "mafft", 
                  aa_aln_path = opt$mafft,
                  codon_aln_tool = "pal2nal", 
                  dnds_est.method = "YN",
                  kaks_calc_path = opt$kaks,
                  comp_cores = 38, 
                  clean_folders = TRUE,
                  delete_corrupt_cds = TRUE)

if (opt$qry=="Coprinopsis_cinerea.longest-gene.cds.fa") {
  IDmatch<-read.delim("Coprinopsis_cinerea.GenematchID", header = F)
  names(IDmatch)<-c("Genes","query_id")
  dnds.result<-merge(dnds.result, IDmatch, by = "query_id", all.x =T)
} else {
  dnds.result$Genes<-dnds.result$query_id
}
write.table(dnds.result, opt$output, 
            row.names = F, sep = "\t", quote = F)

save.image(paste0(opt$output,".RData"))