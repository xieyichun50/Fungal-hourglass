##usage
#Rscript PSfinal.R -i *.diamond.tab.sorted.PS.txt -o *.PSfinal.txt

#Output the following
##Species.taxidN.txt (unique performed to reduce record numbers)

#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))

##Create parameters
option_list <- list(
  make_option(c("-i","--input"), type="character", default=NULL,
              help="Species*.diamond.taxon.sorted.tab generate by diamond.sh' [default %default]",
              dest="input"),
  make_option(c("-o","--output"), type="character", default=NULL,
              help="output filename prefix [default %default]",
              dest="output_filename"),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print out all parameter settings [default]")
)

options(error=traceback)

parser <- OptionParser(usage = "%prog -i orthologs.txt -o out [options]",option_list=option_list)
opt = parse_args(parser)

##manual add files
{
  #opt$input<-"Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.proteins.fa.diamond.tab.sorted.PS.txt"
  #opt$output<-paste0(gsub(".diamond.tab.sorted.PS.txt","",opt$input),".PSfinal.txt")
}

cat(paste0("File will be written to ",opt$output))

PStable<-read.delim(opt$input, header = F)
names(PStable)<-c("qseqid","PS")
PStable<-PStable[PStable$qseqid != "qseqid",]
PStable$PS<-as.numeric(PStable$PS)
if (opt$input == "Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.proteins.fa.tab.sorted.PS.txt") {
  PStable<-separate(PStable, qseqid, c("Gene","Transcript"), sep = "-")
} else if (opt$input == "Fusarium_graminearum.protein.fa.tab.sorted.PS.txt"){
  IDmatch<-read.delim("Fusarium_graminearum.GenematchID", header = F)
  names(IDmatch)<-c("qseqid","Gene")
  PStable<-merge(PStable, IDmatch, by = "qseqid", all.x =T)
} else {
  PStable<-separate(PStable, qseqid, c("Gene","Transcript"), sep = "-")
}

PStable<-PStable %>% group_by(Gene) %>% summarise(PS=min(PS))
cat(nrow(PStable))

write.table(PStable, opt$output, sep = "\t", row.names = F, quote = F)
