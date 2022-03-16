##usage
#Rscript assign_taxon_level_to_proteins.R -i *.diamond.taxon.tab.sorted -p *.simplePS.level.txt

#Output the following
##*Species.PS.txt (each splitted file, unique performed to reduce record numbers)

#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))

##Create parameters
option_list <- list(
  make_option(c("-i","--input"), type="character", default=NULL,
              help="Species*.diamond.taxon.sorted.tab generate by diamond.sh' [default %default]",
              dest="input"),
  make_option(c("-p","--ps"), type="character", default=NULL,
              help="*.simplePS.level.txt generate by 02mapping.sh' [default %default]",
              dest="psfile"),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print out all parameter settings [default]")
)

options(error=traceback)

parser <- OptionParser(usage = "%prog -i orthologs.txt -o out [options]",option_list=option_list)
opt = parse_args(parser)

##manual add files
{
  #opt$input<-"x9183Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.proteins.fa.diamond.taxon.tab.sorted"
  #opt$psfile<-"Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.proteins.fa.simplePS.level.txt"
}

cat(paste0("File will be written to ",opt$input,".PS.txt\n"))

protout<-read.delim(opt$input, header = F)
names(protout)<-c("qseqid", "staxids")
pathways<-protout[grep(";",protout$staxids),]
protout<-protout[-grep(";",protout$staxids),]
qsetaxon.1v1<-matrix(NA, nrow = 1, ncol = 2)
qsetaxon.1v1<-as.data.frame(qsetaxon.1v1)
names(qsetaxon.1v1)<-c("qseqid", "staxids")

for (i in 1:nrow(pathways)) {
  subtable<-pathways[i,]
  rcnames<-list(c(strsplit(subtable$qseqid[1], ',')[[1]]),c(strsplit(subtable$staxids[1], ';')[[1]]))
  pairtable<-matrix(data = NA, nrow = length(rcnames[[1]]), ncol = length(rcnames[[2]]), dimnames = rcnames)
  pairtable<-as.data.frame(pairtable)
  pairtable$qseqid<-rownames(pairtable)
  rownames(pairtable)<-1:nrow(pairtable)
  pairtable<-as.data.frame(pairtable)
  pairtable.new<-pairtable %>% gather(staxids, pair, c(1:ncol(pairtable)-1))
  pairtable.new<-pairtable.new[,c(1:2)]
  qsetaxon.1v1<-rbind(qsetaxon.1v1, pairtable.new)
}

protout<-rbind(protout,qsetaxon.1v1)

psfile<-read.delim(opt$psfile, header = T)
names(psfile)<-c("staxids","PS")
psfile<-psfile[psfile$PS != "0" & psfile$PS != "PS",]
psfile$PS<-as.numeric(psfile$PS)
protmatch<-merge(protout, psfile, by = "staxids", all.x = T)
protmatch<-protmatch[is.na(protmatch$PS)==F,c("qseqid","PS")]
protmatch<-unique(protmatch)
PS<-protmatch %>% group_by(qseqid) %>% summarise(PS=min(PS))
write.table(PS, 
            paste0(opt$input,".PS.txt"), 
            row.names = F, quote = F, sep = "\t")