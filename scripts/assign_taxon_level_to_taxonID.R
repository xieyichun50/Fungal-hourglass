#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(taxonomizr))
suppressPackageStartupMessages(library(optparse))

##Create parameters
option_list <- list(
  make_option(c("-i","--input"), type="character", default=NULL,
              help="Species*.taxonlist.txt generate by '01protein_ortholog_blast.sh' [default %default]",
              dest="taxonlist"),
  make_option(c("-s","--ref-species"), type="character", default=NULL,
              help="taxon id of your target species[default %default]",
              dest="refspecies"),
  make_option(c("-o","--output"), type="character", default="level",
              help="output filename suffix [default %default]",
              dest="output_filename"),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print out all parameter settings [default]")
)

options(error=traceback)

parser <- OptionParser(usage = "%prog -i orthologs.txt -o out [options]",option_list=option_list)
opt = parse_args(parser)

##For the first time, unmask
#prepareDatabase("nameNode.sqlite", getAccessions=FALSE)

##manual add files
{
  #opt$taxonlist<-"Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.proteins.fa.taxonlist.txt"
  #opt$refspecies<-"5346" ##Coprinopsis_cinerea
  #opt$refspecies<-"229533" ##Fusarium graminearum PH-1
}

taxonlist<-read.delim(paste0(opt$taxonlist), header = F)
names(taxonlist)[1]<-"txn"
ref.species<-as.data.frame(t(getTaxonomy(opt$refspecies,
                                         sqlFile = "nameNode.sqlite",
                                         desiredTaxa = c("no rank",
                                                         "superkingdom",
                                                         "kingdom",
                                                         "subkingdom",
                                                         "phylum", 
                                                         "subphylum", 
                                                         "class",
                                                         "subclass",
                                                         "order", 
                                                         "family", 
                                                         "genus",
                                                         "species"))))
print(paste0("Start mapping ", opt$taxonlist, ", target species ", opt$refspecies, "\n"))
print(ref.species)

ref.species<-ref.species[,1]

##add columns to taxonlist table
taxonlist$`no rank`<-NA
taxonlist$superkingdom<-NA
taxonlist$kingdom<-NA
taxonlist$subkingdom<-NA
taxonlist$phylum<-NA
taxonlist$subphylum<-NA
taxonlist$class<-NA
taxonlist$subclass<-NA
taxonlist$order<-NA
taxonlist$family<-NA
taxonlist$genus<-NA
taxonlist$species<-NA
taxonlist$PS<-NA

for (i in 1:nrow(taxonlist)) {
  if (i/ceiling(nrow(taxonlist)/100) %% 1 ==0) {
    print(paste0("Finished ", i/ceiling(nrow(taxonlist)/100),"%"))
  }
  test.species<-as.data.frame(t(getTaxonomy(taxonlist$txn[i],
                                            sqlFile = "/home/yichun/tools/nameNode.sqlite",
                                            desiredTaxa = c("no rank",
                                                            "superkingdom",
                                                            "kingdom",
                                                            "subkingdom",
                                                            "phylum", 
                                                            "subphylum", 
                                                            "class",
                                                            "subclass",
                                                            "order", 
                                                            "family", 
                                                            "genus",
                                                            "species"))))
  taxonlist[i,2:13]<-test.species[,1]
  taxonlist$PS[i]<-length(which(test.species[,1] %in% ref.species)==T)
}

taxonlist<-unique(taxonlist[is.na(taxonlist$`no rank`)==F,])
print(paste0(nrow(taxonlist)," has mapped taxonomy, start writing results to files: \n"))

print(paste0("Writing full mapping result to file '", opt$taxonlist,".full.",opt$output_filename,".txt'\n"))
write.table(taxonlist, 
            paste0(opt$taxonlist,".full.",opt$output_filename,".txt"),
            row.names = F, quote = F, sep = "\t")

print(paste0("Writing simplified PS mapping result to file '", opt$taxonlist,".simplePS.",opt$output_filename,".txt'\n"))
write.table(taxonlist[,c("txn","PS")], 
            paste0(opt$taxonlist,".simplePS.",opt$output_filename,".txt"),
            row.names = F, quote = F, sep = "\t")

cat(paste0("Finished ", opt$taxonlist, "!\n"))