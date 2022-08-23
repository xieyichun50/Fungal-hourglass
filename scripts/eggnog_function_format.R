library(dplyr)
library(tidyr)
library(stringr)

#Whole genome Annotation file

eggnog<-read.delim(paste0("eggnog/eggnog.emapper.annotations"),
                   header = TRUE, skip = 3)
eggnog<-separate(eggnog, X.query_name, c("Genes", "TS"), sep = "-", remove = FALSE)

##KEGG from eggnog
{
#read in kegg2name
kegg2name <- read.delim("D:/3enrichment/kegg2name.txt",
                        sep = "\t", colClasses = "character")

#eggnog<-separate(eggnog, V1, c("Genes", "Transcript"), sep = "-", remove = FALSE)
pathways<-eggnog[,c("Genes","KEGG_Pathway")]
pathways<-separate(pathways, KEGG_Pathway, c("ko","map"), sep = ",map", remove = FALSE)
pathways<-pathways[,c(1,3)]
names(pathways)[1]="Genes"
names(pathways)[2]="KEGG"
pathways<-subset(pathways, KEGG != "-" & KEGG != "" & is.na(KEGG)==FALSE,
                 select = c("Genes", "KEGG"))
pathways<-unique(pathways)
#Format GenesKEGGpair
GenesKEGGpair.1v1<-matrix(NA, nrow = 1, ncol = 2)
GenesKEGGpair.1v1<-as.data.frame(GenesKEGGpair.1v1)
names(GenesKEGGpair.1v1)[1]="Genes"
names(GenesKEGGpair.1v1)[2]="KEGG"

for (i in 1:nrow(pathways)) {
  subtable<-pathways[i,]
  rcnames<-list(c(strsplit(subtable$Genes[1], ',')[[1]]),c(strsplit(subtable$KEGG[1], ',')[[1]]))
  pairtable<-matrix(data = NA, nrow = length(rcnames[[1]]), ncol = length(rcnames[[2]]), dimnames = rcnames)
  pairtable<-as.data.frame(pairtable)
  pairtable$Genes<-rownames(pairtable)
  rownames(pairtable)<-1:nrow(pairtable)
  pairtable<-as.data.frame(pairtable)
  pairtable.new<-pairtable %>% gather(KEGG, pair, c(1:ncol(pairtable)-1))
  pairtable.new<-pairtable.new[,c(1:2)]
  GenesKEGGpair.1v1<-rbind(GenesKEGGpair.1v1, pairtable.new)
}
GenesKEGGpair.1v1<-subset(GenesKEGGpair.1v1,
                          is.na(GenesKEGGpair.1v1$Genes)==FALSE,
                          select = c("KEGG", "Genes"))
GenesKEGGpair.1v1<-unique(GenesKEGGpair.1v1)

write.table(GenesKEGGpair.1v1,
            file = "eggnog/KEGG.1v1.txt",
            sep = '\t', row.names = FALSE,
            quote = FALSE)
rm(pairtable, pairtable.new, rcnames, subtable, pathways)
}
KEGGfreq<-as.data.frame(xtabs(~Genes, GenesKEGGpair.1v1))
##KO from eggnog
{
  #read in kegg2name
  ko2name <- read.delim("D:/3enrichment/ko2name.txt",
                          sep = "\t", colClasses = "character")

  #eggnog<-separate(eggnog, V1, c("Genes", "Transcript"), sep = "-", remove = FALSE)
  pathways<-eggnog[,c("Genes","KEGG_ko")]
  names(pathways)[1]="Genes"
  names(pathways)[2]="ko"
  pathways<-subset(pathways, ko != "-" & ko != "" & is.na(ko)==FALSE,
                   select = c("Genes", "ko"))
  pathways<-unique(pathways)
  #Format GenesKEGGpair
  Geneskopair.1v1<-matrix(NA, nrow = 1, ncol = 2)
  Geneskopair.1v1<-as.data.frame(Geneskopair.1v1)
  names(Geneskopair.1v1)[1]="Genes"
  names(Geneskopair.1v1)[2]="ko"

  for (i in 1:nrow(pathways)) {
    subtable<-pathways[i,]
    rcnames<-list(c(strsplit(subtable$Genes[1], ',')[[1]]),c(strsplit(subtable$ko[1], ',')[[1]]))
    pairtable<-matrix(data = NA, nrow = length(rcnames[[1]]), ncol = length(rcnames[[2]]), dimnames = rcnames)
    pairtable<-as.data.frame(pairtable)
    pairtable$Genes<-rownames(pairtable)
    rownames(pairtable)<-1:nrow(pairtable)
    pairtable<-as.data.frame(pairtable)
    pairtable.new<-pairtable %>% gather(ko, pair, c(1:ncol(pairtable)-1))
    pairtable.new<-pairtable.new[,c(1:2)]
    Geneskopair.1v1<-rbind(Geneskopair.1v1, pairtable.new)
  }
  Geneskopair.1v1<-subset(Geneskopair.1v1,
                            is.na(Geneskopair.1v1$Genes)==FALSE,
                            select = c("ko", "Genes"))
  Geneskopair.1v1<-unique(Geneskopair.1v1)
  Geneskopair.1v1$ko<-gsub("ko:","", Geneskopair.1v1$ko)

  write.table(Geneskopair.1v1,
              file = "eggnog/ko.1v1.txt",
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  rm(pairtable, pairtable.new, rcnames, subtable, pathways)
}
KOfreq<-as.data.frame(xtabs(~Genes, Geneskopair.1v1))
##KOG from eggnog
{
#read in kog2name
kog2name<-read.delim("D:/3enrichment/kog2name.txt",
                     sep = "\t", colClasses = "character")

pathways<-eggnog[,c("Genes","COG.Functional.cat.")]
names(pathways)[1]="Genes"
names(pathways)[2]="KOG"
pathways<-subset(pathways, KOG != "" & KOG != "-" & is.na(KOG)==FALSE)
pathways<-unique(pathways)

#Format GenesKOGpair
GenesKOGpair.1v1<-matrix(NA, nrow = 1, ncol = 2)
GenesKOGpair.1v1<-as.data.frame(GenesKOGpair.1v1)
names(GenesKOGpair.1v1)[1]="Genes"
names(GenesKOGpair.1v1)[2]="KOG"

for (i in 1:nrow(pathways)) {
  subtable<-pathways[i,]
  rcnames<-list(c(strsplit(subtable$Genes[1], ',')[[1]]),c(strsplit(subtable$KOG[1], '')[[1]]))
  pairtable<-matrix(data = NA, nrow = length(rcnames[[1]]), ncol = length(rcnames[[2]]), dimnames = rcnames)
  pairtable<-as.data.frame(pairtable)
  pairtable$Genes<-rownames(pairtable)
  rownames(pairtable)<-1:nrow(pairtable)
  pairtable<-as.data.frame(pairtable)
  pairtable.new<-pairtable %>% gather(KOG, pair, c(1:ncol(pairtable)-1))
  pairtable.new<-pairtable.new[,c(1:2)]
  GenesKOGpair.1v1<-rbind(GenesKOGpair.1v1, pairtable.new)
}
GenesKOGpair.1v1<-subset(GenesKOGpair.1v1,
                         is.na(GenesKOGpair.1v1$Genes)==FALSE,
                         select = c("KOG", "Genes"))
GenesKOGpair.1v1<-unique(GenesKOGpair.1v1)

write.table(GenesKOGpair.1v1,
            file = "eggnog/KOG.1v1.txt",
            sep = '\t', row.names = FALSE,
            quote = FALSE)
rm(pairtable, pairtable.new, rcnames, subtable, pathways)

}

KOGfreq<-as.data.frame(xtabs(~Genes, GenesKOGpair.1v1))
##GOfreq
GenesGOpair.1v1<-read.delim("C:/coprinopsis/genome/annotate_misc/annotations.GO.txt",header = FALSE)
GenesGOpair.1v1<-GenesGOpair.1v1[grep("GO",GenesGOpair.1v1$V3),c(3,1)]
names(GenesGOpair.1v1)<-c("GO", "Genes")
GenesGOpair.1v1<-separate(GenesGOpair.1v1, Genes, c("Genes"), sep = "-", remove = TRUE)
GenesGOpair.1v1<-unique(GenesGOpair.1v1)
write.table(GenesGOpair.1v1,
            file = "eggnog/GO.1v1.txt",
            sep = '\t', row.names = FALSE,
            quote = FALSE)
GOfreq<-as.data.frame(xtabs(~Genes, GenesGOpair.1v1))

