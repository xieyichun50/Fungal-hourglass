setwd("hourglass/Neurospora_myTAI/")
library(myTAI)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotrix)
library(agricolae)
library(eoffice)
library(ape)
library(tidytree)
library(ggtree)
library(flextable)
library(stringr)
library(ggplotify)
library(pheatmap)
library(sva)
library(edgeR)
library(BatchQC)

##Read in TAI features from Neurospora crassa
species="Neurospora_crassa"
genes.PS<-read.delim(paste0("/home/yichun/RNAmodification/hourglass/blast/",
                            species,"/e10-5aa30/", 
                            species,".proteins.fa.PSfinal.txt"), header = T)
names(genes.PS)[names(genes.PS)=="Gene"]<-"Genes"
genes.PS.stat<-as.data.frame(xtabs(~PS, genes.PS))

##PS distribution
basetree<-read.tree(text = '(((((((((((((Species,Spc)PS12,Genus)PS11,Family)PS10,Order)PS9,Subclass)PS8,Class)PS7,Subphylum)PS6,Phylum)PS5,Subkingdom)PS4,Kingdom)PS3,Superkingdom)PS2,Org.)PS1,Life)PS0;')
tip.order<-data.frame(node=1:Ntip(basetree), basetree.label = basetree$tip.label)
tip.order$Name<-NA
tip.order$Taxon<-NA
tip.order$Spc<-NA
node.order<-data.frame(node=1:Nnode(basetree) + Ntip(basetree), basetree.label = basetree$node.label)
genes.PS.stat$PS<-paste0("PS",genes.PS.stat$PS)
names(genes.PS.stat)<-c("basetree.label","Name")
node.order<-merge(node.order,genes.PS.stat, by = "basetree.label", all.x = T)

taxoname<-read.delim(paste0(species,".taxon.level.txt"), header = T)
names(taxoname)<-c("basetree.label", "Taxon", "Spc")
node.order<-merge(node.order,taxoname, by = "basetree.label", all.x = T)

nt.order<-rbind(tip.order, node.order)

PS.tree<-full_join(basetree, nt.order, by = "node")

p<-ggtree(PS.tree, size = 1)+
  geom_nodepoint()+
  geom_tiplab(aes(label=label),
              size = 3)+
  geom_text2(aes(subset=!isTip & is.na(Spc)==FALSE, label=Spc), 
             nudge_x = -0.1, nudge_y = 1.5,  
             size=3, color = "black")+
  geom_text2(aes(subset=!isTip, 
                 label=basetree.label), 
             nudge_x = 0.5, nudge_y = 0.3,
             size=3, color = "black")+
  geom_text2(aes(subset=!isTip & is.na(Name)==FALSE, 
                 label=paste0("(",Name,")")),
             nudge_x = 0.5, nudge_y = -0.3,
             size=3, color = "black")
p

f=paste0(species,".PS.tree.pptx")
topptx(p, f, width = 5.8, height = 5, units = "in")
rm(basetree,node.order,nt.order,p,PS.tree, tip.order,taxoname)

#### Expression levels
rawcount1<-read.csv("expression/ballgown/germ_MSM/gene_count_matrix.csv", header = T)
rawcount2<-read.csv("expression/ballgown/sexW/gene_count_matrix.csv", header = T)
rawcount<-merge(rawcount1,rawcount2, by = "gene_id", all = T)
rm(rawcount1,rawcount2)

row.names(rawcount)<-gsub("gene-","", rawcount$gene_id)
rawcount<-rawcount[,2:ncol(rawcount)]
names(rawcount)<-gsub("_ballgown","",names(rawcount))

samplelist<-read.delim("expression/samples_n_reads_decribed_all.txt", header = F)
names(samplelist)<-c("group","sample","SRR","batch")

## batch correction
rawcount.sub<-rawcount[, names(rawcount) %in% samplelist$sample[samplelist$batch %in% c("II","III")]]
rawcount.sub<-as.matrix(rawcount.sub)
samplelist$group.new<-samplelist$group

normcount<-ComBat_seq(rawcount.sub, batch = samplelist$batch[samplelist$batch %in% c("II","III")], 
                      group = samplelist$group.new[samplelist$batch %in% c("II","III")], full_mod = TRUE)
normcount<-as.data.frame(normcount)
normcount.bk<-normcount
id1<-c("T15.rep1","T15.rep2","T15.rep3",
       "T120.rep1","T120.rep2","T120.rep3",
       "T240.rep1","T240.rep2","T240.rep3",
       "T360.rep1","T360.rep2","T360.rep3")
id2<-c("T0h.rep1","T0h.rep2",
       "T2h.rep1","T2h.rep2", 
       "T24h.rep1","T24h.rep2",
       "T48h.rep1","T48h.rep2",
       "T72h.rep1","T72h.rep2",
       "T96h.rep1","T96h.rep2",
       "T120h.rep1","T120h.rep2",
       "T144h.rep1","T144h.rep2")

newcount<-merge(rawcount[,id1], normcount[,id2], by = 0, all = T)
row.names(newcount)<-newcount$Row.names
newcount<-newcount[,2:ncol(newcount)]

newcount<-rawcount
group <- factor(names(newcount))
y<-DGEList(counts=newcount,group=group)
y <- calcNormFactors(y, method = "TMM")
TMM.norm<-cpm(y)

write.table(TMM.norm, file = "gene_count_matrix.TMM.xls", sep = "\t")

TMM.norm.log2<-log2(TMM.norm+1)
write.table(TMM.norm.log2, file = "gene_count_matrix.TMM_log2.xls", sep = "\t")

##select stages
treat.order<-c("T15","T120","T240","T360",
               "T0h","T2h","T24h",
               "T48h","T72h","T96h","T120h","T144h")

#TMM mean
TMM<-as.data.frame(TMM.norm.log2)

TMM.mean<-as.data.frame(matrix(NA, nrow = nrow(TMM), ncol = length(treat.order)))
row.names(TMM.mean)<-row.names(TMM)
names(TMM.mean)<-treat.order

for (i in 1:4) {
  TMM.mean[,i]<-(TMM[,names(TMM)==paste0(names(TMM.mean)[i],".rep1")]+TMM[,names(TMM)==paste0(names(TMM.mean)[i],".rep2")]+TMM[,names(TMM)==paste0(names(TMM.mean)[i],".rep3")])/3
}
for (i in 5:ncol(TMM.mean)) {
  TMM.mean[,i]<-(TMM[,names(TMM)==paste0(names(TMM.mean)[i],".rep1")]+TMM[,names(TMM)==paste0(names(TMM.mean)[i],".rep2")])/2
}
TMM.mean[TMM.mean<1]<-0

TMM.mean$Genes<-row.names(TMM.mean)
row.names(TMM.mean)<-1:nrow(TMM.mean)
TMM.matrix<-TMM.mean[,c("Genes",treat.order)]
TMM.matrix$Genes<-gsub("gene-","",TMM.matrix$Genes)
rm(i, TMM.mean, TMM)

write.table(TMM.matrix,
            "/home/yichun/RNAmodification/hourglass/Neurospora_myTAI/expression/gene_count_log2TMMmatrixmean.txt",
            sep = "\t", quote = F, row.names = F)
write.table(TMM.matrix,
            "/home/yichun/RNAmodification/hourglass/Neurospora_myTAI/Neurospora_crassa.log2TMMmatrixmean.txt",
            sep = "\t", quote = F, row.names = F)
##########################################################
####TAI feature
##########################################################
treat.order<-c("T15","T120","T240","T360",
               "T0h","T2h","T24h",
               "T48h","T72h","T96h","T120h","T144h")
TMM.matrix<-read.delim("Neurospora_crassa.log2TMMmatrixmean.txt", header = TRUE)
TMM.matrix<-TMM.matrix[,c("Genes", treat.order)]

##Read in TAI features from Neurospora crassa
species="Neurospora_crassa"
genes.PS<-read.delim(paste0("/home/yichun/RNAmodification/hourglass/blast/",
                            species,"/e10-5aa30/", 
                            species,".proteins.fa.PSfinal.txt"), header = T)
names(genes.PS)[names(genes.PS)=="Gene"]<-"Genes"

##TAI analysis
####All
ref="All"
treat.order<-c("Spore","Polar","Elong","Brch",
               "Proto","Mate","PIni","PWall",
               "Para","Asci","YP","MP")
names(TMM.matrix)<-c("Genes", treat.order)
plot.order<-c("MP","YP","Asci","Para",
              "PWall","PIni","Mate","Proto",
              "Brch","Elong","Polar","Spore")
genes.hour<-genes.PS
module.list.germ<-list(early = 1, mid = 2:3, late = 4)
module.list.sex<-list(early = 1, mid = 2:3, late = 4:8)

##Merge PS to each table
TMM.matrix.new<-merge(TMM.matrix, genes.hour, by = "Genes", all = F)
TMM.matrix.new<-TMM.matrix.new[,c("PS","Genes",treat.order)]
TMM.matrix.new[is.na(TMM.matrix.new)==T]<-0
##Calculate TAI
TXI.result<-TAI(TMM.matrix.new)
TXI.result
TXI.result<-as.data.frame(TXI.result)
TXI.result$group<-NA
TXI.result$p.value<-NA
TXI.result$std.dev<-NA
TXI.result$statistics<-NA
TXI.result.FlatLineTest<-TXI.result
TXI.result.ReductiveHourglassTest<-TXI.result
TXI.result.list<-as.data.frame(matrix(NA, ncol = 7, nrow = 0))
names(TXI.result.list)<-c("Stage","TXI.result","group","p.value","std.dev","statistics","Process")

##############################################################################################
##Pvalue/contribution/ratio for each combination
##############################################################################################

FlatLineTest.result<-FlatLineTest(TMM.matrix.new[1:6],
                                  permutations = 1000,
                                  plotHistogram = FALSE,
                                  runs = 100,
                                  parallel = 14,
                                  custom.perm.matrix = NULL)
TXI.result.germ<-TXI.result[row.names(TXI.result) %in% c("Spore","Polar","Elong","Brch"),]
TXI.result.germ$Stage<-row.names(TXI.result.germ)
TXI.result.germ$group<-ref
TXI.result.germ$statistics<-"FlatLineTest"
TXI.result.germ$Process<-"Germ"
TXI.result.germ$p.value<-FlatLineTest.result$p.value
TXI.result.germ$std.dev<-FlatLineTest.result$std.dev
TXI.result.list<-rbind(TXI.result.list, TXI.result.germ)

FlatLineTest.result<-FlatLineTest(TMM.matrix.new[c(1:2,7:14)],
                                  permutations = 1000,
                                  plotHistogram = FALSE,
                                  runs = 100,
                                  parallel = 14,
                                  custom.perm.matrix = NULL)
TXI.result.sex<-TXI.result[row.names(TXI.result) %in% c("Proto","Mate","PIni","PWall",
                                                        "Para","Asci","YP","MP"),]
TXI.result.sex$Stage<-row.names(TXI.result.sex)
TXI.result.sex$group<-ref
TXI.result.sex$statistics<-"FlatLineTest"
TXI.result.sex$Process<-"Sex"
TXI.result.sex$p.value<-FlatLineTest.result$p.value
TXI.result.sex$std.dev<-FlatLineTest.result$std.dev
TXI.result.list<-rbind(TXI.result.list, TXI.result.sex)

TXI.result.sub<-rbind(TXI.result.germ,TXI.result.sex)

TAI<-TXI.result.sub %>%
  mutate(Stage = factor(Stage, levels = plot.order)) %>%
  ggplot(aes(x = Stage, y = TXI.result, group = Process))+
  geom_rect(aes(xmin=nrow(TXI.result.sub)-min(module.list.germ[[2]]-0.5)+1, 
                xmax=nrow(TXI.result.sub)-max(module.list.germ[[2]]-0.5), 
                ymin=-Inf, ymax=Inf),
            fill='dodgerblue',alpha = 0.05)+
  geom_rect(aes(xmin=nrow(TXI.result.sub)-min(module.list.sex[[2]]-0.5+4)+1, 
                xmax=nrow(TXI.result.sub)-max(module.list.sex[[2]]-0.5+4), 
                ymin=-Inf, ymax=Inf),
            fill='dodgerblue',alpha = 0.05)+
  geom_ribbon(aes(ymin=TXI.result-std.dev, ymax=TXI.result+std.dev), fill = "gray50", alpha = 0.5)+
  geom_line(aes(linetype="solid"))+
  scale_y_continuous(limits = c(2.30,1.80),
                     breaks = seq(1.80,2.30,0.1),
                     trans = "reverse", position = "right")+
  scale_x_discrete(position = "top")+
  annotate(geom = "text", 
           label=paste0("Flat Line Test\nP.value=", 
                        signif(unique(TXI.result.sub$p.value[TXI.result.sub$Process=="Germ"]),2)),
           x=nrow(TXI.result.sub)-mean(module.list.germ[[2]])+1, y=2.2, 
           colour="black", size=2)+
  annotate(geom = "text", 
           label=paste0("Flat Line Test\nP.value=", 
                        signif(unique(TXI.result.sub$p.value[TXI.result.sub$Process=="Sex"]),2)),
           x=nrow(TXI.result.sub)-mean(module.list.germ[[2]]+4)+1, y=2.2, 
           colour="black", size=2)+
  labs(title = "", x = "", y = "Transcriptome Age Index", colour = NULL)+
  coord_flip()+
  theme(axis.line = element_line(linetype = "solid"),
        axis.ticks.y = element_line(colour = "black", size = 0.5),
        axis.ticks.x = element_line(colour = "black", size = 0.5),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        plot.title = element_text(size = 8, hjust = 0.5, face = "plain"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 0),
        panel.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        legend.position = "none",
        strip.text = element_text(size = 9))
TAI
ggsave(paste0(species,".TAI.Flatline.",ref,".png"), width = 4, height = 3, units = "in", dpi = 300)
f=paste0(species,".TAI.Flatline.",ref,".pptx")
topptx(TAI, f, width = 4, height = 3, units = "in")
rm(TXI.result.sub)

##ReductiveHourglassTest (High-Low-High)
ReductiveHourglassTest.result<-ReductiveHourglassTest(TMM.matrix.new[1:6],
                                                      modules = module.list.germ,
                                                      permutations = 1000,
                                                      lillie.test = FALSE,
                                                      plotHistogram = FALSE,
                                                      runs = 100,
                                                      parallel = 14,
                                                      gof.warning = FALSE,
                                                      custom.perm.matrix = NULL)
TXI.result.germ<-TXI.result[row.names(TXI.result) %in% c("Spore","Polar","Elong","Brch"),]
TXI.result.germ$Stage<-row.names(TXI.result.germ)
TXI.result.germ$group<-ref
TXI.result.germ$statistics<-"ReductiveHourglassTest"
TXI.result.germ$Process<-"Germ"
TXI.result.germ$p.value<-ReductiveHourglassTest.result$p.value
TXI.result.germ$std.dev<-ReductiveHourglassTest.result$std.dev
TXI.result.list<-rbind(TXI.result.list, TXI.result.germ)

ReductiveHourglassTest.result<-ReductiveHourglassTest(TMM.matrix.new[c(1:2,7:14)],
                                                      modules = module.list.sex,
                                                      permutations = 1000,
                                                      lillie.test = FALSE,
                                                      plotHistogram = FALSE,
                                                      runs = 100,
                                                      parallel = 14,
                                                      gof.warning = FALSE,
                                                      custom.perm.matrix = NULL)

TXI.result.sex<-TXI.result[row.names(TXI.result) %in% c("Proto","Mate","PIni","PWall",
                                                        "Para","Asci","YP","MP"),]
TXI.result.sex$Stage<-row.names(TXI.result.sex)
TXI.result.sex$group<-ref
TXI.result.sex$statistics<-"ReductiveHourglassTest"
TXI.result.sex$Process<-"Sex"
TXI.result.sex$p.value<-ReductiveHourglassTest.result$p.value
TXI.result.sex$std.dev<-ReductiveHourglassTest.result$std.dev
TXI.result.list<-rbind(TXI.result.list, TXI.result.sex)

TXI.result.sub<-rbind(TXI.result.germ, TXI.result.sex)

TAI<-TXI.result.sub %>%
  mutate(Stage = factor(Stage, levels = plot.order)) %>%
  ggplot(aes(x = Stage, y = TXI.result, group = Process))+
  geom_rect(aes(xmin=nrow(TXI.result.sub)-min(module.list.germ[[2]]-0.5)+1, 
                xmax=nrow(TXI.result.sub)-max(module.list.germ[[2]]-0.5), 
                ymin=-Inf, ymax=Inf),
            fill='dodgerblue',alpha = 0.05)+
  geom_rect(aes(xmin=nrow(TXI.result.sub)-min(module.list.sex[[2]]-0.5+4)+1, 
                xmax=nrow(TXI.result.sub)-max(module.list.sex[[2]]-0.5+4), 
                ymin=-Inf, ymax=Inf),
            fill='dodgerblue',alpha = 0.05)+
  geom_ribbon(aes(ymin=TXI.result-std.dev, ymax=TXI.result+std.dev), fill = "gray50", alpha = 0.5)+
  geom_line(aes(linetype="solid"))+
  scale_y_continuous(limits = c(2.30,1.80),
                     breaks = seq(1.80,2.30,0.1),
                     trans = "reverse", position = "right")+
  scale_x_discrete(position = "top")+
  annotate(geom = "text", 
           label=paste0("Reductive Hourglass Test\n(High-Low-High)\nP.value=", 
                        signif(unique(TXI.result.sub$p.value[TXI.result.sub$Process=="Germ"]),2)),
           x=nrow(TXI.result.sub)-mean(module.list.germ[[2]])+1, y=2.2, 
           colour="black", size=2)+
  annotate(geom = "text", 
           label=paste0("Reductive Hourglass Test\n(High-Low-High)\nP.value=", 
                        signif(unique(TXI.result.sub$p.value[TXI.result.sub$Process=="Sex"]),2)),
           x=nrow(TXI.result.sub)-mean(module.list.germ[[2]]+4)+1, y=2.2, 
           colour="black", size=2)+
  labs(title = "", x = "", y = "Transcriptome Age Index", colour = NULL)+
  coord_flip()+
  theme(axis.line = element_line(linetype = "solid"),
        axis.ticks.y = element_line(colour = "black", size = 0.5),
        axis.ticks.x = element_line(colour = "black", size = 0.5),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        plot.title = element_text(size = 8, hjust = 0.5, face = "plain"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 0),
        panel.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        legend.position = "none",
        strip.text = element_text(size = 9))
TAI
ggsave(paste0(species,".TAI.reductive.",ref,".png"), width = 4, height = 3, units = "in", dpi = 300)
f=paste0(species,".TAI.reductive.",ref,".pptx")
topptx(TAI, f, width = 4, height = 3, units = "in")

##Relative expression of all stages by PS
Relative.exp<-as.data.frame(REMatrix(TMM.matrix.new))

##Group ratio
REmatrix <- age.apply(ExpressionSet = TMM.matrix.new, RE)
Groups = list(group_1 = 1:2, group_2 = 3:12)
MeanREClassValues <- matrix(NA_real_,length(Groups),ncol(TMM.matrix.new)-2)
StdErr.RE.ClassValues <- matrix(NA_real_,length(Groups),ncol(TMM.matrix.new)-2)
for(i in 1:length(Groups)){
  MeanREClassValues[i , ] <- colMeans(REmatrix[match(as.character(Groups[[i]]),
                                                     rownames(REmatrix)) , ])
  
  StdErr.RE.ClassValues[i , ] <- apply(REmatrix[match(as.character(Groups[[i]]),
                                                      rownames(REmatrix)) , ],2,std.error)
}   

FoldChangeOfMeanREValues<-apply(MeanREClassValues,2,function(x){return((x[1])/(x[2]))})
FoldChangeOfMeanREValues[FoldChangeOfMeanREValues == "NaN"]<-0
REFoldChangeOfMeanREValues<-(FoldChangeOfMeanREValues-
                               min(FoldChangeOfMeanREValues))/(max(FoldChangeOfMeanREValues)-
                                                                 min(FoldChangeOfMeanREValues))

MeanREClassValues<-rbind(REmatrix,MeanREClassValues,REFoldChangeOfMeanREValues,StdErr.RE.ClassValues)
MeanREClassValues<-as.data.frame(MeanREClassValues)
names(MeanREClassValues)<-names(TMM.matrix.new[3:ncol(TMM.matrix.new)])
row.names(MeanREClassValues)<-c(paste0("PS",1:12),"PS1-2","PS3-12","Ratio","stderr.PS1-2","stderr.PS3-12")
write.table(MeanREClassValues, paste0(species,".PS.RE.",ref,".txt"),
            row.names = T, quote = F, sep ="\t")

##contribution
percentTAI<-as.data.frame(pStrata(TMM.matrix.new))
row.names(percentTAI)<-paste0("PS", 1:12)
write.table(percentTAI, paste0(species,".PS.contribution.",ref,".txt"),
            row.names = T, quote = F, sep ="\t")

PlotContribution( ExpressionSet = TMM.matrix.new,
                  legendName = "PS",
                  xlab = "Ontogeny",
                  ylab = "Transcriptome Age Index",
                  y.ticks = 10)
RE.bar<-PlotBarRE(TMM.matrix.new,
                  Groups = list(group_1 = 1:2, group_2 = 3:12),
                  ratio = T,
                  p.adjust.method = "BH")

##########################################################
####TDI feature
##########################################################
pairspecies<-c("Aspergillus_oryzae",
               "Botryosphaeria_dothidea",
               "Cercophora_caudata",
               "Chaetomium_globosum",
               "Colletotrichum_higginsianum",
               "Daldinia_childiae",
               "Fusarium_graminearum",
               "Lasiosphaeria_miniovina",
               "Morchella_importuna",
               "Neurospora_hispaniola",
               "Neurospora_tetrasperma",
               "Saccharomyces_cerevisiae",
               "Sordaria_macrospora",
               "Sporotrichum_thermophile",
               "Trichoderma_reesei",
               "Tuber_melanosporum")
species="Neurospora_crassa"
treat.order<-c("T15","T120","T240","T360",
               "T0h","T2h","T24h",
               "T48h","T72h","T96h","T120h","T144h")
TMM.matrix<-read.delim("Neurospora_crassa.log2TMMmatrixmean.txt", header = TRUE)
TMM.matrix<-TMM.matrix[,c("Genes", treat.order)]
ref="All"
treat.order<-c("Spore","Polar","Elong","Brch",
               "Proto","Mate","PIni","PWall",
               "Para","Asci","YP","MP")
names(TMM.matrix)<-c("Genes", treat.order)
plot.order<-c("MP","YP","Asci","Para",
              "PWall","PIni","Mate","Proto",
              "Brch","Elong","Polar","Spore")
genes.hour<-genes.PS
module.list.germ<-list(early = 1, mid = 2:3, late = 4)
module.list.sex<-list(early = 1, mid = 2:3, late = 4:8)

for (k in 1:length(pairspecies)) {
  ##Read in TDI features from Neurospora
  refspecies=pairspecies[k]
  genes.NS<-read.delim(paste0("/home/yichun/RNAmodification/hourglass/blast/",
                              species,
                              "/dnds/",species,"_VS_",
                              refspecies,".NSfinal.txt"))
  names(genes.NS)[names(genes.NS)=="dNdS"]<-"NS"
  genes.NS<-separate(genes.NS, "Genes", c("Genes"), sep = "T")
  genes.hour<-genes.NS[is.na(genes.NS$NS)==F & genes.NS$NS<=2,c("Genes","NS")]
  
  p<-ggplot(data = genes.hour, aes(x = NS, y = stat(density)/sum(stat(density))))+
    geom_histogram(color = "black", fill = "black", binwidth = 0.01)+
    #scale_x_continuous(limits = c(0,2), breaks = seq(0,2,0.5))+
    scale_y_continuous(limits = c(0,0.13),
                       breaks = seq(0,0.13, 0.02))+
    annotate(geom = "text", 
             label=paste0(gsub("_"," ",refspecies)),
             x = max(genes.hour$NS)/2, y=0.11, 
             colour="black", size=2, fontface="italic")+
    labs(x = "dN/dS", y = "Ratio")+
    theme(plot.subtitle = element_text(vjust = 1),
          plot.caption = element_text(vjust = 1),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 8, colour = "black"),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          plot.title = element_text(size = 0, hjust = 0.5),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          axis.line = element_line(size = 0.5),
          panel.background = element_rect(fill = NA),
          legend.position = "none")
  
  p
  ggsave(paste0(species,".dNdS.",refspecies,".png"), width = 3, height = 3, units = "in", dpi = 300)
  
  ####All
  ##Merge dNdS to each table
  TMM.matrix.new<-merge(TMM.matrix, genes.hour, by = "Genes", all = F)
  TMM.matrix.new<-TMM.matrix.new[,c("NS","Genes",treat.order)]
  TMM.matrix.new[is.na(TMM.matrix.new)==T]<-0
  
  ##Calculate TDI
  TXI.result<-TDI(TMM.matrix.new)
  TXI.result
  TXI.result<-as.data.frame(TXI.result)
  TXI.result$group<-NA
  TXI.result$p.value<-NA
  TXI.result$std.dev<-NA
  TXI.result$statistics<-NA
  TXI.result.FlatLineTest<-TXI.result
  TXI.result.ReductiveHourglassTest<-TXI.result
  TXI.result.list<-as.data.frame(matrix(NA, ncol = 7, nrow = 0))
  names(TXI.result.list)<-c("Stage","TXI.result","group","p.value","std.dev","statistics","Process")
  
  ##Pvalue for each combination
  FlatLineTest.result<-FlatLineTest(TMM.matrix.new[1:6],
                                    permutations = 1000,
                                    plotHistogram = FALSE,
                                    runs = 100,
                                    parallel = 14,
                                    custom.perm.matrix = NULL)
  TXI.result.germ<-TXI.result[row.names(TXI.result) %in% c("Spore","Polar","Elong","Brch"),]
  TXI.result.germ$Stage<-row.names(TXI.result.germ)
  TXI.result.germ$group<-ref
  TXI.result.germ$statistics<-"FlatLineTest"
  TXI.result.germ$Process<-"Germ"
  TXI.result.germ$p.value<-FlatLineTest.result$p.value
  TXI.result.germ$std.dev<-FlatLineTest.result$std.dev
  TXI.result.list<-rbind(TXI.result.list, TXI.result.germ)
  
  FlatLineTest.result<-FlatLineTest(TMM.matrix.new[c(1:2,7:14)],
                                    permutations = 1000,
                                    plotHistogram = FALSE,
                                    runs = 100,
                                    parallel = 14,
                                    custom.perm.matrix = NULL)
  TXI.result.sex<-TXI.result[row.names(TXI.result) %in% c("Proto","Mate","PIni","PWall",
                                                          "Para","Asci","YP","MP"),]
  TXI.result.sex$Stage<-row.names(TXI.result.sex)
  TXI.result.sex$group<-ref
  TXI.result.sex$statistics<-"FlatLineTest"
  TXI.result.sex$Process<-"Sex"
  TXI.result.sex$p.value<-FlatLineTest.result$p.value
  TXI.result.sex$std.dev<-FlatLineTest.result$std.dev
  TXI.result.list<-rbind(TXI.result.list, TXI.result.sex)
  
  TXI.result.sub<-rbind(TXI.result.germ,TXI.result.sex)
  
  scale.bar<-c(floor(min(TXI.result.sub$TXI.result-TXI.result.sub$std.dev)/0.001-0.5)*0.001,
               ceiling(max(TXI.result.sub$TXI.result+TXI.result.sub$std.dev)/0.001+0.5)*0.001)
  if (max(scale.bar)-min(scale.bar)<=0.005) {
    break.by=0.001
  } else if (max(scale.bar)-min(scale.bar)>0.005 & max(scale.bar)-min(scale.bar)<=0.009) {
    break.by=0.0015
  } else if (max(scale.bar)-min(scale.bar)>0.009 & max(scale.bar)-min(scale.bar)<=0.012) {
    break.by=0.002
  } else if (max(scale.bar)-min(scale.bar)>0.012 & max(scale.bar)-min(scale.bar)<=0.018) {
    break.by=0.003
  } else if (max(scale.bar)-min(scale.bar)>0.018 & max(scale.bar)-min(scale.bar)<=0.024) {
    break.by=0.004
  } else if (max(scale.bar)-min(scale.bar)>0.024) {
    break.by=0.005
  }
  
  TDI<-TXI.result.sub %>%
    mutate(Stage = factor(Stage, levels = plot.order)) %>%
    ggplot(aes(x = Stage, y = TXI.result, group = Process))+
    geom_rect(aes(xmin=nrow(TXI.result.sub)-min(module.list.germ[[2]]-0.5)+1, 
                  xmax=nrow(TXI.result.sub)-max(module.list.germ[[2]]-0.5), 
                  ymin=-Inf, ymax=Inf),
              fill='dodgerblue',alpha = 0.05)+
    geom_rect(aes(xmin=nrow(TXI.result.sub)-min(module.list.sex[[2]]-0.5+4)+1, 
                  xmax=nrow(TXI.result.sub)-max(module.list.sex[[2]]-0.5+4), 
                  ymin=-Inf, ymax=Inf),
              fill='dodgerblue',alpha = 0.05)+
    geom_ribbon(aes(ymin=TXI.result-std.dev, ymax=TXI.result+std.dev), fill = "gray50", alpha = 0.5)+
    geom_line(aes(linetype="solid"))+
    scale_y_continuous(limits = scale.bar,
                       breaks = seq(min(scale.bar),max(scale.bar),break.by),
                       #trans = "reverse", 
                       position = "right")+
    scale_x_discrete(position = "bottom")+
    annotate(geom = "text", 
             label=paste0("Flat Line Test\nP.value=", 
                          signif(unique(TXI.result.germ$p.value),2)),
             x=nrow(TXI.result.sub)-mean(module.list.germ[[2]])+1, y=max(scale.bar)-break.by, 
             colour="black", size=2)+
    annotate(geom = "text", 
             label=paste0("Flat Line Test\nP.value=", 
                          signif(unique(TXI.result.sex$p.value),2)),
             x=nrow(TXI.result.sub)-mean(module.list.sex[[2]]+4)+1, y=max(scale.bar)-break.by, 
             colour="black", size=2)+
    labs(title = "", x = "", y = "Transcriptome Divergence Index", colour = NULL)+
    coord_flip()+
    theme(axis.line = element_line(linetype = "solid"),
          axis.ticks.y = element_line(colour = "black", size = 0.5),
          axis.ticks.x = element_line(colour = "black", size = 0.5),
          axis.text.x = element_text(size = 8, colour = "black"),
          axis.text.y = element_text(size = 8, colour = "black"),
          plot.title = element_text(size = 8, hjust = 0.5, face = "plain"),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 0),
          panel.background = element_rect(fill = NA),
          legend.key = element_rect(fill = NA),
          legend.background = element_rect(fill = NA),
          legend.position = "none",
          strip.text = element_text(size = 9))
  TDI
  ggsave(paste0(species,".TDI.Flatline.",ref,".",refspecies,".png"), width = 4, height = 3, units = "in", dpi = 300)
  f=paste0(species,".TDI.Flatline.",ref,".",refspecies,".pptx")
  topptx(TDI, f, width = 5, height = 3, units = "in")
  
  TDI<-TXI.result.sub %>%
    mutate(Stage = factor(Stage, levels = treat.order)) %>%
    ggplot(aes(x = Stage, y = TXI.result, group = Process))+
    geom_rect(aes(xmin=min(module.list.germ[[2]])-0.5,
                  xmax=max(module.list.germ[[2]])+0.5,
                  ymin=-Inf, ymax=max(scale.bar)-0.5*break.by),
              fill='dodgerblue',alpha = 0.05)+
    geom_rect(aes(xmin=min(module.list.sex[[2]]-0.5+4), 
                  xmax=max(module.list.sex[[2]]+0.5+4), 
                  ymin=-Inf, ymax=max(scale.bar)-0.5*break.by),
              fill='dodgerblue',alpha = 0.05)+
    geom_ribbon(aes(ymin=TXI.result-std.dev, ymax=TXI.result+std.dev), fill = "gray50", alpha = 0.5)+
    geom_line(aes(linetype="solid"))+
    scale_y_continuous(limits = scale.bar,
                       breaks = seq(min(scale.bar),max(scale.bar),break.by),
                       #trans = "reverse", 
                       position = "left")+
    scale_x_discrete(position = "bottom")+
    annotate(geom = "text", 
             label=paste0(gsub("_"," ",refspecies),"\n"),
             x=nrow(TXI.result.sub)/2, y=max(scale.bar), 
             colour="black", size=2, fontface = "italic")+
    annotate(geom = "text", 
             label=paste0("\n(Flat Line Test)"),
             x=nrow(TXI.result.sub)/2, y=max(scale.bar), 
             colour="black", size=2)+
    annotate(geom = "text", 
             label=paste0("P.value=", 
                          signif(unique(TXI.result.germ$p.value),2)),
             x=mean(module.list.germ[[2]]), y=max(scale.bar)-1.25*break.by, 
             colour="black", size=2)+
    annotate(geom = "text", 
             label=paste0("P.value=", 
                          signif(unique(TXI.result.sex$p.value),2)),
             x=mean(module.list.germ[[2]])+4, y=max(scale.bar)-1.25*break.by, 
             colour="black", size=2)+
    labs(title = "", x = "", y = "Transcriptome Divergence Index", colour = NULL)+
    theme(axis.line = element_line(linetype = "solid"),
          axis.ticks.y = element_line(colour = "black", size = 0.5),
          axis.ticks.x = element_line(colour = "black", size = 0.5),
          axis.text.x = element_text(size = 8, colour = "black"),
          axis.text.y = element_text(size = 8, colour = "black"),
          plot.title = element_text(size = 8, hjust = 0.5, face = "plain"),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 0),
          panel.background = element_rect(fill = NA),
          legend.key = element_rect(fill = NA),
          legend.background = element_rect(fill = NA),
          legend.position = "none",
          strip.text = element_text(size = 9))
  TDI
  ggsave(paste0(species,".TDI.Flatline.",ref,".",refspecies,".B.png"), 
         width = 5, height = 3, units = "in", dpi = 300)
  rm(TXI.result.sub)
  
  ReductiveHourglassTest.result<-ReductiveHourglassTest(TMM.matrix.new[1:6],
                                                        modules = module.list.germ,
                                                        permutations = 1000,
                                                        lillie.test = FALSE,
                                                        plotHistogram = FALSE,
                                                        runs = 100,
                                                        parallel = 14,
                                                        gof.warning = FALSE,
                                                        custom.perm.matrix = NULL)
  TXI.result.germ<-TXI.result[row.names(TXI.result) %in% c("Spore","Polar","Elong","Brch"),]
  TXI.result.germ$Stage<-row.names(TXI.result.germ)
  TXI.result.germ$group<-ref
  TXI.result.germ$statistics<-"ReductiveHourglassTest"
  TXI.result.germ$Process<-"Germ"
  TXI.result.germ$p.value<-ReductiveHourglassTest.result$p.value
  TXI.result.germ$std.dev<-ReductiveHourglassTest.result$std.dev
  TXI.result.list<-rbind(TXI.result.list, TXI.result.germ)
  
  ReductiveHourglassTest.result<-ReductiveHourglassTest(TMM.matrix.new[c(1:2,7:14)],
                                                        modules = module.list.sex,
                                                        permutations = 1000,
                                                        lillie.test = FALSE,
                                                        plotHistogram = FALSE,
                                                        runs = 100,
                                                        parallel = 14,
                                                        gof.warning = FALSE,
                                                        custom.perm.matrix = NULL)
  
  TXI.result.sex<-TXI.result[row.names(TXI.result) %in% c("Proto","Mate","PIni","PWall",
                                                          "Para","Asci","YP","MP"),]
  TXI.result.sex$Stage<-row.names(TXI.result.sex)
  TXI.result.sex$group<-ref
  TXI.result.sex$statistics<-"ReductiveHourglassTest"
  TXI.result.sex$Process<-"Sex"
  TXI.result.sex$p.value<-ReductiveHourglassTest.result$p.value
  TXI.result.sex$std.dev<-ReductiveHourglassTest.result$std.dev
  TXI.result.list<-rbind(TXI.result.list, TXI.result.sex)
  
  TXI.result.sub<-rbind(TXI.result.germ, TXI.result.sex)
  
  TDI<-TXI.result.sub %>%
    mutate(Stage = factor(Stage, levels = plot.order)) %>%
    ggplot(aes(x = Stage, y = TXI.result, group = Process))+
    geom_rect(aes(xmin=nrow(TXI.result.sub)-min(module.list.germ[[2]]-0.5)+1, 
                  xmax=nrow(TXI.result.sub)-max(module.list.germ[[2]]-0.5), 
                  ymin=-Inf, ymax=Inf),
              fill='dodgerblue',alpha = 0.05)+
    geom_rect(aes(xmin=nrow(TXI.result.sub)-min(module.list.sex[[2]]-0.5+4)+1, 
                  xmax=nrow(TXI.result.sub)-max(module.list.sex[[2]]-0.5+4), 
                  ymin=-Inf, ymax=Inf),
              fill='dodgerblue',alpha = 0.05)+
    geom_ribbon(aes(ymin=TXI.result-std.dev, ymax=TXI.result+std.dev), fill = "gray50", alpha = 0.5)+
    geom_line(aes(linetype="solid"))+
    scale_y_continuous(limits = scale.bar,
                       breaks = seq(min(scale.bar),max(scale.bar),break.by),
                       #trans = "reverse", 
                       position = "right")+
    scale_x_discrete(position = "bottom")+
    annotate(geom = "text", 
             label=paste0("Reductive Hourglass Test\n(High-Low-High)\nP.value=", 
                          signif(unique(TXI.result.germ$p.value),2)),
             x=nrow(TXI.result.sub)-mean(module.list.germ[[2]])+1, y=max(scale.bar)-break.by, 
             colour="black", size=2)+
    annotate(geom = "text", 
             label=paste0("Reductive Hourglass Test\n(High-Low-High)\nP.value=", 
                          signif(unique(TXI.result.sex$p.value),2)),
             x=nrow(TXI.result.sub)-mean(module.list.sex[[2]]+4)+1, y=max(scale.bar)-break.by, 
             colour="black", size=2)+
    labs(title = "", x = "", y = "Transcriptome Divergence Index", colour = NULL)+
    coord_flip()+
    theme(axis.line = element_line(linetype = "solid"),
          axis.ticks.y = element_line(colour = "black", size = 0.5),
          axis.ticks.x = element_line(colour = "black", size = 0.5),
          axis.text.x = element_text(size = 8, colour = "black"),
          axis.text.y = element_text(size = 8, colour = "black"),
          plot.title = element_text(size = 8, hjust = 0.5, face = "plain"),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 0),
          panel.background = element_rect(fill = NA),
          legend.key = element_rect(fill = NA),
          legend.background = element_rect(fill = NA),
          legend.position = "none",
          strip.text = element_text(size = 9))
  TDI
  ggsave(paste0(species,".TDI.reductive.",ref,".",refspecies,".png"),
         width = 4, height = 3, units = "in", dpi = 300)
  f=paste0(species,".TDI.reductive.",ref,".",refspecies,".pptx")
  topptx(TDI, f, width = 4, height = 3, units = "in")
  
  TDI<-TXI.result.sub %>%
    mutate(Stage = factor(Stage, levels = treat.order)) %>%
    ggplot(aes(x = Stage, y = TXI.result, group = Process))+
    geom_rect(aes(xmin=min(module.list.germ[[2]])-0.5,
                  xmax=max(module.list.germ[[2]])+0.5,
                  ymin=-Inf, ymax=max(scale.bar)-break.by),
              fill='dodgerblue',alpha = 0.05)+
    geom_rect(aes(xmin=min(module.list.sex[[2]]-0.5+4), 
                  xmax=max(module.list.sex[[2]]+0.5+4), 
                  ymin=-Inf, ymax=max(scale.bar)-break.by),
              fill='dodgerblue',alpha = 0.05)+
    geom_ribbon(aes(ymin=TXI.result-std.dev, ymax=TXI.result+std.dev), fill = "gray50", alpha = 0.5)+
    geom_line(aes(linetype="solid"))+
    scale_y_continuous(limits = scale.bar,
                       breaks = seq(min(scale.bar),max(scale.bar),break.by),
                       #trans = "reverse", 
                       position = "left")+
    scale_x_discrete(position = "bottom")+
    annotate(geom = "text", 
             label=paste0(gsub("_"," ",refspecies),"\n"),
             x=nrow(TXI.result.sub)/2, y=max(scale.bar), 
             colour="black", size=2, fontface = "italic")+
    annotate(geom = "text", 
             label=paste0("\n(Reductive Hourglass Test)"),
             x=nrow(TXI.result.sub)/2, y=max(scale.bar), 
             colour="black", size=2)+
    annotate(geom = "text", 
             label=paste0("P.value=", 
                          signif(unique(TXI.result.germ$p.value),2)),
             x=mean(module.list.germ[[2]]), y=max(scale.bar)-1.25*break.by, 
             colour="black", size=2)+
    annotate(geom = "text", 
             label=paste0("P.value=", 
                          signif(unique(TXI.result.sex$p.value),2)),
             x=mean(module.list.germ[[2]])+4, y=max(scale.bar)-1.25*break.by, 
             colour="black", size=2)+
    labs(title = "", x = "", y = "Transcriptome Divergence Index", colour = NULL)+
    theme(axis.line = element_line(linetype = "solid"),
          axis.ticks.y = element_line(colour = "black", size = 0.5),
          axis.ticks.x = element_line(colour = "black", size = 0.5),
          axis.text.x = element_text(size = 8, colour = "black"),
          axis.text.y = element_text(size = 8, colour = "black"),
          plot.title = element_text(size = 8, hjust = 0.5, face = "plain"),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 0),
          panel.background = element_rect(fill = NA),
          legend.key = element_rect(fill = NA),
          legend.background = element_rect(fill = NA),
          legend.position = "none",
          strip.text = element_text(size = 9))
  TDI
  ggsave(paste0(species,".TDI.reductive.",ref,".",refspecies,".B.png"),
         width = 5, height = 3, units = "in", dpi = 300)
  f=paste0(species,".TDI.reductive.",ref,".",refspecies,".B.pptx")
  topptx(TDI, f, width = 5, height = 3, units = "in")
  rm(TXI.result.sub)
  
  write.table(TXI.result.list, 
              paste0(species,".NS.",refspecies,".txt"),
              quote = F, row.names = F, sep = "\t")
}

