setwd("hourglass/Coprinopsis_myTAI/")
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

plot.order<-c("YFB","Pri","Knot",
              "Myc",
              "Brch","Ger","BS",
              "Oidia","Scl")
taxoncolor<-c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
              "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
              "#BC80BD", "#CCEBC5", "#FFED6F", "#D9D9D9",  "white")
####All
ref="All"
treat.order<-c("Scl", "Oidia",
               "BS", "Ger", "Brch",
               "Myc", "Knot", "Pri", "YFB")
module.list<-list(early = 1:3, mid = 4:5, late = 6:9)

##Read in gene expression level
TMM.matrix<-read.delim("Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.log2TMMmatrixmean.txt", header = TRUE)
TMM.matrix<-TMM.matrix[,c("Genes", treat.order)]

##########################################################
####TAI feature
##########################################################

##Read in TAI features from CC2G
species="Coprinopsis_cinerea_A43mutB43mut_pab1-1_326"
genes.PS<-read.delim(paste0("/home/yichun/RNAmodification/hourglass/blast/",
                            species,"/e10-5aa30/",
                            species,".proteins.fa.PSfinal.txt"), header = T)
names(genes.PS)[names(genes.PS)=="Gene"]<-"Genes"
#genes.PS$PS[genes.PS$PS == 11]<-10
#genes.PS$PS[genes.PS$PS == 12]<-10
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

taxoname<-read.delim("Coprinopsis_cinerea.taxon.level.txt", header = T)
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

f=paste0("Coprinopsis_cinerea.PS.tree.pptx")
topptx(p, f, width = 5.8, height = 5, units = "in")
rm(TXI.result.sub)

##TAI analysis
genes.hour<-genes.PS

##Merge PS to each table
TMM.matrix.new<-merge(TMM.matrix, genes.hour, by = "Genes", all = F)
TMM.matrix.new<-TMM.matrix.new[,c("PS","Genes",treat.order)]
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
TXI.result.list<-as.data.frame(matrix(NA, ncol = 6, nrow = 0))
names(TXI.result.list)<-c("Stage","TXI.result","group","p.value","std.dev","statistics")

##############################################################################################
##Pvalue/contribution/ratio for each combination
##############################################################################################
####All
ref="All"
treat.order<-c("Scl", "Oidia",
               "BS", "Ger", "Brch",
               "Myc", "Knot", "Pri", "YFB")
module.list<-list(early = 1:4, mid = 5, late = 6:9)
FlatLineTest.result<-FlatLineTest(TMM.matrix.new,
                                  permutations = 1000,
                                  plotHistogram = FALSE,
                                  runs = 100,
                                  parallel = 14,
                                  custom.perm.matrix = NULL)
TXI.result.sub<-TXI.result[row.names(TXI.result) %in% treat.order,]
TXI.result.sub$Stage<-row.names(TXI.result.sub)
TXI.result.sub$group<-ref
TXI.result.sub$statistics<-"FlatLineTest"
TXI.result.sub$p.value<-FlatLineTest.result$p.value
TXI.result.sub$std.dev<-FlatLineTest.result$std.dev
TXI.result.list<-rbind(TXI.result.list, TXI.result.sub)
TAI<-TXI.result.sub %>%
  mutate(Stage = factor(Stage, levels = plot.order)) %>%
  ggplot(aes(x = Stage, y = TXI.result, group = 1))+
  geom_rect(aes(xmin=nrow(TXI.result.sub)-min(module.list[[2]]-0.5)+1,
                xmax=nrow(TXI.result.sub)-max(module.list[[2]]-0.5),
                ymin=-Inf, ymax=Inf),
            fill='dodgerblue',alpha = 0.05)+
  geom_ribbon(aes(ymin=TXI.result-std.dev, ymax=TXI.result+std.dev), fill = "gray50", alpha = 0.5)+
  geom_line(aes(linetype="solid"))+
  scale_y_continuous(limits = c(3.55,3.05),
                     breaks = seq(3.05,3.55,0.1),
                     trans = "reverse", position = "right")+
  scale_x_discrete(position = "top")+
  annotate(geom = "text",
           label=paste0("Flat Line Test\nP.value=",
                        signif(unique(TXI.result.sub$p.value),2)),
           x=nrow(TXI.result.sub)-min(module.list[[2]])+1, y=3.45,
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
ggsave(paste0("Coprinopsis_cinerea.TAI.Flatline.",ref,".png"), width = 4, height = 3, units = "in", dpi = 300)
f=paste0("Coprinopsis_cinerea.TAI.Flatline.",ref,".pptx")
topptx(TAI, f, width = 4, height = 3, units = "in")
rm(TXI.result.sub)

##ReductiveHourglassTest (High-Low-High)
ReductiveHourglassTest.result<-ReductiveHourglassTest(TMM.matrix.new,
                                                      modules = module.list,
                                                      permutations = 1000,
                                                      lillie.test = FALSE,
                                                      plotHistogram = FALSE,
                                                      runs = 100,
                                                      parallel = 14,
                                                      gof.warning = FALSE,
                                                      custom.perm.matrix = NULL)
TXI.result.sub<-TXI.result[row.names(TXI.result) %in% treat.order,]
TXI.result.sub$Stage<-row.names(TXI.result.sub)
TXI.result.sub$group<-ref
TXI.result.sub$statistics<-"ReductiveHourglassTest"
TXI.result.sub$p.value<-ReductiveHourglassTest.result$p.value
TXI.result.sub$std.dev<-ReductiveHourglassTest.result$std.dev
TXI.result.list<-rbind(TXI.result.list, TXI.result.sub)

TAI<-TXI.result.sub %>%
  mutate(Stage = factor(Stage, levels = plot.order)) %>%
  ggplot(aes(x = Stage, y = TXI.result, group = 1))+
  geom_rect(aes(xmin=nrow(TXI.result.sub)-min(module.list[[2]]-0.5)+1,
                xmax=nrow(TXI.result.sub)-max(module.list[[2]]-0.5),
                ymin=-Inf, ymax=Inf),
            fill='dodgerblue',alpha = 0.05)+
  geom_ribbon(aes(ymin=TXI.result-std.dev, ymax=TXI.result+std.dev), fill = "gray50", alpha = 0.5)+
  geom_line(aes(linetype="solid"))+
  scale_y_continuous(limits = c(3.55,3.05),
                     breaks = seq(3.05,3.55,0.1),
                     trans = "reverse", position = "right")+
  scale_x_discrete(position = "top")+
  annotate(geom = "text",
           label=paste0("Reductive Hourglass Test (High-Low-High)\nP.value=",
                                       signif(unique(TXI.result.sub$p.value),2)),
           x=nrow(TXI.result.sub)-min(module.list[[2]])+1, y=3.45,
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
ggsave(paste0("Coprinopsis_cinerea.TAI.reductive.",ref,".png"), width = 4, height = 3, units = "in", dpi = 300)
f=paste0("Coprinopsis_cinerea.TAI.reductive.",ref,".pptx")
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
REFoldChangeOfMeanREValues<-(FoldChangeOfMeanREValues-
                               min(FoldChangeOfMeanREValues))/(max(FoldChangeOfMeanREValues)-
                                                                 min(FoldChangeOfMeanREValues))

MeanREClassValues<-rbind(REmatrix,MeanREClassValues,REFoldChangeOfMeanREValues,StdErr.RE.ClassValues)
MeanREClassValues<-as.data.frame(MeanREClassValues)
names(MeanREClassValues)<-names(TMM.matrix.new[3:ncol(TMM.matrix.new)])
row.names(MeanREClassValues)<-c(paste0("PS",1:12),"PS1-2","PS3-12","Ratio","stderr.PS1-2","stderr.PS3-12")
write.table(MeanREClassValues, paste0("Coprinopsis_cinerea.PS.RE.",ref,".txt"),
            row.names = T, quote = F, sep ="\t")

##contribution
percentTAI<-as.data.frame(pStrata(TMM.matrix.new))
row.names(percentTAI)<-paste0("PS", 1:12)
write.table(percentTAI, paste0("Coprinopsis_cinerea.PS.contribution.",ref,".txt"),
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

##############################################################################################
####Sexual cycle
##############################################################################################
ref="Sex"
treat.order<-c("BS", "Ger", "Brch",
               "Myc","Knot","Pri", "YFB")
module.list<-list(early = 1:2, mid = 3, late = 4:7)

FlatLineTest.result<-FlatLineTest(TMM.matrix.new[,c("PS","Genes",treat.order)],
                                  permutations = 1000,
                                  plotHistogram = FALSE,
                                  runs = 100,
                                  parallel = 14,
                                  custom.perm.matrix = NULL)
TXI.result.sub<-TXI.result[row.names(TXI.result) %in% treat.order,]
TXI.result.sub$Stage<-row.names(TXI.result.sub)
TXI.result.sub$group<-ref
TXI.result.sub$statistics<-"FlatLineTest"
TXI.result.sub$p.value<-FlatLineTest.result$p.value
TXI.result.sub$std.dev<-FlatLineTest.result$std.dev
TXI.result.list<-rbind(TXI.result.list, TXI.result.sub)
TAI<-TXI.result.sub %>%
  mutate(Stage = factor(Stage, levels = plot.order)) %>%
  ggplot(aes(x = Stage, y = TXI.result, group = 1))+
  geom_rect(aes(xmin=nrow(TXI.result.sub)-min(module.list[[2]]-0.5)+1,
                xmax=nrow(TXI.result.sub)-max(module.list[[2]]-0.5),
                ymin=-Inf, ymax=Inf),
            fill='dodgerblue',alpha = 0.05)+
  geom_ribbon(aes(ymin=TXI.result-std.dev, ymax=TXI.result+std.dev), fill = "gray50", alpha = 0.5)+
  geom_line(aes(linetype="solid"))+
  scale_y_continuous(limits = c(3.55,3.05),
                     breaks = seq(3.05,3.55,0.1),
                     trans = "reverse", position = "right")+
  scale_x_discrete(position = "top")+
  annotate(geom = "text",
           label=paste0("Flat Line Test\nP.value=",
                        signif(unique(TXI.result.sub$p.value),2)),
           x=nrow(TXI.result.sub)-min(module.list[[2]])+1, y=3.45,
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
ggsave(paste0("Coprinopsis_cinerea.TAI.Flatline.",ref,".png"), width = 4, height = 2.5, units = "in", dpi = 300)
f=paste0("Coprinopsis_cinerea.TAI.Flatline.",ref,".pptx")
topptx(TAI, f, width = 4, height = 2.5, units = "in")
rm(TXI.result.sub)

##ReductiveHourglassTest (High-Low-High)
ReductiveHourglassTest.result<-ReductiveHourglassTest(TMM.matrix.new[,c("PS","Genes",treat.order)],
                                                      modules = module.list,
                                                      permutations = 1000,
                                                      lillie.test = FALSE,
                                                      plotHistogram = FALSE,
                                                      runs = 100,
                                                      parallel = 14,
                                                      gof.warning = FALSE,
                                                      custom.perm.matrix = NULL)
TXI.result.sub<-TXI.result[row.names(TXI.result) %in% treat.order,]
TXI.result.sub$Stage<-row.names(TXI.result.sub)
TXI.result.sub$group<-ref
TXI.result.sub$statistics<-"ReductiveHourglassTest"
TXI.result.sub$p.value<-ReductiveHourglassTest.result$p.value
TXI.result.sub$std.dev<-ReductiveHourglassTest.result$std.dev
TXI.result.list<-rbind(TXI.result.list, TXI.result.sub)

TAI<-TXI.result.sub %>%
  mutate(Stage = factor(Stage, levels = plot.order)) %>%
  ggplot(aes(x = Stage, y = TXI.result, group = 1))+
  geom_rect(aes(xmin=nrow(TXI.result.sub)-min(module.list[[2]]-0.5)+1,
                xmax=nrow(TXI.result.sub)-max(module.list[[2]]-0.5),
                ymin=-Inf, ymax=Inf),
            fill='dodgerblue',alpha = 0.05)+
  geom_ribbon(aes(ymin=TXI.result-std.dev, ymax=TXI.result+std.dev), fill = "gray50", alpha = 0.5)+
  geom_line(aes(linetype="solid"))+
  scale_y_continuous(limits = c(3.55,3.05),
                     breaks = seq(3.05,3.55,0.1),
                     trans = "reverse", position = "right")+
  scale_x_discrete(position = "top")+
  annotate(geom = "text",
           label=paste0("Reductive Hourglass Test (High-Low-High)\nP.value=",
                        signif(unique(TXI.result.sub$p.value),2)),
           x=nrow(TXI.result.sub)-min(module.list[[2]])+1, y=3.45,
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
ggsave(paste0("Coprinopsis_cinerea.TAI.reductive.",ref,".png"), width = 4, height = 2.5, units = "in", dpi = 300)
f=paste0("Coprinopsis_cinerea.TAI.reductive.",ref,".pptx")
topptx(TAI, f, width = 4, height = 2.5, units = "in")

##Relative expression of all stages by PS
Relative.exp<-as.data.frame(REMatrix(TMM.matrix.new[,c("PS","Genes",treat.order)]))

##Group ratio
REmatrix <- age.apply(ExpressionSet = TMM.matrix.new[,c("PS","Genes",treat.order)], RE)
Groups = list(group_1 = 1:2, group_2 = 3:12)
MeanREClassValues <- matrix(NA_real_,length(Groups),ncol(TMM.matrix.new[,c("PS","Genes",treat.order)])-2)
StdErr.RE.ClassValues <- matrix(NA_real_,length(Groups),ncol(TMM.matrix.new[,c("PS","Genes",treat.order)])-2)
for(i in 1:length(Groups)){
  MeanREClassValues[i , ] <- colMeans(REmatrix[match(as.character(Groups[[i]]),
                                                     rownames(REmatrix)) , ])

  StdErr.RE.ClassValues[i , ] <- apply(REmatrix[match(as.character(Groups[[i]]),
                                                      rownames(REmatrix)) , ],2,std.error)
}
FoldChangeOfMeanREValues<-apply(MeanREClassValues,2,function(x){return((x[1])/(x[2]))})
REFoldChangeOfMeanREValues<-(FoldChangeOfMeanREValues-
                               min(FoldChangeOfMeanREValues))/(max(FoldChangeOfMeanREValues)-
                                                                 min(FoldChangeOfMeanREValues))

MeanREClassValues<-rbind(REmatrix,MeanREClassValues,REFoldChangeOfMeanREValues,StdErr.RE.ClassValues)
MeanREClassValues<-as.data.frame(MeanREClassValues)
names(MeanREClassValues)<-names(TMM.matrix.new[,(2+3):ncol(TMM.matrix.new)])
row.names(MeanREClassValues)<-c(paste0("PS",1:12),"PS1-2","PS3-12","Ratio","stderr.PS1-2","stderr.PS3-12")
write.table(MeanREClassValues, paste0("Coprinopsis_cinerea.PS.RE.",ref,".txt"),
            row.names = T, quote = F, sep ="\t")

##contribution
percentTAI<-as.data.frame(pStrata(TMM.matrix.new[,c("PS","Genes",treat.order)]))
row.names(percentTAI)<-paste0("PS", 1:12)
write.table(percentTAI, paste0("Coprinopsis_cinerea.PS.contribution.",ref,".txt"),
            row.names = T, quote = F, sep ="\t")

PlotContribution( ExpressionSet = TMM.matrix.new[,c("PS","Genes",treat.order)],
                  legendName = "PS",
                  xlab = "Ontogeny",
                  ylab = "Transcriptome Age Index",
                  y.ticks = 10)
RE.bar<-PlotBarRE(TMM.matrix.new[,c("PS","Genes",treat.order)],
                  Groups = list(group_1 = 1:2, group_2 = 3:12),
                  ratio = T,
                  p.adjust.method = "BH")


##########################################################
####TDI feature
##########################################################
pairspecies<-c("Agaricus_bisporus",
               "Boletus_edulis",
               "Botryobasidium_botryosum",
               "Calocera_viscosa",
               "Coprinellus_pellucidus",
               "Coprinopsis_marcescibilis",
               "Coprinopsis_sclerotiger",
               "Cryptococcus_neoformans",
               "Laccaria_bicolor",
               "Lentinula_edodes",
               "Pleurotus_ostreatus",
               "Puccinia_triticina",
               "Sphaerobolus_stellatus",
               "Ustilago_maydis",
               "Wolfiporia_cocos")
species="Coprinopsis_cinerea_A43mutB43mut_pab1-1_326"
for (k in 1:length(pairspecies)) {
  ##Read in TDI features from CC2G
  refspecies=pairspecies[k]
  genes.NS<-read.delim(paste0("/home/yichun/RNAmodification/hourglass/blast/",
                              species,
                              "/dnds/Coprinopsis_cinerea_VS_",
                              refspecies,".NSfinal.txt"))
  names(genes.NS)[names(genes.NS)=="dNdS"]<-"NS"
  genes.NS<-separate(genes.NS, "Genes", c("Genes"), sep = "T")
  genes.hour<-genes.NS[is.na(genes.NS$NS)==F & genes.NS$NS<=2,c("Genes","NS")]

  p<-ggplot(data = genes.hour, aes(x = NS, y = stat(density)/sum(stat(density))))+
    geom_histogram(color = "black", fill = "black", binwidth = 0.01)+
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2))+
    scale_y_continuous(limits = c(0,0.12),
                       breaks = seq(0,0.12, 0.02))+
    annotate(geom = "text",
             label=paste0(gsub("_"," ",refspecies)),
             x=0.5, y=0.11,
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
  ggsave(paste0("Coprinopsis_cinerea.dNdS.",refspecies,".png"), width = 3, height = 3, units = "in", dpi = 300)

  ####All
  ref="All"
  treat.order<-c("Scl", "Oidia",
                 "BS", "Ger", "Brch",
                 "Myc", "Knot", "Pri", "YFB")
  module.list<-list(early = 1:4, mid = 5, late = 6:9)
  ##Read in gene expression level
  TMM.matrix<-read.delim("Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.log2TMMmatrixmean.txt", header = TRUE)
  TMM.matrix<-TMM.matrix[,c("Genes", treat.order)]
  ##Merge dNdS to each table
  TMM.matrix.new<-merge(TMM.matrix, genes.hour, by = "Genes", all = F)
  TMM.matrix.new<-TMM.matrix.new[,c("NS","Genes",treat.order)]

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
  TXI.result.list<-as.data.frame(matrix(NA, ncol = 6, nrow = 0))
  names(TXI.result.list)<-c("Stage","TXI.result","group","p.value","std.dev","statistics")

  ##Pvalue for each combination
  ####All
  FlatLineTest.result<-FlatLineTest(TMM.matrix.new,
                                    permutations = 1000,
                                    plotHistogram = FALSE,
                                    runs = 100,
                                    parallel = 14,
                                    custom.perm.matrix = NULL)
  TXI.result.sub<-TXI.result[row.names(TXI.result) %in% treat.order,]
  TXI.result.sub$Stage<-row.names(TXI.result.sub)
  TXI.result.sub$group<-ref
  TXI.result.sub$statistics<-"FlatLineTest"
  TXI.result.sub$p.value<-FlatLineTest.result$p.value
  TXI.result.sub$std.dev<-FlatLineTest.result$std.dev
  TXI.result.list<-rbind(TXI.result.list, TXI.result.sub)

  scale.bar<-c(floor(min(TXI.result.sub$TXI.result+TXI.result.sub$std.dev)/0.001-0.5)*0.001,
               ceiling(max(TXI.result.sub$TXI.result+TXI.result.sub$std.dev)/0.001+0.5)*0.001)

  TDI<-TXI.result.sub %>%
    mutate(Stage = factor(Stage, levels = plot.order)) %>%
    ggplot(aes(x = Stage, y = TXI.result, group = 1))+
    geom_rect(aes(xmin=nrow(TXI.result.sub)-min(module.list[[2]]-0.5)+1,
                  xmax=nrow(TXI.result.sub)-max(module.list[[2]]-0.5),
                  ymin=-Inf, ymax=Inf),
              fill='dodgerblue',alpha = 0.05)+
    geom_ribbon(aes(ymin=TXI.result-std.dev, ymax=TXI.result+std.dev), fill = "gray50", alpha = 0.5)+
    geom_line(aes(linetype="solid"))+
    scale_y_continuous(limits = scale.bar,
                       breaks = seq(min(scale.bar),max(scale.bar),0.001),
                       #trans = "reverse",
                       position = "right")+
    scale_x_discrete(position = "bottom")+
    annotate(geom = "text",
             label=paste0("Flat Line Test\nP.value=",
                          signif(unique(TXI.result.sub$p.value),2)),
             x=nrow(TXI.result.sub)-min(module.list[[2]])+1, y=max(scale.bar)-0.001,
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
  ggsave(paste0("Coprinopsis_cinerea.TDI.Flatline.",ref,".",refspecies,".png"), width = 4, height = 3, units = "in", dpi = 300)
  f=paste0("Coprinopsis_cinerea.TDI.Flatline.",ref,".",refspecies,".pptx")
  topptx(TDI, f, width = 4, height = 3, units = "in")

  TDI<-TXI.result.sub %>%
    mutate(Stage = factor(Stage, levels = treat.order)) %>%
    ggplot(aes(x = Stage, y = TXI.result, group = 1))+
    geom_rect(aes(xmin=min(module.list[[2]])-0.5,
                  xmax=max(module.list[[2]])+0.5,
                  ymin=-Inf, ymax=max(scale.bar)-0.001),
              fill='dodgerblue',alpha = 0.05)+
    geom_ribbon(aes(ymin=TXI.result-std.dev, ymax=TXI.result+std.dev), fill = "gray50", alpha = 0.5)+
    geom_line(aes(linetype="solid"))+
    scale_y_continuous(limits = scale.bar,
                       breaks = seq(min(scale.bar),max(scale.bar),0.001),
                       #trans = "reverse",
                       position = "left")+
    scale_x_discrete(position = "bottom")+
    annotate(geom = "text",
             label=paste0(gsub("_"," ",refspecies),"\n"),
             x=nrow(TXI.result.sub)/2, y=max(scale.bar),
             colour="black", size=2, fontface = "italic")+
    annotate(geom = "text",
             label=paste0("\n(Flat Line Test: P.value=",
                          signif(unique(TXI.result.sub$p.value),2),")"),
             x=nrow(TXI.result.sub)/2, y=max(scale.bar),
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
  ggsave(paste0("Coprinopsis_cinerea.TDI.Flatline.",ref,".",refspecies,".B.png"),
         width = 4, height = 3, units = "in", dpi = 300)
  rm(TXI.result.sub)

  ReductiveHourglassTest.result<-ReductiveHourglassTest(TMM.matrix.new,
                                                        modules = module.list,
                                                        permutations = 1000,
                                                        lillie.test = FALSE,
                                                        plotHistogram = FALSE,
                                                        runs = 100,
                                                        parallel = 14,
                                                        gof.warning = FALSE,
                                                        custom.perm.matrix = NULL)
  TXI.result.sub<-TXI.result[row.names(TXI.result) %in% treat.order,]
  TXI.result.sub$Stage<-row.names(TXI.result.sub)
  TXI.result.sub$group<-ref
  TXI.result.sub$statistics<-"ReductiveHourglassTest"
  TXI.result.sub$p.value<-ReductiveHourglassTest.result$p.value
  TXI.result.sub$std.dev<-ReductiveHourglassTest.result$std.dev
  TXI.result.list<-rbind(TXI.result.list, TXI.result.sub)
  TDI<-TXI.result.sub %>%
    mutate(Stage = factor(Stage, levels = plot.order)) %>%
    ggplot(aes(x = Stage, y = TXI.result, group = 1))+
    geom_rect(aes(xmin=nrow(TXI.result.sub)-min(module.list[[2]]-0.5)+1,
                  xmax=nrow(TXI.result.sub)-max(module.list[[2]]-0.5),
                  ymin=-Inf, ymax=Inf),
              fill='dodgerblue',alpha = 0.05)+
    geom_ribbon(aes(ymin=TXI.result-std.dev, ymax=TXI.result+std.dev), fill = "gray50", alpha = 0.5)+
    geom_line(aes(linetype="solid"))+
    scale_y_continuous(limits = scale.bar,
                       breaks = seq(min(scale.bar),max(scale.bar),0.001),
                       #trans = "reverse",
                       position = "right")+
    scale_x_discrete(position = "bottom")+
    annotate(geom = "text",
             label=paste0("Reductive Hourglass Test (High-Low-High)\nP.value=",
                          signif(unique(TXI.result.sub$p.value),2)),
             x=nrow(TXI.result.sub)-min(module.list[[2]])+1, y=max(scale.bar)-0.002,
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
  ggsave(paste0("Coprinopsis_cinerea.TDI.reductive.",ref,".",refspecies,".png"),
         width = 4, height = 3, units = "in", dpi = 300)
  f=paste0("Coprinopsis_cinerea.TDI.reductive.",ref,".",refspecies,".pptx")
  topptx(TDI, f, width = 4, height = 3, units = "in")

  TDI<-TXI.result.sub %>%
    mutate(Stage = factor(Stage, levels = treat.order)) %>%
    ggplot(aes(x = Stage, y = TXI.result, group = 1))+
    geom_rect(aes(xmin=min(module.list[[2]])-0.5,
                  xmax=max(module.list[[2]])+0.5,
                  ymin=-Inf, ymax=max(scale.bar)-0.001),
              fill='dodgerblue',alpha = 0.05)+
    geom_ribbon(aes(ymin=TXI.result-std.dev, ymax=TXI.result+std.dev), fill = "gray50", alpha = 0.5)+
    geom_line(aes(linetype="solid"))+
    scale_y_continuous(limits = scale.bar,
                       breaks = seq(min(scale.bar),max(scale.bar),0.001),
                       #trans = "reverse",
                       position = "left")+
    scale_x_discrete(position = "bottom")+
    annotate(geom = "text",
             label=paste0(gsub("_"," ",refspecies),"\n"),
             x=nrow(TXI.result.sub)/2, y=max(scale.bar),
             colour="black", size=2, fontface = "italic")+
    annotate(geom = "text",
             label=paste0("\n",
                          "(Reductive Hourglass Test: P.value=",
                          signif(unique(TXI.result.sub$p.value),2),")"),
             x=nrow(TXI.result.sub)/2, y=max(scale.bar),
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
  ggsave(paste0("Coprinopsis_cinerea.TDI.reductive.",ref,".",refspecies,".B.png"),
         width = 4, height = 3, units = "in", dpi = 300)
  f=paste0("Coprinopsis_cinerea.TDI.reductive.",ref,".",refspecies,".B.pptx")
  topptx(TDI, f, width = 4, height = 3, units = "in")
  rm(TXI.result.sub)

  ####Sexual cycle
  ref="Sex"
  treat.order<-c("BS", "Ger", "Brch",
                 "Myc", "Knot", "Pri", "YFB")
  module.list<-list(early = 1:2, mid = 3, late = 4:7)

  FlatLineTest.result<-FlatLineTest(TMM.matrix.new[,c("NS","Genes",treat.order)],
                                    permutations = 1000,
                                    plotHistogram = FALSE,
                                    runs = 100,
                                    parallel = 14,
                                    custom.perm.matrix = NULL)
  TXI.result.sub<-TXI.result[row.names(TXI.result) %in% treat.order,]
  TXI.result.sub$Stage<-row.names(TXI.result.sub)
  TXI.result.sub$group<-ref
  TXI.result.sub$statistics<-"FlatLineTest"
  TXI.result.sub$p.value<-FlatLineTest.result$p.value
  TXI.result.sub$std.dev<-FlatLineTest.result$std.dev
  TXI.result.list<-rbind(TXI.result.list, TXI.result.sub)
  scale.bar<-c(floor(min(TXI.result.sub$TXI.result+TXI.result.sub$std.dev)/0.001-0.5)*0.001,
               ceiling(max(TXI.result.sub$TXI.result+TXI.result.sub$std.dev)/0.001+0.5)*0.001)

  TDI<-TXI.result.sub %>%
    mutate(Stage = factor(Stage, levels = plot.order)) %>%
    ggplot(aes(x = Stage, y = TXI.result, group = 1))+
    geom_rect(aes(xmin=nrow(TXI.result.sub)-min(module.list[[2]]-0.5)+1,
                  xmax=nrow(TXI.result.sub)-max(module.list[[2]]-0.5),
                  ymin=-Inf, ymax=Inf),
              fill='dodgerblue',alpha = 0.05)+
    geom_ribbon(aes(ymin=TXI.result-std.dev, ymax=TXI.result+std.dev), fill = "gray50", alpha = 0.5)+
    geom_line(aes(linetype="solid"))+
    scale_y_continuous(limits = scale.bar,
                       breaks = seq(min(scale.bar),max(scale.bar),0.001),
                       #trans = "reverse",
                       position = "right")+
    scale_x_discrete(position = "bottom")+
    annotate(geom = "text",
             label=paste0("Flat Line Test\nP.value=",
                          signif(unique(TXI.result.sub$p.value),2)),
             x=nrow(TXI.result.sub)-min(module.list[[2]])+1, y=max(scale.bar)-0.001,
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
  ggsave(paste0("Coprinopsis_cinerea.TDI.Flatline.",ref,".",refspecies,".png"),
         width = 4, height = 3, units = "in", dpi = 300)
  f=paste0("Coprinopsis_cinerea.TDI.Flatline.",ref,".",refspecies,".pptx")
  topptx(TDI, f, width = 4, height = 3, units = "in")
  rm(TXI.result.sub)

  ReductiveHourglassTest.result<-ReductiveHourglassTest(TMM.matrix.new[,c("NS","Genes",treat.order)],
                                                        modules = module.list,
                                                        permutations = 1000,
                                                        lillie.test = FALSE,
                                                        plotHistogram = FALSE,
                                                        runs = 100,
                                                        parallel = 14,
                                                        gof.warning = FALSE,
                                                        custom.perm.matrix = NULL)
  TXI.result.sub<-TXI.result[row.names(TXI.result) %in% treat.order,]
  TXI.result.sub$Stage<-row.names(TXI.result.sub)
  TXI.result.sub$group<-ref
  TXI.result.sub$statistics<-"ReductiveHourglassTest"
  TXI.result.sub$p.value<-ReductiveHourglassTest.result$p.value
  TXI.result.sub$std.dev<-ReductiveHourglassTest.result$std.dev
  TXI.result.list<-rbind(TXI.result.list, TXI.result.sub)

  scale.bar<-c(floor(min(TXI.result.sub$TXI.result+TXI.result.sub$std.dev)/0.001-0.5)*0.001,
               ceiling(max(TXI.result.sub$TXI.result+TXI.result.sub$std.dev)/0.001+0.5)*0.001)

  TDI<-TXI.result.sub %>%
    mutate(Stage = factor(Stage, levels = plot.order)) %>%
    ggplot(aes(x = Stage, y = TXI.result, group = 1))+
    geom_rect(aes(xmin=nrow(TXI.result.sub)-min(module.list[[2]]-0.5)+1,
                  xmax=nrow(TXI.result.sub)-max(module.list[[2]]-0.5),
                  ymin=-Inf, ymax=Inf),
              fill='dodgerblue',alpha = 0.05)+
    geom_ribbon(aes(ymin=TXI.result-std.dev, ymax=TXI.result+std.dev), fill = "gray50", alpha = 0.5)+
    geom_line(aes(linetype="solid"))+
    scale_y_continuous(limits = scale.bar,
                       breaks = seq(min(scale.bar),max(scale.bar),0.001),
                       #trans = "reverse",
                       position = "right")+
    scale_x_discrete(position = "bottom")+
    annotate(geom = "text",
             label=paste0("Reductive Hourglass Test (High-Low-High)\nP.value=",
                          signif(unique(TXI.result.sub$p.value),2)),
             x=nrow(TXI.result.sub)-min(module.list[[2]])+1, y=max(scale.bar)-0.001,
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
  ggsave(paste0("Coprinopsis_cinerea.TDI.reductive.",ref,".",refspecies,".png"),
         width = 4, height = 2.5, units = "in", dpi = 300)
  f=paste0("Coprinopsis_cinerea.TDI.reductive.",ref,".",refspecies,".pptx")
  topptx(TDI, f, width = 4, height = 2.5, units = "in")
  rm(TXI.result.sub)

  write.table(TXI.result.list,
              paste0("Coprinopsis_cinerea.NS.",refspecies,".txt"),
              quote = F, row.names = F, sep = "\t")
}

