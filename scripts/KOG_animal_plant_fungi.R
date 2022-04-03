library(dplyr)
library(tidyr)
library(myTAI)
library(ggplot2)
library(pheatmap)
setwd("hourglass/enrich/")

##Drosophila
################################################################
## KOG ~ all stages 
################################################################
##Read in gene expression level
species="Drosophila_melanogaster"
IDmatch<-read.delim(paste0("/home/yichun/RNAmodification/hourglass/eggnog/",
                           species,".GenematchID"), header = F)
names(IDmatch)<-c("Protein","Gene")

TMM.matrix<-read.delim(paste0("expression.dmel.txt"),header = T)
names(TMM.matrix)[2]<-"Genes"
TMM.matrix$Genes<-gsub("fbpp","FBpp",TMM.matrix$Genes)
names(TMM.matrix)<-gsub("X","",names(TMM.matrix))
names(TMM.matrix)<-gsub("\\.embryo","",names(TMM.matrix))
treat.order<-names(TMM.matrix)[3:ncol(TMM.matrix)]
TMM.matrix<-TMM.matrix[,c("Genes", treat.order)]

module.list<-list(early = 1:3, mid = 4:5, late = 6:12)

##RE by kogclass
kog2name<-read.delim("/home/yichun/RNAmodification/hourglass/eggnog/kog2name.txt", header = TRUE)
kog2name<-kog2name[order(kog2name$ONTOLOGY,kog2name$kogClass),]
kog2name$KS<-1:nrow(kog2name)
kog2name$ONTOLOGY<-gsub("INFORMATION STORAGE AND PROCESSING", "Information storage and processing", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("CELLULAR PROCESSES AND SIGNALING", "Cellular processes and signaling", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("POORLY CHARACTERIZED", "Poor", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("METABOLISM", "Metabolism", kog2name$ONTOLOGY)

GenesKOGpair.1v1<-read.delim(paste0("/home/yichun/RNAmodification/hourglass/eggnog/",
                                    species,".proteins.fa.KOG.1v1.txt"), header = TRUE)
GenesKOGpair.1v1<-unique(GenesKOGpair.1v1[,c("KOG","Genes")])
names(GenesKOGpair.1v1)[1]<-"kogClass"
GenesKOGpair.1v1<-merge(GenesKOGpair.1v1, kog2name[,c("KS","kogClass")], by = "kogClass", all.x = T)

TMM.matrix.new<-merge(unique(GenesKOGpair.1v1[,c("KS","Genes")]), TMM.matrix, by = "Genes", all.y = T)
TMM.matrix.new<-TMM.matrix.new[is.na(TMM.matrix.new$KS)==F, c("KS","Genes",treat.order)]

##Relative expression of all stages by PS
TMM.matrix.new[TMM.matrix.new<1]<-NA
Relative.exp<-as.data.frame(age.apply(ExpressionSet = TMM.matrix.new, RE))

Relative.exp$KS<-row.names(Relative.exp)
Relative.exp<-merge(Relative.exp, kog2name, by = "KS", all.x = T)
Relative.exp<-Relative.exp[order(as.numeric(Relative.exp$KS)),]
Relative.exp$Description<-paste0(" [",Relative.exp$kogClass,"] ", Relative.exp$kogName)

write.table(Relative.exp, paste0(species,".KOG.RE.txt"),
            quote = F, row.names = F, sep = "\t")
heatmap.input<-Relative.exp[,treat.order]
row.names(heatmap.input)<-Relative.exp$Description

row.anno<-data.frame(ONTOLOGY = factor(Relative.exp[,c("ONTOLOGY")]))
row.names(row.anno)<-Relative.exp$Description

col.anno<-data.frame(ONTOGENY = factor(c(rep("Early", length(module.list[[1]])),
                                         rep("Mid", length(module.list[[2]])),
                                         rep("Late", length(module.list[[3]])))))
row.names(col.anno)<-names(heatmap.input)

pheatmap(heatmap.input,
         legend = TRUE,
         scale = "none",
         #scale = "row",
         cluster_rows = F, cluster_cols = F,
         gaps_row = c(nrow(Relative.exp[Relative.exp$ONTOLOGY %in% c("Cellular processes and signaling"),]),
                      nrow(Relative.exp[Relative.exp$ONTOLOGY %in% c("Cellular processes and signaling",
                                                                     "Information storage and processing"),]),
                      nrow(Relative.exp)-1),
         annotation_row = row.anno,
         annotation_col = col.anno,
         annotation_colors = list(ONTOLOGY = c(`Cellular processes and signaling` = "#fdbf6f",
                                               `Information storage and processing` = "#a6cee3", 
                                               `Metabolism` = "#b2df8a",
                                               `Poor` = "gray75"),
                                  ONTOGENY = c(Early = "#33a02c",
                                               Mid = "#ff7f00",
                                               Late = "#6a3d9a")),
         cellwidth = 30, cellheight = 30, fontsize = 24,
         angle_col = c("90"),
         width = 24, height = 16,
         show_colnames = TRUE,
         show_rownames = T,
         filename = paste0(species,".KOG.RE.heatmap.png"))

##RE by ONTOLOGY
kog2name<-read.delim("/home/yichun/RNAmodification/hourglass/eggnog/kog2name.txt", header = TRUE)
kog2name<-kog2name[order(kog2name$ONTOLOGY,kog2name$kogClass),]
kog2name$KS[kog2name$ONTOLOGY=="CELLULAR PROCESSES AND SIGNALING"]<- 1
kog2name$KS[kog2name$ONTOLOGY=="INFORMATION STORAGE AND PROCESSING"]<- 2
kog2name$KS[kog2name$ONTOLOGY=="METABOLISM"]<- 3
kog2name$KS[kog2name$ONTOLOGY=="POORLY CHARACTERIZED"]<- 4

kog2name$ONTOLOGY<-gsub("INFORMATION STORAGE AND PROCESSING", "Information storage and processing", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("CELLULAR PROCESSES AND SIGNALING", "Cellular processes and signaling", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("POORLY CHARACTERIZED", "Poor", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("METABOLISM", "Metabolism", kog2name$ONTOLOGY)

GenesKOGpair.1v1<-read.delim(paste0("/home/yichun/RNAmodification/hourglass/eggnog/",
                                    species,".proteins.fa.KOG.1v1.txt"), header = TRUE)
GenesKOGpair.1v1<-unique(GenesKOGpair.1v1[,c("KOG","Genes")])
names(GenesKOGpair.1v1)[1]<-"kogClass"
GenesKOGpair.1v1<-merge(GenesKOGpair.1v1, kog2name[,c("KS","kogClass")], by = "kogClass", all.x = T)

TMM.matrix.new<-merge(unique(GenesKOGpair.1v1[,c("KS","Genes")]), TMM.matrix, by = "Genes", all.y = T)
TMM.matrix.new<-TMM.matrix.new[is.na(TMM.matrix.new$KS)==F, c("KS","Genes",treat.order)]
##Relative expression of all stages by PS
TMM.matrix.new[TMM.matrix.new<1]<-NA
Relative.exp<-as.data.frame(age.apply(ExpressionSet = TMM.matrix.new, RE))
Relative.exp$KS<-row.names(Relative.exp)
Relative.exp<-merge(Relative.exp, unique(kog2name[,c("KS","ONTOLOGY")]), by = "KS", all.x = T)
Relative.exp<-Relative.exp[order(as.numeric(Relative.exp$KS)),]

write.table(Relative.exp, paste0(species,".KOG.RE.ONTOLOGY.txt"),
            quote = F, row.names = F, sep = "\t")
heatmap.input<-Relative.exp[,treat.order]
row.names(heatmap.input)<-Relative.exp$ONTOLOGY

col.anno<-data.frame(ONTOGENY = factor(c(rep("Early", length(module.list[[1]])),
                                         rep("Mid", length(module.list[[2]])),
                                         rep("Late", length(module.list[[3]])))))
row.names(col.anno)<-names(heatmap.input)

pheatmap(heatmap.input,
         legend = TRUE,
         scale = "none",
         #scale = "row",
         cluster_rows = F, cluster_cols = F,
         gaps_row = c(nrow(Relative.exp[Relative.exp$ONTOLOGY %in% c("Cellular processes and signaling"),]),
                      nrow(Relative.exp[Relative.exp$ONTOLOGY %in% c("Cellular processes and signaling",
                                                                     "Information storage and processing"),]),
                      nrow(Relative.exp)-1),
         annotation_col = col.anno,
         annotation_colors = list(ONTOGENY = c(Early = "#33a02c",
                                               Mid = "#ff7f00",
                                               Late = "#6a3d9a")),
         cellwidth = 30, cellheight = 30, fontsize = 24,
         angle_col = c("90"),
         width = 16, height = 6,
         show_colnames = TRUE,
         show_rownames = T,
         filename = paste0(species,".KOG.RE.ONTOLOGY.heatmap.png"))

line.input<-as.data.frame(matrix(NA, nrow = 0,ncol = 3))
names(line.input)<-c("Stage","RE","ONTOLOGY")
for (i in 2:(1+length(treat.order))) {
  line.input.sub<-Relative.exp[,c(i,ncol(Relative.exp))]
  line.input.sub$Stage<-names(Relative.exp[i])
  names(line.input.sub)[1]<-"RE"
  line.input<-rbind(line.input,line.input.sub)
}
p<-line.input %>%
  mutate(Stage = factor(Stage, levels = treat.order)) %>%
  ggplot(aes(x = Stage, y = RE, group = ONTOLOGY))+
  geom_vline(xintercept = min(module.list[[2]])-0.5, color = "grey75", linetype = "dashed")+
  geom_vline(xintercept = max(module.list[[2]])+0.5, color = "grey75", linetype = "dashed")+
  geom_line(aes(color = ONTOLOGY))+
  scale_color_manual(values =  c("#fdbf6f","#a6cee3", "#b2df8a","gray75"))+
  scale_y_continuous(limits = c(0,1.1),
                     breaks = seq(0,1,0.2),
                     position = "left")+
  scale_x_discrete(position = "bottom")+
  annotate(geom = "text", 
           label="Early",
           x=mean(module.list[[1]]), y=1.1,
           colour="black", size=2, fontface = "italic")+
  annotate(geom = "text", 
           label="Mid",
           x=mean(module.list[[2]]), y=1.1,
           colour="black", size=2, fontface = "italic")+
  annotate(geom = "text", 
           label="Late",
           x=mean(module.list[[3]]), y=1.1,
           colour="black", size=2, fontface = "italic")+
  labs(title = "", x = "", y = "Relative expression", colour = NULL)+
  theme(axis.line = element_line(linetype = "solid"),
        axis.ticks.y = element_line(colour = "black", size = 0.5),
        axis.ticks.x = element_line(colour = "black", size = 0.5),
        axis.text.x = element_text(size = 8, colour = "black", angle =45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 8, colour = "black"),
        plot.title = element_text(size = 8, hjust = 0.5, face = "plain"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 0),
        panel.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        legend.position = "right",
        strip.text = element_text(size = 9))
p
ggsave(paste0(species,".KOG.RE.ONTOLOGY.line.png"), 
       width = 6, height = 2.5, units = "in", dpi = 300)
f=paste0(species,".KOG.RE.ONTOLOGY.line.pptx")
topptx(p, f, width = 6, height = 2.5, units = "in")

################################################################
## KOG ~ Early-Mid-Late 
################################################################
TMM.matrix.EML<-data.frame(Genes = TMM.matrix$Genes)
if (length(module.list[[1]])<2) {
  TMM.matrix.EML$Early<-TMM.matrix[,1+module.list[[1]]]
} else {
  TMM.matrix.EML$Early<-apply(TMM.matrix[,1+module.list[[1]]], 1, mean)
}
TMM.matrix.EML$Mid<-apply(TMM.matrix[,1+module.list[[2]]], 1, mean)
TMM.matrix.EML$Late<-apply(TMM.matrix[,1+module.list[[3]]], 1, mean)
treat.order<-c("Early","Mid","Late")

##RE by kogclass
kog2name<-read.delim("/home/yichun/RNAmodification/hourglass/eggnog/kog2name.txt", header = TRUE)
kog2name<-kog2name[order(kog2name$ONTOLOGY,kog2name$kogClass),]
kog2name$KS<-1:nrow(kog2name)
kog2name$ONTOLOGY<-gsub("INFORMATION STORAGE AND PROCESSING", "Information storage and processing", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("CELLULAR PROCESSES AND SIGNALING", "Cellular processes and signaling", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("POORLY CHARACTERIZED", "Poor", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("METABOLISM", "Metabolism", kog2name$ONTOLOGY)

GenesKOGpair.1v1<-read.delim(paste0("/home/yichun/RNAmodification/hourglass/eggnog/",
                                    species,".proteins.fa.KOG.1v1.txt"), header = TRUE)
GenesKOGpair.1v1<-unique(GenesKOGpair.1v1[,c("KOG","Genes")])
names(GenesKOGpair.1v1)[1]<-"kogClass"
GenesKOGpair.1v1<-merge(GenesKOGpair.1v1, kog2name[,c("KS","kogClass")], by = "kogClass", all.x = T)

TMM.matrix.new<-merge(unique(GenesKOGpair.1v1[,c("KS","Genes")]), TMM.matrix.EML, by = "Genes", all.y = T)
TMM.matrix.new<-TMM.matrix.new[is.na(TMM.matrix.new$KS)==F, c("KS","Genes",treat.order)]
##Relative expression of all stages by PS
TMM.matrix.new[TMM.matrix.new<1]<-NA
Relative.exp<-as.data.frame(age.apply(ExpressionSet = TMM.matrix.new, RE))

Relative.exp$KS<-row.names(Relative.exp)
Relative.exp<-merge(Relative.exp, kog2name, by = "KS", all.x = T)
Relative.exp<-Relative.exp[order(as.numeric(Relative.exp$KS)),]
Relative.exp$Description<-paste0(" [",Relative.exp$kogClass,"] ", Relative.exp$kogName)

write.table(Relative.exp, paste0(species,".EML.KOG.RE.txt"),
            quote = F, row.names = F, sep = "\t")
heatmap.input<-Relative.exp[,treat.order]
row.names(heatmap.input)<-Relative.exp$Description

row.anno<-data.frame(ONTOLOGY = factor(Relative.exp[,c("ONTOLOGY")]))
row.names(row.anno)<-Relative.exp$Description

pheatmap(heatmap.input,
         legend = TRUE,
         scale = "none",
         #scale = "row",
         cluster_rows = F, cluster_cols = F,
         gaps_row = c(nrow(Relative.exp[Relative.exp$ONTOLOGY %in% c("Cellular processes and signaling"),]),
                      nrow(Relative.exp[Relative.exp$ONTOLOGY %in% c("Cellular processes and signaling",
                                                                     "Information storage and processing"),]),
                      nrow(Relative.exp)-1),
         annotation_row = row.anno,
         annotation_colors = list(ONTOLOGY = c(`Cellular processes and signaling` = "#fdbf6f",
                                               `Information storage and processing` = "#a6cee3", 
                                               `Metabolism` = "#b2df8a",
                                               `Poor` = "gray75"),
                                  ONTOGENY = c(Early = "#33a02c",
                                               Mid = "#ff7f00",
                                               Late = "#6a3d9a")),
         cellwidth = 30, cellheight = 30, fontsize = 24,
         angle_col = c("90"),
         width = 24, height = 16,
         show_colnames = TRUE,
         show_rownames = T,
         filename = paste0(species,".EML.KOG.RE.heatmap.png"))

##RE by ONTOLOGY
kog2name<-read.delim("/home/yichun/RNAmodification/hourglass/eggnog/kog2name.txt", header = TRUE)
kog2name<-kog2name[order(kog2name$ONTOLOGY,kog2name$kogClass),]
kog2name$KS[kog2name$ONTOLOGY=="CELLULAR PROCESSES AND SIGNALING"]<- 1
kog2name$KS[kog2name$ONTOLOGY=="INFORMATION STORAGE AND PROCESSING"]<- 2
kog2name$KS[kog2name$ONTOLOGY=="METABOLISM"]<- 3
kog2name$KS[kog2name$ONTOLOGY=="POORLY CHARACTERIZED"]<- 4

kog2name$ONTOLOGY<-gsub("INFORMATION STORAGE AND PROCESSING", "Information storage and processing", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("CELLULAR PROCESSES AND SIGNALING", "Cellular processes and signaling", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("POORLY CHARACTERIZED", "Poor", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("METABOLISM", "Metabolism", kog2name$ONTOLOGY)

GenesKOGpair.1v1<-read.delim(paste0("/home/yichun/RNAmodification/hourglass/eggnog/",
                                    species,".proteins.fa.KOG.1v1.txt"), header = TRUE)
GenesKOGpair.1v1<-unique(GenesKOGpair.1v1[,c("KOG","Genes")])
names(GenesKOGpair.1v1)[1]<-"kogClass"
GenesKOGpair.1v1<-merge(GenesKOGpair.1v1, kog2name[,c("KS","kogClass")], by = "kogClass", all.x = T)

TMM.matrix.new<-merge(GenesKOGpair.1v1[,c("KS","Genes")], TMM.matrix.EML, by = "Genes", all.y = T)
TMM.matrix.new<-TMM.matrix.new[is.na(TMM.matrix.new$KS)==F, c("KS","Genes",treat.order)]
##Relative expression of all stages by PS
TMM.matrix.new[TMM.matrix.new<1]<-NA
Relative.exp<-as.data.frame(age.apply(ExpressionSet = TMM.matrix.new, RE))

Relative.exp$KS<-row.names(Relative.exp)
Relative.exp<-merge(Relative.exp, unique(kog2name[,c("KS","ONTOLOGY")]), by = "KS", all.x = T)
Relative.exp<-Relative.exp[order(as.numeric(Relative.exp$KS)),]

write.table(Relative.exp, paste0(species,".EML.KOG.RE.ONTOLOGY.txt"),
            quote = F, row.names = F, sep = "\t")
heatmap.input<-Relative.exp[,treat.order]
row.names(heatmap.input)<-Relative.exp$ONTOLOGY

pheatmap(heatmap.input,
         legend = TRUE,
         scale = "none",
         #scale = "row",
         cluster_rows = F, cluster_cols = F,
         gaps_row = c(nrow(Relative.exp[Relative.exp$ONTOLOGY %in% c("Cellular processes and signaling"),]),
                      nrow(Relative.exp[Relative.exp$ONTOLOGY %in% c("Cellular processes and signaling",
                                                                     "Information storage and processing"),]),
                      nrow(Relative.exp)-1),
         cellwidth = 30, cellheight = 30, fontsize = 24,
         angle_col = c("90"),
         width = 8, height = 6,
         show_colnames = TRUE,
         show_rownames = T,
         filename = paste0(species,".EML.KOG.RE.ONTOLOGY.heatmap.png"))

line.input<-as.data.frame(matrix(NA, nrow = 0,ncol = 3))
names(line.input)<-c("Stage","RE","ONTOLOGY")
for (i in 2:(1+length(treat.order))) {
  line.input.sub<-Relative.exp[,c(i,ncol(Relative.exp))]
  line.input.sub$Stage<-names(Relative.exp[i])
  names(line.input.sub)[1]<-"RE"
  line.input<-rbind(line.input,line.input.sub)
}
p<-line.input %>%
  mutate(Stage = factor(Stage, levels = treat.order)) %>%
  ggplot(aes(x = Stage, y = RE, group = ONTOLOGY))+
  geom_vline(xintercept = 1.5, color = "grey75", linetype = "dashed")+
  geom_vline(xintercept = 2.5, color = "grey75", linetype = "dashed")+
  geom_line(aes(color = ONTOLOGY))+
  scale_color_manual(values =  c("#fdbf6f","#a6cee3", "#b2df8a","gray75"))+
  scale_y_continuous(limits = c(0,1.1),
                     breaks = seq(0,1,0.2),
                     position = "left")+
  scale_x_discrete(position = "bottom")+
  annotate(geom = "text", 
           label="Early",
           x=1, y=1.1,
           colour="black", size=2, fontface = "italic")+
  annotate(geom = "text", 
           label="Mid",
           x=2, y=1.1,
           colour="black", size=2, fontface = "italic")+
  annotate(geom = "text", 
           label="Late",
           x=3, y=1.1,
           colour="black", size=2, fontface = "italic")+
  labs(title = "", x = "", y = "Relative expression", colour = NULL)+
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
        legend.position = "right",
        strip.text = element_text(size = 9))
p
ggsave(paste0(species,".EML.KOG.RE.ONTOLOGY.line.png"), 
       width = 5, height = 2.5, units = "in", dpi = 300)
f=paste0(species,".EML.KOG.RE.ONTOLOGY.line.pptx")
topptx(p, f, width = 5, height = 2.5, units = "in")


##Danio
################################################################
## KOG ~ all stages 
################################################################
##Read in gene expression level
species="Danio_rerio"
IDmatch<-read.delim(paste0("/home/yichun/RNAmodification/hourglass/eggnog/",
                           species,".GenematchID"), header = F)
names(IDmatch)<-c("Protein","Gene")

TMM.matrix<-read.delim(paste0("expression.drer.txt"),header = T)
names(TMM.matrix)[2]<-"Genes"
names(TMM.matrix)<-gsub("X","",names(TMM.matrix))
treat.order<-names(TMM.matrix)[3:ncol(TMM.matrix)]
TMM.matrix<-TMM.matrix[,c("Genes", treat.order)]

module.list<-list(early = 1:18, mid = 19:36, late = 37:40)

##RE by kogclass
kog2name<-read.delim("/home/yichun/RNAmodification/hourglass/eggnog/kog2name.txt", header = TRUE)
kog2name<-kog2name[order(kog2name$ONTOLOGY,kog2name$kogClass),]
kog2name$KS<-1:nrow(kog2name)
kog2name$ONTOLOGY<-gsub("INFORMATION STORAGE AND PROCESSING", "Information storage and processing", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("CELLULAR PROCESSES AND SIGNALING", "Cellular processes and signaling", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("POORLY CHARACTERIZED", "Poor", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("METABOLISM", "Metabolism", kog2name$ONTOLOGY)

IDmatch<-read.delim(paste0("/home/yichun/RNAmodification/hourglass/eggnog/",
                           species,".GenematchID"), header = F)
names(IDmatch)<-c("Protein","Genes")
GenesKOGpair.1v1<-read.delim(paste0("/home/yichun/RNAmodification/hourglass/eggnog/",
                                    species,".proteins.fa.KOG.1v1.txt"), header = TRUE)
names(GenesKOGpair.1v1)<-c("KOG","Protein")
GenesKOGpair.1v1<-merge(GenesKOGpair.1v1, IDmatch, by = "Protein", all.x = T)
GenesKOGpair.1v1<-unique(GenesKOGpair.1v1[,c("KOG","Genes")])
names(GenesKOGpair.1v1)[1]<-"kogClass"
GenesKOGpair.1v1<-merge(GenesKOGpair.1v1, kog2name[,c("KS","kogClass")], by = "kogClass", all.x = T)

TMM.matrix.new<-merge(unique(GenesKOGpair.1v1[,c("KS","Genes")]), TMM.matrix, by = "Genes", all.y = T)
TMM.matrix.new<-TMM.matrix.new[is.na(TMM.matrix.new$KS)==F, c("KS","Genes",treat.order)]

##Relative expression of all stages by PS
TMM.matrix.new[TMM.matrix.new<1]<-NA
Relative.exp<-as.data.frame(age.apply(ExpressionSet = TMM.matrix.new, RE))

Relative.exp$KS<-row.names(Relative.exp)
Relative.exp<-merge(Relative.exp, kog2name, by = "KS", all.x = T)
Relative.exp<-Relative.exp[order(as.numeric(Relative.exp$KS)),]
Relative.exp$Description<-paste0(" [",Relative.exp$kogClass,"] ", Relative.exp$kogName)

write.table(Relative.exp, paste0(species,".KOG.RE.txt"),
            quote = F, row.names = F, sep = "\t")
heatmap.input<-Relative.exp[,treat.order]
row.names(heatmap.input)<-Relative.exp$Description

row.anno<-data.frame(ONTOLOGY = factor(Relative.exp[,c("ONTOLOGY")]))
row.names(row.anno)<-Relative.exp$Description

col.anno<-data.frame(ONTOGENY = factor(c(rep("Early", length(module.list[[1]])),
                                         rep("Mid", length(module.list[[2]])),
                                         rep("Late", length(module.list[[3]])))))
row.names(col.anno)<-names(heatmap.input)

pheatmap(heatmap.input,
         legend = TRUE,
         scale = "none",
         #scale = "row",
         cluster_rows = F, cluster_cols = F,
         gaps_row = c(nrow(Relative.exp[Relative.exp$ONTOLOGY %in% c("Cellular processes and signaling"),]),
                      nrow(Relative.exp[Relative.exp$ONTOLOGY %in% c("Cellular processes and signaling",
                                                                     "Information storage and processing"),]),
                      nrow(Relative.exp)-1),
         annotation_row = row.anno,
         annotation_col = col.anno,
         annotation_colors = list(ONTOLOGY = c(`Cellular processes and signaling` = "#fdbf6f",
                                               `Information storage and processing` = "#a6cee3", 
                                               `Metabolism` = "#b2df8a",
                                               `Poor` = "gray75"),
                                  ONTOGENY = c(Early = "#33a02c",
                                               Mid = "#ff7f00",
                                               Late = "#6a3d9a")),
         cellwidth = 30, cellheight = 30, fontsize = 24,
         angle_col = c("90"),
         width = 36, height = 16,
         show_colnames = TRUE,
         show_rownames = T,
         filename = paste0(species,".KOG.RE.heatmap.png"))

##RE by ONTOLOGY
kog2name<-read.delim("/home/yichun/RNAmodification/hourglass/eggnog/kog2name.txt", header = TRUE)
kog2name<-kog2name[order(kog2name$ONTOLOGY,kog2name$kogClass),]
kog2name$KS[kog2name$ONTOLOGY=="CELLULAR PROCESSES AND SIGNALING"]<- 1
kog2name$KS[kog2name$ONTOLOGY=="INFORMATION STORAGE AND PROCESSING"]<- 2
kog2name$KS[kog2name$ONTOLOGY=="METABOLISM"]<- 3
kog2name$KS[kog2name$ONTOLOGY=="POORLY CHARACTERIZED"]<- 4

kog2name$ONTOLOGY<-gsub("INFORMATION STORAGE AND PROCESSING", "Information storage and processing", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("CELLULAR PROCESSES AND SIGNALING", "Cellular processes and signaling", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("POORLY CHARACTERIZED", "Poor", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("METABOLISM", "Metabolism", kog2name$ONTOLOGY)

IDmatch<-read.delim(paste0("/home/yichun/RNAmodification/hourglass/eggnog/",
                           species,".GenematchID"), header = F)
names(IDmatch)<-c("Protein","Genes")
GenesKOGpair.1v1<-read.delim(paste0("/home/yichun/RNAmodification/hourglass/eggnog/",
                                    species,".proteins.fa.KOG.1v1.txt"), header = TRUE)
names(GenesKOGpair.1v1)<-c("KOG","Protein")
GenesKOGpair.1v1<-merge(GenesKOGpair.1v1, IDmatch, by = "Protein", all.x = T)
GenesKOGpair.1v1<-unique(GenesKOGpair.1v1[,c("KOG","Genes")])
names(GenesKOGpair.1v1)[1]<-"kogClass"
GenesKOGpair.1v1<-merge(GenesKOGpair.1v1, kog2name[,c("KS","kogClass")], by = "kogClass", all.x = T)

TMM.matrix.new<-merge(unique(GenesKOGpair.1v1[,c("KS","Genes")]), TMM.matrix, by = "Genes", all.y = T)
TMM.matrix.new<-TMM.matrix.new[is.na(TMM.matrix.new$KS)==F, c("KS","Genes",treat.order)]
##Relative expression of all stages by PS
TMM.matrix.new[TMM.matrix.new<1]<-NA
Relative.exp<-as.data.frame(age.apply(ExpressionSet = TMM.matrix.new, RE))

Relative.exp$KS<-row.names(Relative.exp)
Relative.exp<-merge(Relative.exp, unique(kog2name[,c("KS","ONTOLOGY")]), by = "KS", all.x = T)
Relative.exp<-Relative.exp[order(as.numeric(Relative.exp$KS)),]

write.table(Relative.exp, paste0(species,".KOG.RE.ONTOLOGY.txt"),
            quote = F, row.names = F, sep = "\t")
heatmap.input<-Relative.exp[,treat.order]
row.names(heatmap.input)<-Relative.exp$ONTOLOGY

col.anno<-data.frame(ONTOGENY = factor(c(rep("Early", length(module.list[[1]])),
                                         rep("Mid", length(module.list[[2]])),
                                         rep("Late", length(module.list[[3]])))))
row.names(col.anno)<-names(heatmap.input)

pheatmap(heatmap.input,
         legend = TRUE,
         scale = "none",
         #scale = "row",
         cluster_rows = F, cluster_cols = F,
         gaps_row = c(nrow(Relative.exp[Relative.exp$ONTOLOGY %in% c("Cellular processes and signaling"),]),
                      nrow(Relative.exp[Relative.exp$ONTOLOGY %in% c("Cellular processes and signaling",
                                                                     "Information storage and processing"),]),
                      nrow(Relative.exp)-1),
         annotation_col = col.anno,
         annotation_colors = list(ONTOGENY = c(Early = "#33a02c",
                                               Mid = "#ff7f00",
                                               Late = "#6a3d9a")),
         cellwidth = 30, cellheight = 30, fontsize = 24,
         angle_col = c("90"),
         width = 36, height = 6,
         show_colnames = TRUE,
         show_rownames = T,
         filename = paste0(species,".KOG.RE.ONTOLOGY.heatmap.png"))

line.input<-as.data.frame(matrix(NA, nrow = 0,ncol = 3))
names(line.input)<-c("Stage","RE","ONTOLOGY")
for (i in 2:(1+length(treat.order))) {
  line.input.sub<-Relative.exp[,c(i,ncol(Relative.exp))]
  line.input.sub$Stage<-names(Relative.exp[i])
  names(line.input.sub)[1]<-"RE"
  line.input<-rbind(line.input,line.input.sub)
}
p<-line.input %>%
  mutate(Stage = factor(Stage, levels = treat.order)) %>%
  ggplot(aes(x = Stage, y = RE, group = ONTOLOGY))+
  geom_vline(xintercept = min(module.list[[2]])-0.5, color = "grey75", linetype = "dashed")+
  geom_vline(xintercept = max(module.list[[2]])+0.5, color = "grey75", linetype = "dashed")+
  geom_line(aes(color = ONTOLOGY))+
  scale_color_manual(values =  c("#fdbf6f","#a6cee3", "#b2df8a","gray75"))+
  scale_y_continuous(limits = c(0,1.1),
                     breaks = seq(0,1,0.2),
                     position = "left")+
  scale_x_discrete(position = "bottom")+
  annotate(geom = "text", 
           label="Early",
           x=mean(module.list[[1]]), y=1.1,
           colour="black", size=2, fontface = "italic")+
  annotate(geom = "text", 
           label="Mid",
           x=mean(module.list[[2]]), y=1.1,
           colour="black", size=2, fontface = "italic")+
  annotate(geom = "text", 
           label="Late",
           x=mean(module.list[[3]]), y=1.1,
           colour="black", size=2, fontface = "italic")+
  labs(title = "", x = "", y = "Relative expression", colour = NULL)+
  theme(axis.line = element_line(linetype = "solid"),
        axis.ticks.y = element_line(colour = "black", size = 0.5),
        axis.ticks.x = element_line(colour = "black", size = 0.5),
        axis.text.x = element_text(size = 8, colour = "black", angle =45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 8, colour = "black"),
        plot.title = element_text(size = 8, hjust = 0.5, face = "plain"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 0),
        panel.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        legend.position = "right",
        strip.text = element_text(size = 9))
p
ggsave(paste0(species,".KOG.RE.ONTOLOGY.line.png"), 
       width = 8, height = 2.5, units = "in", dpi = 300)
f=paste0(species,".KOG.RE.ONTOLOGY.line.pptx")
topptx(p, f, width = 8, height = 2.5, units = "in")

################################################################
## KOG ~ Early-Mid-Late 
################################################################
TMM.matrix.EML<-data.frame(Genes = TMM.matrix$Genes)
if (length(module.list[[1]])<2) {
  TMM.matrix.EML$Early<-TMM.matrix[,1+module.list[[1]]]
} else {
  TMM.matrix.EML$Early<-apply(TMM.matrix[,1+module.list[[1]]], 1, mean)
}
TMM.matrix.EML$Mid<-apply(TMM.matrix[,1+module.list[[2]]], 1, mean)
TMM.matrix.EML$Late<-apply(TMM.matrix[,1+module.list[[3]]], 1, mean)
treat.order<-c("Early","Mid","Late")

##RE by kogclass
kog2name<-read.delim("/home/yichun/RNAmodification/hourglass/eggnog/kog2name.txt", header = TRUE)
kog2name<-kog2name[order(kog2name$ONTOLOGY,kog2name$kogClass),]
kog2name$KS<-1:nrow(kog2name)
kog2name$ONTOLOGY<-gsub("INFORMATION STORAGE AND PROCESSING", "Information storage and processing", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("CELLULAR PROCESSES AND SIGNALING", "Cellular processes and signaling", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("POORLY CHARACTERIZED", "Poor", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("METABOLISM", "Metabolism", kog2name$ONTOLOGY)

IDmatch<-read.delim(paste0("/home/yichun/RNAmodification/hourglass/eggnog/",
                           species,".GenematchID"), header = F)
names(IDmatch)<-c("Protein","Genes")
GenesKOGpair.1v1<-read.delim(paste0("/home/yichun/RNAmodification/hourglass/eggnog/",
                                    species,".proteins.fa.KOG.1v1.txt"), header = TRUE)
names(GenesKOGpair.1v1)<-c("KOG","Protein")
GenesKOGpair.1v1<-merge(GenesKOGpair.1v1, IDmatch, by = "Protein", all.x = T)
GenesKOGpair.1v1<-unique(GenesKOGpair.1v1[,c("KOG","Genes")])
names(GenesKOGpair.1v1)[1]<-"kogClass"
GenesKOGpair.1v1<-merge(GenesKOGpair.1v1, kog2name[,c("KS","kogClass")], by = "kogClass", all.x = T)

TMM.matrix.new<-merge(unique(GenesKOGpair.1v1[,c("KS","Genes")]), TMM.matrix.EML, by = "Genes", all.y = T)
TMM.matrix.new<-TMM.matrix.new[is.na(TMM.matrix.new$KS)==F, c("KS","Genes",treat.order)]
##Relative expression of all stages by PS
TMM.matrix.new[TMM.matrix.new<1]<-NA
Relative.exp<-as.data.frame(age.apply(ExpressionSet = TMM.matrix.new, RE))

Relative.exp$KS<-row.names(Relative.exp)
Relative.exp<-merge(Relative.exp, kog2name, by = "KS", all.x = T)
Relative.exp<-Relative.exp[order(as.numeric(Relative.exp$KS)),]
Relative.exp$Description<-paste0(" [",Relative.exp$kogClass,"] ", Relative.exp$kogName)

write.table(Relative.exp, paste0(species,".EML.KOG.RE.txt"),
            quote = F, row.names = F, sep = "\t")
heatmap.input<-Relative.exp[,treat.order]
row.names(heatmap.input)<-Relative.exp$Description

row.anno<-data.frame(ONTOLOGY = factor(Relative.exp[,c("ONTOLOGY")]))
row.names(row.anno)<-Relative.exp$Description

pheatmap(heatmap.input,
         legend = TRUE,
         scale = "none",
         #scale = "row",
         cluster_rows = F, cluster_cols = F,
         gaps_row = c(nrow(Relative.exp[Relative.exp$ONTOLOGY %in% c("Cellular processes and signaling"),]),
                      nrow(Relative.exp[Relative.exp$ONTOLOGY %in% c("Cellular processes and signaling",
                                                                     "Information storage and processing"),]),
                      nrow(Relative.exp)-1),
         annotation_row = row.anno,
         annotation_colors = list(ONTOLOGY = c(`Cellular processes and signaling` = "#fdbf6f",
                                               `Information storage and processing` = "#a6cee3", 
                                               `Metabolism` = "#b2df8a",
                                               `Poor` = "gray75"),
                                  ONTOGENY = c(Early = "#33a02c",
                                               Mid = "#ff7f00",
                                               Late = "#6a3d9a")),
         cellwidth = 30, cellheight = 30, fontsize = 24,
         angle_col = c("90"),
         width = 24, height = 16,
         show_colnames = TRUE,
         show_rownames = T,
         filename = paste0(species,".EML.KOG.RE.heatmap.png"))

##RE by ONTOLOGY
kog2name<-read.delim("/home/yichun/RNAmodification/hourglass/eggnog/kog2name.txt", header = TRUE)
kog2name<-kog2name[order(kog2name$ONTOLOGY,kog2name$kogClass),]
kog2name$KS[kog2name$ONTOLOGY=="CELLULAR PROCESSES AND SIGNALING"]<- 1
kog2name$KS[kog2name$ONTOLOGY=="INFORMATION STORAGE AND PROCESSING"]<- 2
kog2name$KS[kog2name$ONTOLOGY=="METABOLISM"]<- 3
kog2name$KS[kog2name$ONTOLOGY=="POORLY CHARACTERIZED"]<- 4

kog2name$ONTOLOGY<-gsub("INFORMATION STORAGE AND PROCESSING", "Information storage and processing", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("CELLULAR PROCESSES AND SIGNALING", "Cellular processes and signaling", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("POORLY CHARACTERIZED", "Poor", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("METABOLISM", "Metabolism", kog2name$ONTOLOGY)

IDmatch<-read.delim(paste0("/home/yichun/RNAmodification/hourglass/eggnog/",
                           species,".GenematchID"), header = F)
names(IDmatch)<-c("Protein","Genes")
GenesKOGpair.1v1<-read.delim(paste0("/home/yichun/RNAmodification/hourglass/eggnog/",
                                    species,".proteins.fa.KOG.1v1.txt"), header = TRUE)
names(GenesKOGpair.1v1)<-c("KOG","Protein")
GenesKOGpair.1v1<-merge(GenesKOGpair.1v1, IDmatch, by = "Protein", all.x = T)
GenesKOGpair.1v1<-unique(GenesKOGpair.1v1[,c("KOG","Genes")])
names(GenesKOGpair.1v1)[1]<-"kogClass"
GenesKOGpair.1v1<-merge(GenesKOGpair.1v1, kog2name[,c("KS","kogClass")], by = "kogClass", all.x = T)

TMM.matrix.new<-merge(GenesKOGpair.1v1[,c("KS","Genes")], TMM.matrix.EML, by = "Genes", all.y = T)
TMM.matrix.new<-TMM.matrix.new[is.na(TMM.matrix.new$KS)==F, c("KS","Genes",treat.order)]
##Relative expression of all stages by PS
TMM.matrix.new[TMM.matrix.new<1]<-NA
Relative.exp<-as.data.frame(age.apply(ExpressionSet = TMM.matrix.new, RE))

Relative.exp$KS<-row.names(Relative.exp)
Relative.exp<-merge(Relative.exp, unique(kog2name[,c("KS","ONTOLOGY")]), by = "KS", all.x = T)
Relative.exp<-Relative.exp[order(as.numeric(Relative.exp$KS)),]

write.table(Relative.exp, paste0(species,".EML.KOG.RE.ONTOLOGY.txt"),
            quote = F, row.names = F, sep = "\t")
heatmap.input<-Relative.exp[,treat.order]
row.names(heatmap.input)<-Relative.exp$ONTOLOGY

pheatmap(heatmap.input,
         legend = TRUE,
         scale = "none",
         #scale = "row",
         cluster_rows = F, cluster_cols = F,
         gaps_row = c(nrow(Relative.exp[Relative.exp$ONTOLOGY %in% c("Cellular processes and signaling"),]),
                      nrow(Relative.exp[Relative.exp$ONTOLOGY %in% c("Cellular processes and signaling",
                                                                     "Information storage and processing"),]),
                      nrow(Relative.exp)-1),
         cellwidth = 30, cellheight = 30, fontsize = 24,
         angle_col = c("90"),
         width = 8, height = 6,
         show_colnames = TRUE,
         show_rownames = T,
         filename = paste0(species,".EML.KOG.RE.ONTOLOGY.heatmap.png"))

line.input<-as.data.frame(matrix(NA, nrow = 0,ncol = 3))
names(line.input)<-c("Stage","RE","ONTOLOGY")
for (i in 2:(1+length(treat.order))) {
  line.input.sub<-Relative.exp[,c(i,ncol(Relative.exp))]
  line.input.sub$Stage<-names(Relative.exp[i])
  names(line.input.sub)[1]<-"RE"
  line.input<-rbind(line.input,line.input.sub)
}
p<-line.input %>%
  mutate(Stage = factor(Stage, levels = treat.order)) %>%
  ggplot(aes(x = Stage, y = RE, group = ONTOLOGY))+
  geom_vline(xintercept = 1.5, color = "grey75", linetype = "dashed")+
  geom_vline(xintercept = 2.5, color = "grey75", linetype = "dashed")+
  geom_line(aes(color = ONTOLOGY))+
  scale_color_manual(values =  c("#fdbf6f","#a6cee3", "#b2df8a","gray75"))+
  scale_y_continuous(limits = c(0,1.1),
                     breaks = seq(0,1,0.2),
                     position = "left")+
  scale_x_discrete(position = "bottom")+
  annotate(geom = "text", 
           label="Early",
           x=1, y=1.1,
           colour="black", size=2, fontface = "italic")+
  annotate(geom = "text", 
           label="Mid",
           x=2, y=1.1,
           colour="black", size=2, fontface = "italic")+
  annotate(geom = "text", 
           label="Late",
           x=3, y=1.1,
           colour="black", size=2, fontface = "italic")+
  labs(title = "", x = "", y = "Relative expression", colour = NULL)+
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
        legend.position = "right",
        strip.text = element_text(size = 9))
p
ggsave(paste0(species,".EML.KOG.RE.ONTOLOGY.line.png"), 
       width = 5, height = 2.5, units = "in", dpi = 300)
f=paste0(species,".EML.KOG.RE.ONTOLOGY.line.pptx")
topptx(p, f, width = 5, height = 2.5, units = "in")


##Arabidopsis
################################################################
## KOG ~ all stages 
################################################################
##Read in gene expression level
species="Arabidopsis_thaliana"

TMM.matrix<-read.delim(paste0("expression.atha.txt"),header = T)
names(TMM.matrix)[2]<-"Genes"
TMM.matrix$Genes<-gsub("at1g","AT1G",TMM.matrix$Genes)
TMM.matrix$Genes<-gsub("at2g","AT2G",TMM.matrix$Genes)
TMM.matrix$Genes<-gsub("at3g","AT3G",TMM.matrix$Genes)
TMM.matrix$Genes<-gsub("at4g","AT4G",TMM.matrix$Genes)
TMM.matrix$Genes<-gsub("at5g","AT5G",TMM.matrix$Genes)
TMM.matrix$Genes<-gsub("atcg","ATCG",TMM.matrix$Genes)
TMM.matrix$Genes<-gsub("atmg","ATMG",TMM.matrix$Genes)
treat.order<-names(TMM.matrix)[3:ncol(TMM.matrix)]
TMM.matrix<-TMM.matrix[,c("Genes", treat.order)]

module.list<-list(early = 1:2, mid = 3:5, late = 6:7)

##RE by kogclass
kog2name<-read.delim("/home/yichun/RNAmodification/hourglass/eggnog/kog2name.txt", header = TRUE)
kog2name<-kog2name[order(kog2name$ONTOLOGY,kog2name$kogClass),]
kog2name$KS<-1:nrow(kog2name)
kog2name$ONTOLOGY<-gsub("INFORMATION STORAGE AND PROCESSING", "Information storage and processing", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("CELLULAR PROCESSES AND SIGNALING", "Cellular processes and signaling", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("POORLY CHARACTERIZED", "Poor", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("METABOLISM", "Metabolism", kog2name$ONTOLOGY)

GenesKOGpair.1v1<-read.delim(paste0("/home/yichun/RNAmodification/hourglass/eggnog/",
                                    species,".protein.fa.KOG.1v1.txt"), header = TRUE)
names(GenesKOGpair.1v1)<-c("KOG","Protein")
GenesKOGpair.1v1<-separate(GenesKOGpair.1v1, Protein, c("Genes"), sep = "\\|")
GenesKOGpair.1v1<-unique(GenesKOGpair.1v1[,c("KOG","Genes")])
names(GenesKOGpair.1v1)[1]<-"kogClass"
GenesKOGpair.1v1<-merge(GenesKOGpair.1v1, kog2name[,c("KS","kogClass")], by = "kogClass", all.x = T)

TMM.matrix.new<-merge(unique(GenesKOGpair.1v1[,c("KS","Genes")]), TMM.matrix, by = "Genes", all.y = T)
TMM.matrix.new<-TMM.matrix.new[is.na(TMM.matrix.new$KS)==F, c("KS","Genes",treat.order)]

##Relative expression of all stages by PS
TMM.matrix.new[TMM.matrix.new<1]<-NA
Relative.exp<-as.data.frame(age.apply(ExpressionSet = TMM.matrix.new, RE))

Relative.exp$KS<-row.names(Relative.exp)
Relative.exp<-merge(Relative.exp, kog2name, by = "KS", all.x = T)
Relative.exp<-Relative.exp[order(as.numeric(Relative.exp$KS)),]
Relative.exp$Description<-paste0(" [",Relative.exp$kogClass,"] ", Relative.exp$kogName)

write.table(Relative.exp, paste0(species,".KOG.RE.txt"),
            quote = F, row.names = F, sep = "\t")
heatmap.input<-Relative.exp[,treat.order]
row.names(heatmap.input)<-Relative.exp$Description

row.anno<-data.frame(ONTOLOGY = factor(Relative.exp[,c("ONTOLOGY")]))
row.names(row.anno)<-Relative.exp$Description

col.anno<-data.frame(ONTOGENY = factor(c(rep("Early", length(module.list[[1]])),
                                         rep("Mid", length(module.list[[2]])),
                                         rep("Late", length(module.list[[3]])))))
row.names(col.anno)<-names(heatmap.input)

pheatmap(heatmap.input,
         legend = TRUE,
         scale = "none",
         #scale = "row",
         cluster_rows = F, cluster_cols = F,
         gaps_row = c(nrow(Relative.exp[Relative.exp$ONTOLOGY %in% c("Cellular processes and signaling"),]),
                      nrow(Relative.exp[Relative.exp$ONTOLOGY %in% c("Cellular processes and signaling",
                                                                     "Information storage and processing"),]),
                      nrow(Relative.exp)-1),
         annotation_row = row.anno,
         annotation_col = col.anno,
         annotation_colors = list(ONTOLOGY = c(`Cellular processes and signaling` = "#fdbf6f",
                                               `Information storage and processing` = "#a6cee3", 
                                               `Metabolism` = "#b2df8a",
                                               `Poor` = "gray75"),
                                  ONTOGENY = c(Early = "#33a02c",
                                               Mid = "#ff7f00",
                                               Late = "#6a3d9a")),
         cellwidth = 30, cellheight = 30, fontsize = 24,
         angle_col = c("90"),
         width = 24, height = 16,
         show_colnames = TRUE,
         show_rownames = T,
         filename = paste0(species,".KOG.RE.heatmap.png"))

##RE by ONTOLOGY
kog2name<-read.delim("/home/yichun/RNAmodification/hourglass/eggnog/kog2name.txt", header = TRUE)
kog2name<-kog2name[order(kog2name$ONTOLOGY,kog2name$kogClass),]
kog2name$KS[kog2name$ONTOLOGY=="CELLULAR PROCESSES AND SIGNALING"]<- 1
kog2name$KS[kog2name$ONTOLOGY=="INFORMATION STORAGE AND PROCESSING"]<- 2
kog2name$KS[kog2name$ONTOLOGY=="METABOLISM"]<- 3
kog2name$KS[kog2name$ONTOLOGY=="POORLY CHARACTERIZED"]<- 4

kog2name$ONTOLOGY<-gsub("INFORMATION STORAGE AND PROCESSING", "Information storage and processing", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("CELLULAR PROCESSES AND SIGNALING", "Cellular processes and signaling", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("POORLY CHARACTERIZED", "Poor", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("METABOLISM", "Metabolism", kog2name$ONTOLOGY)

GenesKOGpair.1v1<-read.delim(paste0("/home/yichun/RNAmodification/hourglass/eggnog/",
                                    species,".protein.fa.KOG.1v1.txt"), header = TRUE)
names(GenesKOGpair.1v1)<-c("KOG","Protein")
GenesKOGpair.1v1<-separate(GenesKOGpair.1v1, Protein, c("Genes"), sep = "\\|")
GenesKOGpair.1v1<-unique(GenesKOGpair.1v1[,c("KOG","Genes")])
names(GenesKOGpair.1v1)[1]<-"kogClass"
GenesKOGpair.1v1<-merge(GenesKOGpair.1v1, kog2name[,c("KS","kogClass")], by = "kogClass", all.x = T)

TMM.matrix.new<-merge(unique(GenesKOGpair.1v1[,c("KS","Genes")]), TMM.matrix, by = "Genes", all.y = T)
TMM.matrix.new<-TMM.matrix.new[is.na(TMM.matrix.new$KS)==F, c("KS","Genes",treat.order)]
##Relative expression of all stages by PS
TMM.matrix.new[TMM.matrix.new<1]<-NA
Relative.exp<-as.data.frame(age.apply(ExpressionSet = TMM.matrix.new, RE))

Relative.exp$KS<-row.names(Relative.exp)
Relative.exp<-merge(Relative.exp, unique(kog2name[,c("KS","ONTOLOGY")]), by = "KS", all.x = T)
Relative.exp<-Relative.exp[order(as.numeric(Relative.exp$KS)),]

write.table(Relative.exp, paste0(species,".KOG.RE.ONTOLOGY.txt"),
            quote = F, row.names = F, sep = "\t")
heatmap.input<-Relative.exp[,treat.order]
row.names(heatmap.input)<-Relative.exp$ONTOLOGY

col.anno<-data.frame(ONTOGENY = factor(c(rep("Early", length(module.list[[1]])),
                                         rep("Mid", length(module.list[[2]])),
                                         rep("Late", length(module.list[[3]])))))
row.names(col.anno)<-names(heatmap.input)

pheatmap(heatmap.input,
         legend = TRUE,
         scale = "none",
         #scale = "row",
         cluster_rows = F, cluster_cols = F,
         gaps_row = c(nrow(Relative.exp[Relative.exp$ONTOLOGY %in% c("Cellular processes and signaling"),]),
                      nrow(Relative.exp[Relative.exp$ONTOLOGY %in% c("Cellular processes and signaling",
                                                                     "Information storage and processing"),]),
                      nrow(Relative.exp)-1),
         annotation_col = col.anno,
         annotation_colors = list(ONTOGENY = c(Early = "#33a02c",
                                               Mid = "#ff7f00",
                                               Late = "#6a3d9a")),
         cellwidth = 30, cellheight = 30, fontsize = 24,
         angle_col = c("90"),
         width = 12, height = 6,
         show_colnames = TRUE,
         show_rownames = T,
         filename = paste0(species,".KOG.RE.ONTOLOGY.heatmap.png"))

line.input<-as.data.frame(matrix(NA, nrow = 0,ncol = 3))
names(line.input)<-c("Stage","RE","ONTOLOGY")
for (i in 2:(1+length(treat.order))) {
  line.input.sub<-Relative.exp[,c(i,ncol(Relative.exp))]
  line.input.sub$Stage<-names(Relative.exp[i])
  names(line.input.sub)[1]<-"RE"
  line.input<-rbind(line.input,line.input.sub)
}
p<-line.input %>%
  mutate(Stage = factor(Stage, levels = treat.order)) %>%
  ggplot(aes(x = Stage, y = RE, group = ONTOLOGY))+
  geom_vline(xintercept = min(module.list[[2]])-0.5, color = "grey75", linetype = "dashed")+
  geom_vline(xintercept = max(module.list[[2]])+0.5, color = "grey75", linetype = "dashed")+
  geom_line(aes(color = ONTOLOGY))+
  scale_color_manual(values =  c("#fdbf6f","#a6cee3", "#b2df8a","gray75"))+
  scale_y_continuous(limits = c(0,1.1),
                     breaks = seq(0,1,0.2),
                     position = "left")+
  scale_x_discrete(position = "bottom")+
  annotate(geom = "text", 
           label="Early",
           x=mean(module.list[[1]]), y=1.1,
           colour="black", size=2, fontface = "italic")+
  annotate(geom = "text", 
           label="Mid",
           x=mean(module.list[[2]]), y=1.1,
           colour="black", size=2, fontface = "italic")+
  annotate(geom = "text", 
           label="Late",
           x=mean(module.list[[3]]), y=1.1,
           colour="black", size=2, fontface = "italic")+
  labs(title = "", x = "", y = "Relative expression", colour = NULL)+
  theme(axis.line = element_line(linetype = "solid"),
        axis.ticks.y = element_line(colour = "black", size = 0.5),
        axis.ticks.x = element_line(colour = "black", size = 0.5),
        axis.text.x = element_text(size = 8, colour = "black", angle =45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 8, colour = "black"),
        plot.title = element_text(size = 8, hjust = 0.5, face = "plain"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 0),
        panel.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        legend.position = "right",
        strip.text = element_text(size = 9))
p
ggsave(paste0(species,".KOG.RE.ONTOLOGY.line.png"), 
       width = 6, height = 2.5, units = "in", dpi = 300)
f=paste0(species,".KOG.RE.ONTOLOGY.line.pptx")
topptx(p, f, width = 6, height = 2.5, units = "in")

################################################################
## KOG ~ Early-Mid-Late 
################################################################
TMM.matrix.EML<-data.frame(Genes = TMM.matrix$Genes)
if (length(module.list[[1]])<2) {
  TMM.matrix.EML$Early<-TMM.matrix[,1+module.list[[1]]]
} else {
  TMM.matrix.EML$Early<-apply(TMM.matrix[,1+module.list[[1]]], 1, mean)
}
TMM.matrix.EML$Mid<-apply(TMM.matrix[,1+module.list[[2]]], 1, mean)
TMM.matrix.EML$Late<-apply(TMM.matrix[,1+module.list[[3]]], 1, mean)
treat.order<-c("Early","Mid","Late")

##RE by kogclass
kog2name<-read.delim("/home/yichun/RNAmodification/hourglass/eggnog/kog2name.txt", header = TRUE)
kog2name<-kog2name[order(kog2name$ONTOLOGY,kog2name$kogClass),]
kog2name$KS<-1:nrow(kog2name)
kog2name$ONTOLOGY<-gsub("INFORMATION STORAGE AND PROCESSING", "Information storage and processing", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("CELLULAR PROCESSES AND SIGNALING", "Cellular processes and signaling", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("POORLY CHARACTERIZED", "Poor", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("METABOLISM", "Metabolism", kog2name$ONTOLOGY)

GenesKOGpair.1v1<-read.delim(paste0("/home/yichun/RNAmodification/hourglass/eggnog/",
                                    species,".protein.fa.KOG.1v1.txt"), header = TRUE)
names(GenesKOGpair.1v1)<-c("KOG","Protein")
GenesKOGpair.1v1<-separate(GenesKOGpair.1v1, Protein, c("Genes"), sep = "\\|")
GenesKOGpair.1v1<-unique(GenesKOGpair.1v1[,c("KOG","Genes")])
names(GenesKOGpair.1v1)[1]<-"kogClass"
GenesKOGpair.1v1<-merge(GenesKOGpair.1v1, kog2name[,c("KS","kogClass")], by = "kogClass", all.x = T)

TMM.matrix.new<-merge(unique(GenesKOGpair.1v1[,c("KS","Genes")]), TMM.matrix.EML, by = "Genes", all.y = T)
TMM.matrix.new<-TMM.matrix.new[is.na(TMM.matrix.new$KS)==F, c("KS","Genes",treat.order)]
##Relative expression of all stages by PS
TMM.matrix.new[TMM.matrix.new<1]<-NA
Relative.exp<-as.data.frame(age.apply(ExpressionSet = TMM.matrix.new, RE))

Relative.exp$KS<-row.names(Relative.exp)
Relative.exp<-merge(Relative.exp, kog2name, by = "KS", all.x = T)
Relative.exp<-Relative.exp[order(as.numeric(Relative.exp$KS)),]
Relative.exp$Description<-paste0(" [",Relative.exp$kogClass,"] ", Relative.exp$kogName)

write.table(Relative.exp, paste0(species,".EML.KOG.RE.txt"),
            quote = F, row.names = F, sep = "\t")
heatmap.input<-Relative.exp[,treat.order]
row.names(heatmap.input)<-Relative.exp$Description

row.anno<-data.frame(ONTOLOGY = factor(Relative.exp[,c("ONTOLOGY")]))
row.names(row.anno)<-Relative.exp$Description

pheatmap(heatmap.input,
         legend = TRUE,
         scale = "none",
         #scale = "row",
         cluster_rows = F, cluster_cols = F,
         gaps_row = c(nrow(Relative.exp[Relative.exp$ONTOLOGY %in% c("Cellular processes and signaling"),]),
                      nrow(Relative.exp[Relative.exp$ONTOLOGY %in% c("Cellular processes and signaling",
                                                                     "Information storage and processing"),]),
                      nrow(Relative.exp)-1),
         annotation_row = row.anno,
         annotation_colors = list(ONTOLOGY = c(`Cellular processes and signaling` = "#fdbf6f",
                                               `Information storage and processing` = "#a6cee3", 
                                               `Metabolism` = "#b2df8a",
                                               `Poor` = "gray75"),
                                  ONTOGENY = c(Early = "#33a02c",
                                               Mid = "#ff7f00",
                                               Late = "#6a3d9a")),
         cellwidth = 30, cellheight = 30, fontsize = 24,
         angle_col = c("90"),
         width = 24, height = 16,
         show_colnames = TRUE,
         show_rownames = T,
         filename = paste0(species,".EML.KOG.RE.heatmap.png"))

##RE by ONTOLOGY
kog2name<-read.delim("/home/yichun/RNAmodification/hourglass/eggnog/kog2name.txt", header = TRUE)
kog2name<-kog2name[order(kog2name$ONTOLOGY,kog2name$kogClass),]
kog2name$KS[kog2name$ONTOLOGY=="CELLULAR PROCESSES AND SIGNALING"]<- 1
kog2name$KS[kog2name$ONTOLOGY=="INFORMATION STORAGE AND PROCESSING"]<- 2
kog2name$KS[kog2name$ONTOLOGY=="METABOLISM"]<- 3
kog2name$KS[kog2name$ONTOLOGY=="POORLY CHARACTERIZED"]<- 4

kog2name$ONTOLOGY<-gsub("INFORMATION STORAGE AND PROCESSING", "Information storage and processing", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("CELLULAR PROCESSES AND SIGNALING", "Cellular processes and signaling", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("POORLY CHARACTERIZED", "Poor", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("METABOLISM", "Metabolism", kog2name$ONTOLOGY)

GenesKOGpair.1v1<-read.delim(paste0("/home/yichun/RNAmodification/hourglass/eggnog/",
                                    species,".protein.fa.KOG.1v1.txt"), header = TRUE)
names(GenesKOGpair.1v1)<-c("KOG","Protein")
GenesKOGpair.1v1<-separate(GenesKOGpair.1v1, Protein, c("Genes"), sep = "\\|")
GenesKOGpair.1v1<-unique(GenesKOGpair.1v1[,c("KOG","Genes")])
names(GenesKOGpair.1v1)[1]<-"kogClass"
GenesKOGpair.1v1<-merge(GenesKOGpair.1v1, kog2name[,c("KS","kogClass")], by = "kogClass", all.x = T)

TMM.matrix.new<-merge(GenesKOGpair.1v1[,c("KS","Genes")], TMM.matrix.EML, by = "Genes", all.y = T)
TMM.matrix.new<-TMM.matrix.new[is.na(TMM.matrix.new$KS)==F, c("KS","Genes",treat.order)]
##Relative expression of all stages by PS
TMM.matrix.new[TMM.matrix.new<1]<-NA
Relative.exp<-as.data.frame(age.apply(ExpressionSet = TMM.matrix.new, RE))

Relative.exp$KS<-row.names(Relative.exp)
Relative.exp<-merge(Relative.exp, unique(kog2name[,c("KS","ONTOLOGY")]), by = "KS", all.x = T)
Relative.exp<-Relative.exp[order(as.numeric(Relative.exp$KS)),]

write.table(Relative.exp, paste0(species,".EML.KOG.RE.ONTOLOGY.txt"),
            quote = F, row.names = F, sep = "\t")
heatmap.input<-Relative.exp[,treat.order]
row.names(heatmap.input)<-Relative.exp$ONTOLOGY

pheatmap(heatmap.input,
         legend = TRUE,
         scale = "none",
         #scale = "row",
         cluster_rows = F, cluster_cols = F,
         gaps_row = c(nrow(Relative.exp[Relative.exp$ONTOLOGY %in% c("Cellular processes and signaling"),]),
                      nrow(Relative.exp[Relative.exp$ONTOLOGY %in% c("Cellular processes and signaling",
                                                                     "Information storage and processing"),]),
                      nrow(Relative.exp)-1),
         cellwidth = 30, cellheight = 30, fontsize = 24,
         angle_col = c("90"),
         width = 8, height = 6,
         show_colnames = TRUE,
         show_rownames = T,
         filename = paste0(species,".EML.KOG.RE.ONTOLOGY.heatmap.png"))

line.input<-as.data.frame(matrix(NA, nrow = 0,ncol = 3))
names(line.input)<-c("Stage","RE","ONTOLOGY")
for (i in 2:(1+length(treat.order))) {
  line.input.sub<-Relative.exp[,c(i,ncol(Relative.exp))]
  line.input.sub$Stage<-names(Relative.exp[i])
  names(line.input.sub)[1]<-"RE"
  line.input<-rbind(line.input,line.input.sub)
}
p<-line.input %>%
  mutate(Stage = factor(Stage, levels = treat.order)) %>%
  ggplot(aes(x = Stage, y = RE, group = ONTOLOGY))+
  geom_vline(xintercept = 1.5, color = "grey75", linetype = "dashed")+
  geom_vline(xintercept = 2.5, color = "grey75", linetype = "dashed")+
  geom_line(aes(color = ONTOLOGY))+
  scale_color_manual(values =  c("#fdbf6f","#a6cee3", "#b2df8a","gray75"))+
  scale_y_continuous(limits = c(0,1.1),
                     breaks = seq(0,1,0.2),
                     position = "left")+
  scale_x_discrete(position = "bottom")+
  annotate(geom = "text", 
           label="Early",
           x=1, y=1.1,
           colour="black", size=2, fontface = "italic")+
  annotate(geom = "text", 
           label="Mid",
           x=2, y=1.1,
           colour="black", size=2, fontface = "italic")+
  annotate(geom = "text", 
           label="Late",
           x=3, y=1.1,
           colour="black", size=2, fontface = "italic")+
  labs(title = "", x = "", y = "Relative expression", colour = NULL)+
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
        legend.position = "right",
        strip.text = element_text(size = 9))
p
ggsave(paste0(species,".EML.KOG.RE.ONTOLOGY.line.png"), 
       width = 5, height = 2.5, units = "in", dpi = 300)
f=paste0(species,".EML.KOG.RE.ONTOLOGY.line.pptx")
topptx(p, f, width = 5, height = 2.5, units = "in")

################################################################
##Plot all species
################################################################
kog2name<-read.delim("/home/yichun/RNAmodification/hourglass/eggnog/kog2name.txt", header = TRUE)
kog2name<-kog2name[order(kog2name$ONTOLOGY,kog2name$kogClass),]
kog2name$KS<-1:nrow(kog2name)
kog2name$ONTOLOGY<-gsub("INFORMATION STORAGE AND PROCESSING", "Information storage and processing", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("CELLULAR PROCESSES AND SIGNALING", "Cellular processes and signaling", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("POORLY CHARACTERIZED", "Poor", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("METABOLISM", "Metabolism", kog2name$ONTOLOGY)
specieslist<-c("Rhizopus_delemar",
               "Fusarium_graminearum",
               "Coprinopsis_cinerea",
               "Drosophila_melanogaster",
               "Danio_rerio",
               "Arabidopsis_thaliana")
####pheatmap
KOG.all<-data.frame(kogClass = kog2name$kogClass[kog2name$kogClass!="R"])

for (i in 1:length(specieslist)) {
  species<-specieslist[i]
  KOGtable.sub<-read.delim(paste0(species,".EML.KOG.RE.txt"))
  KOGtable.sub<-KOGtable.sub[,c("Early","Mid","Late","kogClass")]
  names(KOGtable.sub)[1:3]<-paste0(species,"  ",names(KOGtable.sub)[1:3])
  KOG.all<-merge(KOG.all,KOGtable.sub, by = c("kogClass"), all = T)
}
rm(KOGtable.sub)
KOG.all<-merge(KOG.all, kog2name, by = "kogClass", all.x = T)
KOG.all<-KOG.all[order(KOG.all$ONTOLOGY,KOG.all$kogClass),]
KOG.all$Description<-paste0("[",KOG.all$kogClass,"] ",KOG.all$kogName)
row.names(KOG.all)<-KOG.all$Description

heatmap.input<-KOG.all[,2:19]

row.anno<-data.frame(ONTOLOGY = factor(KOG.all[,c("ONTOLOGY")]))
row.names(row.anno)<-KOG.all$Description

heatmap.input<-heatmap.input[,c(paste0(specieslist,"  ","Early"),
                                paste0(specieslist,"  ","Mid"),
                                paste0(specieslist,"  ","Late"))]
col.anno<-data.frame(ONTOGENY = factor(c(rep("Early", length(specieslist)),
                                         rep("Mid", length(specieslist)),
                                         rep("Late", length(specieslist)))))
row.names(col.anno)<-names(heatmap.input)

pheatmap(heatmap.input,
         legend = TRUE,
         scale = "none",
         #scale = "row",
         cluster_rows = F, cluster_cols = F,
         gaps_row = c(nrow(KOG.all[KOG.all$ONTOLOGY %in% c("Cellular processes and signaling"),]),
                      nrow(KOG.all[KOG.all$ONTOLOGY %in% c("Cellular processes and signaling",
                                                           "Information storage and processing"),]),
                      nrow(KOG.all)-1),
         gaps_col = c(length(specieslist),2*length(specieslist)),
         annotation_row = row.anno,
         annotation_col = col.anno,
         annotation_colors = list(ONTOLOGY = c(`Cellular processes and signaling` = "#fdbf6f",
                                               `Information storage and processing` = "#a6cee3", 
                                               `Metabolism` = "#b2df8a",
                                               `Poor` = "gray75"),
                                  ONTOGENY = c(Early = "#33a02c",
                                               Mid = "#ff7f00",
                                               Late = "#6a3d9a")),
         cellwidth = 30, cellheight = 30, fontsize = 24,
         angle_col = c("90"),
         width = 28, height = 16,
         show_colnames = TRUE,
         show_rownames = T,
         filename = paste0("All.KOG.RE.EML.heatmap.png"))

####ggplot2
KOG.all.list.bk<-as.data.frame(matrix(NA, nrow = 0, ncol = 4))
names(KOG.all.list.bk)<-c("kogClass","species","stage","RE")
KOG.all.list<-KOG.all.list.bk

for (i in 1:length(specieslist)) {
  KOGtable.species<-KOG.all.list.bk
  species<-specieslist[i]
  KOGtable.sub<-read.delim(paste0(species,".EML.KOG.RE.txt"))
  
  KOGtable.sub<-KOGtable.sub[,c("kogClass","Early","Mid","Late")]
  for (j in 2:ncol(KOGtable.sub)) {
    KOGtable.sub.list<-KOGtable.sub[,c(1,j)]
    names(KOGtable.sub.list)[2]<-"RE"
    KOGtable.sub.list$stage<-names(KOGtable.sub)[j]
    KOGtable.sub.list$species<-species
    KOGtable.species<-rbind(KOGtable.species,KOGtable.sub.list)
  }
  KOG.all.list<-rbind(KOG.all.list,KOGtable.species)
}

kog2name<-kog2name[kog2name$kogClass!="R",]
kog2name$Description<-paste0(kog2name$kogName," [",kog2name$kogClass,"]")
kog2name<-kog2name[order(kog2name$ONTOLOGY,kog2name$kogClass),]
KOG.all.list<-merge(KOG.all.list,kog2name, by = "kogClass", all.x = T)
kog2name.order<-kog2name$Description

KOG.all.list$ONTOLOGY<-gsub("Cellular processes and signaling",
                            "Cellular processes\nand signaling",
                            KOG.all.list$ONTOLOGY)
KOG.all.list$ONTOLOGY<-gsub("Information storage and processing",
                            "Information storage\nand processing",
                            KOG.all.list$ONTOLOGY)

KOG.all.list$species<-gsub("Rhizopus_delemar","R._delemar",KOG.all.list$species)
KOG.all.list$species<-gsub("Fusarium_graminearum","F._graminearum",KOG.all.list$species)
KOG.all.list$species<-gsub("Coprinopsis_cinerea","C._cinerea",KOG.all.list$species)
KOG.all.list$species<-gsub("Drosophila_melanogaster","D._melanogaster",KOG.all.list$species)
KOG.all.list$species<-gsub("Danio_rerio","Da._rerio",KOG.all.list$species)
KOG.all.list$species<-gsub("Arabidopsis_thaliana","A._thaliana",KOG.all.list$species)
KOG.all.list$species<-gsub("_"," ",KOG.all.list$species)
specieslist.new<-c("R._delemar",
               "F._graminearum",
               "C._cinerea",
               "D._melanogaster",
               "Da._rerio",
               "A._thaliana")

p<-KOG.all.list %>%
  mutate(species=factor(species, levels = c(gsub("_"," ",specieslist.new))),
         stage = factor(stage, levels = c("Early","Mid","Late")),
         Description = factor(Description, levels = kog2name.order)) %>% 
  ggplot(aes(x=species, y=Description))+
  geom_tile(aes(fill = RE-0.5))+
  #scale_fill_gradient2(low = "#91BFDB", mid = "#FFFFBF", high = "#FC8D59", 
  scale_fill_gradient2(low = "#4575B4", mid = "#FFFFBF", high = "#D73027",
                       limits = c(-0.5,0.5), n.breaks = 3, labels = c("Low","Medium","High"))+
  labs(y = "", x = "", fill = "Relative\nexpression", title = "")+
  theme(plot.subtitle = element_text(vjust = 1),
        plot.caption = element_text(vjust = 1),
        axis.text.x = element_text(size = 12, colour = "black", angle = 45, face = "italic", hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "left",
        axis.line = element_line(size = 0.5, colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.background = element_rect(fill = NA),
        panel.grid.minor = element_line(colour = "gray50", size = 0.5),
        strip.text = element_text(size = 10))+
  facet_grid(ONTOLOGY~stage, scales = "free_y", space = "free")
p
ggsave(paste0("All.KOG.RE.EML.heatmapgg.png"), width = 12, height = 8, units = "in", dpi = 300)
ggsave(paste0("All.KOG.RE.EML.heatmapgg.tiff"), width = 12, height = 8, units = "in", dpi = 300)
f=paste0("All.KOG.RE.EML.heatmapgg.pptx")
topptx(p, f, width = 12, height = 8, units = "in")
