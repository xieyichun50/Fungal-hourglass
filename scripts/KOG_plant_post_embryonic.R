library(dplyr)
library(tidyr)
library(myTAI)
library(ggplot2)
library(pheatmap)
library(eoffice)
setwd("hourglass/enrich/")


##Arabidopsis
################################################################
## KOG ~ all stages 
################################################################
##Read in gene expression level
species="Arabidopsis_thaliana"
treat="germ"
module.list<-list(early = 1:2, mid = 3:5, late = 6:7)
treat="flor"
module.list<-list(early = 1:3, mid = 4:6, late = 7:9)

TMM.matrix<-read.delim(paste0("expression.atha.",treat,".txt"),header = T)
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
GenesKOGpair.1v1<-separate(GenesKOGpair.1v1, Genes, c("Genes"), sep = "\\.")
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

write.table(Relative.exp, paste0(species,".",treat,".KOG.RE.txt"),
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
         filename = paste0(species,".",treat,".KOG.RE.heatmap.png"))

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
GenesKOGpair.1v1<-separate(GenesKOGpair.1v1, Genes, c("Genes"), sep = "\\.")
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

write.table(Relative.exp, paste0(species,".",treat,".KOG.RE.ONTOLOGY.txt"),
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
         filename = paste0(species,".",treat,".KOG.RE.ONTOLOGY.heatmap.png"))

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
f=paste0(species,".",treat,".KOG.RE.ONTOLOGY.line.pptx")
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
GenesKOGpair.1v1<-separate(GenesKOGpair.1v1, Genes, c("Genes"), sep = "\\.")
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

write.table(Relative.exp, paste0(species,".",treat,".EML.KOG.RE.txt"),
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
         filename = paste0(species,".",treat,".EML.KOG.RE.heatmap.png"))

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
GenesKOGpair.1v1<-separate(GenesKOGpair.1v1, Genes, c("Genes"), sep = "\\.")
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

write.table(Relative.exp, paste0(species,".",treat,".EML.KOG.RE.ONTOLOGY.txt"),
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
         filename = paste0(species,".",treat,".EML.KOG.RE.ONTOLOGY.heatmap.png"))

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
ggsave(paste0(species,".",treat,".EML.KOG.RE.ONTOLOGY.line.png"), 
       width = 5, height = 2.5, units = "in", dpi = 300)
f=paste0(species,".",treat,".EML.KOG.RE.ONTOLOGY.line.pptx")
topptx(p, f, width = 5, height = 2.5, units = "in")

