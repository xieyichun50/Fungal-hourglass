library(dplyr)
library(tidyr)
library(myTAI)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
setwd("hourglass/enrich/")

normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}
kog2name<-read.delim("/home/yichun/3enrichment/kog2name.txt", header = TRUE)

##########################################################
####Drosophila
##########################################################
species="Drosophila_melanogaster"
IDmatch<-read.delim(paste0("/home/yichun/RNAmodification/hourglass/eggnog/",
                           species,".GenematchID"), header = F)
names(IDmatch)<-c("Protein","Gene")
exp<-read.delim(paste0("expression.dmel.txt"),header = T)
names(exp)[2]<-"Protein"
exp$Protein<-gsub("fbpp","FBpp",exp$Protein)
exp[,3:ncol(exp)]<-log2(exp[,3:ncol(exp)]+1)
##calling with diffgene
exp.diff <-DiffGenes(exp,
                     method = "log-foldchange",
                     p.adjust.method = "BH",
                     nrep = c(3,2,7),
                     stage.names = c("Early", "Mid", "Late"))
exp.diff$Groups[exp.diff$`Early->Mid` < -1 & exp.diff$`Mid->Late` > 1]<-"LHL"
exp.diff$Groups[exp.diff$`Early->Mid` > 1 & exp.diff$`Mid->Late` < -1]<-"HLH"
exp.diff<-merge(IDmatch,exp.diff,by="Protein", all.y = T)

heatmap.input<-exp[exp$Protein %in% exp.diff$Protein[exp.diff$Groups %in% c("LHL","HLH")], ]
row.names(heatmap.input)<-heatmap.input$Gene
heatmap.input<-heatmap.input[,4:ncol(heatmap.input)]
names(heatmap.input)<-gsub("X","",names(heatmap.input))
names(heatmap.input)<-gsub("\\.embryo","",names(heatmap.input))

heatmap.input<-as.data.frame(t(apply(heatmap.input, 1, normalize)), row.names = row.names(heatmap.input))
heatp<-pheatmap(heatmap.input,
                legend = TRUE,
                scale = "none",
                #scale = "row",
                cluster_rows = TRUE, cluster_cols = FALSE,
                cutree_rows = NA, treeheight_row = 100,
                border_color = NA,
                cellwidth = 30, cellheight = 2, fontsize = 24,
                angle_col = c("90"),
                width = ncol(heatmap.input)*0.6+1, height = nrow(heatmap.input)*0.04,
                show_colnames = TRUE,
                show_rownames = F,
                filename = paste0(species,".heatmap.LHL.png"))

## Enrich
GenesKOGpair.1v1<-read.delim(paste0("/home/yichun/RNAmodification/hourglass/eggnog/",
                                    species,".proteins.fa.KOG.1v1.txt"), header = TRUE)
names(GenesKOGpair.1v1)<-c("KOG","Protein")
GenesKOGpair.1v1<-merge(GenesKOGpair.1v1, IDmatch, by = "Protein", all.x = T)
GenesKOGpair.1v1<-unique(GenesKOGpair.1v1[,c("KOG","Gene")])
waist<-exp.diff[exp.diff$Groups %in% c("LHL","HLH"),c("Gene","Groups")]
View(as.data.frame(xtabs(~Gene, waist)))

KOG.all.1<-compareCluster(Gene ~ Groups,
                          data = waist,
                          fun = 'enricher',
                          TERM2GENE = GenesKOGpair.1v1,
                          TERM2NAME = kog2name,
                          pvalueCutoff = 1,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 1,
                          minGSSize = 1,
                          maxGSSize = nrow(GenesKOGpair.1v1))
KOG.all.1<-as.data.frame(KOG.all.1)
KOG.all.1<-subset(KOG.all.1, is.na(Description)==FALSE)
write.table(KOG.all.1, paste0(species,".KOG.txt"),
            row.names = F, quote = F, sep = "\t")

##########################################################
####Danio
##########################################################
species="Danio_rerio"
IDmatch<-read.delim(paste0("/home/yichun/RNAmodification/hourglass/eggnog/",
                           species,".GenematchID"), header = F)
names(IDmatch)<-c("Protein","Gene")
exp<-read.delim(paste0("expression.drer.txt"),header = T)
names(exp)[2]<-"Gene"
exp[,3:ncol(exp)]<-log2(exp[,3:ncol(exp)]+1)
##calling with diffgene
exp.diff <-DiffGenes(exp,
                     method = "log-foldchange",
                     p.adjust.method = "BH",
                     nrep = c(18,18,4),
                     stage.names = c("Early", "Mid", "Late"))
exp.diff$Groups[exp.diff$`Early->Mid` < -1 & exp.diff$`Mid->Late` > 1]<-"LHL"
exp.diff$Groups[exp.diff$`Early->Mid` > 1 & exp.diff$`Mid->Late` < -1]<-"HLH"

heatmap.input<-exp[exp$Gene %in% exp.diff$Gene[exp.diff$Groups %in% c("LHL","HLH")], ]
row.names(heatmap.input)<-heatmap.input$Gene
heatmap.input<-heatmap.input[,4:ncol(heatmap.input)]
names(heatmap.input)<-gsub("X","",names(heatmap.input))
names(heatmap.input)<-gsub("\\.embryo","",names(heatmap.input))

heatmap.input<-as.data.frame(t(apply(heatmap.input, 1, normalize)), row.names = row.names(heatmap.input))
heatp<-pheatmap(heatmap.input,
                legend = TRUE,
                scale = "none",
                #scale = "row",
                cluster_rows = TRUE, cluster_cols = FALSE,
                cutree_rows = NA, treeheight_row = 100,
                border_color = NA,
                cellwidth = 30, cellheight = 2, fontsize = 24,
                angle_col = c("90"),
                width = ncol(heatmap.input)*0.6+1, height = nrow(heatmap.input)*0.04,
                show_colnames = TRUE,
                show_rownames = F,
                filename = paste0(species,".heatmap.LHL.png"))

## Enrich
GenesKOGpair.1v1<-read.delim(paste0("/home/yichun/RNAmodification/hourglass/eggnog/",
                                    species,".proteins.fa.KOG.1v1.txt"), header = TRUE)
names(GenesKOGpair.1v1)<-c("KOG","Protein")
GenesKOGpair.1v1<-merge(GenesKOGpair.1v1, IDmatch, by = "Protein", all.x = T)
GenesKOGpair.1v1<-unique(GenesKOGpair.1v1[,c("KOG","Gene")])
waist<-exp.diff[exp.diff$Groups %in% c("LHL","HLH"),c("Gene","Groups")]
View(as.data.frame(xtabs(~Gene, waist)))

KOG.all.1<-compareCluster(Gene ~ Groups,
                          data = waist,
                          fun = 'enricher',
                          TERM2GENE = GenesKOGpair.1v1,
                          TERM2NAME = kog2name,
                          pvalueCutoff = 1,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 1,
                          minGSSize = 1,
                          maxGSSize = nrow(GenesKOGpair.1v1))
KOG.all.1<-as.data.frame(KOG.all.1)
KOG.all.1<-subset(KOG.all.1, is.na(Description)==FALSE)
write.table(KOG.all.1, paste0(species,".KOG.txt"),
            row.names = F, quote = F, sep = "\t")

##########################################################
####Arabidopsis
##########################################################
species="Arabidopsis_thaliana"
exp<-read.delim(paste0("expression.atha.txt"),header = T)
names(exp)[2]<-"Gene"

exp$Gene<-gsub("at1g","AT1G",exp$Gene)
exp$Gene<-gsub("at2g","AT2G",exp$Gene)
exp$Gene<-gsub("at3g","AT3G",exp$Gene)
exp$Gene<-gsub("at4g","AT4G",exp$Gene)
exp$Gene<-gsub("at5g","AT5G",exp$Gene)
exp$Gene<-gsub("atcg","ATCG",exp$Gene)
exp$Gene<-gsub("atmg","ATMG",exp$Gene)

exp[,3:ncol(exp)]<-log2(exp[,3:ncol(exp)]+1)
##calling with diffgene
exp.diff <-DiffGenes(exp,
                     method = "log-foldchange",
                     p.adjust.method = "BH",
                     nrep = c(2,3,2),
                     stage.names = c("Early", "Mid", "Late"))
exp.diff$Groups[exp.diff$`Early->Mid` < -1 & exp.diff$`Mid->Late` > 1]<-"LHL"
exp.diff$Groups[exp.diff$`Early->Mid` > 1 & exp.diff$`Mid->Late` < -1]<-"HLH"

heatmap.input<-exp[exp$Gene %in% exp.diff$Gene[exp.diff$Groups %in% c("LHL","HLH")], ]
row.names(heatmap.input)<-heatmap.input$Gene
heatmap.input<-heatmap.input[,4:ncol(heatmap.input)]

heatmap.input<-as.data.frame(t(apply(heatmap.input, 1, normalize)), row.names = row.names(heatmap.input))
heatp<-pheatmap(heatmap.input,
                legend = TRUE,
                scale = "none",
                #scale = "row",
                cluster_rows = TRUE, cluster_cols = FALSE,
                cutree_rows = NA, treeheight_row = 100,
                border_color = NA,
                cellwidth = 30, cellheight = 2, fontsize = 24,
                angle_col = c("90"),
                width = 8, height = nrow(heatmap.input)*0.04,
                show_colnames = TRUE,
                show_rownames = F,
                filename = paste0(species,".heatmap.LHL.png"))

## Enrich
GenesKOGpair.1v1<-read.delim(paste0("/home/yichun/RNAmodification/hourglass/eggnog/",
                                    species,".protein.fa.KOG.1v1.txt"), header = TRUE)
names(GenesKOGpair.1v1)<-c("KOG","Protein")
GenesKOGpair.1v1<-separate(GenesKOGpair.1v1, Protein, c("Gene"), sep = "\\|")
GenesKOGpair.1v1<-unique(GenesKOGpair.1v1[,c("KOG","Gene")])
waist<-exp.diff[exp.diff$Groups %in% c("LHL","HLH"),c("Gene","Groups")]
View(as.data.frame(xtabs(~Gene, waist)))

KOG.all.1<-compareCluster(Gene ~ Groups,
                          data = waist,
                          fun = 'enricher',
                          TERM2GENE = GenesKOGpair.1v1,
                          TERM2NAME = kog2name,
                          pvalueCutoff = 1,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 1,
                          minGSSize = 1,
                          maxGSSize = nrow(GenesKOGpair.1v1))
KOG.all.1<-as.data.frame(KOG.all.1)
KOG.all.1<-subset(KOG.all.1, is.na(Description)==FALSE)
write.table(KOG.all.1, paste0(species,".KOG.txt"),
            row.names = F, quote = F, sep = "\t")

##########################################################
####Rhizopus
##########################################################
species="Rhizopus_delemar"
IDmatch<-read.delim(paste0("/home/yichun/RNAmodification/hourglass/eggnog/",
                           species,".GenematchID"), header = F)
names(IDmatch)<-c("Protein","Gene")
exp<-read.delim(paste0("Rhizopus_delemar.log2TMMmatrixmean.txt"),header = T)
exp$Groups<-1:nrow(exp)
exp<-exp[,c(ncol(exp),1:2,2:(ncol(exp)-1))]
##calling with diffgene
exp.diff <-DiffGenes(exp,
                     method = "log-foldchange",
                     p.adjust.method = "BH",
                     nrep = c(2,4,5),
                     stage.names = c("Early", "Mid", "Late"))
exp.diff$Groups[exp.diff$`Early->Mid` < -1 & exp.diff$`Mid->Late` > 1]<-"LHL"
exp.diff$Groups[exp.diff$`Early->Mid` > 1 & exp.diff$`Mid->Late` < -1]<-"HLH"

heatmap.input<-exp[exp$Gene %in% exp.diff$Gene[exp.diff$Groups %in% c("LHL","HLH")], ]
row.names(heatmap.input)<-heatmap.input$Gene
heatmap.input<-heatmap.input[,c(3,5:ncol(heatmap.input))]

heatmap.input<-as.data.frame(t(apply(heatmap.input, 1, normalize)), row.names = row.names(heatmap.input))
heatp<-pheatmap(heatmap.input,
                legend = TRUE,
                scale = "none",
                #scale = "row",
                cluster_rows = TRUE, cluster_cols = FALSE,
                cutree_rows = NA, treeheight_row = 100,
                border_color = NA,
                cellwidth = 30, cellheight = 0.5, fontsize = 24,
                angle_col = c("90"),
                width = ncol(heatmap.input)*0.6+1, height = nrow(heatmap.input)*0.01,
                show_colnames = TRUE,
                show_rownames = F,
                filename = paste0(species,".heatmap.LHL.png"))

## Enrich
GenesKOGpair.1v1<-read.delim(paste0("/home/yichun/RNAmodification/hourglass/eggnog/",
                                    species,".protein.fa.KOG.1v1.txt"), header = TRUE)
names(GenesKOGpair.1v1)<-c("KOG","Protein")
GenesKOGpair.1v1<-merge(GenesKOGpair.1v1, IDmatch, by = "Protein", all.x = T)
GenesKOGpair.1v1<-unique(GenesKOGpair.1v1[,c("KOG","Gene")])
waist<-exp.diff[exp.diff$Groups %in% c("LHL","HLH"),c("Genes","Groups")]
View(as.data.frame(xtabs(~Genes, waist)))

KOG.all.1<-compareCluster(Genes ~ Groups,
                          data = waist,
                          fun = 'enricher',
                          TERM2GENE = GenesKOGpair.1v1,
                          TERM2NAME = kog2name,
                          pvalueCutoff = 1,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 1,
                          minGSSize = 1,
                          maxGSSize = nrow(GenesKOGpair.1v1))
KOG.all.1<-as.data.frame(KOG.all.1)
KOG.all.1<-subset(KOG.all.1, is.na(Description)==FALSE)
write.table(KOG.all.1, paste0(species,".KOG.txt"),
            row.names = F, quote = F, sep = "\t")

##########################################################
####Fusarium
##########################################################
species="Fusarium_graminearum"
IDmatch<-read.delim(paste0("/home/yichun/RNAmodification/hourglass/eggnog/",
                           species,".GenematchID"), header = F)
names(IDmatch)<-c("Protein","Gene")
exp<-read.delim(paste0("Fusarium_graminearum.log2TMMmatrixmean.txt"),header = T)
exp$Groups<-1:nrow(exp)
exp<-exp[,c(ncol(exp),1:2,2:(ncol(exp)-1))]
##calling with diffgene
exp.diff <-DiffGenes(exp,
                     method = "log-foldchange",
                     p.adjust.method = "BH",
                     nrep = c(2,3,6),
                     stage.names = c("Early", "Mid", "Late"))
exp.diff$Groups[exp.diff$`Early->Mid` < -1 & exp.diff$`Mid->Late` > 1]<-"LHL"
exp.diff$Groups[exp.diff$`Early->Mid` > 1 & exp.diff$`Mid->Late` < -1]<-"HLH"

heatmap.input<-exp[exp$Gene %in% exp.diff$Gene[exp.diff$Groups %in% c("LHL","HLH")], ]
row.names(heatmap.input)<-heatmap.input$Gene
heatmap.input<-heatmap.input[,c(3,5:ncol(heatmap.input))]

heatmap.input<-as.data.frame(t(apply(heatmap.input, 1, normalize)), row.names = row.names(heatmap.input))
heatp<-pheatmap(heatmap.input,
                legend = TRUE,
                scale = "none",
                #scale = "row",
                cluster_rows = TRUE, cluster_cols = FALSE,
                cutree_rows = NA, treeheight_row = 100,
                border_color = NA,
                cellwidth = 30, cellheight = 0.5, fontsize = 24,
                angle_col = c("90"),
                width = ncol(heatmap.input)*0.6+1, height = nrow(heatmap.input)*0.01,
                show_colnames = TRUE,
                show_rownames = F,
                filename = paste0(species,".heatmap.LHL.png"))

## Enrich
GenesKOGpair.1v1<-read.delim(paste0("/home/yichun/RNAmodification/hourglass/eggnog/",
                                    species,".protein.fa.KOG.1v1.txt"), header = TRUE)
names(GenesKOGpair.1v1)<-c("KOG","Protein")
GenesKOGpair.1v1<-merge(GenesKOGpair.1v1, IDmatch, by = "Protein", all.x = T)
GenesKOGpair.1v1<-unique(GenesKOGpair.1v1[,c("KOG","Gene")])
waist<-exp.diff[exp.diff$Groups %in% c("LHL","HLH"),c("Genes","Groups")]
View(as.data.frame(xtabs(~Genes, waist)))

KOG.all.1<-compareCluster(Genes ~ Groups,
                          data = waist,
                          fun = 'enricher',
                          TERM2GENE = GenesKOGpair.1v1,
                          TERM2NAME = kog2name,
                          pvalueCutoff = 1,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 1,
                          minGSSize = 1,
                          maxGSSize = nrow(GenesKOGpair.1v1))
KOG.all.1<-as.data.frame(KOG.all.1)
KOG.all.1<-subset(KOG.all.1, is.na(Description)==FALSE)
write.table(KOG.all.1, paste0(species,".KOG.txt"),
            row.names = F, quote = F, sep = "\t")

##########################################################
####Coprinopsis
##########################################################
species="Coprinopsis_cinerea"

exp<-read.delim(paste0("Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.log2TMMmatrixmean.txt"),header = T)
exp<-exp[,c("Genes","Spore", "Ger", "Brch",
            "Myc", "Knot", "Pri", "YFB")]
exp$Groups<-1:nrow(exp)
exp<-exp[,c(ncol(exp),1:2,2:(ncol(exp)-1))]
##calling with diffgene
exp.diff <-DiffGenes(exp,
                     method = "log-foldchange",
                     p.adjust.method = "BH",
                     nrep = c(2,2,4),
                     stage.names = c("Early", "Mid", "Late"))
exp.diff$Groups[exp.diff$`Early->Mid` < -1 & exp.diff$`Mid->Late` > 1]<-"LHL"
exp.diff$Groups[exp.diff$`Early->Mid` > 1 & exp.diff$`Mid->Late` < -1]<-"HLH"

heatmap.input<-exp[exp$Gene %in% exp.diff$Gene[exp.diff$Groups %in% c("LHL","HLH")], ]
row.names(heatmap.input)<-heatmap.input$Gene
heatmap.input<-heatmap.input[,c(3,5:ncol(heatmap.input))]

heatmap.input<-as.data.frame(t(apply(heatmap.input, 1, normalize)), row.names = row.names(heatmap.input))
heatp<-pheatmap(heatmap.input,
                legend = TRUE,
                scale = "none",
                #scale = "row",
                cluster_rows = TRUE, cluster_cols = FALSE,
                cutree_rows = NA, treeheight_row = 100,
                border_color = NA,
                cellwidth = 30, cellheight = 0.5, fontsize = 24,
                angle_col = c("90"),
                width = ncol(heatmap.input)*0.6+1, height = nrow(heatmap.input)*0.01,
                show_colnames = TRUE,
                show_rownames = F,
                filename = paste0(species,".heatmap.LHL.png"))

## Enrich
GenesKOGpair.1v1<-read.delim(paste0("/home/yichun/RNAmodification/hourglass/eggnog/",
                                    species,".proteins.fa.KOG.1v1.txt"), header = TRUE)
GenesKOGpair.1v1<-unique(GenesKOGpair.1v1[,c("KOG","Genes")])
waist<-exp.diff[exp.diff$Groups %in% c("LHL","HLH"),c("Genes","Groups")]
View(as.data.frame(xtabs(~Genes, waist)))

KOG.all.1<-compareCluster(Genes ~ Groups,
                          data = waist,
                          fun = 'enricher',
                          TERM2GENE = GenesKOGpair.1v1,
                          TERM2NAME = kog2name,
                          pvalueCutoff = 1,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 1,
                          minGSSize = 1,
                          maxGSSize = nrow(GenesKOGpair.1v1))
KOG.all.1<-as.data.frame(KOG.all.1)
KOG.all.1<-subset(KOG.all.1, is.na(Description)==FALSE)
write.table(KOG.all.1, paste0(species,".KOG.txt"),
            row.names = F, quote = F, sep = "\t")

##########################################################
##Plot
##########################################################
kog2name<-read.delim("/home/yichun/3enrichment/kog2name.txt", header = TRUE)
koglist<-rbind(kog2name[,c("kogClass", "kogName")],kog2name[,c("kogClass", "kogName")])
names(koglist)<-c("ID","Description")
koglist$Groups[1:25]<-"HLH"
koglist$Groups[26:50]<-"LHL"
specieslist<-c("Coprinopsis_cinerea",
               "Fusarium_graminearum",
               "Rhizopus_delemar",
               "Drosophila_melanogaster",
               "Danio_rerio",
               "Arabidopsis_thaliana")
KOG.all.1<-as.data.frame(matrix(NA, nrow = 0, ncol = 12))
names(KOG.all.1)<-c("Cluster", "Groups",
                    "ID", "Description",
                    "GeneRatio", "BgRatio",
                    "pvalue", "p.adjust", "qvalue",
                    "geneID", "Count", "species")
for (i in 1:length(specieslist)) {
  species<-specieslist[i]
  KOGtable.sub<-read.delim(paste0(species,".KOG.txt"))
  KOGtable.sub<-merge(KOGtable.sub, koglist, by = c("ID","Description","Groups"), all = T)
  KOGtable.sub$species<-species
  KOG.all.1<-rbind(KOG.all.1,KOGtable.sub)
}
plotin<-as.data.frame(KOG.all.1)
plotin<-subset(plotin, is.na(Description)==FALSE)

#### Plot
kog2name<-read.delim("/home/yichun/3enrichment/kog2name.txt", header = TRUE)
kog2name$ONTOLOGY<-gsub("INFORMATION STORAGE AND PROCESSING", "Information storage\nand processing", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("CELLULAR PROCESSES AND SIGNALING", "Cellular processes\nand signaling", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("POORLY CHARACTERIZED", "Poor", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("METABOLISM", "Metabolism", kog2name$ONTOLOGY)

plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)

colnames(plotinsep)[names(plotinsep)=="ID"]<-"kogClass"

plotdata<-merge(plotinsep, kog2name, by = c("kogClass"), all.x = TRUE)

plotdata$ratio=plotdata$Genenumerator/plotdata$Genedenominator
plotdata$enrich.factor=(plotdata$Genenumerator*plotdata$BGdenominator)/(plotdata$Genedenominator*plotdata$BGnumerator)
plotdata$Description<-paste0(plotdata$Description," [",plotdata$kogClass,"]")
plotdata$species<-gsub("_"," ", plotdata$species)
plotdata$ratio[is.na(plotdata$ratio)]<-0
plotdata$enrich.factor[is.na(plotdata$enrich.factor)]<-1
plotdata<-plotdata[plotdata$kogClass != "R",]

plotdata1<-plotdata[plotdata$Groups=="LHL",]

a<-plotdata1 %>%
  mutate(species=factor(species, levels = c(gsub("_"," ",specieslist)))) %>%
  ggplot(aes(x = log2(enrich.factor), y = Description, size = enrich.factor, colour = p.adjust))+
  labs(title = "", size = "Enrich factor", colour = "P value", x = "Enrich factor", y = "")+
  geom_point(shape = 19)+
  scale_size_area(limits = c(-3,3))+
  scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
  scale_x_continuous(limits = c(-3,3), breaks = seq(-3,3,1))+
  guides(size = guide_legend(order = 1))+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        axis.title.x = element_text(colour = "black", size = 11),
        axis.title.y = element_text(colour = "white"),
        axis.text = element_text(colour = "black", size = 11),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 11),
        plot.title = element_text(size = 11),
        strip.text.x = element_text(size = 11, face = "italic"),
        strip.text.y = element_text(size = 11))+
  facet_grid(ONTOLOGY~species, scales = "free_y", space = "free")

a

ggsave(paste0("KOG.LHL.png"), width = 20, height = 8, units = "in", dpi = 300)
ggsave(paste0("KOG.LHL.tiff"), width = 20, height = 8, units = "in", dpi = 300)

plotdata1<-plotdata[plotdata$Groups=="HLH",]

a<-plotdata1 %>%
  mutate(species=factor(species, levels = c(gsub("_"," ",specieslist)))) %>%
  ggplot(aes(x = log2(enrich.factor), y = Description, size = ratio, colour = p.adjust))+
  labs(title = "", size = "Gene ratio", colour = "Adjust P value", x = "Enrich factor", y = "")+
  geom_point(shape = 19)+
  scale_size_area(limits = c(0.01,0.4))+
  scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
  scale_x_continuous(limits = c(-3.5,3.5), breaks = seq(-3,3,1))+
  guides(size = guide_legend(order = 1))+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        axis.title.x = element_text(colour = "black", size = 11),
        axis.title.y = element_text(colour = "white"),
        axis.text = element_text(colour = "black", size = 11),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 11),
        plot.title = element_text(size = 11),
        strip.text.x = element_text(size = 11, face = "italic"),
        strip.text.y = element_text(size = 11))+
  facet_grid(ONTOLOGY~species, scales = "free_y", space = "free")

a

ggsave(paste0("KOG.HLH.png"), width = 20, height = 8, units = "in", dpi = 300)
ggsave(paste0("KOG.HLH.tiff"), width = 20, height = 8, units = "in", dpi = 300)

plotdata$Groups<-gsub("LHL","Low-High-Low",plotdata$Groups)
plotdata$Groups<-gsub("HLH","High-Low-High",plotdata$Groups)
p<-plotdata %>%
  mutate(species=factor(species, levels = c(gsub("_"," ",specieslist))),
         Groups = factor(Groups, levels = c("Low-High-Low","High-Low-High"))) %>%
  ggplot(aes(x=species, y=Description))+
  geom_tile(aes(fill = log10(enrich.factor)))+
  #scale_fill_gradient2(low = "#91BFDB", mid = "#FFFFBF", high = "#FC8D59",
  scale_fill_gradient2(low = "#4575B4", mid = "#FFFFBF", high = "#D73027",
                       limits = c(-1,1), n.breaks = 3, labels = c("Under-\nrepresent","Even","Over-\nrepresent"))+
  labs(y = "", x = "", fill = "Enrich factor\n\n", title = "")+
  theme(plot.subtitle = element_text(vjust = 1),
        plot.caption = element_text(vjust = 1),
        axis.text.x = element_text(size = 12, colour = "black", angle = 45, face = "italic", hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.line = element_line(size = 0.5, colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.background = element_rect(fill = NA),
        panel.grid.minor = element_line(colour = "gray50", size = 0.5))+
  facet_grid(ONTOLOGY~Groups, scales = "free_y", space = "free")
p
ggsave(paste0("KOG.heatmap.png"), width = 12, height = 8, units = "in", dpi = 300)
ggsave(paste0("KOG.heatmap.tiff"), width = 12, height = 8, units = "in", dpi = 300)
