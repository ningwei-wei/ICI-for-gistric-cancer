####### load R package
library(Scissor)
library(Seurat)
library(pheatmap)
library(ComplexHeatmap)
library(ggsci)
library(data.table)
library(stringr)
library(morandipal)

colors <- c(get_palette("complete_17"),get_palette("gradient_full"))
###############
############################# OMIX001073

OMIX001073<- Read10X(data.dir = "nw/scRNA/gist cancer/data/OMIX001073/")
OMIX001073 <- CreateSeuratObject(counts = OMIX001073, project = "OMIX001073",min.cells = 10, min.features = 200)

meta_data<-read.table("nw/scRNA/gist cancer/data/OMIX001073/OMIX001073-20-04.tsv",sep = "\t",header = T,fill = T)

OMIX001073 <- NormalizeData(OMIX001073)
OMIX001073 <- FindVariableFeatures(OMIX001073, selection.method = "vst", nfeatures = 2000)
OMIX001073 <- ScaleData(OMIX001073, features = rownames(OMIX001073))
OMIX001073 <- RunPCA(OMIX001073, features = VariableFeatures(OMIX001073),reduction.name = "pca")
ElbowPlot(OMIX001073,ndims = 50,reduction = "pca")
OMIX001073 <- RunUMAP(OMIX001073,reduction = "pca", dims = 1:40)
OMIX001073 <- FindNeighbors(OMIX001073, reduction = "pca", dims = 1:40)
OMIX001073 <- FindClusters(OMIX001073, resolution = 0.1,random.seed = 123)

OMIX001073$celltype<-meta_data$cluster
OMIX001073$patient<-meta_data$patient
OMIX001073$tiss<-meta_data$tissue

OMIX001073<-subset(OMIX001073,subset = tiss=="Tumor")

################ Figure 1A

DimPlot(OMIX001073, reduction = 'umap', label = T,cols = colors)

scissor_maker<-FindAllMarkers(OMIX001073)

################################ Figure 1B
library(scRNAtoolVis)
library(tidyverse)
library(RColorBrewer)
colour_palette <- brewer.pal(12, "Set3")

sce_maker<-scissor_maker
sce_maker<-sce_maker%>%filter(p_val_adj<0.05)
data<-sce_maker

jjVolcano(diffData = data,
          log2FC.cutoff = 1, 
          pSize  = 3, 
          fontface = 'italic', 
          aesCol = c('purple','orange'), 
          #tile.col = colour, 
          #col.type = "adjustP", 
          topGeneN = 5, 
          tile.col = colour_palette
)

############################### Figure 1C

Tcell_data<-subset(OMIX001073,subset = celltype=="T cells & NK cells")

Tcell_data <- NormalizeData(Tcell_data)
Tcell_data <- FindVariableFeatures(Tcell_data, selection.method = "vst", nfeatures = 2000)
Tcell_data <- ScaleData(Tcell_data, features = rownames(Tcell_data))
Tcell_data <- RunPCA(Tcell_data, features = VariableFeatures(Tcell_data),reduction.name = "pca")
ElbowPlot(Tcell_data,ndims = 50,reduction = "pca")
Tcell_data <- RunUMAP(Tcell_data,reduction = "pca", dims = 1:40)
Tcell_data <- FindNeighbors(Tcell_data, reduction = "pca", dims = 1:40)
Tcell_data <- FindClusters(Tcell_data, resolution = 0.2,random.seed = 123)

library(RColorBrewer)
colour_palette <- brewer.pal(12, "Set3")

color_use<-c("#e69501","#60a3e1","#017552","#f7da48","#01549d","#b45300","#bd7a9d",
             "#646464","#956402","#0977d1",colour_palette[-2])
DimPlot(Tcell_data,cols = color_use)

p<-DimPlot(Tcell_data, reduction = 'umap',  order = T,cols =color_use,label = T,label.box = T)+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) 

############################# Figure 1D

bulk_data<-readRDS("nw/scRNA/gist cancer/data/STAD-PRJEB25780.Response.Rds")
bulk_data<-as.data.frame(bulk_data)
rownames(bulk_data)<-bulk_data$GENE_SYMBOL
bulk_data<-bulk_data[,-1]
bulk_data<-log2(bulk_data+1)

bulk_pheno<-read.table("nw/scRNA/gist cancer/data/STAD-PRJEB25780.Response.tsv",sep = "\t")
bulk_pheno<-bulk_pheno[,c(1,5,7)]
bulk_pheno<-bulk_pheno[-which(bulk_pheno$Treatment=="Normal"),]

bulk_data<-bulk_data[,bulk_pheno$sample_id]
bulk_data<-as.matrix(bulk_data)

############################
library(Scissor)
tag <- c('N', 'R')
pheno<-ifelse(bulk_pheno$response_NR=="R",1,0)

infos1 <- Scissor(bulk_data, Tcell_data,pheno,tag = tag, alpha = 0.05,
                  family = "binomial", Save_file = 'nw/scRNA/gist cancer/data/Scissor_STAD_11_17.RData')
#############################
Scissor_select <- rep("Non_defined", ncol(Tcell_data))
names(Scissor_select) <- colnames(Tcell_data)
Scissor_select[infos1$Scissor_pos] <- "PD-L1_sen"
Scissor_select[infos1$Scissor_neg] <- "PD-L1_res"

Tcell_data <- AddMetaData(Tcell_data, metadata = Scissor_select, col.name = "scissor")

plot_data <- Embeddings(Tcell_data, "umap") %>% 
  as.data.frame() %>%
  cbind(Tcell_data@meta.data$scissor) 
colnames(plot_data) <- c("UMAP1", "UMAP2", "Cluster")

p<-ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = Cluster, alpha = Cluster)) +
  geom_point(size = 0.8) +
  scale_color_manual(values = c('grey','royalblue','indianred1')) +
  scale_alpha_manual(values = c(0.1,0.5,0.5))+
  theme_classic()

############################# Figure 1E
library(AUCell)
library(GSEABase)
library(tidyverse)
genes<-c("CD274","PDCD1","PIK3CD","PIK3R1","PIK3R2","PTPN6","PTPN11","PDCD1LG2","PIK3R3")
genes1<-c("AKT3","GNB5","PIK3R6","AKT1","AKT2","PIK3R5","MTOR","GNB1","GNB2","GNB3","GNG3","GNG4","GNG5","GNG7","GNG10","GNG11","GNGT1","GNGT2","CXCL8","CXCR2","GNG13","PIK3CG","GNG2","GNG12","GNB4","GNG8")
genes2<-c("SMAD2","SMAD3","SMAD4","TGFB1","TGFB2","TGFB3","TGFBR1","TGFBR2")

# 创建关联表
df1 <- data.frame(pathway = "PD-L1/PD1", gene = genes)
df2 <- data.frame(pathway = "PI3K/AKT/mTOR", gene = genes1)
df3 <- data.frame(pathway = "TGF-β/Smad", gene = genes2)

# 合并
data_combined <- rbind(df1, df2, df3)

geneSets <- c(GeneSet(genes, setName="PD1"),GeneSet(genes1, setName="AKT"),
              GeneSet(genes, setName="TGFβ"))
geneSets<-GeneSetCollection(geneSets)


exprMatrix<-Tcell_data@assays$RNA@data

cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=T,nCores = 50)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)

Tcell_data$PD1<-getAUC(cells_AUC)["PD1",]
Tcell_data$AKT<-getAUC(cells_AUC)["AKT",]
Tcell_data$TGFβ<-getAUC(cells_AUC)["TGFβ",]

FeaturePlot(Tcell_data,features = "PD1",cols = colorRampPalette(c("white","#af434a", "red"))(100),order = T)
DotPlot(Tcell_data,features = c("PD1","AKT","TGFβ"),group.by = "scissor")

p1 <- DotPlot(object = Tcell_data,features = c("PD1","AKT","TGFβ"), group.by="scissor", scale.by ="size", scale=5, col.max = 1.3) +
  RotatedAxis() +
  scale_color_gradientn(values = seq(0,1,0.1),colors = c("#5bacdb","#62b8e3","#c5e1f0","white","#f0b4c9","#ea87aa","#ea7ea3"))

################################ Figure 1F

library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)

GDSC2_Expr = readRDS('D:\\DataFiles\\DataFiles\\DataFiles\\Training Data/GDSC2_Expr (RMA Normalized and Log Transformed).rds')
GDSC2_Res = readRDS("D:\\DataFiles\\DataFiles\\DataFiles\\Training Data/GDSC2_Res.rds")
GDSC2_Res <- exp(GDSC2_Res) 
colnames(GDSC2_Res)<-str_split_fixed(colnames(GDSC2_Res),"_",2)[,1]

a<-read.csv("g:/gastric cancer/Drug_listTue Nov 18 03_01_14 2025.csv")
a<-a%>%filter(Target.pathway =="PI3K/MTOR signaling",Datasets=="GDSC2")

GDSC2_Res<-GDSC2_Res[,intersect(a$Name,colnames(GDSC2_Res))]

testexpr<-read.csv("g:/gastric cancer/tmpdata.csv")
rownames(testexpr)<-testexpr$Cell_barcode
testexpr<-testexpr[,-c(24568:24569)]
a<-rownames(testexpr)
testexpr<-t(testexpr)

setwd("g:/gastric cancer/")
calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = testexpr,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )

#####################################
library(ggplot2)
library(ggpubr)
library(ggprism)

data<-read.csv("g:/gastric cancer/calcPhenotype_Output/DrugPredictions.csv")
colnames(data)[1]<-"Cell_barcode"
testexpr<-read.csv("g:/gastric cancer/tmpdata.csv")
a<-testexpr[,c(24568:24569)]

data<-inner_join(a,data)

a<-colnames(data)[-c(1:2)]

p1<-ggplot(data,aes(Condition,Rapamycin,fill = Condition))+geom_boxplot()+stat_compare_means()+
  scale_fill_manual(values = c("#5bacdb","#ea7ea3"))+
  theme_prism()

############################## Figure 1G
Idents(Tcell_data)<-Tcell_data$scissor
deg_res<-FindMarkers(Tcell_data,ident.1 = "PD-L1_res",logfc.threshold = 0,min.pct = 0)

deg_res$difference <- deg_res$pct.1 - deg_res$pct.2
deg_res_sig <- deg_res[which(deg_res$p_val_adj<0.05 & abs(deg_res$avg_log2FC) >0.3),]
deg_res_sig$label <- rownames(deg_res_sig)


library(ggplot2)
library(ggrepel)
library(ggsci)

p<-ggplot(deg_res, aes(x=difference, y=avg_log2FC)) + 
  geom_point(size=2, color="grey60") + 
  geom_text_repel(data=deg_res_sig, aes(label=label),
                  color="black",fontface="italic")+
  geom_point(data=deg_res[which(deg_res$p_val_adj<0.05 & deg_res$avg_log2FC>0.3),],
             aes(x=difference, y=avg_log2FC),
             size=2, color="#DC0000FF")+
  geom_point(data=deg_res[which(deg_res$p_val_adj<0.05 & deg_res$avg_log2FC< -0.1),],
             aes(x=difference, y=avg_log2FC),
             size=2, color= "#3C5488FF")+
  theme_classic()+
  theme(axis.text.x = element_text(colour = 'black',size = 12),
        axis.text.y = element_text(colour = 'black',size = 12),
        axis.title = element_text(colour = 'black',size = 15),
        axis.line = element_line(color = 'black', size = 1))+
  geom_hline(yintercept = 0,lty=2,lwd = 1)+
  geom_vline(xintercept = 0,lty=2,lwd = 1)+
  ylab("Log-fold Change")+
  xlab("Delta Percent")

################################ Figure 1H
###################################
library(ggplot2)
library(tidyverse)
library(ggsci)

color<-c(pal_atlassian()(8))
color<-color[c(1,4)]

data<-read.csv("d:/gastric cancer/result/gProfiler_hsapiens_2025-11-21 15-40-21__intersections.csv")
data<-data[c(2,3,7),]
data$source<-"Up"

data1<-read.csv("d:/gastric cancer/result/gProfiler_hsapiens_2025-11-21 15-43-02__intersections——down.csv")
data1<-data1[c(1,6,13,16,17,50,96),]
data1$source<-"Down"

data<-rbind(data,data1)

data<-data%>%group_by(source)%>%arrange(intersection_size)
data$term_name<-factor(data$term_name,levels = data$term_name)

p1 <- ggplot(data,
             aes(x=term_name,y=intersection_size, fill=source)) + 
  geom_bar(stat="identity", width=0.8) +  
  scale_fill_manual(values = color[1:length(unique(data$source))] ) +  
  coord_flip() +  #
  xlab("term name") +  
  ylab("Gene Number") +  
  labs(title = "Terms Enrich")+  
  theme_bw()+facet_grid(source~., scale = 'free_y', space = 'free_y')+theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = "g:/gastric cancer/plot/Figure1/kegg.pdf",p1,width = 10,height = 6)

                               







