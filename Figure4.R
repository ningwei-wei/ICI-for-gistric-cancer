library(Seurat)
library(dorothea)
library(tidyverse)
library(viper)

load("nw/scRNA/gist cancer/data/Tcell_data.Rdata")

## We read Dorothea Regulons for Human:
dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))

## We obtain the regulons based on interactions with confidence level A, B and C
regulon <- dorothea_regulon_human %>% dplyr::filter(confidence %in% c("A","B","C"))

pbmc <- run_viper(Tcell_data, regulon,options = list(method = "scale", minsize = 4, 
                                                     eset.filter = FALSE, cores = 10, verbose = FALSE))
Assays(pbmc)

save(pbmc,file = "nw/scRNA/gist cancer/data/pbmc.Rdata")
load("nw/scRNA/gist cancer/data/pbmc.Rdata")

DotPlot(pbmc,features = c("IRF1","SOX9","ERG"))

DefaultAssay(pbmc) <- "dorothea"
Idents(pbmc)<-pbmc$scissor
table(Idents(pbmc))
pbmc<-subset(pbmc,subset= scissor%in%c("PD-L1_res","PD-L1_sen"))

pbmc.markers <- FindAllMarkers(object = pbmc,logfc.threshold = 0.3)   
pbmc <- ScaleData(pbmc)

DoHeatmap(pbmc,rownames(pbmc.markers),size=3,slot = 'scale.data')

############################
viper_scores_df <- GetAssayData(pbmc, slot = "scale.data", assay = "dorothea") %>%
  data.frame(check.names = F) %>% t()
viper_scores_df[1:4,1:4]

## We create a data frame containing the cells and their clusters
CellsClusters <- data.frame(cell = names(Idents(pbmc)), 
                            cell_type = as.character(Idents(pbmc)),
                            check.names = F)
head(CellsClusters)

## We create a data frame with the Viper score per cell and its clusters
viper_scores_clusters <- viper_scores_df  %>%
  data.frame() %>% 
  rownames_to_column("cell") %>%
  gather(tf, activity, -cell) %>%
  inner_join(CellsClusters)


## We summarize the Viper scores by cellpopulation
summarized_viper_scores <- viper_scores_clusters %>% 
  group_by(tf, cell_type) %>%
  summarise(avg = mean(activity),
            std = sd(activity))

# For visualization purposes, we select the 20 most variable TFs across clusters according to our scores.
head(summarized_viper_scores)


## We select the 20 most variable TFs. (20*9 populations = 180)
highly_variable_tfs <- summarized_viper_scores %>%
  group_by(tf) %>%
  mutate(var = var(avg))  %>%
  ungroup() %>%
  top_n(40, var) %>%
  distinct(tf)
highly_variable_tfs
write.csv(highly_variable_tfs,file = "nw/scRNA/gist cancer/data/TF_gene.csv",quote = F,row.names = F)

highly_variable_tfs[21,1]<-"NFKB1"

## We prepare the data for the plot
summarized_viper_scores_df <- summarized_viper_scores %>%
  semi_join(highly_variable_tfs, by = "tf") %>%
  dplyr::select(-std) %>%   
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE) 

palette_length = 100
my_color = colorRampPalette(c("#1e466e", "#376795","#528fad","#72bcd5","#aadce0",
                              "#ffe6b7","#ffd06f","#f7aa58","#ef8a47","#e76254"))(palette_length)


data_min <- min(summarized_viper_scores_df)
data_max <- max(summarized_viper_scores_df)


neg_breaks <- seq(data_min, 0, length.out = ceiling(palette_length/2) + 1)

pos_breaks <- seq(0, data_max, length.out = ceiling(palette_length/2) + 1)[-1] # 去掉第一个0

my_breaks <- c(neg_breaks, pos_breaks)


library(pheatmap)
viper_hmap <- pheatmap(t(summarized_viper_scores_df),fontsize=14, 
                       fontsize_row = 10, 
                       color=my_color, breaks = my_breaks, 
                       main = "DoRothEA (ABC)", angle_col = 45,
                       treeheight_col = 0,  border_color = NA) 

####################### Figure 4A
pdf("nw/scRNA/gist cancer/2025_11_11/Figure5/heatmap_new.pdf",width = 8,height = 6)
print(viper_hmap)
dev.off()


a<-colnames(summarized_viper_scores_df)
write.csv(a,file = "nw/scRNA/gist cancer/data/gene.csv",quote = F,row.names = F)

####################################### Figure 4C,4D

library(scTenifoldKnk)
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggrepel)  
library(pbapply)
options(mc.cores = 50)

load("nw/scRNA/gist cancer/data/Tcell_data.Rdata")
highVarGenes <- FindVariableFeatures(Tcell_data,nfeatures = 20000)

countMatrix <- GetAssayData(Tcell_data, slot = "counts")
countMatrix <- countMatrix[highVarGenes@assays[["RNA"]]@var.features, ]

result <- scTenifoldKnk(countMatrix = countMatrix, 
                        gKO = 'IRF1', 
                        qc_mtThreshold = 0.1,
                        qc_minLSize = 200,
                        nc_nNet = 10, 
                        nc_nCells = 500, 
                        nc_nComp = 3 )
top_genes <- head(result$diffRegulation[order(-result$diffRegulation$FC), ], 21)
top_genes<-top_genes[-1,]

a<-result$diffRegulation

p<-ggplot(top_genes, aes(x=reorder(gene, FC), y=FC)) +
  geom_bar(stat='identity', fill='steelblue') +
  coord_flip() +
  labs(title="Top 20 Differentially Regulated Genes",
       x="Gene", y="FC") +
  theme_minimal()
ggsave("nw/scRNA/gist cancer/2025_11_11/Figure5/nfkb1.pdf",p,width = 8,height = 6)


df <- result$diffRegulation

df$log_pval <- -log10(df$p.adj)

label_genes <- subset(df, abs(Z) > 2 & p.adj < 0.01)

ggplot(df, aes(x=Z, y=log_pval)) +
  geom_point(alpha=0.5) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="red") +
  geom_vline(xintercept=c(2), linetype="dashed", color="blue") +
  geom_text_repel(data=label_genes, aes(label=gene),
                  size=3, max.overlaps=50) +
  labs(title="Z vs -log10(p-value)",
       x="Z-score", y="-log10(p-value)") +
  theme_classic()

############################## Figure 4E
library(scMetabolism)
library(tidyverse)
library(rsvd)
library(Seurat)
library(pheatmap)
library(ComplexHeatmap)
library(ggsci)
library(data.table)
library(stringr)

load("nw/scRNA/gist cancer/data/Tcell_data.Rdata")
DotPlot(Tcell_data,features = c("IDO1","ARG1"))

countexp.Seurat <- sc.metabolism.Seurat(obj = Tcell_data,  
                                        method = "AUCell", 
                                        imputation = F, 
                                        ncores = 10,
                                        metabolism.type = c("KEGG","REACTOME")) ## KEGG; REACTOME


score <- countexp.Seurat@assays$METABOLISM$score
score[1:4,1:4]

colnames(score)<-rownames(countexp.Seurat@meta.data)

identical(colnames(score), rownames(countexp.Seurat@meta.data))

countexp.Seurat@meta.data <- cbind(countexp.Seurat@meta.data,t(score) )


p1 <- FeaturePlot(countexp.Seurat,features = "Glycolysis / Gluconeogenesis",order = T)
p2 <- VlnPlot(countexp.Seurat,features = "Glycolysis / Gluconeogenesis")
p1 + p2


df <- countexp.Seurat@meta.data

avg_df = aggregate(df[,96:ncol(df)],
                   list(df$scissor),
                   mean)


avg_df <- avg_df %>% 
  #dplyr::select(1:20) %>% 
  column_to_rownames("Group.1") 

which(colnames(avg_df)=="Glutathione metabolism")

avg_df<-avg_df[,c(1,2,3,6,7,10,20,46,42,37,28,53,43,83,73)]

h_state <- pheatmap(t(avg_df),
                    scale = "row",
                    column_title = "",
                    col = colorRampPalette(c("#1e466e", "#376795","#528fad","#72bcd5","#aadce0",
                                             "#ffe6b7","#ffd06f","#f7aa58","#ef8a47","#e76254"))(100),
                    angle_col = "45")

pdf("nw/scRNA/gist cancer/2025_11_11/Figure5/Meta_heatmap.pdf",width = 8,height = 6)
print(h_state)
dev.off()














