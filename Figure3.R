library(monocle)
library(Seurat)
library(tidyverse)
library(patchwork)

load("nw/scRNA/gist cancer/data/Tcell_data.Rdata")
Idents(Tcell_data)<-Tcell_data$scissor

DimPlot(Tcell_data,label = T)
expr_matrix <- as(as.matrix(Tcell_data@assays$RNA@counts), 'sparseMatrix')
p_data <- Tcell_data@meta.data 

f_data <- data.frame(gene_short_name = row.names(Tcell_data),row.names = row.names(Tcell_data))


pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)

cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)


disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
cds <- setOrderingFilter(cds, disp.genes)
plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')

cds <- orderCells(cds)
cds <- orderCells(cds, root_state = 2)

save(cds,file="nw/scRNA/gist cancer/data/cds.Rdata")
load("nw/scRNA/gist cancer/data/cds.Rdata")

################# Figure 3A
pdf("nw/scRNA/gist cancer/2025_11_11/Figure4/train.monocle.pseudotime.pdf",width = 8,height = 6)
plot_cell_trajectory(cds,color_by="Pseudotime", cell_size=0.1,show_backbone=TRUE)+
  scale_color_viridis_c(option="C")
dev.off()

################# Figure 3B
pdf("nw/scRNA/gist cancer/2025_11_11/Figure4/train.monocle.celltype.pdf",width = 8,height = 6)
plot_cell_trajectory(cds,color_by="scissor", cell_size =0.2,show_backbone=TRUE)+
  scale_color_manual(values=c("PD-L1_sen"="#E41A1C", "PD-L1_res"="#4DAF4A", "Non_defined"="grey"))
dev.off()

pdf("nw/scRNA/gist cancer/2025_11_11/Figure4/trajectory.pdf",width = 8,height = 10)
plot_complex_cell_trajectory(cds, x = 1, y = 2,
                             color_by = "scissor")+
  scale_color_manual(values=c("PD-L1_sen"="#E41A1C", "PD-L1_res"="#4DAF4A", "Non_defined"="grey"))+
  theme(legend.title = element_blank()) 
dev.off()


############################ Figure 3C
library(ggpubr)
library(ggsci)
df <- pData(cds) 
View(df)
ClusterName_color_panel <- c(pal_aaas()(10),pal_npg()(10))
p<-ggplot(df, aes(Pseudotime, colour = scissor, fill=scissor)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2()+
  scale_color_manual(values=c("PD-L1_sen"="#E41A1C", "PD-L1_res"="#4DAF4A", "Non_defined"="grey"))+
  scale_fill_manual(values=c("PD-L1_sen"="#E41A1C", "PD-L1_res"="#4DAF4A", "Non_defined"="grey"))
scale_fill_manual(name = "", values = ClusterName_color_panel)+
  scale_color_manual(name = "", values = ClusterName_color_panel)
ggsave("nw/scRNA/gist cancer/2025_11_11/Figure4/density.pdf",p,width = 10,height = 6)

############################## Figure 3D
#
Time_diff <- differentialGeneTest(cds[disp.genes,], cores = 10, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_diff <- Time_diff[,c(5,2,3,4,1,6)] 
write.csv(Time_diff, "nw/scRNA/sun/data/Time_diff_all.csv", row.names = F)
Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
p=plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=6, show_rownames=T, return_heatmap=T)

hclust_obj <- as.hclust(p$tree_row)

gene_clusters <- cutree(hclust_obj, k = 6)
cluster_genes <- split(names(gene_clusters), gene_clusters)
Tmp_gene<-c(cluster_genes[[1]][1:20],cluster_genes[[3]][1:20])

p=plot_pseudotime_heatmap(cds[Tmp_gene,], num_clusters=2, show_rownames=T,
                          return_heatmap=T)

ggsave("nw/scRNA/gist cancer/2025_11_11/Figure4/Time_heatmapAll.pdf", p, 
       width = 8, height = 6)

######################## Figure 3F
library(RColorBrewer)
BEAM_res <- BEAM(cds[c(disp.genes,"CD3D","CD3E","CD3G","CD2","CD247"),],branch_point = 4, cores = 20)


BEAM_res <-BEAM_res[order(BEAM_res$qval),]

BEAM_res <-BEAM_res[,c("gene_short_name", "pval", "qval")]

head(BEAM_res)

write.csv(BEAM_res,"BEAM_res.csv", row.names = F)

pdf("nw/scRNA/gist cancer/12_27/plot/BEAM.pdf",width = 8,height = 15)
plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,qval < 1e-4))[1:100],],branch_point = 1, num_clusters = 4, cores = 10,use_gene_short_name= T,show_rownames = T)
dev.off()

tmp1=plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,qval < 1e-4)),],
                                 branch_point = 1,
                                 num_clusters = 6, 
                                 cores = 1,
                                 #branch_labels = c("Cell fate 1", "Cell fate 2"),
                                 #hmcols = NULL, 
                                 hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62),
                                 branch_colors = c("#979797", "#F05662", "#7990C8"), 
                                 use_gene_short_name = T,
                                 show_rownames = T,
                                 return_heatmap= T 
)


pdf("nw/scRNA/gist cancer/12_27/plot/branched_heatmap_1.pdf",width = 8,height = 35)
tmp1$ph_res
dev.off()

a<-tmp1$annotation_row
a<-a%>%filter(Cluster==c(1,5))

tmp1=plot_genes_branched_heatmap(cds[row.names(a),],
                                 branch_point = 1,
                                 num_clusters = 2, 
                                 cores = 1,
                                 branch_labels = c("Cell fate 1", "Cell fate 2"),
                                 #hmcols = NULL, 
                                 hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62),
                                 branch_colors = c("#979797", "#F05662", "#7990C8"), 
                                 use_gene_short_name = T,
                                 show_rownames = T,
                                 return_heatmap= T 
)
pdf("nw/scRNA/gist cancer/2025_11_11/Figure4/branched_heatmap_1.pdf",width = 8,height = 12)
tmp1$ph_res
dev.off()

a<-tmp1$annotation_row
a<-a%>%filter(Cluster==c(2))
write.csv(rownames(a),file="nw/scRNA/gist cancer/result/fate2.csv",quote = F,row.names = F)

write.csv(a,file = "nw/scRNA/gist cancer/12_27/data/gene_cluster.csv")
