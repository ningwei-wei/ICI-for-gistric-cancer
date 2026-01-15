library(CellChat)
library(Seurat)

load("nw/scRNA/gist cancer/data/OMIX001073.Rdata")
cellchat <- createCellChat(object=OMIX001073,meta = OMIX001073@meta.data,
                           group.by = "scissor_celltype")

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)

cellchat@DB <- CellChatDB

# This step is necessary even if using the whole database
cellchat <- subsetData(cellchat) 
# do parallel
future::plan("multicore", workers = 40) 

####  Identify overexpressed genes
cellchat <- identifyOverExpressedGenes(cellchat)
### Identify overexpressed ligand-receptor pairs
cellchat <- identifyOverExpressedInteractions(cellchat)

#project gene expression data onto PPI 
#(Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cellchat <- projectData(cellchat, PPI.human)
cellchat@data.project[1:4,1:4]

options(future.globals.maxSize = 2 * 1024^3)
cellchat <- computeCommunProb(cellchat, raw.use = F, population.size = TRUE) 
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "nw/scRNA/sun/data/cell-cell_communications.all.csv")

### Calculate the communication results of all ligand-receptor interactions related to each signaling pathway
cellchat <- computeCommunProbPathway(cellchat)
### Calculate the communication results between integrated cell types
cellchat <- aggregateNet(cellchat)

save(cellchat,file = "nw/scRNA/gist cancer/data/cellchat.Rdata")
load("nw/scRNA/gist cancer/data/cellchat.Rdata")

groupSize <- as.numeric(table(cellchat@idents))

##################################  Figure 2A
pdf("nw/scRNA/gist cancer/2025_11_11/Figure3/interaction_weight.pdf",width = 8,height = 6)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 remove.isolate = F,arrow.size = 0.1,arrow.width = 0.1,
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

############### Figure 2B
mat <- cellchat@net$weight
pdf("nw/scRNA/gist cancer/2025_11_11/Figure3/celltype.pdf",width = 8,height = 6)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 10:11) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames
                 = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, 
                   arrow.size = 0.1,
                   weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()
##################
cellchat@netP$pathways
levels(cellchat@idents)  

pathways.show <- c("MHC-I")  
vertex.receiver = c(5,10,11) 

################ Figure 2C
pdf("nw/scRNA/gist cancer/2025_11_11/Figure3/MHC_hier.pdf",width = 10,height = 6)
netVisual_aggregate(cellchat, signaling = pathways.show, arrow.size=0.1,
                    vertex.weight = 0.1,
                    vertex.receiver = vertex.receiver,layout = "hierarchy")
dev.off()
# Circle plot
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

############# Figure 2D
pathways.show <- c("LAMININ")  
pdf("nw/scRNA/gist cancer/2025_11_11/Figure3/LAMININ_heatmap.pdf",width = 10,height = 6)
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
dev.off()


cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
pdf("nw/scRNA/gist cancer/12_27/plot/out.pdf",width = 15,height = 30)
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",width = 15,height = 30)
dev.off()

pdf("nw/scRNA/gist cancer/12_27/plot/incom.pdf",width = 15,height = 30)
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",width = 15,height = 30)
dev.off()


# Chord diagram
par(mfrow=c(1,1))
pdf(file ="nw/scRNA/sun/plot/cellchat/pathway.pdf", width = 18, height =16)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
dev.off()

############################## Figure 2E
immune_checkpoint_pathways <- c("MHC-I", "PVR", "CD80", "CD86", "LIGHT", "BTLA","CXCL")

pdf("nw/scRNA/gist cancer/2025_11_11/Figure3/bulle.pdf",width = 8,height = 6)
netVisual_bubble(cellchat, sources.use = c(5), targets.use = c(10,11), angle.x = 45,
                 signaling = immune_checkpoint_pathways, remove.isolate = T)
dev.off()

netVisual_bubble(cellchat, sources.use = c(5), 
                 targets.use = c(10,11), remove.isolate = T)


############################## Figure 2H
pdf("nw/scRNA/gist cancer/2025_11_11/Figure3/gene.pdf",width = 8,height = 6)
plotGeneExpression(cellchat, signaling = "MHC-I", enriched.only = TRUE, type = "violin"
)
dev.off()


