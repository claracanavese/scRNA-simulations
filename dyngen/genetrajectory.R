library(devtools)
require(Seurat)
require(scales)
require(ggplot2)
require(viridis)
require(dplyr)
require(GeneTrajectory)
require(Matrix)
require(plot3D)

data_S <- readRDS("./GeneTrajectory/dataset_seurat_cycle_1000.rds")
data_dyn <- readRDS("./GeneTrajectory/dataset_seurat_dis_100600.rds")

# Assign cell types
data_dyn@meta.data$celltype = data_dyn@misc$traj_progressions$from
DimPlot(data_dyn, reduction = "MDS", group.by = "celltype")

# Seurat standard preprocessing
data_dyn <- NormalizeData(data_dyn)
data_dyn <- FindVariableFeatures(data_dyn)

# Remove burner genes that Dyngen uses to initialise the simulation
all_genes <- data_dyn@assays[["RNA"]]@var.features
genes <- all_genes[!all_genes %in% grep('Burn', all_genes, value=TRUE)]
data_dyn@assays[["RNA"]]@var.features <- genes
# Scale the data
data_dyn <- ScaleData(data_dyn)
# PCA
data_dyn <- RunPCA(data_dyn)
DimPlot(data_dyn, reduction = "pca", group.by = "celltype") + labs(title = element_blank())
ElbowPlot(data_dyn, ndims = 40)
# UMAP
data_dyn <- FindNeighbors(data_dyn, dims = 1:40)
# data_dyn <- FindClusters(data_dyn, resolution = 0.5)
data_dyn <- RunUMAP(data_dyn, dims = 1:40, min.dist = 0.2, spread = 0.3)
# plot 2D embedding
DimPlot(data_dyn, reduction = "umap", group.by = "celltype") + labs(title = element_blank())

# Select genes expressed by 1% to 50% of cells among the top 500 variable genes
all_genes <- data_dyn@assays[["RNA"]]@var.features
expr_percent <- apply(as.matrix(data_dyn[["RNA"]]@data[all_genes, ]) > 0, 1, sum)/ncol(data_dyn)
genes <- all_genes[which(expr_percent > 0.01 & expr_percent < 0.5)]
length(genes)

# Compute the Diffusion Map cell embedding
data_dyn <- GeneTrajectory::RunDM(data_dyn, reduction = "pca", dims = 1:40)
# Calculate cell-cell graph distances over a cell-cell kNN graph
cell.graph.dist <- GetGraphDistance(data_dyn, dims = 1:30, K = 25)
# Coarse-grain the cell graph by grouping cells into `N`=500 "meta-cells"
cg_output <- CoarseGrain(data_dyn, cell.graph.dist, genes, N=900)

# Create a virtualenv using reticulate
if(!reticulate::virtualenv_exists('gene_trajectory')){
  reticulate::virtualenv_create('gene_trajectory', packages=c('gene_trajectory'))
}
reticulate::use_virtualenv('gene_trajectory')
# Import the function to compute gene-gene distances
cal_ot_mat_from_numpy <- reticulate::import('gene_trajectory.compute_gene_distance_cmd')$cal_ot_mat_from_numpy
# Compute gene-gene distances 
# Coarse-grain
gene.dist.mat <- cal_ot_mat_from_numpy(ot_cost = cg_output[["graph.dist"]], gene_expr = cg_output[["gene.expression"]], num_iter_max = 50000, show_progress_bar = TRUE)
rownames(gene.dist.mat) <- cg_output[["features"]]
colnames(gene.dist.mat) <- cg_output[["features"]]
# NO Coarse-grain
gene.dist.mat <- cal_ot_mat_from_numpy(ot_cost = cell.graph.dist, gene_expr = t(as.matrix(data_dyn[["RNA"]]@data[genes,])), num_iter_max = 50000, show_progress_bar = TRUE)
rownames(gene.dist.mat) <- data_dyn[["RNA"]]@var.features
colnames(gene.dist.mat) <- data_dyn[["RNA"]]@var.features
dim(gene.dist.mat)

# construct the gene embedding by employing Diffusion Map
gene_embedding <- GetGeneEmbedding(gene.dist.mat, K = 15)$diffu.emb


# Extract 3 gene trajectories
gene_trajectory <- ExtractGeneTrajectory(gene_embedding, gene.dist.mat, N = 3, t.list = c(4,4,4), K = 5)
table(gene_trajectory$selected)

# Extract the ordered list of genes along each gene trajectory
gene_list <- list()
for (i in 1:3){
  gene_trajectory_sub <- gene_trajectory[which(gene_trajectory$selected == paste0("Trajectory-", i)),]
  genes <- rownames(gene_trajectory_sub)[order(gene_trajectory_sub[, paste0("Pseudoorder", i)])]
  gene_list[[i]] <- genes
}

# Select transcription factors
regs <- data_dyn@assays[["RNA"]]@var.features
regs <- regs[grepl("TF",data_dyn@assays[["RNA"]]@var.features)]

regs2 <- data_dyn@assays[["RNA"]]@var.features
regs2 <- regs2[!grepl("Target",data_dyn@assays[["RNA"]]@var.features)]

regs_left <- regs[grepl("left",regs)]
regs_left_A <- regs_left[grepl("A",regs_left)]
regs_left_B <- regs_left[grepl("B",regs_left)]
regs_left_C <- regs_left[grepl("C",regs_left)]
regs_right <- regs[grepl("right",regs)]
regs_all <- regs[!((grepl("left",regs)) | (grepl("right",regs)))]

genes_left <- which(row.names(gene_trajectory) %in% regs_left)
genes_left_A <- which(row.names(gene_trajectory) %in% regs_left_A)
genes_left_B <- which(row.names(gene_trajectory) %in% regs_left_B)
genes_left_C <- which(row.names(gene_trajectory) %in% regs_left_C)
genes_right <- which(row.names(gene_trajectory) %in% regs_right)
genes_all <- which(row.names(gene_trajectory) %in% regs_all)

regs_A <- regs[grepl("^A",regs)]
regs_B <- regs[grepl("^B",regs)]
regs_C <- regs[grepl("^C",regs)]
regs_D <- regs[grepl("^D",regs)]

genes_A <- which(row.names(gene_trajectory) %in% regs_A)
genes_B <- which(row.names(gene_trajectory) %in% regs_B)
genes_C <- which(row.names(gene_trajectory) %in% regs_C)
genes_D <- which(row.names(gene_trajectory) %in% regs_D)

# Visualize gene trajectories
par(mar = c(2.5,2.5,2.5,2.5))

# scatter3D(gene_embedding[,1],
#           gene_embedding[,2],
#           gene_embedding[,3],
#           bty = "b2", colvar = as.integer(as.factor(gene_trajectory$selected)),
#           main = "trajectory", pch = 19, cex = 1, theta = 30, phi = 45,
#           col = ramp.col(c(hue_pal()(3))))

scatter2D(gene_embedding[,1],
          gene_embedding[,2],
          colvar = as.integer(as.factor(gene_trajectory$selected)),
          main = "Trajectories", pch = 19, cex.main = 1.8, cex.axis = 1.2, cex = 1.3,
          col = c(hue_pal()(2)), colkey = FALSE, type = "p")
legend("topright",
       legend = c("Trajectory 1", "Trajectory 2"),
       col = c(hue_pal()(2)), # Fill color for inverted triangle
       pch = 19,      # Shapes for each item
       pt.cex = 1.2,
       cex = 1)


# Highlight TFs
points(gene_embedding[genes_left, 1], 
       gene_embedding[genes_left, 2],
       col = "black", pch = 15, cex = 1.5)
points(gene_embedding[genes_right, 1], 
       gene_embedding[genes_right, 2],
       col = "black", pch = 17, cex = 1.5)
points(gene_embedding[genes_all, 1], 
       gene_embedding[genes_all, 2],
       col = "black", pch = 25, bg = "black", cex = 1.5)

points(gene_embedding[genes_A, 1], 
       gene_embedding[genes_A, 2],
       col = "black", pch = 15, cex = 1.5)
points(gene_embedding[genes_B, 1], 
       gene_embedding[genes_B, 2],
       col = "black", pch = 17, cex = 1.5)
points(gene_embedding[genes_C, 1], 
       gene_embedding[genes_C, 2],
       col = "black", pch = 25, bg = "black", cex = 1.5)

# Plot legend
# legend("topright",
#        legend = c("Trajectory 1", "Trajectory 2", "Trajectory 3", "Group A TFs", "Group C TFs", "Group D TFs"),
#        col = c(hue_pal()(3), "black", "black", "black"),
#        pt.bg = c(NA, NA, NA, NA, NA, "black"), # Fill color for inverted triangle
#        pch = c(19, 19, 19, 15, 17, 25),        # Shapes for each item
#        pt.cex = 1.5,
#        cex = 1.2)

legend("bottomright",
       legend = c("A TFs", "B TFs", "C TFs"),
       col = c("black", "black", "black"),
       pt.bg = c(NA,NA,"black"), # Fill color for inverted triangle
       pch = c(15,17,25),        # Shapes for each item
       pt.cex = 1.2,
       cex = 1)

legend("topright",
       legend = c("Trajectory 1", "Trajectory 2", "Left TFs", "Right TFs", "Common TFs"),
       col = c(hue_pal()(2), "black", "black", "black"),
       pt.bg = c(NA,NA,NA,NA,"black"), # Fill color for inverted triangle
       pch = c(19,19,15,17,25),        # Shapes for each item
       pt.cex = 1.2,
       cex = 1)

          


# Highlight specific genes
gene_A1 <- which(row.names(gene_trajectory) == "A1-TF1")  # Find the index of "A1-TF1"
gene_A2 <- which(row.names(gene_trajectory) == "A2-TF1")
points(gene_embedding[gene_A2, 1], 
       gene_embedding[highlight_gene, 2],
       col = "black", pch = 19, cex = 1.5)



library(SeuratWrappers)
data_S <- RunALRA(data_S)
data_dyn <- RunALRA(data_dyn)

data_dyn <- AddGeneBinScore(data_dyn, gene_trajectory, N.bin = 4, trajectories = 1:2, assay = "RNA", reverse = c(F))

# Visualize gene bin plots for each gene trajectory
FeaturePlot(data_dyn, pt.size = 0.5, features = paste0("Trajectory",1,"_genes", 1:4), ncol = 2, order = T, reduction = "umap") &
  scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10)))  & NoLegend() & NoAxes() & theme(title = element_text(size = 10))

FeaturePlot(data_dyn, pt.size = 0.5, features = paste0("Trajectory",2,"_genes", 1:4), ncol = 2, order = T, reduction = "umap") &
  scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoLegend() & NoAxes() & theme(title = element_text(size = 10))

FeaturePlot(data_dyn, pt.size = 0.5, features = paste0("Trajectory",3,"_genes", 1:4), ncol = 2, order = T, reduction = "umap") &
  scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoLegend() & NoAxes() & theme(title = element_text(size = 10))

FeaturePlot(data_dyn, pt.size = 0.5, features = paste0("Trajectory",4,"_genes", 1:4), ncol = 2, order = T, reduction = "umap") &
  scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoLegend() & NoAxes() & theme(title = element_text(size = 10))



FeaturePlot(data_dyn, pt.size = 1, features = c("A1-TF1","A3-TF1","A5-TF1"), ncol = 3, order = T, reduction = "umap") &
  scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoLegend() & theme(title = element_text(size = 10))
FeaturePlot(data_dyn, pt.size = 1, features = c("B1-TF1","B3-TF1","B5-TF1"), ncol = 3, order = T, reduction = "umap") &
  scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoLegend() & theme(title = element_text(size = 10))
FeaturePlot(data_dyn, pt.size = 1, features = c("C1-TF1","C3-TF1","C5-TF1"), ncol = 3, order = T, reduction = "umap") &
  scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoLegend() & theme(title = element_text(size = 10))



SaveSeuratRds(data_dyn, file = "./cycle_1000/genetraj_cycle.rds")
data_dyn <- readRDS("./disconnected/genetraj_dis_cycle_linear.rds")
