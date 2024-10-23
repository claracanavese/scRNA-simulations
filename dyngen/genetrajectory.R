library(devtools)
require(Seurat)
require(scales)
require(ggplot2)
require(viridis)
require(dplyr)
require(GeneTrajectory)
require(Matrix)
require(plot3D)

data_S <- readRDS("human_myeloid_seurat_obj.rds")
data_dyn <- readRDS("./GeneTrajectory/dataset_seurat_cycle_1000.rds")

# Assign cell types
data_dyn@meta.data$celltype = data_dyn@misc$traj_progressions$from
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
DimPlot(data_dyn, reduction = "pca", group.by = "celltype")
ElbowPlot(data_dyn)
# UMAP
data_dyn <- FindNeighbors(data_dyn, dims = 1:10)
# data_dyn <- FindClusters(data_dyn, resolution = 0.5)
data_dyn <- RunUMAP(data_dyn, dims = 1:10)
# plot 2D embedding
DimPlot(data_dyn, reduction = "umap", group.by = "celltype")

# Select genes expressed by 1% to 50% of cells among the top 500 variable genes
all_genes <- data_dyn@assays[["RNA"]]@var.features
expr_percent <- apply(as.matrix(data_dyn[["RNA"]]@data[all_genes, ]) > 0, 1, sum)/ncol(data_dyn)
genes <- all_genes[which(expr_percent > 0.01 & expr_percent < 0.5)]
length(genes)

# Compute the Diffusion Map cell embedding
data_dyn <- GeneTrajectory::RunDM(data_dyn, reduction = "pca", dims = 1:10)
# Calculate cell-cell graph distances over a cell-cell kNN graph
cell.graph.dist <- GetGraphDistance(data_dyn, dims = 1:10, K = 15)
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
gene.dist.mat <- cal_ot_mat_from_numpy(ot_cost = cg_output[["graph.dist"]], gene_expr = cg_output[["gene.expression"]], num_iter_max = 50000, show_progress_bar = TRUE)
gene.dist.mat <- cal_ot_mat_from_numpy(ot_cost = cell.graph.dist, gene_expr = t(as.matrix(data_dyn[["RNA"]]@data[all_genes,])), num_iter_max = 50000, show_progress_bar = TRUE)
rownames(gene.dist.mat) <- cg_output[["features"]]
colnames(gene.dist.mat) <- cg_output[["features"]]
rownames(gene.dist.mat) <- data_dyn[["RNA"]]@var.features
colnames(gene.dist.mat) <- data_dyn[["RNA"]]@var.features
dim(gene.dist.mat)

# construct the gene embedding by employing Diffusion Map
gene_embedding <- GetGeneEmbedding(gene.dist.mat, K = 15)$diffu.emb

# Extract 3 gene trajectories
gene_trajectory <- ExtractGeneTrajectory(gene_embedding, gene.dist.mat, N = 1, t.list = c(80), K = 15)
table(gene_trajectory$selected)

# Visualize gene trajectories
par(mar = c(1.5,1.5,1.5,1.5))

scatter3D(gene_embedding[,1],
          gene_embedding[,2],
          gene_embedding[,3],
          bty = "b2", colvar = as.integer(as.factor(gene_trajectory$selected))-1,
          main = "trajectory", pch = 19, cex = 1, theta = 45, phi = 0,
          col = ramp.col(c(hue_pal()(3))))

scatter2D(gene_embedding[,1],
          gene_embedding[,2],
          colvar = as.integer(as.factor(gene_trajectory$selected))-1,
          main = "trajectory", pch = 19, cex = 1,
          col = ramp.col(c(hue_pal()(3))))



library(SeuratWrappers)
data_dyn <- RunALRA(data_dyn)


data_dyn <- AddGeneBinScore(data_dyn, gene_trajectory, N.bin = 5, trajectories = 1:2, assay = "alra", reverse = c(F, F, T))

# Visualize gene bin plots for each gene trajectory
FeaturePlot(data_dyn, pt.size = 0.05, features = paste0("Trajectory",1,"_genes", 1:5), ncol = 5, order = T) &
  scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoLegend() & NoAxes() & theme(title = element_text(size = 10))

FeaturePlot(data_dyn, pt.size = 0.05, features = paste0("Trajectory",2,"_genes", 1:5), ncol = 5, order = T) &
  scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoLegend() & NoAxes() & theme(title = element_text(size = 10))
