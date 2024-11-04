library(tidyverse)
library(dyngen)
library(dynwrap)
library(dyno)
library(tidygraph)

set.seed(1)

# Define backbone and other parameters
names(list_backbones())
backbone <- backbone_disconnected(left_backbone = backbone_cycle(), 
                                  right_backbone = backbone_linear(),
                                  num_common_modules = 0)
backbone$expression_patterns 
backbone$module_network
backbone$module_info
write.csv(backbone$expression_patterns,"dis_1000/expression_patterns_disconnected_cycle_linear.csv")

config <- initialise_model(
  num_tfs = nrow(backbone$module_info),
  num_targets = 100,
  num_hks = 0,
  backbone = backbone,
  num_cells = 300,
  simulation_params = simulation_default(
    burn_time = simtime_from_backbone(backbone, burn = TRUE) * 1.5,
    compute_rna_velocity = TRUE,
    store_reaction_propensities = TRUE,
    compute_cellwise_grn = TRUE,
    experiment_params = simulation_type_wild_type(
    num_simulations = 100
    )
  ),
  verbose = FALSE
)

backbone$module_info

plot_backbone_statenet(config)
plot_backbone_modulenet(config)

# Generate transcription factors (TFs)
model <- generate_tf_network(config)
plot_feature_network(model, show_targets = FALSE)

# Sample target genes and housekeeping genes (HKs)
model <- generate_feature_network(model)
plot_feature_network(model)

plot_feature_network(model, show_hks = TRUE)

#  Generate kinetics
model <- generate_kinetics(model)
plot_feature_network(model)

plot_feature_network(model, show_hks = TRUE)

# Simulate gold standard
model <- generate_gold_standard(model)
plot_gold_simulations(model) + scale_colour_brewer(palette = "Dark2")
# expression of the modules (average of TFs)
plot_gold_expression(model, what = "mol_mrna") # mrna
plot_gold_expression(model, label_changing = FALSE) # premrna, mrna, and protein

# Simulate cells
model <- generate_cells(model)
plot_simulations(model)
# gold standard and simulations overlayed
plot_gold_simulations(model) + scale_colour_brewer(palette = "Dark2")
# simulations mapped to gold standard
plot_gold_mappings(model, do_facet = FALSE) + scale_colour_brewer(palette = "Dark2")
# expression of the modules (average of TFs)
plot_simulation_expression(model_cycle, 1:4, what = "mol_mrna")

# Experiment emulation
model <- generate_experiment(model)

# Convert to a dyno object
dataset <- as_dyno(model)
plot_dimred(dataset)

plot_graph(dataset)


# alternative: Convert to an anndata/SCE/Seurat object
library(anndata)
ad <- as_anndata(model)
ad$write_h5ad("dyngen/dataset_grn.h5ad")

library(SingleCellExperiment)
sce <- as_sce(model)
write_rds(sce, "dataset_sce.rds")

library(Seurat)
seurat <- as_seurat(model)
write_rds(seurat, "GeneTrajectory/dataset_seurat_dis_cycle_linear_1000.rds")

# One-shot function
out <- generate_dataset(
  config,
  format = "dyno",
  make_plots = TRUE
)

dataset <- out$dataset
model <- out$model
print(out$plot)
# Save
write_rds(model, "dis_1000/model_dis_1000.rds", compress = "gz")
write_rds(dataset, "dis_1000/dataset_dis_1000.rds", compress = "gz")

out$plot[[9]]

out$plot


# Load
data_cycle <- readRDS("dis_1000/dataset_dis_1000.rds")
model_cycle <- readRDS("dis_1000/model_dis_1000.rds")



# GRN

feature_info <- dataset$feature_info
feature_network <-
  dataset$regulatory_network %>%
  arrange(as.character(regulator) == as.character(target))

# add extra edges invisible between regulators from the same module
feature_network_tmp <-
  bind_rows(
    feature_network,
    feature_info %>%
      filter(is_tf) %>%
      select(module_id, feature_id) %>%
      group_by(module_id) %>%
      do({
        crossing(regulator = .$feature_id, target = .$feature_id) %>%
          mutate(effect = -2)
      }) %>%
      ungroup() %>%
      filter(regulator < target)
  )

gr <- tbl_graph(nodes = feature_info, edges = feature_network_tmp)
gr_df <- as.data.frame(gr)
tf_df <- gr_df %>% filter(is_tf == TRUE)


layout <- igraph::layout_with_fr(gr) %>%
  dynutils::scale_minmax() %>%
  magrittr::set_rownames(feature_info$feature_id) %>%
  magrittr::set_colnames(c("x", "y")) %>%
  as.data.frame()

cell_wps <- do.call(dynwrap::select_waypoint_cells, c(dataset[c("milestone_ids", "milestone_network", "milestone_percentages", "progressions", "divergence_regions")], list(num_cells_selected = 1)))
cells <- c("cell223","cell494","cell263","cell1")

node_df <- bind_cols(feature_info %>% select(-mol_premrna:-mol_protein), layout)

cap <- .015
capend <- .02

edge_df <-
  bind_rows(
    feature_network %>% mutate(group = "Global GRN"),
    dataset$regulatory_network_sc %>%
      filter(cell_id %in% cell_wps) %>%
      arrange(as.character(regulator) == as.character(target)) %>%
      left_join(feature_network %>% select(regulator, target, effect), by = c("regulator", "target")) %>%
      mutate(effect = strength, group = paste0("GRN of cell ", match(cell_id, cell_wps)))
  ) %>%
  mutate(
    regulator = as.character(regulator),
    target = as.character(target),
    x_ = layout[regulator, "x"],
    y_ = layout[regulator, "y"],
    xend_ = layout[target, "x"],
    yend_ = layout[target, "y"],
    len = sqrt( (x_ - xend_)^2 + (y_ - yend_)^2 ),
    al = (len - cap) / len,
    alend = (len - capend) / len,
    x = al * x_ + (1 - al) * xend_,
    y = al * y_ + (1 - al) * yend_,
    xend = alend * xend_ + (1 - alend) * x_,
    yend = alend * yend_ + (1 - alend) * y_
  )

arrow_up <- grid::arrow(type = "closed", angle = 30, length = grid::unit(1.3, "mm"))
arrow_down <- grid::arrow(type = "closed", angle = 89, length = grid::unit(1.3, "mm"))
g <- ggplot() +
  geom_point(aes(x, y), colour = "gray", node_df) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend, colour = effect), edge_df %>% filter(effect >= 0), arrow = arrow_up) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend, colour = effect), edge_df %>% filter(effect < 0), arrow = arrow_down) +
  #theme_graph(base_family = "Helvetica") +
  coord_equal() +
  labs(colour = "Regulatory\nactivity") +
  scale_colour_gradientn(colours = c(rev(RColorBrewer::brewer.pal(5, "Reds")), RColorBrewer::brewer.pal(5, "Greens")), limits = c(-1, 1)) +
  facet_wrap(~group) +
  theme(
    plot.margin = margin(0, 0, 0, 0, "cm")
  )
g
