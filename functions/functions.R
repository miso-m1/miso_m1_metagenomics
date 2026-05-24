# =========================================================
# Global setup
# =========================================================
# Check and install devtools if necessary
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Download and install SpiecEasi from GitHub if necessary
if (!requireNamespace("SpiecEasi", quietly = TRUE)) {
  devtools::install_github("zdk123/SpiecEasi")
}

# Download and install pairwiseAdonis from GitHub if necessary
if (!requireNamespace("pairwiseAdonis", quietly = TRUE)) {
  # Prevent the installation from failing due to minor warnings
  Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = TRUE)
  # The installation must point to the specific sub-directory within the repository
  devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
}

library(phyloseq)
library(ggplot2)
library(viridis)
library(dplyr)
library(tidyr)
library(tibble)
library(readxl)
library(vegan)
library(igraph)
library(SpiecEasi)
library(here)
library(knitr)
library(tidyverse)
library(broom)
library(FSA)
library(broom)
library(pairwiseAdonis)
library(gridExtra)
library(patchwork)


options(lifecycle_verbosity = "quiet")

theme_set(theme_minimal(base_size = 13))

options(
  ggplot2.discrete.fill = function(...)
    scale_fill_viridis_d(option = "cividis", ...),
  
  ggplot2.discrete.colour = function(...)
    scale_color_viridis_d(option = "cividis", ...)
)

# =========================================================
# Alpha diversity functions
# =========================================================

# Compute Shannon diversity from relative abundances (vector of proportions)
shannon_from_rel <- function(x) {
  x <- x[x > 0]
  p <- x / sum(x)          # ensure they sum to 1 (if not already)
  -sum(p * log(p))
}

# Compute Simpson (1 - D) from relative abundances
simpson_from_rel <- function(x) {
  x <- x[x > 0]
  p <- x / sum(x)
  1 - sum(p^2)
}

# Compute Inverse Simpson from relative abundances
invsimpson_from_rel <- function(x) {
  x <- x[x > 0]
  p <- x / sum(x)
  1 / sum(p^2)
}

# Compute Observed richness (number of taxa with abundance > 0)
observed_from_rel <- function(x) {
  sum(x > 0)
}

# Apply to a phyloseq object (assumes OTU table contains relative abundances)
compute_alpha_from_rel <- function(ps) {
  otu <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) {
    otu <- t(otu)   # now rows = samples, columns = taxa
  }
  
  data.frame(
    Sample_ID = sample_names(ps),
    Shannon = apply(otu, 1, shannon_from_rel),
    Simpson = apply(otu, 1, simpson_from_rel),
    InvSimpson = apply(otu, 1, invsimpson_from_rel),
    Observed = apply(otu, 1, observed_from_rel),
    Disease = sample_data(ps)$disease,   # adjust column name as needed
    Gender = sample_data(ps)$gender,
    Country = sample_data(ps)$country
  )
}

normality_tests <- function(alpha_df, metrics) {
  
  bind_rows(lapply(metrics, function(metric) {
    
    test <- shapiro.test(alpha_df[[metric]])
    
    tidy(test) %>%
      mutate(Metric = metric)
    
  }))
}

kruskal_tests <- function(alpha_df, metrics) {
  
  bind_rows(lapply(metrics, function(metric) {
    
    test <- kruskal.test(alpha_df[[metric]] ~ Disease,
                         data = alpha_df)
    
    tidy(test) %>%
      mutate(Metric = metric)
    
  }))
}

pairwise_wilcox_tests <- function(alpha_df, metrics) {
  
  bind_rows(lapply(metrics, function(metric) {
    
    test <- pairwise.wilcox.test(
      alpha_df[[metric]],
      alpha_df$Disease,
      p.adjust.method = "BH"
    )
    
    tidy(test) %>%
      mutate(Metric = metric)
    
  }))
}

# =========================================================
# Network analysis functions
# =========================================================

# Function to build the multi-method consensus network
build_consensus_network <- function(ps_obj, group_name, group_var, tax_level, prevalence_threshold) {
  
  # cat("  -> Inferring sub-network for group:", group_name, "...\n")
  
  # A. Subsampling and phylogenetic taxonomic agglomeration
  sdata <- sample_data(ps_obj)
  ps_sub <- prune_samples(sdata[[group_var]] == group_name, ps_obj)
  ps_tax <- tax_glom(ps_sub, taxrank = tax_level)
  ps_tax <- prune_taxa(taxa_sums(ps_tax) > 0, ps_tax)
  
  # B. Count matrix extraction and orientation
  otu <- as(otu_table(ps_tax), "matrix")
  if(!taxa_are_rows(ps_tax)) otu <- t(otu)
  
  # C. Prevalence filtering to control zero-inflation noise
  prevalence <- rowMeans(otu > 0)
  otu_filt <- otu[prevalence >= prevalence_threshold, ]
  
  if(nrow(otu_filt) < 3) return(NULL)
  
  # Transposition: SpiecEasi and SparCC require samples in rows (n x p)
  otu_input <- t(otu_filt)
  
  # D. Sequential execution of inference algorithms (Silent Mode)
  se_mb <- suppressWarnings(spiec.easi(otu_input, method = 'mb', lambda.min.ratio = 1e-2, nlambda = 20, pulsar.params = list(rep.num = 20)))
  se_gl <- suppressWarnings(spiec.easi(otu_input, method = 'glasso', lambda.min.ratio = 1e-2, nlambda = 20, pulsar.params = list(rep.num = 20)))
  sparcc_res <- sparcc(otu_input)
  
  # E. Adjacency matrix extraction and binarization
  adj_mb <- as.matrix(getRefit(se_mb))
  adj_gl <- as.matrix(getRefit(se_gl))
  adj_sparcc <- ifelse(abs(sparcc_res$Cor) >= 0.4, 1, 0)
  diag(adj_sparcc) <- 0 # Removing self-loops
  
  # F. Computing Majority Vote Consensus (OneNet-Mean)
  adj_mean <- (adj_mb + adj_gl + adj_sparcc) / 3
  adj_consensus <- ifelse(adj_mean >= 0.66, 1, 0)
  
  # Retaining weights only for consensus-validated edges
  adj_weighted <- adj_mean
  adj_weighted[adj_consensus == 0] <- 0
  
  # G. Creating a weighted, undirected igraph object
  net <- graph_from_adjacency_matrix(adj_weighted, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  # Mapping formal taxonomy to node names
  tax_info <- as.data.frame(tax_table(ps_tax))
  V(net)$name  <- tax_info[rownames(otu_filt), tax_level]
  V(net)$label <- V(net)$name
  
  return(net)
}

# Function to compile and summarize topological metrics
print_network_stats <- function(network_list, tax_level) {
  stats_df <- do.call(rbind, lapply(names(network_list), function(g) {
    net <- network_list[[g]]
    if(is.null(net)) return(NULL)
    
    data.frame(
      Tax_Level       = tax_level,
      Disease_Cohort  = g,
      Nodes_Count     = vcount(net),
      Edges_Count     = ecount(net),
      Graph_Density   = round(edge_density(net), 3),
      Average_Degree  = round(mean(degree(net)), 2)
    )
  }))
  return(stats_df)
}

plot_global_network <- function(network_list, tax_level) {
  all_edges <- data.frame()
  
  # Compiling edge lists across all individual sub-networks
  for(g in names(network_list)) {
    net <- network_list[[g]]
    if(!is.null(net) && ecount(net) > 0) {
      edgelist <- igraph::as_data_frame(net, what = "edges")
      edgelist$Disease <- g
      all_edges <- rbind(all_edges, edgelist)
    }
  }
  
  if(nrow(all_edges) == 0) return(NULL)
  
  # Instantiating the master multi-edge graph
  master_net <- graph_from_data_frame(all_edges, directed = FALSE)
  master_net <- delete.vertices(master_net, degree(master_net) == 0) # Pruning isolated nodes
  
  # Mapping disease-specific color codes
  unique_diseases <- unique(E(master_net)$Disease)
  color_palette <- setNames(viridis(length(unique_diseases), option = "cividis"), unique_diseases)
  E(master_net)$color <- color_palette[E(master_net)$Disease]
  
  # Computing geometric curvatures for parallel shared edges
  edge_curvatures <- curve_multiple(master_net)
  
  # Adaptive sizing based on taxonomic rank resolution
  if(tax_level == "Phylum") { v_size <- 8; v_cex <- 0.7 }
  else if(tax_level == "Family") { v_size <- 5; v_cex <- 0.6 }
  else { v_size <- 3; v_cex <- 0.45 }
  
  # Inline HTML Quarto Report Rendering (L'affichage en direct)
  plot(master_net, layout = layout_with_fr(master_net), vertex.size = v_size,
       vertex.color = "grey90", vertex.frame.color = "white",
       vertex.label.color = "black", vertex.label.cex = v_cex,
       edge.color = E(master_net)$color, edge.width = 1.5, edge.curved = edge_curvatures,
       main = paste("Global Consensus Multi-Disease Network -", tax_level))
  
  legend("bottomright", title = "Diseases", legend = names(color_palette), 
         col = color_palette, lwd = 3, bty = "n", cex = 0.8)
  
  return(master_net)
}

# Session info for reproducibility
print_session_info <- function() {
  sessionInfo()
}