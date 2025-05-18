###################################
## Phylogenetic Tree Plot 
###################################

# ==== 1. Load Libraries ====
library(ape)          # Phylogenetic tree handling
library(maps)
library(geiger)
library(phytools)     # Advanced tree plotting with plotTree()
library(RColorBrewer) # Color palettes
library(tidyverse)    # Data manipulation
library(png)          # For silhouette images
library(grid)         # Image rasterization



# ==== 2. Load Data ====
# Load the consensus tree and the trait data
load("Chapters/Consensus_Tree.Rda")      # Loads "phylogeny"
phylo <- phylogeny
trait_data <- read_csv("Chapters/Bird_data_clean.csv")

# Check tree properties
cat("Is the tree ultrametric? ", is.ultrametric(phylo), "\n")
cat("Is the tree rooted? ", is.rooted(phylo), "\n")

# ==== 3. Map Traits to Tree Tips ====
# Initialize tip colors: default to grey
tip_colors <- setNames(rep("grey", length(phylo$tip.label)), phylo$tip.label)

# Find matching species and assign colors
matched <- intersect(trait_data$Species, phylo$tip.label)
tip_traits <- trait_data$CRB_Final[match(matched, trait_data$Species)]
tip_colors[matched] <- ifelse(tip_traits == 1, "red", "black")

# ==== 4. Plot the Tree ====
plot(
  phylo,
  type = "fan",
  tip.color = tip_colors,
  edge.color = tip_colors,
  cex = 0.1

)

# ==== 5. Add Clade Labels for Each Order ====
# Define the list of orders to label
actual_orders <- unique(trait_data$Order)
colors_order <- colorRampPalette(brewer.pal(8, "Dark2"))(length(actual_orders))

# Loop through each order and find its MRCA (Most Recent Common Ancestor)
for (i in seq_along(actual_orders)) {
  order <- actual_orders[i]
  
  # Get species in the order
  species_in_order <- trait_data$Species[trait_data$Order == order]
  species_in_tree <- species_in_order[species_in_order %in% phylo$tip.label]
  
  # Only proceed if there are two or more species in the tree
  if (length(species_in_tree) >= 2) {
    mrca_node <- getMRCA(phylo, species_in_tree)
    
    if (!is.null(mrca_node)) {
      # Add the clade label with an arc
      arc.cladelabels(
        phy = phylo,
        node = mrca_node,
        text = order,
        orientation = "horizontal",
        cex = 0.6,
        col = colors_order[i]
      )
    }
  }
}


# ==== 5. Add Order-Level Labels ====
# Define the list of orders to label
actual_orders <- unique(trait_data$Order) #  "Accipitriformes", "Bucerotiformes", "Cariamiformes", "Cathartiformes","Coliiformes", "Coraciiformes", "Falconiformes", "Passeriformes","Piciformes", "Psittaciformes", "Strigiformes"
# Generate color palette for order labels
colors_order <- colorRampPalette(brewer.pal(8, "Dark2"))(length(actual_orders))

# Plot order labels
for (i in seq_along(actual_orders)) {
  order <- actual_orders[i]
  
  # Get species in the order
  species_in_order <- trait_data$Species[trait_data$Order == order]
  species_in_tree <- species_in_order[species_in_order %in% phylo$tip.label]
  
  if (length(species_in_tree) > 0) {
    if (length(species_in_tree) == 1) {
      # If only one species, place the label at the tip
      tip_num <- which(phylo$tip.label == species_in_tree)
      tiplabels(text = order, tip = tip_num, adj = -0.2, frame = "none", cex = 0.6, col = colors_order[i])
    } else {
      # If multiple species, place the label around the MRCA
      mrca_node <- getMRCA(phylo, species_in_tree)
      if (!is.null(mrca_node)) {
        arc.cladelabels(
          phylo, node = mrca_node, text = order,
          orientation = "horizontal",
          cex = 0.7, col = colors_order[i]
        )
      }
    }
  }
}

# ==== 6. Add Silhouette Images for Each Order ====
image_folder <- "Order Images"
for (i in seq_along(actual_orders)) {
  # Load the image dynamically based on order number
  image_path <- file.path(image_folder, paste0("Order ", i, ".png"))
  
  # Only load if the image exists
  if (file.exists(image_path)) {
    silhouette <- readPNG(image_path)
    
    # Calculate positioning around the fan
    angle <- (i / length(actual_orders)) * 2 * pi
    x_img <- 0.5 + 0.4 * cos(angle)
    y_img <- 0.5 + 0.4 * sin(angle)
    
    # Add image to plot
    grid.raster(silhouette, x = x_img, y = y_img, width = 0.12)
  } else {
    warning(paste("Image not found for:", image_path))
  }
}

# ==== 7. Final Touches ====
title(main = "Phylogenetic Tree with Order Labels and Traits", cex.main = 1.2)
cat("Tree plotted successfully!\n")


# to test fan plots
#library(phytools)
#test_tree <- rtree(20)
#plotTree(test_tree, type = "fan")


# this was the old way of plotting
#plot_tree <- plot.phylo(
#  phylo, type = "fan", cex = 0.1, 
#  tip.color = tip_colors, label.offset = 0.5,
#  edge.color = edge_colors,
# #show.tip.label = FALSE,  # hide tip labels to focus on branch colors
#  no.margin = TRUE
#)
