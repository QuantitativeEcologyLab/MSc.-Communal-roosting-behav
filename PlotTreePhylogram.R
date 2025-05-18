###################################
## Phylogenetic Tree Plot - Phylogram Layout
###################################

# ==== 1. Load Libraries ====
library(ape)          # Phylogenetic tree handling
library(maps)         # For geographic mapping (not used here but included in your original script)
library(geiger)       # Phylogenetic analysis
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

# ==== 3. Reduce the Trait Data to 10% ====
set.seed(42)  # For reproducibility
sampled_species <- sample(trait_data$Species, size = ceiling(0.1 * nrow(trait_data)))

# Only keep species that exist in the tree
sampled_species <- sampled_species[sampled_species %in% phylo$tip.label]

# Prune the trait data as well
trait_data_small <- trait_data %>% filter(Species %in% sampled_species)

# ==== 4. Prune the Phylogenetic Tree ====
phylo_small <- keep.tip(phylo, sampled_species)

phylo <- phylo_small
trait_data <- trait_data_small

# Check the new size
cat("Original tree size:", length(phylo$tip.label), "\n")
cat("Reduced tree size:", length(phylo_small$tip.label), "\n")
# Check tree properties
cat("Is the tree ultrametric? ", is.ultrametric(phylo), "\n")
cat("Is the tree rooted? ", is.rooted(phylo), "\n")

# ==== 3. Map Traits to Tree Tips ====
# Initialize tip colors: default to grey
tip_colors <- setNames(rep("grey", length(phylo$tip.label)), phylo$tip.label)

# Find matching species and assign colors
matched <- intersect(trait_data$Species, phylo$tip.label)
tip_traits <- trait_data$CRB_Final[match(matched, trait_data$Species)]
tip_colors[matched] <- ifelse(tip_traits == 1, "red", "blue")

# edge colors
edge_colors <- rep("grey", nrow(phylo$edge))
tip_indices <- match(phylo$tip.label, trait_data$Species)
species_colors <- ifelse(trait_data$CRB_Final[tip_indices] == 1, "red", "blue")
edge_colors[match(1:length(phylo$tip.label), phylo$edge[, 2])] <- species_colors

# ==== 4. Plot the Tree - Phylogram Layout ====
plot(
  phylo,
  type = "phylogram",
  #ftype = "off",       # Turn off default labels
  colors = species_colors,
  #color = tip_colors,
  lwd = 1.2,
  fsize = 0.6
)



# ==== 5. Add Order-Level Labels ====
# Define the list of orders to label
actual_orders <- unique(trait_data$Order)
colors_order <- colorRampPalette(brewer.pal(8, "Dark2"))(length(actual_orders))

# Get the vertical positions for each order
y_positions <- seq(0.9, 0.1, length.out = length(actual_orders))

# Plot each order label at the right
for (i in seq_along(actual_orders)) {
  order <- actual_orders[i]
  
  # Plot the label at the right of the tree
  grid.text(
    label = order,
    x = 1.15,            # 15% further to the right of the tree
    y = y_positions[i],   # Evenly spaced vertically
    just = "left",        # Left-aligned text
    gp = gpar(col = colors_order[i], fontsize = 10)
  )
}

# # ==== 6. Add Silhouette Images for Each Order ====
# image_folder <- "Chapters/Order Images"
# 
# # Adjust the image placement for phylogram
# y_positions <- seq(0.1, 0.9, length.out = length(actual_orders))
# 
# for (i in seq_along(actual_orders)) {
#   # Load the image dynamically based on order number
#   image_path <- file.path(image_folder, paste0("Order ", i, ".png"))
# 
#   # Only load if the image exists
#   if (file.exists(image_path)) {
#     silhouette <- readPNG(image_path)[,,1]  # Only take the grayscale layer
# 
#     # Calculate positioning vertically along the side
#     x_img <- 1.05   # Slightly off the tree for visibility
#     y_img <- y_positions[i]  # Evenly spaced down the plot
# 
#     # Add image to plot (use grayscale as black)
#     grid.raster(silhouette, x = x_img, y = y_img, width = 0.1, interpolate = FALSE)
#   } else {
#     warning(paste("Image not found for:", image_path))
#   }
# }


# ==== 7. Final Touches ====
title(main = "Phylogenetic Tree (Phylogram) with Order Labels and Traits", cex.main = 1.2)
cat("Tree plotted successfully!\n")
