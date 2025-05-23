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
library(rphylopic)    # for loading phylo images



# ==== 2. Load Data ====
# Load the consensus tree and the trait data
load("C:/Users/Sandra/OneDrive - UBC/PhD proposal/Chapter 2/phylo_data/Consensus_Tree.Rda")
phylo <- phylogeny

trait_data <- read_csv("C:/Users/Sandra/OneDrive - UBC/PhD proposal/Chapter 2/Updated database/Bird_data_clean.csv")


####THIS IS THE ONE#####
# ==== 1. Load Libraries ====
library(ape)
library(phytools)
library(RColorBrewer)

# ==== 2. Prepare the Data ====
# Use all species in the phylogeny
sampled_species <- phylogeny$tip.label

# Prune the phylo tree to match the full species list
phylo <- keep.tip(phylogeny, sampled_species)

# ==== 3. Map Traits to Tree Edges ====
# Find matching species and their edges
matched <- intersect(trait_data$Species, phylo$tip.label)
tip_traits <- trait_data$CRB_Final[match(matched, trait_data$Species)]

# Create a named vector for states (0 or 1)
tip_states <- ifelse(tip_traits == 1, "Presence", "Absence")
names(tip_states) <- matched

# Initialize all edges as "Absence"
edge_states <- setNames(rep("Absence", length(phylo$tip.label)), phylo$tip.label)

# Apply the trait states to the edges
edge_states[names(tip_states)] <- tip_states

# ==== 4. Identify Tip Edges ====
tip_edges <- which(phylo$edge[, 2] <= length(phylo$tip.label))
edge_colors <- rep("black", nrow(phylo$edge))

for (i in seq_along(tip_edges)) {
  tip_index <- phylo$edge[tip_edges[i], 2]
  species_name <- phylo$tip.label[tip_index]
  
  if (species_name %in% names(tip_states)) {
    edge_colors[tip_edges[i]] <- ifelse(tip_states[species_name] == "Presence", "red", "blue")
  }
}

# ==== 5. Plot the Full Tree with Tip Edge Coloring ====
plot(
  phylo,
  type = "fan",
  edge.color = edge_colors,
  cex = 0.1,
  label.offset = 0.7,
  no.margin = TRUE,
  show.node.label = FALSE
)

# ==== 6. Add Legend with Custom Labels ====
legend(
  "topright",
  legend = c("0 = Absence", "1 = Presence"),
  title= "Communal Roosting Behaviour",
  col = c("blue", "red"),
  lwd = 2,
  bty = "n",
  cex=0.7
)

# ==== 7. Plot All Clade Labels with Abbreviation for <30 and Hide <5 ====
actual_orders <- unique(trait_data$Order)
colors_order <- colorRampPalette(brewer.pal(8, "Dark2"))(length(actual_orders))
clade_positions <- data.frame()

# Function to abbreviate clade names
abbreviate_label <- function(label) {
  words <- unlist(strsplit(label, "_| "))
  if (length(words) > 1) {
    paste0(substr(words, 1, 3), collapse = ".")
  } else {
    substr(label, 1, 4)
  }
}

# Loop through each order and find its MRCA (Most Recent Common Ancestor)
for (i in seq_along(actual_orders)) {
  order <- actual_orders[i]
  species_in_order <- trait_data$Species[trait_data$Order == order]
  species_in_tree <- species_in_order[species_in_order %in% phylo$tip.label]
  
  # Proceed only if there are at least 5 species
  if (length(species_in_tree) >= 10) {
    mrca_node <- getMRCA(phylo, species_in_tree)
    
    if (!is.null(mrca_node)) {
      # Abbreviate the name if the clade has fewer than 30 species
      clade_label <- if (length(species_in_tree) < 25) {
        abbreviate_label(order)
      } else {
        order
      }
      
      # Set font size based on the clade size
      font_size <- ifelse(length(species_in_tree) > 15, 0.7, ifelse(length(species_in_tree) < 15, 0.4, 0.5))
      
      # Add the clade label
      arc.cladelabels(
        phy = phylo,
        node = mrca_node,
        text = clade_label,
        cex = font_size,
        col = colors_order[i],
        lwd = 2,
        ln.offset = 1.1,
        lab.offset = 1.2,
        mark.node=FALSE
      )
      
      # Store for reference
      clade_positions <- rbind(clade_positions, data.frame(
        label = order,
        color = colors_order[i],
        node = mrca_node
      ))
    }
  }
}

# ==== adding the phylo images ====
unique_orders <- unique(trait_data$Order)

# Function to compute x and y coordinates based on radius and angle
calculate_xy_positions <- function(df) {
  if (!all(c("radius", "angle") %in% colnames(df))) {
    stop("DataFrame must contain 'radius' and 'angle' columns.")
  }
  
  df$angle <- df$angle * (pi / 180)
  # Calculate x and y based on radius and angle
  df$x <- df$radius * cos(df$angle)
  df$y <- df$radius * sin(df$angle)
  
  # Return the updated dataframe
  return(df)
}

# First entry creates the dataframe
image_df <- data.frame(
  species = c("Corvus_corax", 
              "Acanthagenys_rufogularis",
              "Accipiter_gentilis",
              "Amazona_albifrons",
              "Bubo_bubo",
              "Cathartes_aura",
              "Contopus_sordidulus",
              "Dryocopus_martius",
              "Falco_peregrinus",
              "Merops_apiaster",
              "Ramphastos_sulfuratus",
              "Troglodytes_hiemalis"),
  radius = c(100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100),
  angle = c(190, 160, 320,75, 290, 345,
            140, 300, 10,310, 280, 200)  # in degrees
)
# calculate x and y based on angle and radius
image_df <- calculate_xy_positions(image_df)

# Try to get PhyloPic UUIDs for the species names
image_df$uuid <- sapply(image_df$species, function(x) {
  tryCatch(
    get_uuid(x), 
    error = function(e) "95c59456-77ac-489a-af08-b01001831727"
  )
})

#find images for each uuid
image_df$svg <- lapply(image_df$uuid, get_phylopic)

# apply the images
add_phylopic_base(img = image_df$svg, x = image_df$x, y = image_df$y, height = 4)


# ==== 10. Final Touches ====
title(main = "Phylogenetic Tree with Order Labels and Traits", cex.main = 1.2)
cat("Tree plotted successfully!\n")


