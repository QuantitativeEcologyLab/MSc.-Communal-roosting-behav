###################################
## Phylogenetic Tree Plot 

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
library(tribble)
library(ragg)
#Color Scheme
primary_colors <- c(
  "#0072B2", # Aqua/Teal
  "#E69F00", # Salmon/Orange
  "#8DA0CB", # Light Blue
  "#E78AC3"  # Pink/Lilac
)

# Lighter Versions
lighter_colors <- sapply(primary_colors, function(col) adjustcolor(col, red.f=1.4, green.f=1.4, blue.f=1.4))
lighter_colors <- as.character(lighter_colors)

# ==== 2. Load Data ====
# Load the consensus tree and the trait data
load("Models/Consensus_Tree.Rda")
phylo <- phylogeny

trait_data <- read_csv("Data/Bird_data_clean.csv")

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
    edge_colors[tip_edges[i]] <- ifelse(tip_states[species_name] == "Presence", primary_colors[1], primary_colors[2])
  }
}

# ==== 5. Plot the Full Tree with Tip Edge Coloring and Save====
png(
  filename = "Figures/Figure 2.2 Evolution of Communal Roosting Behaviour in Core Land Birds.png",
  width    = 10,     # inches
  height   = 8,     # inches
  units    = "in",
  res      = 300    # DPI
)

par(mar = c(6, 5, 5, 5))

plot(
  phylo,
  type = "fan",
  edge.color = edge_colors,
  cex = 0.1,
  label.offset = 0.7,
  #no.margin = TRUE,
  main = "",
  show.node.label = FALSE
)

# ==== 6. Add Legend with Custom Labels ====
legend(
  "topright",
  legend = c("Absence", "Presence"),
  title= "CRB",
  col = c(primary_colors[1], primary_colors[2]),
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

# Method 1 Manually Position the Images
image_df <- tribble(
  ~species,                ~radius, ~angle,
  "Falco_peregrinus",         110,       5,
  "Amazona_aestiva",          110,      35,  #psittaciformes
  "Contopus cooperi",         110,     125,
  "Philemon",                 110,     150,    #changed name
  "Corvus_corax",             110,     170, #passeriformes
  "Aethopyga_siparaja",       110,     195,
  "Troglodytes_hiemalis",     110,     230,
  "Turdus_pilaris",           110,     260,
  "Bubo_bubo",                110,     295,
  "Dryocopus_martius",        110,     305,
  "Merops",                   110,     310,     #changed name
  "Accipiter_gentilis",       110,     320,
  "Cathartes_aura",           110,     350
)

#function to update this dataframe
update_image_df <- function(df,
                            default_uuid = "95c59456-77ac-489a-af08-b01001831727") {
  # sanity check
  if (!all(c("species","radius","angle") %in% names(df))) {
    stop("Data frame must contain columns: species, radius, angle")
  }
  
  # compute radians & coords
  df$angle_rad <- df$angle * pi/180
  df$x         <- df$radius * cos(df$angle_rad)
  df$y         <- df$radius * sin(df$angle_rad)
  
  # lookup UUIDs (fallback on default)
  df$uuid <- sapply(df$species, function(sp) {
    tryCatch(get_uuid(sp), error = function(e) default_uuid)
  })
  
  # fetch SVGs (NULL if missing)
  df$svg <- lapply(df$uuid, function(id) {
    tryCatch(get_phylopic(uuid = id), error = function(e) NULL)
  })
  
  df
}

#call function to update the df
image_df <- update_image_df(image_df)

# STEP 5: add images to the fan chart
add_phylopic_base(
  img    = image_df$svg,
  x      = image_df$x,
  y      = image_df$y,
  height = 8
)


dev.off()


#Version for printing
png(
  filename = "Figures/Figure 2.2 Evolution of Communal Roosting Behaviour.png",
  width    = 30,     # inches
  height   = 24,     # inches
  units    = "in",
  res      = 300    # DPI
)


par(mar = c(6, 5, 5, 5))

plot(
  phylo,
  type = "fan",
  edge.color = edge_colors,
  cex = 0.3,
  label.offset = 0.7,
  #no.margin = TRUE,
  main = "",
  cex.main=4,
  show.node.label = FALSE
)

# ==== 6. Add Legend with Custom Labels ====
legend(
  "topright",
  legend = c("Absence", "Presence"),
  title= "CRB",
  col = c(primary_colors[1], primary_colors[2]),
  lwd = 2,
  bty = "n",
  cex=1.5
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
      font_size <- ifelse(length(species_in_tree) > 15, 1.5, ifelse(length(species_in_tree) < 15, 1.3, 1.4))
      
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

# Method 1 Manually Position the Images
image_df <- tribble(
  ~species,                ~radius, ~angle,
  "Falco_peregrinus",         110,       5,
  "Amazona_aestiva",          110,      35,  #psittaciformes
  "Contopus cooperi",         110,     125,
  "Philemon",                 110,     150,    #changed name
  "Corvus_corax",             110,     170, #passeriformes
  "Aethopyga_siparaja",       110,     195,
  "Troglodytes_hiemalis",     110,     230,
  "Turdus_pilaris",           110,     260,
  "Bubo_bubo",                110,     295,
  "Dryocopus_martius",        110,     305,
  "Merops",                   110,     310,     #changed name
  "Accipiter_gentilis",       110,     320,
  "Cathartes_aura",           110,     350
)

#function to update this dataframe
update_image_df <- function(df,
                            default_uuid = "95c59456-77ac-489a-af08-b01001831727") {
  # sanity check
  if (!all(c("species","radius","angle") %in% names(df))) {
    stop("Data frame must contain columns: species, radius, angle")
  }
  
  # compute radians & coords
  df$angle_rad <- df$angle * pi/180
  df$x         <- df$radius * cos(df$angle_rad)
  df$y         <- df$radius * sin(df$angle_rad)
  
  # lookup UUIDs (fallback on default)
  df$uuid <- sapply(df$species, function(sp) {
    tryCatch(get_uuid(sp), error = function(e) default_uuid)
  })
  
  # fetch SVGs (NULL if missing)
  df$svg <- lapply(df$uuid, function(id) {
    tryCatch(get_phylopic(uuid = id), error = function(e) NULL)
  })
  
  df
}

#call function to update the df
image_df <- update_image_df(image_df)

# STEP 5: add images to the fan chart
add_phylopic_base(
  img    = image_df$svg,
  x      = image_df$x,
  y      = image_df$y,
  height = 8
)


dev.off()


































## Alternate, Method two lookup position based on their order in the fan chart
n_images   <- 10
tip_labels <- phylo$tip.label
N          <- length(tip_labels)
radius     <- 110

# STEP 1: generate a buffer of evenly spaced tip indices
even_idx <- round(seq(1, N, length.out = n_images * 2))

# STEP 2: pick the first n_images with valid UUIDs
selected <- data.frame(species = character(),
                       index   = integer(),
                       uuid    = character(),
                       stringsAsFactors = FALSE)

i <- 1
while(nrow(selected) < n_images && i <= length(even_idx)) {
  idx <- even_idx[i]
  sp  <- tip_labels[idx]
  
  uuid <- tryCatch(get_uuid(sp), error = function(e) NA)
  if(!is.na(uuid)) {
    selected <- rbind(selected,
                      data.frame(species = sp,
                                 index   = idx,
                                 uuid    = uuid,
                                 stringsAsFactors = FALSE))
  }
  i <- i + 1
}

# STEP 3: compute angles and Cartesian coords
selected$angle    <- (selected$index - 1) * (360 / N)
selected$angle_r  <- selected$angle * pi/180
selected$x        <- radius * cos(selected$angle_r)
selected$y        <- radius * sin(selected$angle_r)

# STEP 4: fetch SVGs
selected$svg <- lapply(selected$uuid, function(id) {
  tryCatch(get_phylopic(uuid = id), error = function(e) NULL)
})

# STEP 5: add images to the fan chart
add_phylopic_base(
  img    = selected$svg,
  x      = selected$x,
  y      = selected$y,
  height = 8
)
################################################################################
