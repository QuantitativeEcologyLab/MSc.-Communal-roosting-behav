###################################
## Phylogenetic Tree Plot 
###################################

# ==== 1. Load Libraries ====
library(ape)          # Phylogenetic tree handling
library(RColorBrewer) # Color palettes
library(tidyverse)    # Data manipulation
library(rphylopic)


# ==== 1. Load Data ====
load("Chapters/Consensus_Tree.Rda")
trait_data <- read_csv("Chapters/Bird_data_clean.csv")

# ==== 2. Filter Only Strigiformes and Prune Tree ====
trait_data <- trait_data %>% filter(Order == "Strigiformes")
phylogeny <- keep.tip(phylogeny, trait_data$Species[trait_data$Species %in% phylogeny$tip.label])


# ==== 3. Map Traits to Tree Edges (for Strigiformes only) ====
tip_states <- setNames(
  ifelse(trait_data$CRB_Final == 1, "Presence", "Absence"),
  trait_data$Species
)

edge_states <- tip_states

# ==== 4. Identify Tip Edges ====
tip_edges <- which(phylogeny$edge[, 2] <= length(phylogeny$tip.label))
tip_indices <- phylogeny$edge[tip_edges, 2]
species_names <- phylogeny$tip.label[tip_indices]
edge_colors <- rep("black", nrow(phylogeny$edge))
edge_colors[tip_edges] <- ifelse(tip_states[species_names] == "Presence", "red", "blue")

# ==== Phylo List of Species =====

# Try to get PhyloPic UUIDs for the species names
trait_data$uuid <- sapply(trait_data$Species, function(x) {
  tryCatch(
    get_uuid(x), 
    error = function(e) "95c59456-77ac-489a-af08-b01001831727"
  )
})


trait_data$svg <- lapply(trait_data$uuid, get_phylopic)

# ==== 5. Plot  ====
par( oma = c(0, 0, 3, 0) ) # Bottom, left, top, right (top is now 6 lines)

# Plot with main title
plot(
  phylogeny,
  type = "phylogram",
  edge.color = edge_colors,
  cex = 0.5,
  label.offset = 0.7,
  no.margin = TRUE,
  show.node.label = FALSE,
)
mtext("Strigiforms and roosting behaviour", side = 3, outer = TRUE, line = 1, cex = 1.2)

# add_phylopic_base(img = trait_data$svg,
#                   x = max(nodeHeights(phylogeny)), y = 1:26, height = 0.5)

add_phylopic_base(img = trait_data$svg,
                   x = max(nodeHeights(phylogeny))-1, y = 1:Ntip(phylogeny), height = 0.5)


legend(
  "topleft",
  legend = c("Absence", "Presence"),
  title = "Communal Roosting Behaviour",
  col = c("blue", "red"),
  lwd = 2,
  bty = "n",
  cex = 0.7
)




