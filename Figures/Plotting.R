# Install Packages and Load Libraries
install.packages("ape")
install.packages("ggtree")
install.packages("tidyverse")
install.packages("ggmosaic")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtree")
BiocManager::install("ggtreeExtra")

library(ape)
library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(viridis)
library(stringr)
library(dplyr)
library(ggplot2)
library(scales)
library(ggnewscale)
library(ggmosaic)
library(scales)

# Load data and tree #---------------------------------------------------------#
Bird_data_clean <- read_csv("Chapters/Bird_data_clean.csv")
load("Chapters/Consensus_Tree.Rda") #loads "phylogeny"

# Keep only matching species
tree <- keep.tip(phylogeny, phylogeny$tip.label[phylogeny$tip.label %in% Bird_data_clean$Species])
traits <- Bird_data_clean %>% filter(Species %in% tree$tip.label)

#---------------------------------------------------------#
# Primary Set2 Colors
primary_colors <- c(
  "#66C2A5", # Aqua/Teal
  "#FC8D62", # Salmon/Orange
  "#8DA0CB", # Light Blue
  "#E78AC3"  # Pink/Lilac
)

# Lighter Versions
lighter_colors <- sapply(primary_colors, function(col) adjustcolor(col, red.f=1.4, green.f=1.4, blue.f=1.4))
lighter_colors <- as.character(lighter_colors)

# 4 Shades of Aqua (first primary color)
aqua_shades <- c(
  adjustcolor(primary_colors[1], red.f=1, green.f=1, blue.f=1),
  adjustcolor(primary_colors[1], red.f=1.3, green.f=1.3, blue.f=1.3),
  adjustcolor(primary_colors[1], red.f=1.6, green.f=1.6, blue.f=1.6),
  adjustcolor(primary_colors[1], red.f=1.8, green.f=1.8, blue.f=1.8)
)

# Four Matching Grays
grays <- c("#F2F2F2", "#CCCCCC", "#888888", "#222222")

# Combine all for visualization
all_colors <- c(primary_colors, lighter_colors, aqua_shades, grays)
color_labels <- c(
  paste0("Primary ", c("Aqua","Salmon","Lt Blue","Pink")),
  paste0("Lighter ", c("Aqua","Salmon","Lt Blue","Pink")),
  paste0("Aqua Shade ", 1:4),
  paste0("Gray ", 1:4)
)

# Visualize
par(mar=c(1,10,2,1))
barplot(
  rep(1, length(all_colors)), horiz=TRUE, col=all_colors, border=NA,
  names.arg=color_labels, las=1, main="Set2 Expanded Palette"
)

#---------------------------------------------------------#
# Plot phylo tree with Mass bar charts #
p <- ggtree(tree, layout = "fan") %<+% traits +
  geom_fruit(
    geom = geom_bar,
    mapping = aes(y = label, x = Mass, fill = Order),
    orientation = "y",
    stat = "identity",
    width = 3
  ) +
  geom_tippoint(aes(color = as.factor(CRB_Final)), size = 1.2) +
  scale_color_manual(values = c("0" = primary_colors[1], "1" = primary_colors[2]))

print(p)

#Plot phylo tree with heat map #---------------------------------------------------------#
p <- ggtree(tree, layout = "fan") %<+% traits +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = label, fill = Order),
    stat = "identity",
    width = 1
  ) +
  geom_tree(aes(color = as.factor(CRB_Final)), size = 0.8) +
  geom_tippoint(aes(color = as.factor(CRB_Final)), size = 1.2) +
  scale_color_manual(values = c("0" = primary_colors[1], "1" = primary_colors[2])) +
  scale_fill_viridis_d(option = "plasma") +  # or scale_fill_brewer(palette = "Set3")
  theme(legend.position = "right")

print(p)

#phylo plot trial #---------------------------------------------------------#
p <- ggtree(tree, layout = "rectangular") %<+% traits +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = label, fill = Order),
    stat = "identity",
    width = 1
  ) +
  geom_tippoint(aes(color = as.factor(CRB_Final)), size = 1.2) +
  geom_tree(aes(color = as.factor(CRB_Final)), size = 0.8) +
  scale_color_viridis_d(option = "plasma") +
  scale_fill_brewer(palette = "Set3") +  # optional
  theme(legend.position = "right")

print(p)



#---------------------------------------------------------#

# GENERIC PLOT HWI MASS AND CRB #
# x y chart of HWI vs Mass for CRB
cols <- c(primary_colors[1], primary_colors[2])[as.factor(traits$CRB_Final)]
plot(
  traits$HWI, traits$mass_kg,
  col = cols,
  pch = 20,
  xlab = "HWI",
  ylab = "Mass (kg)",
  main = "CRB prevalence in landbirds"
)

legend(
  "topright",
  legend = c("No CRB", "CRB"),
  col = c(primary_colors[1], primary_colors[2]),
  pch = 20
)


# x and y HWI and mass for CRB
trophic_levels <- unique(traits$Trophic_level)
n_levels <- length(trophic_levels)
base_colors <- c(primary_colors[1], primary_colors[2], primary_colors[3], primary_colors[4])[1:n_levels]
names(base_colors) <- trophic_levels

# Assign color by trophic level
traits$color <- base_colors[traits$Trophic_level]

# Assign pch: 19 = solid, 21 = empty circle with border color
traits$pch <- ifelse(traits$CRB_Final == 1, 19, 21)
plot(
  traits$HWI, traits$mass_kg,
  type = "n",
  xlab = "HWI",
  ylab = "Mass (kg)",
  main = "CRB prevalence in landbirds per trophic level"
)

# Plot empty circles (CRB 0)
with(traits[traits$CRB_Final == 0, ],
     points(HWI, mass_kg, col = color, pch = 21, bg = "white", lwd = 2))

# Plot solid circles (CRB 1)
with(traits[traits$CRB_Final == 1, ],
     points(HWI, mass_kg, col = color, pch = 19))
# Make legend items for each trophic level
legend(
  "topleft",
  legend = paste(rep(trophic_levels, each = 2),
                 rep(c("No CRB", "CRB"), times = n_levels)),
  col = rep(base_colors, each = 2),
  pch = rep(c(21, 19), times = n_levels),
  pt.bg = c(rep("white", n_levels), rep(base_colors, n_levels)),
  pt.lwd = 2,
  cex = 0.6
)

#---------------------------------------------------------#
########### PERCENTAGE CRB PER TROPHIC LEVEL #############
# Calculate percentage of species by Trophic_level and CRB_Final
df_summary <- Bird_data_clean %>%
  group_by(Trophic_level, CRB_Final) %>%
  summarise(count = n(), .groups = 'drop') 

fill_colors <- c(primary_colors, lighter_colors)
names(fill_colors) <- c(
  paste0("primary_colors[", seq_along(primary_colors), "]"),
  paste0("lighter_colors[", seq_along(lighter_colors), "]")
)


df_summary$fill_group <- mapply(function(troph, crb) {
  idx <- as.integer(factor(troph, levels = unique(df_summary$Trophic_level)))
  if (crb == 1) {
    paste0("primary_colors[", idx, "]")
  } else {
    paste0("lighter_colors[", idx, "]")
  }
}, df_summary$Trophic_level, df_summary$CRB_Final)

n <- length(primary_colors)
legend_labels <- c(rep("CRB Absent", n), rep("CRB Present", n))  # or "Lighter"/"Darker"

ggplot(df_summary, aes(x = Trophic_level, y = count, fill = fill_group)) +
  geom_col(position = "fill") +  # fill makes it 100% stacked
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(
    values = fill_colors,
    labels = legend_labels,
    name = "Communal Roosting"   # or "CRB Status"
  ) +
  labs(
    title = "Proportion of CRB by Trophic Level",
    x = "Trophic Level",
    y = "Percentage",
    fill = "CRB_Final") +
  theme(
    plot.title = element_text(hjust = 0.5)  # Centers the title
  )

#---------------------------------------------------------#
# Count of species by Trophic_level and CRB_Final

tab <- table(traits$Trophic_level, traits$CRB_Final)

# Optional: Give better column names
colnames(tab) <- c("No CRB", "CRB")

# Plot stacked bar chart
par(mfrow = c(1,1))         # Single plot layout (default)
par(mar = c(5, 4, 4, 2) + 0.1)  # Default plot margins
trophic_labels <- rownames(tab)
barplot(
  t(tab),
  beside = FALSE,
  legend = TRUE,
  col = c(primary_colors[1], lighter_colors[1]),
  xlab = "Trophic Level",
  ylab = "Number of Species",
  main = "Species with/without CRB by Trophic Level",
  names.arg = trophic_labels
)

# 2 side by side bar charts
par(mfrow = c(1,2))

# --- Left: Stacked bar chart (counts) ---
barplot(
  t(tab),                  # transpose so bars stack correctly
  beside = FALSE,          # stacked
  legend = TRUE,
  col = c(primary_colors[1], lighter_colors[1]),
  xlab = "Trophic Level",
  ylab = "Number of Species",
  main = "Count by Trophic Level"
)
# --- Right: 100% stacked bar chart (proportions) ---
tab_prop <- prop.table(tab, margin = 1)  # proportions by row (trophic level)

barplot(
  t(tab_prop),
  beside = FALSE,
  legend = TRUE,
  col = c(primary_colors[1], lighter_colors[1]),
  xlab = "Trophic Level",
  ylab = "Proportion of Species",
  main = "Proportion by Trophic Level",
  ylim = c(0, 1)
)

# reset to one plot per page
par(mfrow = c(1,1))         # Single plot layout (default)
par(mar = c(5, 4, 4, 2) + 0.1)  # Default plot margins

#---------------------------------------------------------#
########### PERCENTAGE CRB PER HWI #############

## predicting does HWI influence CRB?

# "histogram" style of CRB by handwing
breaks <- seq(0, 60, by = 10)
traits$HWI_bin <- cut(traits$HWI, breaks = breaks, include.lowest = TRUE, right = FALSE)
tab <- table(traits$HWI_bin, traits$CRB_Final)
colnames(tab) <- c("No CRB", "CRB")
barplot(
  t(tab),
  beside = FALSE,   # stacked
  col = c(primary_colors[1], lighter_colors[1]),
  legend = TRUE,
  xlab = "HWI Bin",
  ylab = "Number of Species",
  main = "Number of Species with/without CRB by HWI Bin"
)

# 100% stacked by chart by HWI Bin
breaks <- seq(0, 60, by = 10)
traits$HWI_bin <- cut(traits$HWI, breaks = breaks, include.lowest = TRUE, right = FALSE)
tab <- table(traits$HWI_bin, traits$CRB_Final)
colnames(tab) <- c("No CRB", "CRB")
# Compute proportions by HWI bin (row)
tab_prop <- prop.table(tab, margin = 1)
barplot(
  t(tab_prop),
  beside = FALSE,         # stacked bars
  col = c(primary_colors[1], lighter_colors[1]),
  legend = TRUE,
  xlab = "HWI Bin",
  ylab = "Proportion of Species",
  main = "Proportion by HWI",
  ylim = c(0, 1)
)

#Explore data
top10_mass <- traits[order(-traits$mass_kg), ][1:10, "Species"]
top10_mass$Species  # Print species names (or just use top10_mass)
top10_hwi <- traits[order(-traits$HWI), ][1:10, "Species"]
top10_hwi$Species  # Print species names (or just use top10_hwi)


######################################################################
## Phylogenetic Tree Plot, Subset (Subset of Accipitriforms)

# ==== 1. Load Libraries ====
library(ape)          # Phylogenetic tree handling
library(maps)
library(geiger)
library(ggtree)
library(RColorBrewer) # Color palettes
library(tidyverse)    # Data manipulation
library(png)          # For silhouette images
library(grid)         # Image rasterization
library(rphylopic)    # for loading phylo images
library(tribble)
library(ragg)
#Color Scheme
primary_colors <- c(
  "#66C2A5", # Aqua/Teal
  "#FC8D62", # Salmon/Orange
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
subset_tips <- c("Gyps_indicus", "Harpyopsis_novaeguineae")
#subset_tips <- c("Ninox_ochracea", "Phoeniculus_damarensis") # strigiformes
#subset_tips <- c("Eclectus_roratus", "Alisterus_chloropterus" ) # parrots
#subset_tips <- c("Toxostoma_rufum", "Myadestes_townsendi" ) # parrots

subset_tips %in% phylo$tip.label
mrca_node <- getMRCA(phylo, subset_tips)
sampled_species <- extract.clade(phylo, mrca_node)
#sampled_species <- phylogeny$tip.label

# Prune the phylo tree to the extent of the two species
phylo <- extract.clade(phylo, mrca_node)

phylo$tip.label
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
  filename = "Figures/Figure 2.X Evolution of Communal Roosting Behaviour in Subset of Raptors.png",
  width    = 10,     # inches
  height   = 8,     # inches
  units    = "in",
  res      = 300    # DPI
)

#par(mar = c(6, 5, 5, 5))

plot(
  phylo,
  type = "phylo",
  edge.color = edge_colors,
  cex = 0.5,
  label.offset = 0.7,
  #no.margin = TRUE,
  main = "Evolution of Communal Roosting Behavior in Subset of Raptors",
  show.node.label = FALSE
)

# ==== 6. Add Legend with Custom Labels ====
legend(
  "bottomleft",
  legend = c("Absence", "Presence"),
  title= "CRB",
  col = c(primary_colors[1], primary_colors[2]),
  lwd = 2,
  bty = "n",
  cex=0.7
)

dev.off()


###################
# Plotting a subset with other variables

# ==== 0. Install (if needed) ====
# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ggtreeExtra")

# ==== 1. Libraries ====
library(ape)            # reading & manipulating trees
library(ggtree)         # ggplot2 grammar for trees
library(ggtreeExtra)    # geom_fruit()
library(tidyverse)      # read_csv + dplyr
library(RColorBrewer)   # nice palettes


# ==== 2. Load Data ====
load("Models/Consensus_Tree.Rda")            # gives object `phylogeny`
trait_data <- read_csv("Data/Bird_data_clean.csv")

# ==== 3. Subset to Your Clade ====
subset_tips <- c("Gyps_indicus", "Harpyopsis_novaeguineae")
mrca_node   <- getMRCA(phylogeny, subset_tips)
phylo       <- extract.clade(phylogeny, mrca_node)

# ==== 4. Prepare Tip‐Data ====
tip_df <- trait_data %>%
  filter(Species %in% phylo$tip.label) %>%
  mutate(
    State   = if_else(CRB_Final == 1, "Presence", "Absence"),
    Mass    = Mass               # ensure you have this numeric column
  )

# ==== 5. Build Base Tree Plot ====
base_plot <- ggtree(phylo) %<+% tip_df +
  aes(color = State) +
  geom_tiplab(size = 1, offset = 0.1) +
#  geom_facet(panel = "Mass (g)", data=tip_df, geom = geom_col, 
#             aes(x = trip_df$Mass, color = tip_df$State), orientation = 'y', width = .6) +
  geom_facet(panel = "Mass (g)", geom = geom_col, 
             aes(x = trip_df$Mass, color = tip_df$State), orientation = 'y', width = .6) +
  geom_facet(panel = "HWI", data=tip_df, geom = geom_col, 
             aes(x = tip_df$HWI, color = tip_df$State), orientation = 'y', width = .6) +
labs(
  title = "Raptor Clade: CRB",
  subtitle = "with Mass, HWI"
) +
theme_tree2(legend.position=c(.05, .85))

# ==== 7. Preview ====
base_plot


# ==== 5. Save as SVG ====
svglite::svglite(
  filename = "Figures/Figure_2X_CRB_Raptors.svg",
  width    = 5,    # inches
  height   = 5      # inches
)

print(base_plot)
dev.off()

