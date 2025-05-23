# Install the package if not already installed
install.packages("ape")
install.packages("ggtree")
install.packages("tidyverse")
install.packages("ggmosaic")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtree")
BiocManager::install("ggtreeExtra")

# Load required libraries
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

# Load data and tree
Bird_data_clean <- read_csv("Chapters/Bird_data_clean.csv")
load("Chapters/Consensus_Tree.Rda") #loads "phylogeny"

# Keep only matching species
tree <- keep.tip(phylogeny, phylogeny$tip.label[phylogeny$tip.label %in% Bird_data_clean$Species])
traits <- Bird_data_clean %>% filter(Species %in% tree$tip.label)

# # Base tree + trait data
# p <- ggtree(tree, layout = "fan") %<+% traits
# 
# # Color tips by CRB_Final
# p1 <- p +
#   geom_tippoint(aes(color = CRB_Final), size = 1.2) +
#   scale_color_viridis(option = "inferno")
# 
# # Show the plot, most plain version
# print(p1)

################PHYLO TREE CHARTS TRIALS###############################
#Plot with Order Labels to outside on an angle
# Example base plot:
# p <- ggtree(tree, layout = "fan") %<+% traits + 
#   geom_tippoint(aes(color = CRB_Final), size = 1.2) +
#   scale_color_viridis(option = "inferno")
# 
# # Join traits data with the tree tip data
# tip_data <- data.frame(label = tree$tip.label) %>%
#   left_join(traits, by = c("label" = "Species"))
# 
# # Extract the ggtree plot data for tips
# tip_positions <- p$data %>% 
#   filter(isTip) %>%
#   left_join(tip_data, by = c("label")) 
# 
# # Each tip has an angle, we want the mean angle per Order
# order_positions <- tip_positions %>%
#   group_by(Order.x) %>%
#   summarise(angle = mean(angle), # average angle for label placement
#             x = mean(x), # average x (radius)
#             y = mean(y)) %>%
#   ungroup()
# 
# order_positions <- order_positions %>%
#   mutate(
#     # Adjust label angle for readability
#     angle_adj = ifelse(angle > 90 & angle < 270, angle + 180, angle),
#     
#     # Adjust horizontal justification
#     hjust = ifelse(angle > 90 & angle < 270, 1, 0)
#   )
# 
# 
# 
# # Multiply x and y by a factor to push labels outside the tree radius
# label_radius_factor <- 1.05
# p + 
#   geom_text(data = order_positions,
#             aes(x = x * label_radius_factor, 
#                 y = y * label_radius_factor, 
#                 label = Order.x, 
#                 angle = angle_adj,
#                 hjust = hjust),
#             size = 3)

######################
#Plot phylo tree with Mass bar charts
p <- ggtree(tree, layout = "fan") %<+% traits +
  geom_fruit(
    geom = geom_bar,
    mapping = aes(y = label, x = Mass, fill = Order),
    orientation = "y",
    stat = "identity",
    width = 3
  ) +
  geom_tippoint(aes(color = as.factor(CRB_Final)), size = 1.2) +
  scale_color_manual(values = c("0" = "grey60", "1" = "firebrick"))

print(p)


######################
#Plot phylo tree with heat map
p <- ggtree(tree, layout = "fan") %<+% traits +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = label, fill = Order),
    stat = "identity",
    width = 1
  ) +
  geom_tree(aes(color = as.factor(CRB_Final)), size = 0.8) +
  geom_tippoint(aes(color = as.factor(CRB_Final)), size = 1.2) +
  scale_color_manual(values = c("0" = "grey60", "1" = "firebrick")) +
  scale_fill_viridis_d(option = "plasma") +  # or scale_fill_brewer(palette = "Set3")
  theme(legend.position = "right")

print(p)


######################
#phylo plot trial
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



########### GENERIC PLOT HWI MASS AND CRB #############
# x y chart of HWI vs Mass for CRB
cols <- c("blue", "red")[as.factor(trait_data$CRB_Final)]
plot(
  trait_data$HWI, trait_data$mass_kg,
  col = cols,
  pch = 20,
  xlab = "HWI",
  ylab = "Mass (kg)",
  main = "CRB prevalence in landbirds"
)

legend(
  "topright",
  legend = c("No CRB", "CRB"),
  col = c("blue", "red"),
  pch = 20
)


# x and y HWI and mass for CRB
trophic_levels <- unique(trait_data$Trophic_level)
n_levels <- length(trophic_levels)
base_colors <- c("red", "green", "blue", "orange")[1:n_levels]
names(base_colors) <- trophic_levels

# Assign color by trophic level
trait_data$color <- base_colors[trait_data$Trophic_level]

# Assign pch: 19 = solid, 21 = empty circle with border color
trait_data$pch <- ifelse(trait_data$CRB_Final == 1, 19, 21)
plot(
  trait_data$HWI, trait_data$mass_kg,
  type = "n",
  xlab = "HWI",
  ylab = "Mass (kg)",
  main = "CRB prevalence in landbirds per trophic level"
)

# Plot empty circles (CRB 0)
with(trait_data[trait_data$CRB_Final == 0, ],
     points(HWI, mass_kg, col = color, pch = 21, bg = "white", lwd = 2))

# Plot solid circles (CRB 1)
with(trait_data[trait_data$CRB_Final == 1, ],
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


########### PERCENTAGE CRB PER TROPHIC LEVEL #############
# Option 1 Calculate percentage of species by Trophic_level and CRB_Final
df_summary <- Bird_data_clean %>%
  group_by(Trophic_level, CRB_Final) %>%
  summarise(count = n(), .groups = 'drop') 

#
ggplot(df_summary, aes(x = Trophic_level, y = count, fill = factor(CRB_Final))) +
  geom_col(position = "fill") +  # fill makes it 100% stacked
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "Proportion of CRB by Trophic Level",
    x = "Trophic Level",
    y = "Percentage",
    fill = "CRB_Final"
  ) +
  theme_minimal()


# Option 2 Calculate percentage of species by Trophic_level and CRB_Final
#Load data
trait_data <- read_csv("Chapters/Bird_data_clean.csv")

## CRB by species and trophic 
# Make a table of counts
tab <- table(trait_data$Trophic_level, trait_data$CRB_Final)

# Optional: Give better column names
colnames(tab) <- c("No CRB", "CRB")

# Plot stacked bar chart
barplot(
  t(tab),
  beside = FALSE,
  legend = TRUE,
  col = c("blue", "red"),
  xlab = "Trophic Level",
  ylab = "Number of Species",
  main = "Species with/without CRB by Trophic Level"
)

# 2 side by side bar charts
par(mfrow = c(1,2))

# --- Left: Stacked bar chart (counts) ---
barplot(
  t(tab),                  # transpose so bars stack correctly
  beside = FALSE,          # stacked
  legend = TRUE,
  col = c("blue", "red"),
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
  col = c("blue", "red"),
  xlab = "Trophic Level",
  ylab = "Proportion of Species",
  main = "Proportion by Trophic Level",
  ylim = c(0, 1)
)

# Optional: reset to one plot per page
par(mfrow = c(1,1))





########### PERCENTAGE CRB PER HWI #############

## predicting does HWI influence CRB?

# "histogram" style of CRB by handwing
breaks <- seq(0, 60, by = 10)
trait_data$HWI_bin <- cut(trait_data$HWI, breaks = breaks, include.lowest = TRUE, right = FALSE)
tab <- table(trait_data$HWI_bin, trait_data$CRB_Final)
colnames(tab) <- c("No CRB", "CRB")
barplot(
  t(tab),
  beside = FALSE,   # stacked
  col = c("blue", "red"),
  legend = TRUE,
  xlab = "HWI Bin",
  ylab = "Number of Species",
  main = "Number of Species with/without CRB by HWI Bin"
)

# 100% stacked by chart by HWI Bin
breaks <- seq(0, 60, by = 10)
trait_data$HWI_bin <- cut(trait_data$HWI, breaks = breaks, include.lowest = TRUE, right = FALSE)
tab <- table(trait_data$HWI_bin, trait_data$CRB_Final)
colnames(tab) <- c("No CRB", "CRB")
# Compute proportions by HWI bin (row)
tab_prop <- prop.table(tab, margin = 1)
barplot(
  t(tab_prop),
  beside = FALSE,         # stacked bars
  col = c("blue", "red"),
  legend = TRUE,
  xlab = "HWI Bin",
  ylab = "Proportion of Species",
  main = "Proportion by HWI",
  ylim = c(0, 1)
)

#Explore data
top10_mass <- trait_data[order(-trait_data$mass_kg), ][1:10, "Species"]
top10_mass$Species  # Print species names (or just use top10_mass)
top10_hwi <- trait_data[order(-trait_data$HWI), ][1:10, "Species"]
top10_hwi$Species  # Print species names (or just use top10_hwi)







