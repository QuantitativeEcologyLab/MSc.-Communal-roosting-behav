###################################
## Phylogenetic Tree Plot 
###################################

# ==== 1. Load Libraries ====
library(ape)          # Phylogenetic tree handling
library(RColorBrewer) # Color palettes
library(tidyverse)    # Data manipulation


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
mtext("Phylogenetic Tree with Order Labels and Traits", side = 3, outer = TRUE, line = 1, cex = 1.2)

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
  col = c("grey", "blue"),
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
  col = c("grey", "blue"),
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
  col = c("grey", "blue"),
  xlab = "Trophic Level",
  ylab = "Proportion of Species",
  main = "Proportion by Trophic Level",
  ylim = c(0, 1)
)

# Optional: reset to one plot per page
par(mfrow = c(1,1))


# x y chart of HWI vs Mass for CRB
cols <- c("grey", "blue")[as.factor(trait_data$CRB_Final)]
plot(
  trait_data$HWI, trait_data$mass_kg,
  col = cols,
  pch = 19,
  xlab = "HWI",
  ylab = "Mass (kg)",
  main = "Species HWI vs Mass, Colored by CRB"
)

legend(
  "topright",
  legend = c("No CRB", "CRB"),
  col = c("grey", "blue"),
  pch = 19
)
# end ##########


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
  main = "Species HWI vs Mass\n(Color: Trophic, Shape: CRB)"
)

# Plot empty circles (CRB 0)
with(trait_data[trait_data$CRB_Final == 0, ],
     points(HWI, mass_kg, col = color, pch = 21, bg = "white", lwd = 2))

# Plot solid circles (CRB 1)
with(trait_data[trait_data$CRB_Final == 1, ],
     points(HWI, mass_kg, col = color, pch = 19))
# Make legend items for each trophic level
legend(
  "topright",
  legend = paste(rep(trophic_levels, each = 2),
                 rep(c("No CRB", "CRB"), times = n_levels)),
  col = rep(base_colors, each = 2),
  pch = rep(c(21, 19), times = n_levels),
  pt.bg = c(rep("white", n_levels), rep(base_colors, n_levels)),
  pt.lwd = 2,
  cex = 0.8
)

## predicting does HWI influence CRB?
# Basic logistic regression: CRB ~ HWI
model <- glm(CRB_Final ~ HWI, data = trait_data, family = binomial)
summary(model)

# Make a sequence of HWI values for prediction
HWI_seq <- seq(min(trait_data$HWI, na.rm = TRUE), max(trait_data$HWI, na.rm = TRUE), length.out = 100)

# Predict probability
pred <- predict(model, newdata = data.frame(HWI = HWI_seq), type = "response")

# Scatter plot + probability curve
plot(trait_data$HWI, trait_data$CRB_Final,
     xlab = "HWI",
     ylab = "Probability of CRB",
     pch = 19,
     col = "grey",
     main = "CRB vs HWI with Logistic Regression Fit")
lines(HWI_seq, pred, col = "blue", lwd = 2)

# "histogram" style of CRB by handwing
breaks <- seq(0, 60, by = 10)
trait_data$HWI_bin <- cut(trait_data$HWI, breaks = breaks, include.lowest = TRUE, right = FALSE)
tab <- table(trait_data$HWI_bin, trait_data$CRB_Final)
colnames(tab) <- c("No CRB", "CRB")
barplot(
  t(tab),
  beside = FALSE,   # stacked
  col = c("grey", "blue"),
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
  col = c("grey", "blue"),
  legend = TRUE,
  xlab = "HWI Bin",
  ylab = "Proportion of Species",
  main = "Proportion with/without CRB by HWI Bin",
  ylim = c(0, 1)
)

top10_mass <- trait_data[order(-trait_data$mass_kg), ][1:10, "Species"]
top10_mass$Species  # Print species names (or just use top10_mass)
top10_hwi <- trait_data[order(-trait_data$HWI), ][1:10, "Species"]
top10_hwi$Species  # Print species names (or just use top10_hwi)

# STacked chart by mass
fine_breaks <- seq(0, 1, by = 0.1)
coarse_breaks <- seq(2, 12, by = 1)   # 1 to 2 is a single bin, so start at 2
breaks <- c(fine_breaks, coarse_breaks)

trait_data$mass_bin <- cut(trait_data$mass_kg, breaks = breaks, include.lowest = TRUE, right = FALSE)
tab_mass <- table(trait_data$mass_bin, trait_data$CRB_Final)
colnames(tab_mass) <- c("No CRB", "CRB")

# Stacked bar chart (counts)
barplot(
  t(tab_mass),
  beside = FALSE,
  col = c("grey", "blue"),
  legend = TRUE,
  xlab = "Mass Bin (kg)",
  ylab = "Number of Species",
  main = "Species with/without CRB by Custom Mass Bin"
)

# 100% stacked (proportion)
tab_mass_prop <- prop.table(tab_mass, margin = 1)
barplot(
  t(tab_mass_prop),
  beside = FALSE,
  col = c("grey", "blue"),
  legend = TRUE,
  xlab = "Mass Bin (kg)",
  ylab = "Proportion of Species",
  main = "Proportion with/without CRB by Custom Mass Bin",
  ylim = c(0, 1)
)

