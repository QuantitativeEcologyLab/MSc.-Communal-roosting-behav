library(ape)
library(phytools)
library(RColorBrewer)
library(tidyverse)
library(viridis)
library(png)

# --- Helper function: Converts grayscale PNGs to RGB arrays ---
as_rgb_array <- function(img) {
  if (length(dim(img)) == 2) {
    # Grayscale: replicate to get 3 RGB channels
    array(rep(img, 3), dim = c(dim(img), 3))
  } else if (length(dim(img)) == 3 && dim(img)[3] %in% c(3, 4)) {
    # Already RGB or RGBA
    img
  } else {
    stop("Unknown image format")
  }
}

# --- 1. Load and prep data ---
load("Chapters/Consensus_Tree.Rda")       # loads phylogeny
phylo <- phylogeny
trait_data <- read_csv("Chapters/Bird_data_clean.csv")
trait_data$Species <- gsub(" ", "_", trait_data$Species)

# --- 2. Subsample 10% of species present in both data and tree ---
set.seed(42)
matched_species <- intersect(trait_data$Species, phylo$tip.label)
subsample_size <- round(0.1 * length(matched_species))
subsampled_species <- sample(matched_species, subsample_size)

# --- 3. Prune tree and data to the 10% sample ---
phylo_sub <- keep.tip(phylo, subsampled_species)
trait_sub <- trait_data %>% filter(Species %in% subsampled_species)

# --- 4. Set tip colors for binary trait (e.g., CRB_Final) ---
tip_traits_sub <- setNames(rep(NA, length(phylo_sub$tip.label)), phylo_sub$tip.label)
tip_traits_sub[trait_sub$Species] <- trait_sub$CRB_Final
tip_colors_sub <- ifelse(tip_traits_sub == 1, "#377eb8", "#e41a1c")
tip_colors_sub[is.na(tip_colors_sub)] <- "grey"

# --- 5. Plot fan tree ---
plot(phylo_sub, type="fan", tip.color=tip_colors_sub, cex=0.7, 
     label.offset=0.7, lwd=1)
legend("topright", legend=c("Trait = 1", "Trait = 0", "Missing"), pch=19, 
       col=c("#377eb8", "#e41a1c", "grey"), bty="n")

# --- 6. Add arc labels for Orders in the sample ---
orders.info <- trait_sub %>%
  group_by(Order) %>%
  summarise(num_in_order = n(),
            node = ifelse(n() > 1, findMRCA(phylo_sub, Species), which(phylo_sub$tip.label == Species[1])))

order_cols <- brewer.pal(n = nrow(orders.info), name = "Set3")
order_cex <- 0.4 + orders.info$num_in_order / 20

for (i in 1:nrow(orders.info)) {
  arc.cladelabels(
    text = orders.info$Order[i], 
    node = orders.info$node[i], 
    col = order_cols[i],
    ln.offset = 1.1, 
    lab.offset = 1.18, 
    cex = order_cex[i],
    mark.node = FALSE
  )
}

# --- 7. Overlay one or more phylopic images at specific tips ---

# List of species for which you have silhouette images (filenames must match exactly)
pics <- c("Acanthagenys_rufogularis") # Add more if desired

# Load and convert images to RGB arrays (so they always work with rasterImage)
phylopics <- lapply(pics, function(name) {
  as_rgb_array(readPNG(paste0("Data/Phylopic/", name, ".png")))
})
phylopics.dim.ratio <- sapply(phylopics, function(pic) dim(pic)[2]/dim(pic)[1])

# Find the index and angle of each pic's tip on the fan plot
phylopics.indices <- match(pics, phylo_sub$tip.label)
n_tips <- length(phylo_sub$tip.label)
phylopics.angles <- 2 * pi * (phylopics.indices - 1) / n_tips
radius <- 1.13  # just outside the tips

phylopics.x <- radius * cos(phylopics.angles)
phylopics.y <- radius * sin(phylopics.angles)

# Overlay the images, scaled so their area is similar regardless of aspect ratio
for (i in seq_along(phylopics)) {
  width <- 0.13 * phylopics.dim.ratio[i]
  height <- 0.13
  rasterImage(
    phylopics[[i]], 
    phylopics.x[i] - width/2, phylopics.y[i] - height/2,
    phylopics.x[i] + width/2, phylopics.y[i] + height/2
  )
}

