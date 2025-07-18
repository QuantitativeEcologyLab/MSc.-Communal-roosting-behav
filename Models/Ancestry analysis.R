#install.packages("furrr")  
library(furrr)
library(future)
plan(multisession, workers = 4)  # Adjust cores as needed
library(phytools)
library(dplyr)
library(future)
library(geiger)   # for fitMk
library(ape)
library(tidyr)
#install.packages("scico")
library(scico) #C"scico"library(scico) #Color palette


############PREPARE THE DATA FOR ANALYSIS###############
#call the data
#Load csv dataset
Bird_data <- read.csv("C:/Users/sandracd/OneDrive - UBC (1)/PhD proposal/Chapter 2/Updated database/Chapter_2_PhD_NEWdata_updated.csv")

#Remove underscores that were introduced by mistake
Bird_data$Species <- gsub("_+$", "", Bird_data$Species)

#Subset the variables that will go into the model
Bird_data_clean<- Bird_data %>%
  select(Species, Mass, Trophic_level, HWI, brain_volume_mm3, CRB_Final)

#Remove NA in CRB
Bird_data_clean <- Bird_data %>%
  select(Species, Mass, Trophic_level, HWI, brain_volume_mm3, CRB_Final) %>%
  drop_na(CRB_Final)

#Read tree
load("C:/Users/sandracd/OneDrive - UBC (1)/PhD proposal/Chapter 2/phylo_data/Consensus_Tree.Rda")

#subset objects
sub <- Bird_data_clean[Bird_data_clean$Species %in% phylogeny$tip.label,] 
sub$phylogeny <- sub$Species

#rename stuff
data(sub)
x<-setNames(sub$CRB_Final,
            rownames(sub))
tree<- phylogeny

#rename rows and make sure there are no repeated names
# Remove duplicates by selecting unique species
sub_clean <- sub %>%
  distinct(Species, .keep_all = TRUE)

# Assign row names after cleaning duplicates
rownames(sub_clean) <- gsub(" ", "_", sub_clean$Species)

# Now, proceed with the summarizing
sub_summary <- sub_clean %>%
  group_by(Species) %>%
  summarize(CRB_Final = mean(CRB_Final, na.rm = TRUE)) %>%
  ungroup()

# Set names for trait vector
x <- setNames(sub_summary$CRB_Final, gsub(" ", "_", sub_summary$Species))

# Match to tree
common_species <- intersect(tree$tip.label, names(x))
tree <- drop.tip(tree, setdiff(tree$tip.label, common_species))
x <- x[common_species]

# Convert to factor
x <- factor(x)
x <- droplevels(x)

# Optional: check
levels(x)
table(x)

############FIT THE GLOBAL MODELS###############
#Fit the global models
#Equal rates model assumes the rate of evolution 0-1 and 1-0 are equal
er <- fitMk(tree, x, model = "ER")
#All rates different implies each transition can differ in rates
ard<-fitMk(tree,x,model="ARD")
#Model 0 to 1 CRB
dirr01<-fitMk(tree,x,model=matrix(c(0,0,1,0),2,2))
#Model 1 to 0 CRB
dirr10<-fitMk(tree,x,model=matrix(c(0,1,0,0),2,2))
#best fit model assessment
anova(er,dirr01,dirr10,ard)

#Best fit model supports ard. Lets explore it
## marginal reconstruction under ARD model
anc_ard<-ancr(ard)
anc_ard


############FIT THE LOCAL MODEL###############
#Fit the local model
#Now lets estimate local node whioch means Q will not be constant but
#instead allows for changes at each node
## local estimation (under ARD model)

# Function to compute ancestral state probabilities for a single node
compute_node <- function(i, tree, x, model = "ARD") {
  node <- i + Ntip(tree)
  logL <- setNames(rep(NA, length(levels(x))), levels(x))
  
  for (j in levels(x)) {
    tip_name <- paste0("node_", node)
    tt <- try(bind.tip(tree, tip.label = tip_name, edge.length = 0, where = node), silent = TRUE)
    if (inherits(tt, "try-error")) next
    
    xx <- setNames(c(as.character(x), as.character(j)), c(names(x), tip_name))
    xx <- factor(xx, levels = levels(x))
    
    if (length(unique(xx)) < 2) next
    
    fit <- tryCatch({
      fitMk(tt, xx, model = model)
    }, error = function(e) NULL)
    
    if (!is.null(fit)) {
      logL[j] <- logLik(fit)
    }
  }
  
  if (all(!is.na(logL))) {
    return(exp(logL - max(logL)) / sum(exp(logL - max(logL))))  # Softmax
  } else {
    return(rep(NA, length(levels(x))))
  }
}

# Run the ancestral reconstructions in parallel
anc_list <- future_map(1:tree$Nnode, compute_node, tree = tree, x = x, .progress = TRUE)

#Save object
saveRDS(anc_list, file = "C:/Users/sandracd/OneDrive - UBC (1)/PhD proposal/Chapter 2/Models/anc_list.rds")
anc_list <- readRDS("Models/anc_list.rds")


# Get the root node number (i.e., node number, not index)
root_node_number <- getMRCA(tree, tree$tip.label)

# Convert to index used in anc_list (subtract number of tips)
root_index <- root_node_number - Ntip(tree)

# Retrieve the ancestral state probabilities
root_probs <- anc_list[[root_index]]

# Name the states
names(root_probs) <- levels(x)

# Print the result
print(round(root_probs, 3))

# Convert to matrix
anc_local <- do.call(rbind, anc_list)
rownames(anc_local) <- 1:tree$Nnode + Ntip(tree)
colnames(anc_local) <- levels(x)

############### PLOT ALL DATA AND ANCESTRAL STATE ###############
#Plot version 
#Plot and save
png("Figures/Ancestry analysis.png", width = 2000, height = 2000, res = 300)

# Custom state colors
state_colors <- c("#66C2A5", "#FC8D62")  # for 0 and 1

# Plot the tree with smaller tip label font
plotTree(tree,
         type = "fan",
         ftype = "i", # italic font 
         fsize = 0.1, # font size tip label
         offset = 0.5, 
         mar = rep(0, 4), 
         lwd = 1)

# Transparent frame for clean pie overlays
par(fg = "transparent")

# Node pie charts with new colors
nodelabels(pie = anc_local,
           cex = 0.1,
           piecol = state_colors)

# Smaller tip pie charts
tiplabels(pie = to.matrix(x[tree$tip.label], levels(x)),
          piecol = state_colors,
          cex = 0.1)                           

# Reset frame
par(fg = "black")

# Add legend
legend("topright",            
       legend = c("Absence", "Presence"),
       fill = state_colors,
       border = "black",
       bty = "n",               
       cex = 0.8,
       title = "CRB")

# Add panel label
mtext("b) ARD model (\"local\")", line = 0, adj = 0)
dev.off()





############### PLOT ACCIPITER DATA AND ANCESTRAL STATE ###############
#Plot subset on Accipiters

# Step 1: Get species in Accipitridae
accip_species <- Bird_data$Species[Bird_data$Family == "Accipitridae"]

# Step 2: Check which of these are in the tree
accip_species <- accip_species[accip_species %in% tree$tip.label]

# Step 3: Get MRCA of these species
mrca_node <- getMRCA(tree, accip_species)

# Step 4: Extract the clade
accip_tree <- extract.clade(tree, node = mrca_node)

# Get all internal nodes in accip_tree
accip_node_ids <- (Ntip(accip_tree) + 1):(Ntip(accip_tree) + accip_tree$Nnode)

# Map those to the original full tree's node numbers
accip_node_labels <- matchNodes(accip_tree, tree)[,2]  # columns: accip_tree node, original tree node

# Keep only rows that exist in anc_local
valid_nodes <- accip_node_labels[accip_node_labels %in% rownames(anc_local)]

# Subset and reorder the ancestral reconstructions
anc_accip <- anc_local[as.character(valid_nodes), , drop = FALSE]
rownames(anc_accip) <- accip_node_ids[accip_node_labels %in% rownames(anc_local)]  # match to accip_tree




#Plot and save
png("Figures/Ancestry analysis Accipiters.png", width = 2000, height = 2000, res = 300)
# Plot
plotTree(accip_tree,
         #type = "fan",
         ftype = "i",
         fsize = 0.4,
         lwd = 1,
         offset = 0.5,
         mar = rep(0, 4))
title("")

# Use your CRB state colors
state_colors <- c("#66C2A5", "#FC8D62")  # for 0 and 1

# Transparent plotting for clean overlay
par(fg = "transparent")

# Node pie charts
nodelabels(pie = anc_accip,
           piecol = state_colors,
           cex = 0.2)

tiplabels(
  pie = to.matrix(x[accip_tree$tip.label], levels(x)),
  piecol = state_colors,
  cex = 0.2,
  offset = 0.5
)

# Reset frame color
par(fg = "black")
dev.off()



#Number of independent transitions under parsimony assumption
# STEP 1: Ensure your binary trait is a named vector of tip states
# For example, from your dataset:
x <- setNames(sub_clean$CRB_Final, sub_clean$Species)

# Ensure the species names in x match those in your tree
x <- x[tree$tip.label]

# STEP 2: Run stochastic character mapping (simulate transitions)
set.seed(123)  # for reproducibility
simmap_trees <- make.simmap(tree, x, model = "ARD", nsim = 100)

# STEP 3: Summarize transitions
summary_simmap <- describe.simmap(simmap_trees)

# Print the number of transitions
summary_simmap$counts


#Identify transition clades
# Step 1: Identify most probable state at each node
node_states <- apply(anc_local, 1, which.max) - 1  # -1 to map 1=0, 2=1

# Step 2: Build full edge list with ancestral and descendant node states
# Edge matrix from phylo object
edge_df <- as.data.frame(tree$edge)
colnames(edge_df) <- c("ancestor", "descendant")

# Add state for ancestor
edge_df$descendant_state <- NA

# For internal nodes (descendants > number of tips)
internal_idx <- edge_df$descendant > Ntip(tree)
edge_df$descendant_state[internal_idx] <- node_states[edge_df$descendant[internal_idx] - Ntip(tree)]

# For tips
tip_idx <- !internal_idx
tip_labels <- tree$tip.label[edge_df$descendant[tip_idx]]
edge_df$descendant_state[tip_idx] <- x[tip_labels]

# Add state for descendant
edge_df$descendant_state <- ifelse(edge_df$descendant > Ntip(tree),
                                   node_states[edge_df$descendant - Ntip(tree)],
                                   x[tree$tip.label[edge_df$descendant]])

# Step 3: Filter for transitions
transitions <- edge_df %>%
  filter(ancestor_state != descendant_state)

# Step 4: Extract the descendant nodes that are transition origins
transition_nodes <- unique(transitions$descendant)

# Step 5: For each transition node, extract clade tip labels
clade_changes <- lapply(transition_nodes, function(node) {
  clade <- extract.clade(tree, node)
  clade$tip.label
})

names(clade_changes) <- paste0("Node_", transition_nodes)

# Optional: show example
clade_changes[[1]]


############### ANCESTRAL STATE FOR PREDICTOR VARIABLES ###############
#Ancestral state reconstruction for HWI mass and trophic level
# Clean and name the trait vectors
hwi_vec <- setNames(sub_clean$HWI, sub_clean$Species)
mass_vec <- setNames(sub_clean$Mass / 1000, sub_clean$Species)  # convert to kg if needed

# Match to tree
hwi_vec <- hwi_vec[tree$tip.label]
mass_vec <- mass_vec[tree$tip.label]

# Estimate ancestral states
hwi_anc <- fastAnc(tree, hwi_vec)
mass_anc <- fastAnc(tree, mass_vec)

# Get root node (first ancestral node is always the first index)
hwi_root <- hwi_anc[1]
mass_root <- mass_anc[1]

cat("Root HWI estimate:", round(hwi_root, 3), "\n")
cat("Root Mass estimate (kg):", round(mass_root, 3), "\n")


# === Visualize HWI ===
# Step 1: Prune the tree
accip_tree <- drop.tip(tree, setdiff(tree$tip.label, accip_species))

# Step 2: Subset the trait vector (HWI)
accip_hwi <- hwi_vec[names(hwi_vec) %in% accip_species]

# Step 3: Remap using contMap without plotting
accip_hwi_map <- contMap(accip_tree, accip_hwi, plot = FALSE, res = 3)

# Step 4: Apply the same green-to-orange custom color gradient
# custom_palette <- colorRampPalette(c("#66C2A5", "#FC8D62"))
n_hwi <- length(accip_hwi_map$cols)
accip_hwi_map$cols[] <- scico(n_hwi, palette = "imola")

# Step 5: Plot fan tree with custom colors
plot(accip_hwi_map,
     #type = "fan",
     fsize = 0.4,
     lwd = 2,
     outline = FALSE,
     legend = FALSE,
     main = "Accipitridae HWI — Custom Gradient")

# Manually add a color bar in bottom left
add.color.bar(
  leg = 0.2 * max(nodeHeights(accip_tree)),   # length of legend bar
  cols = accip_hwi_map$cols,
  title = "HWI",
  lims = round(range(accip_hwi, na.rm = TRUE), 2),
  digits = 2,
  prompt = FALSE,
  x = 0,
  y = -0.03 * max(nodeHeights(accip_tree)),   # position below tree
  lwd = 10,
  fsize = 0.6,       # <--- font size
  subtitle = ""
)


# === Visualize mass ===
# Step 1: Prune the tree to only Accipitridae
accip_tree <- drop.tip(tree, setdiff(tree$tip.label, accip_species))

# Step 2: Create a trait vector for mass
accip_mass <- mass_vec[names(mass_vec) %in% accip_species]

# Step 3: transform mass vector
accip_mass_log <- log10(accip_mass)  # log-transform mass


# Generate contMap object for mass
accip_mass_map <- contMap(accip_tree, accip_mass_log, plot = FALSE, res = 3)

# Apply custom color palette (e.g., from scico)
n_mass <- length(accip_mass_map$cols)
accip_mass_map$cols[] <- (scico(n_mass, palette = "imola"))


# # Custom color gradient: green (low) to orange (high)
# custom_palette <- colorRampPalette(c("#66C2A5", "#FC8D62"))

# # Apply the custom color scale correctly
# n_colors <- length(accip_mass_map$cols)
# accip_mass_map$cols[] <- custom_palette(n_colors)

# Plot using the fan layout
plot(accip_mass_map,
     #type = "fan",
     fsize = 0.4,
     lwd = 2,
     outline = FALSE,
     legend = 0.7 * max(nodeHeights(accip_tree)),
     main = "Accipitridae Mass (kg) — Custom Gradient")



# Plot the fan tree without the default legend
plot(accip_mass_map,
     #type = "fan",
     fsize = 0.5,
     lwd = 2,
     outline = FALSE,
     legend = FALSE,
     main = "Accipitridae Mass (kg) — Custom Gradient")

# Manually add a color bar in bottom left
add.color.bar(
  leg = 0.2 * max(nodeHeights(accip_tree)),   # length of legend bar
  cols = accip_mass_map$cols,
  title = "Mass (kg)",
  lims = round(range(accip_mass, na.rm = TRUE), 2),
  digits = 2,
  prompt = FALSE,
  x = 0,
  y = -0.03 * max(nodeHeights(accip_tree)),   # position below tree
  lwd = 10,
  fsize = 0.6,       # <--- font size
  subtitle = ""
)



# Save to file
png("Figures/Ancestral state Accipitridae_HWI_Mass.png", width = 3000, height = 1500, res = 300)

# Set up side-by-side panels
par(mfrow = c(1, 2), mar = c(1, 6, 4, 1))  # more space on left for labels

# --- Panel A: HWI ---
plot(accip_hwi_map,
     #type = "fan",
     fsize = 0.3,
     lwd = 1.5,
     outline = FALSE,
     legend = FALSE,
     main = "A. Accipitridae — HWI")

mtext("A", side = 3, line = -1, adj = 0.05, cex = 1, font = 2)  # moved down

add.color.bar(
  leg = 0.3 * max(nodeHeights(accip_tree)),
  cols = accip_hwi_map$cols,
  title = "HWI",
  lims = round(range(accip_hwi, na.rm = TRUE), 2),
  digits = 2,
  prompt = FALSE,
  x = 0,
  y = -0.03 * max(nodeHeights(accip_tree)),
  lwd = 10,
  fsize = 0.6
)

# --- Panel B: Mass ---
plot(accip_mass_map,
     #type = "fan",
     fsize = 0.3,
     lwd = 1.5,
     outline = FALSE,
     legend = FALSE,
     main = "B. Accipitridae — Mass (kg)")

mtext("B", side = 3, line = -1, adj = 0.05, cex = 1, font = 2)  # moved down

add.color.bar(
  leg = 0.4 * max(nodeHeights(accip_tree)),
  cols = accip_mass_map$cols,
  title = "log Mass (kg)",
  lims = round(range(accip_mass_log, na.rm = TRUE), 2),
  digits = 2,
  prompt = FALSE,
  x = 0,
  y = -0.03 * max(nodeHeights(accip_tree)),
  lwd = 10,
  fsize = 0.6
)
# Finish saving
dev.off()
