install.packages("furrr")  # if not already installed
library(furrr)
library(future)
plan(multisession, workers = 4)  # Adjust cores as needed
library(phytools)
library(dplyr)
library(furrr)
library(future)
library(geiger)   # for fitMk
library(ape)
library(tidyr)


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
anc_list <- readRDS("anc_list.rds")

# Convert to matrix
anc_local <- do.call(rbind, anc_list)
rownames(anc_local) <- 1:tree$Nnode + Ntip(tree)
colnames(anc_local) <- levels(x)


#Plot version 1
plotTree(tree, type = "fan", ftype = "off", mar = c(1,1,1,1), lwd = 0.5)
par(fg = "transparent")
nodelabels(pie = anc_local, cex = 0.3, piecol = hcl.colors(n = 2))
par(fg = "black")
mtext("ARD model – circular layout", line = -1, adj = 0)


#Plot version 2
#Plot and save
png("C:/Users/sandracd/OneDrive - UBC (1)/PhD proposal/Chapter 2/Models/Ancestry analysis.png", width = 2000, height = 2000, res = 300)

# Custom state colors
state_colors <- c("cornsilk1", "deeppink4")  # for 0 and 1

# Plot the tree with smaller tip label font
plotTree(tree,
         type = "fan",
         ftype = "i", # italic font 
         fsize = 0.2, # font size tip label
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

# Add panel label
mtext("b) ARD model (\"local\")", line = 0, adj = 0)
dev.off()







