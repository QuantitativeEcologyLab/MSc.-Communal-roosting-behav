###################################
## Phylogenetic Tree Plot (10% Sample) - Using treedataverse + ggtree
## Displaying Clade Labels Only (No Images)
###################################

# ==== 1. Load Libraries ====  
library(treedataverse)    # loads ggtree, treeio, tidytree, ggtreeExtra  
library(ape)              # phylogenetic tree handling  
library(tidyverse)        # data manipulation  

# ==== 2. Load & Subset Data ====  
load("Chapters/Consensus_Tree.Rda")      # loads object `phylogeny`  
set.seed(42)  
sampled <- sample(phylogeny$tip.label, size = round(0.1 * length(phylogeny$tip.label)))  
phylo_subset <- keep.tip(phylogeny, sampled)  
trait_subset <- read_csv("Chapters/Bird_data_clean.csv") %>% 
  filter(Species %in% sampled)  

# ==== 3. Prepare Tip Traits ====  
tip_traits <- tibble(label = phylo_subset$tip.label) %>%  
  left_join(trait_subset %>% select(Species, Order, CRB_Final), by = c("label" = "Species")) %>%  
  mutate(  
    CRB_Final = replace_na(CRB_Final, 0),                  # default 0 if missing  
    tip_color = if_else(CRB_Final == 1, "red", "black"),  
    tip_shape = if_else(CRB_Final == 1, 16, 17)            # 16=circle, 17=triangle  
  )  

# ==== 4. Identify MRCA Nodes for Each Order ====  
order_mrcas <- tip_traits %>% 
  filter(!is.na(Order)) %>% 
  group_by(Order) %>% 
  filter(n() >= 2) %>% 
  summarise(node = getMRCA(phylo_subset, label), .groups = "drop")  

# ==== 5. Base Fan Plot ====  
plt <- ggtree(phylo_subset, layout = "fan", open.angle = 0) %<+% tip_traits +
  geom_tiplab(aes(label = label), size = 1.5, align = TRUE) +  
  geom_tippoint(aes(color = I(tip_color), shape = factor(CRB_Final)), size = 2) +  
  scale_shape_manual(values = c("0" = 17, "1" = 16), name = "Communal Roosting")  

# ==== 6. Add Clade Labels (Text Only) ====  
for (i in seq_len(nrow(order_mrcas))) {
  plt <- plt + 
    geom_cladelabel(
      node        = order_mrcas$node[i],
      label       = order_mrcas$Order[i],
      align       = FALSE,
      offset      = 60,
      offset.text = 30,
      barsize     = 2,
      angle       = 0,
      fontsize    = 3
    )
}

# ==== 7. Final Styling ====  
plt <- plt +  
  theme_tree2() +  
  ggtitle("Communal Roosting Behaviour") +  
  theme(
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank(),
    legend.position = "top"
  )  

# ==== 8. Display ====  
print(plt)  
cat("Half-circle fan plot with CRB coloring, clade labels (text only), and no images.
")
