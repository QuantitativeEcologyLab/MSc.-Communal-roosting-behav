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

# Base tree + trait data
p <- ggtree(tree, layout = "fan") %<+% traits

# Color tips by CRB_Final
p1 <- p +
  geom_tippoint(aes(color = CRB_Final), size = 1.2) +
  scale_color_viridis(option = "inferno")

# Show the plot, most plain version
print(p1)

###############################################
#Plot with Order Labels to outside on an angle
# Example base plot:
p <- ggtree(tree, layout = "fan") %<+% traits + 
  geom_tippoint(aes(color = CRB_Final), size = 1.2) +
  scale_color_viridis(option = "inferno")

# Join traits data with the tree tip data
tip_data <- data.frame(label = tree$tip.label) %>%
  left_join(traits, by = c("label" = "Species"))

# Extract the ggtree plot data for tips
tip_positions <- p$data %>% 
  filter(isTip) %>%
  left_join(tip_data, by = c("label")) 

# Each tip has an angle, we want the mean angle per Order
order_positions <- tip_positions %>%
  group_by(Order.x) %>%
  summarise(angle = mean(angle), # average angle for label placement
            x = mean(x), # average x (radius)
            y = mean(y)) %>%
  ungroup()

order_positions <- order_positions %>%
  mutate(
    # Adjust label angle for readability
    angle_adj = ifelse(angle > 90 & angle < 270, angle + 180, angle),
    
    # Adjust horizontal justification
    hjust = ifelse(angle > 90 & angle < 270, 1, 0)
  )



# Multiply x and y by a factor to push labels outside the tree radius
label_radius_factor <- 1.05
p + 
  geom_text(data = order_positions,
            aes(x = x * label_radius_factor, 
                y = y * label_radius_factor, 
                label = Order.x, 
                angle = angle_adj,
                hjust = hjust),
            size = 3)

######################
#Version with Mass bar charts
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
#Version with Heatmap
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


###########################
#Rectangular Plot
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

#############################
# Calculate percentage of species by Trophic_level and CRB_Final


df_summary <- Bird_data_clean %>%
  group_by(Trophic_level, CRB_Final) %>%
  summarise(count = n(), .groups = 'drop') 

ggplot(df_summary, aes(x = Trophic_level, y = count, fill = factor(CRB_Final))) +
  geom_col(position = "fill") +  # fill makes it 100% stacked
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "Proportion of CRB_Final by Trophic Level",
    x = "Trophic Level",
    y = "Percentage",
    fill = "CRB_Final"
  ) +
  theme_minimal()

ggplot(Bird_data_clean, aes(x = HWI, y = Mass)) +
  geom_hex() + facet_wrap(~CRB_Final)


###############################
## Hex Bin Plot of CRB by HWI and Mass

# Split the data by CRB_Final group
df_crb0 <- subset(Bird_data_clean, CRB_Final == 0)
df_crb1 <- subset(Bird_data_clean, CRB_Final == 1)

ggplot() +
  # CRB 0 group: Red → Grey
  geom_hex(data = df_crb0, aes(x = HWI, y = Mass, fill = ..count..), bins = 30, alpha = 0.6) +
  scale_fill_gradient(
    low = "yellow",
    high = "red",
    name = "Not Communal Rooster"
  ) +
  
  # Reset the fill scale
  ggnewscale::new_scale_fill() +
  
  # CRB 1 group: Blue → White
  geom_hex(data = df_crb1, aes(x = HWI, y = Mass, fill = ..count..), bins = 30, alpha = 0.6) +
  scale_fill_gradient(
    low = "lightblue",
    high = "blue",
    name = "Communal Rooster"
  ) +
  
  # Labels and theme
  labs(
    title = "HWI vs Mass with CRB Groups Overlaid",
    x = "Hand Wing Index (HWI)",
    y = "Mass"
  ) +
  theme_minimal()

##################################
## Bubble Chart

ggplot(Bird_data_clean, aes(x = HWI, y = Mass, size = brain_volume_mm3, color = as.factor(CRB_Final))) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(name = "Brain Volume (mm³)") +
  scale_color_manual(values = c("0" = "red", "1" = "blue"), labels = c("Not Communal", "Communal")) +
  labs(x = "Hand Wing Index", y = "Mass", color = "CRB Status") +
  theme_minimal()


## Bubble Chart with Tropic Level by Color
ggplot(Bird_data_clean, aes(x = HWI, y = Mass,
                            size = brain_volume_mm3,
                            color = Trophic_level,
                            shape = as.factor(CRB_Final))) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(name = "Brain Volume (mm³)") +
  scale_shape_manual(
    name = "CRB Status",
    values = c("0" = 16, "1" = 17),    # 16 = circle, 17 = triangle
    labels = c("0" = "Not Communal", "1" = "Communal")
  ) +
  labs(
    x = "Hand Wing Index (HWI)",
    y = "Mass",
    color = "Trophic Level"
  ) +
  theme_minimal()


##################################




