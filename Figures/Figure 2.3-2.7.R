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
load("Models/Consensus_Tree.Rda") #loads "phylogeny"

# Keep only matching species
tree <- keep.tip(phylogeny, phylogeny$tip.label[phylogeny$tip.label %in% Bird_data_clean$Species])
traits <- Bird_data_clean %>% filter(Species %in% tree$tip.label)


# Color palette for figures
primary_colors <- c(
  "#66C2A5", # No CRB
  "#FC8D62", # CRB
  "#8DA0CB", # Light Blue
  "#E78AC3"  # Pink/Lilac
)


########### HWI AND MASS/KG AND CRB IN ALL SPECIES ###########
png(
  filename = "Figures/Figure 2.3 Relationship between Mass and HWI and CRB prevalence.png",
  width    = 10,     # inches
  height   = 8,     # inches
  units    = "in",
  res      = 300    # DPI
)

# x y chart of HWI vs Mass for CRB
cols <- c(primary_colors[1], primary_colors[2])[as.factor(traits$CRB_Final)]
par(mar=c(5,5,5,5))

#plot
plot(
  traits$HWI, traits$mass_kg,
  col = cols,
  pch = 20,
  xlab = "HWI",
  ylab = "Mass (kg)",
  main = "Relationship between Mass and HWI and CRB prevalence"
)
legend(
  "topright",
  legend = c("Absence", "Presence"),
  col = c(primary_colors[1], primary_colors[2]),
  pch = 20
)

#close to save
dev.off()




########### HWI AND MASS/KG AND CRB IN ALL SPECIES PER TROPHIC LEVEL ###########
# x and y HWI and mass for CRB
trophic_levels <- unique(traits$Trophic_level)
n_levels <- length(trophic_levels)
base_colors <- c(primary_colors[1], primary_colors[2], primary_colors[3], primary_colors[4])[1:n_levels]
names(base_colors) <- trophic_levels

# Assign color by trophic level
traits$color <- base_colors[traits$Trophic_level]

# Assign pch: 19 = solid, 21 = empty circle with border color
traits$pch <- ifelse(traits$CRB_Final == 1, 19, 21)

png(
  filename = "Figures/Figure 2.4 Relationship between Mass and HWI per trophic level and CRB prevalence.png",
  width    = 10,     # inches
  height   = 8,     # inches
  units    = "in",
  res      = 300    # DPI
)
plot(
  traits$HWI, traits$mass_kg,
  type = "n",
  xlab = "HWI",
  ylab = "Mass (kg)",
  main = "Relationship between Mass and HWI per trophic level and CRB prevalence"
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

#close to save
dev.off()




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
legend_labels <- c(rep("Absence", n), rep("Presence", n))  # or "Lighter"/"Darker"

Fig2.5 <- ggplot(df_summary, aes(x = Trophic_level, y = count, fill = fill_group)) +
  geom_col(position = "fill") +  # fill makes it 100% stacked
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(
    values = fill_colors,
    labels = legend_labels,
    name = ""   # or "CRB Status"
  ) +
  labs(
    title = "CRB Prevalence by Trophic Level",
    x = "",
    y = "Percentage",
    fill = "CRB_Final") +
  theme_classic()+
  theme(
    plot.title = element_text(hjust = 0.5)  # Centers the title
  )

ggsave(
  filename ="Figures/Figure 2.5 Prevalence of CRB per Trophic Level.png",
  plot     = Fig2.5,
  width    = 10,
  height   = 8,
  units    = "in",
  dpi      = 300
)





########### PERCENTAGE CRB PER HWI #############
# 100% stacked by chart by HWI 
#Set limits of values
breaks <- seq(0, 60, by = 10)
traits$HWI_bin <- cut(traits$HWI, breaks = breaks, include.lowest = TRUE, right = FALSE)
tab <- table(traits$HWI_bin, traits$CRB_Final)
colnames(tab) <- c("Absence", "Presence")

# Compute proportions by HWI bin (row)
tab_prop <- prop.table(tab, margin = 1)

par(mar=c(5,5,5,5))

#Open to save
png(
  filename = "Figures/Figure 2.6 Prevalence of CRB per HWI.png",
  width    = 10,     # inches
  height   = 8,     # inches
  units    = "in",
  res      = 300    # DPI
)
barplot(
  t(tab_prop),
  beside = FALSE,         # stacked bars
  col = c(primary_colors[1], primary_colors[2]),
  xlab = "HWI",
  ylab = "Percentage",
  main = "CRB prevalence in relation to HWI",
  ylim = c(0, 1)
)
legend(
  "topright",
  legend = c("Absence", "Presence"),
  fill   = c(primary_colors[1], primary_colors[2]),
  bty    = "n",
  xpd    = TRUE,         # allow drawing outside plot region
  inset  = c(-0.1, 0)
)

#close to save
dev.off()





########### PERCENTAGE CRB PER MASS #############
# Group mass every 2 kg
breaks_mass <- seq(
  from = 0,
  to   = ceiling(max(traits$mass_kg, na.rm = TRUE)),
  by   = 2
)

# Cut into bins and tabulate vs. CRB_Final
traits$mass_bin <- cut(
  traits$mass_kg,
  breaks         = breaks_mass,
  include.lowest = TRUE,
  right          = FALSE
)

tab_mass <- table(traits$mass_bin, traits$CRB_Final)
colnames(tab_mass) <- c("Absence", "Presence")

# Convert to  proportions
tab_mass_prop <- prop.table(tab_mass, margin = 1)

#Plot 
par(mar = c(5, 5, 5, 5))
#Open to save
png(
  filename = "Figures/Figure 2.7 Prevalence of CRB per MASS.png",
  width    = 10,     # inches
  height   = 8,     # inches
  units    = "in",
  res      = 300    # DPI
)

barplot(
  t(tab_mass_prop),
  beside = FALSE,
  col    = c(primary_colors[1], primary_colors[2]),
  xlab   = "Mass (kg)",
  ylab   = "Percentage",
  main   = "CRB prevalence in relation to body mass",
  ylim   = c(0, 1)
)

#Push legend further right
legend(
  "topright",
  legend = c("Absence", "Presence"),
  fill   = c(primary_colors[1], primary_colors[2]),
  bty    = "n",
  xpd    = TRUE,
  inset  = c(-0.1, 0)
)


#close to save
dev.off()

















