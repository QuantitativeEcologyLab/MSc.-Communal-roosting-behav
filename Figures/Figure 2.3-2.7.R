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
library(brms)

# Load data and tree #---------------------------------------------------------#
Bird_data_clean <- read_csv("Chapters/Bird_data_clean.csv")
load("Models/Consensus_Tree.Rda") #loads "phylogeny"

# Keep only matching species
tree <- keep.tip(phylogeny, phylogeny$tip.label[phylogeny$tip.label %in% Bird_data_clean$Species])
traits <- Bird_data_clean %>% filter(Species %in% tree$tip.label)


# Color palette for figures
primary_colors <- c(
  "#0072B2", # No CRB
  "#E69F00", # CRB
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
library(scales)
df_summary <- Bird_data_clean %>%
  group_by(Trophic_level, CRB_Final) %>%
  summarise(count = n(), .groups = 'drop') %>%
  mutate(
    CRB_Final = 
  )


Fig2.5 <- ggplot( df_summary, 
  aes(x = Trophic_level,
      y = count, 
      fill = factor(CRB_Final)
      )) +
  geom_col(position = "fill") +  # fill makes it 100% stacked
  scale_fill_manual(
    values = c(primary_colors[1], 
               primary_colors[2]),
    labels = c("Absence", "Presence"),
    name = "" ) +
  scale_y_continuous(labels = label_percent() )+
  labs(
    title = "CRB Prevalence by Trophic Level",
    x = "",
    y = "",
    fill = "CRB_Final") +
  theme_classic() +
  theme(
        plot.title    = element_text(hjust = 0.5, size=16),
        axis.text.x  = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size=16),
        axis.text.y = element_text(size=16),
        legend.position = "right",
        legend.text = element_text(size=16))  # Centers the title

ggsave(
  filename ="Figures/Figure 2.5 Prevalence of CRB per Trophic Level.png",
  plot     = Fig2.5,
  width    = 10,
  height   = 8,
  units    = "in",
  dpi      = 300
)




# Ensure CRB_Final is a factor or numeric with 0 and 1
Bird_data_clean %>%
  group_by(Trophic_level, CRB_Final) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Trophic_level) %>%
  mutate(percent = round(100 * n / sum(n), 1)) %>%
  arrange(Trophic_level, desc(CRB_Final))


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


####### FIGURES WITH MODEL AND UNIQUE VARIABLE ######
#HWI model and CRB
#load model
test_model_phyl_subset_40_PRIORS_trophic <- readRDS("Models/test_model_phyl_subset_40_PRIORS_trophic.rds")

# Create new data for prediction (fix other covariates, vary HWI)
newdata <- data.frame(
  HWI = seq(min(traits$HWI, na.rm = TRUE),
            max(traits$HWI, na.rm = TRUE),
            length.out = 100),
  mass_kg = mean(traits$mass_kg, na.rm = TRUE),           # Fix mass_kg at mean
  Trophic_level = "Carnivore",                            # Choose one trophic level
  phylogeny = NA                                           # Needed if phylogeny is a group-level effect
)

# Predict from the model
preds <- posterior_epred(test_model_phyl_subset_40_PRIORS_trophic, newdata = newdata, re_formula = NA)  # exclude random effects

# Compute summary statistics across iterations for each data point
pred_summary <- apply(preds, 2, function(x) {
  c(mean = mean(x),
    lower = quantile(x, 0.025),
    upper = quantile(x, 0.975))
})

# Transpose so each row is one prediction with its summary
pred_summary <- t(pred_summary)

# Format as data frame for plotting or saving
pred_df <- data.frame(
  HWI = newdata$HWI,
  Estimate = pred_summary[, "mean"],
  Lower = pred_summary[, "lower.2.5%"],
  Upper = pred_summary[, "upper.97.5%"]
)



#plot together
ggplot(Bird_data_clean, aes(x = HWI, y = CRB_Final)) +
  # Jittered observed values
  geom_jitter(aes(color = factor(CRB_Final), shape = factor(CRB_Final)),
              width = 0.2, height = 0.02, alpha = 0.5) +
  
  # Model predicted line and ribbon
  # geom_ribbon(data = pred_df,
  #             inherit.aes = FALSE,
  #             aes(x = HWI, ymin = Lower, ymax = Upper),
  #             fill = "grey70", alpha = 0.4) +
  geom_line(data = pred_df,
            inherit.aes = FALSE,
            aes(x = HWI, y = Estimate),
            color = "black", size = 1.2, linetype = "dashed") +
  # Color and shape
  scale_color_manual(values = c("#66C2A5", "#FC8D62"),
                     labels = c("CRB absent", "CRB present")) +
  scale_shape_manual(values = c(1, 4),
                     labels = c("CRB absent", "CRB present")) +
  
  # Labels
  labs(
    x = "Hand-wing Index (HWI)",
    y = "Probability of communal roosting",
    color = "CRB Status",
    shape = "CRB Status"
  ) +
  
  # Theme
  theme_classic(base_size = 14) +
  theme(
    legend.position = c(0.8, 0.8),
    legend.box.background = element_rect(color = "black"),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = "Figures/CRB_HWI_model_prediction.png",
  width = 10,          # width in inches
  height = 6,          # height in inches
  dpi = 300            # high resolution
)



# Create sequence of mass values across observed range
mass_seq <- seq(min(traits$mass_kg, na.rm = TRUE),
                max(traits$mass_kg, na.rm = TRUE), length.out = 100)

# Create new data frame for prediction
newdata_mass <- data.frame(
  mass_kg = mass_seq,
  HWI = median(traits$HWI, na.rm = TRUE),
  Trophic_level = "Carnivore"  # or any reference level you prefer
)


# Generate predictions without random effects
preds_mass <- posterior_epred(test_model_phyl_subset_40_PRIORS_trophic,
                              newdata = newdata_mass, re_formula = NA)

# Summarize across posterior draws
# Summarize across posterior draws and transpose
pred_mass_summary <- t(apply(preds_mass, 2, function(x) c(
  mean = mean(x),
  lower = quantile(x, 0.025),
  upper = quantile(x, 0.975)
)))

pred_mass_df <- data.frame(
  mass_kg = newdata_mass$mass_kg,
  Estimate = pred_summary[, "mean"],
  Lower = pred_summary[, "lower.2.5%"],
  Upper = pred_summary[, "upper.97.5%"]
)

#Plot
ggplot(traits, aes(x = mass_kg, y = CRB_Final)) +
  geom_jitter(aes(color = factor(CRB_Final), shape = factor(CRB_Final)),
              width = 0.2, height = 0.02, alpha = 0.5) +
  scale_color_manual(values = c("#66C2A5", "#FC8D62"),
                     labels = c("CRB absent", "CRB present")) +
  scale_shape_manual(values = c(1, 4),
                     labels = c("CRB absent", "CRB present")) +
  geom_line(data = pred_mass_df, aes(x = mass_kg, y = Estimate),
            color = "black", size = 1.2, linetype = "dashed") +
  # geom_ribbon(data = pred_mass_df,
  #             aes(x = mass_kg, ymin = Lower, ymax = Upper),
  #             alpha = 0.2, fill = "black") +
  labs(x = "Body Mass (kg)",
       y = "Estimated probability of communal roosting",
       color = "CRB Status", shape = "CRB Status") +
  theme_classic(base_size = 14) +
  theme(legend.position = c(0.8, 0.8),
        legend.box.background = element_rect(color = "black"),
        panel.grid.minor = element_blank())

ggsave(
  filename = "Figures/CRB_mass_model_prediction.png",
  width = 10,          # width in inches
  height = 6,          # height in inches
  dpi = 300            # high resolution
)
