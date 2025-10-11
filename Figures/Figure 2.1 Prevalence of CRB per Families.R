library(dplyr)
library(ggplot2)
library(scales)
library(forcats)
library(readr)
library(rphylopic)

#Read data
trait_data <- read_csv("Data/Bird_data_clean.csv")

# Data arrangement - drop NA, calculate percentages and lump small families
plot_data <- trait_data %>%
  filter(!is.na(CRB_Final)) %>% 
  mutate(
    Family2 = fct_lump_min(Family, min = 25, other_level = "Other")
  ) %>% 
  count(Family2, CRB_Final) %>% 
  group_by(Family2) %>% 
  mutate(pct = n / sum(n) * 100) %>% 
  ungroup()


# Reorder Family2 by decreasing Presence (%)
plot_data <- plot_data %>%
  group_by(Family2) %>%
  mutate(presence_pct = pct[CRB_Final == 1]) %>%
  ungroup() %>%
  mutate(Family2 = fct_reorder(Family2, presence_pct, .desc = TRUE))

# Plot
Fig2.1 <- ggplot(plot_data, aes(x = Family2, y = pct, fill = factor(CRB_Final))) +
  geom_col() +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  scale_fill_manual(
    name   = "",
    values = c("0" = "#0072B2", "1" = "#E69F00"),  # Blue = Absence, Orange = Presence
    labels = c("0" = "Absence", "1" = "Presence")
  ) +
  labs(
    x        = "",
    y        = "",
    title    = ""
  ) +
  theme_classic() +
  theme(
    plot.title    = element_text(hjust = 0.5, size = 16),
    axis.text.x   = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 16),
    axis.text.y   = element_text(size = 16),
    legend.position = "right",
    legend.text   = element_text(size = 16)
  ) +
  add_phylopic(name = "Ploceidae", x = 1, y = 10, height = 10) +
  add_phylopic(name = "Sturnidae", x = 2, y = 10, height = 10) +
  add_phylopic(name = "Corvidae", x = 3, y = 10, height = 10) +
  add_phylopic(name = "Accipitridae", x = 4, y = 10, height = 10) +
  add_phylopic(name = "Psittacidae", x = 5, y = 10, height = 10, verbose = TRUE) +
  add_phylopic(name = "Other", x = 6, y = 10, height = 10)

Fig2.1

# save it
ggsave(
  filename ="Figures/Figure 2.1 Prevalence of CRB per Families.png",
  plot     = Fig2.1,
  width    = 10,
  height   = 8,
  units    = "in",
  dpi      = 300
)



