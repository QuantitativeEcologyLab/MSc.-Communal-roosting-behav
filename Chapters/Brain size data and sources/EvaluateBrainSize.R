# Load necessary libraries
library(readxl)
library(dplyr)
library(writexl)
library(openxlsx)

#new comment

# Define file paths
#final data source output, and all inputs
brain_mass_file <- file.path(getwd(), "/Chapters/Brain size data and sources/AllDataSourcesCombined.xlsx")

target_file <- file.path(getwd(), "/Chapters/main_database.csv")   #data source 1
output_file <- file.path(getwd(), "BirdDatabase_wBrainMass.csv")   #data source 

# Read brain mass data from Excel file, "Combined" sheet
brain_data <- read_excel(brain_mass_file, sheet = "Combined") %>%
  select(Species, brain.mass.g)

# Read and process target file
target_data <- read.csv(target_file)
target_data$brain_volume_mm3 <- as.numeric(gsub(",", "", target_data$brain_volume_mm3))

# Merge data
updated_data <- target_data %>%
  left_join(brain_data, by = "Species")

# Count and display successful matches
matched_count <- sum(!is.na(updated_data$brain.mass.g))
total_species <- nrow(target_data)
cat("Successfully matched species:", matched_count, "out of", total_species, "\n")

# Calculate average brain density from existing data
data <- read.csv(output_file)
data$brain_volume_mm3 <- as.numeric(gsub(",", "", data$brain_volume_mm3))
clean_data <- data[!is.na(data$brain.mass.g) & !is.na(data$brain_volume_mm3), ]
avg_density <- mean((clean_data$brain.mass.g / clean_data$brain_volume_mm3) * 1000, na.rm = TRUE) #grams/cm3

# Add calculated brain mass
updated_data$brain_mass_calculated <- ifelse(!is.na(updated_data$brain_volume_mm3),
                                             round(updated_data$brain_volume_mm3 * avg_density / 1000, 2),
                                             NA)

# compare new data from BirdlifeAustralia
file_path <- file.path(getwd(), "/Chapters/Brain size data and sources/AllDataSourcesCombined.xlsx")
# Read the bird name columns from each sheet
birds_sheet1 <- read_excel(file_path, sheet = "Combined")$Species
#birds_sheet2 <- read_excel(file_path, sheet = "BirdlifeAustralia")$Species
birds_sheet2 <- read_excel(file_path, sheet = "BFncomms")$species

#birds_sheet2 <- gsub(" ", "_", birds_sheet2)
matching_species <- birds_sheet2[birds_sheet2 %in% birds_sheet1]
num_matches <- length(matching_species)

cat("Number of matches:", num_matches, "\n")
cat("First 5 matches:\n")
print(head(matching_species, 5))


#add data in from BFncomms
df <- read_excel(file_path, sheet = "BFncomms")
density <- 1.09497

df <- df %>%
  mutate(
    calculated_brain_volume = exp(LogBrain),           # back-calculate from ln(volume)
    calculated_brain_mass = calculated_brain_volume / density  # get brain mass in grams
  )
head(df)

write.csv(df, "BFncomms_with_calculated_brain_data.csv", row.names = FALSE)

#search for only non matching rows in BFncomms, and add to combined sheet.

# Read species from 'Combined' sheet
combined_df <- read_excel(file_path, sheet = "Combined")
combined_species <- combined_df$Species
# Read data from 'BFncomms' sheet
df_bfncomms <- read_excel(file_path, sheet = "BFncomms")
# Find non-matching species
non_matching <- df_bfncomms[!df_bfncomms$species %in% combined_species, ]
# Prepare output: species and calculated_brain_mass renamed as brain.mass.g
output <- non_matching %>%
  select(species, calculated_brain_mass) %>%
  rename(Species = species, brain.mass.g = calculated_brain_mass)
# Add empty columns to match Combined structure
missing_cols <- setdiff(names(combined_df), names(output))
output[missing_cols] <- ""
# Reorder columns to match Combined sheet
output <- output[, names(combined_df)]
# Save to CSV
write.csv(output, "non_matching_species.csv", row.names = FALSE)

# Calculate brain density for plotting
data <- read.csv(output_file)
data$brain_volume_mm3 <- as.numeric(gsub(",", "", data$brain_volume_mm3))
clean_data <- data[!is.na(data$brain.mass.g) & !is.na(data$brain_volume_mm3), ]
clean_data$brain_density_gcm3 <- (clean_data$brain.mass.g / clean_data$brain_volume_mm3) * 1000

# Generate plots
plot(clean_data$brain_volume_mm3, clean_data$brain.mass.g)
plot(clean_data$brain_density_gcm3, main = "Brain Density Distribution", xlab = "Observation Number", ylab = "Brain Density (g/cm³)", pch = 16, col = "blue")
abline(h = mean(clean_data$brain_density_gcm3), col = "red", lty = 2)

hist(clean_data$brain_density_gcm3, main = "Distribution of Brain Density", xlab = "Brain Density (g/cm³)", ylab = "Frequency", col = "lightblue", border = "black", breaks = 20)
abline(v = mean(clean_data$brain_density_gcm3), col = "red", lwd = 2, lty = 2)

plot(clean_data$brain.mass.g, clean_data$brain_density_gcm3, main = "Brain Mass vs Density", xlab = "Brain Mass (g)", ylab = "Brain Density (g/cm³)", pch = 16, col = "blue")
abline(lm(brain_density_gcm3 ~ brain.mass.g, data = clean_data), col = "red")

plot(clean_data$brain_volume_mm3, clean_data$brain_density_gcm3, main = "Volume vs Density", xlab = "Brain Volume (mm³)", ylab = "Brain Density (g/cm³)", pch = 16, col = "green")
abline(lm(brain_density_gcm3 ~ brain_volume_mm3, data = clean_data), col = "red")

plot(clean_data$brain.mass.g, clean_data$brain_volume_mm3, main = "Mass vs Volume", xlab = "Brain Mass (g)", ylab = "Brain Volume (mm³)", pch = 16, col = "purple")
abline(lm(brain_volume_mm3 ~ brain.mass.g, data = clean_data), col = "red")

hist(clean_data$brain_density_gcm3, main = "Density Distribution", xlab = "Brain Density (g/cm³)", col = "lightblue", breaks = 20)

plot(clean_data$Mass, clean_data$brain.mass.g, main = "Body Mass vs Brain Mass", xlab = "Body Mass (g)", ylab = "Brain Mass (g)", pch = 16, col = "blue")

