# Load necessary libraries
library(readxl)
library(dplyr)
library(writexl)

# Define file paths
brain_mass_file <- file.path(getwd(), "/Chapters/Brain size data and sources/Bird_Brain_size.xlsx")
target_file <- file.path(getwd(), "/Chapters/main_database.csv")
output_file <- file.path(getwd(), "BirdDatabase_wBrainMass.csv")

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
avg_density <- mean((clean_data$brain.mass.g / clean_data$brain_volume_mm3) * 1000, na.rm = TRUE)

# Add calculated brain mass
updated_data$brain_mass_calculated <- ifelse(!is.na(updated_data$brain_volume_mm3),
                                             round(updated_data$brain_volume_mm3 * avg_density / 1000, 2),
                                             NA)

# Save updated data
write.csv(updated_data, output_file, row.names = FALSE)

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

