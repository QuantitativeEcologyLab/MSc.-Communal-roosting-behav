# Load necessary libraries
library(readxl)
library(dplyr)
library(writexl)  # Or use 'openxlsx' for more control

# Define file paths
brain_mass_file <- file.path(getwd(), "/Chapters/Brain size data and sources/BirdBrainData.csv")
target_file <- file.path(getwd(), "/Chapters/main_database.csv")
output_file <- file.path(getwd(), "BirdDatabase_wBrainMass.csv")  # The new file with added brain mass

# Read brain mass data
brain_data <- read.csv(brain_mass_file)

# Check column names
colnames(brain_data)  

brain_data <- brain_data %>%
  select(Species, brain.mass.g)

# Read the target file
target_data <- read.csv(target_file)

colnames(brain_data)
colnames(target_data)

# Merge (lookup) based on Species
updated_data <- target_data %>%
  left_join(brain_data, by = "Species")  # Merging on 'Species'

#Count successful matches (non-missing BrainMass values)
matched_count <- sum(!is.na(updated_data$brain.mass.g))
total_species <- nrow(target_data)
cat("Successfully matched species:", matched_count, "out of", total_species, "\n")


# Save the updated data to a new Excel file
write.csv(updated_data, output_file, row.names = FALSE)

# Print a success message
cat("Updated data saved to:", output_file)

