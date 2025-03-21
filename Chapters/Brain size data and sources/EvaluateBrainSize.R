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


#now calculate the brain density for all available bird data
# Read the CSV file (replace 'yourfile.csv' with your actual file name)
data <- read.csv(file.path(getwd(), "BirdDatabase_wBrainMass.csv"))

# Remove commas from brain_volume and convert to numeric
data$brain_volume_mm3 <- as.numeric(gsub(",", "", data$brain_volume_mm3))

# Filter rows that have both brain mass and brain volume (removing NA values)
clean_data <- data[!is.na(data$brain.mass.g) & !is.na(data$brain_volume_mm3), ]

# Calculate brain density (assuming mass is in grams and volume is in mm³)
# Density = mass/volume
clean_data$brain_density_gcm3 <- clean_data$brain.mass.g / clean_data$brain_volume_mm3 * 1000

# Create a simple scatter plot of brain density
plot(clean_data$brain_volume_mm3,clean_data$brain.mass.g)

#lm(data=clean_data, brain_volume_mm3 ~ brain.mass.g)

plot(clean_data$brain_density_gcm3,
     main = "Brain Density Distribution",
     xlab = "Observation Number",
     ylab = "Brain Density (g/mm³)",
     pch = 16,  # Solid circle points
     col = "blue")

# Optional: Add a horizontal line at the mean density
abline(h = mean(clean_data$brain_density_gcm3), col = "red", lty = 2)

# Calculate mean and sd
mean_density <- mean(clean_data$brain_density_gcm3)
sd_density <- sd(clean_data$brain_density_gcm3)

# Create histogram
hist(clean_data$brain_density_gcm3,
     main = "Distribution of Brain Density",
     xlab = "Brain Density (g/cm³)",
     ylab = "Frequency",
     col = "lightblue",
     border = "black",
     breaks = 20)  # Adjust breaks for desired bin size

# Add mean line
abline(v = mean_density, col = "red", lwd = 2, lty = 2)

# Add text box with statistics
stats_text <- paste("Mean =", round(mean_density, 4),
                    "\nSD =", round(sd_density, 4),
                    "\nN =", nrow(clean_data))
text(x = max(clean_data$brain_density_gcm3) * 0.85, 
     y = max(hist(clean_data$brain_density_gcm3, plot = FALSE)$counts) * 0.5,
     labels = stats_text,
     pos = 4)  # pos = 4 puts text to right of coordinates

# Optional: Add normal curve overlay
x <- seq(min(clean_data$brain_density_gcm3), max(clean_data$brain_density_gcm3), length = 100)
curve(dnorm(x, mean = mean_density, sd = sd_density) * 
        length(clean_data$brain_density_gcm3) * diff(hist(clean_data$brain_density_gcm3, plot = FALSE)$breaks)[1],
      add = TRUE, col = "blue", lwd = 2)
