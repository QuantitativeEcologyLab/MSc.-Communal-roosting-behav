# Load necessary libraries
library(readxl)
library(dplyr)
library(writexl)  # Or use 'openxlsx' for more control

#Average Brain Density is 1.1046g/cm3 based on the data where both mass and volume exist

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

# Remove commas from brain_volume and convert to numeric
target_data$brain_volume_mm3 <- as.numeric(gsub(",", "", target_data$brain_volume_mm3))

colnames(brain_data)
colnames(target_data)

# Count total and non-NA values for brain.mass.g and brain_volume_mm3
total_brain_mass <- length(data$brain.mass.g)
non_na_brain_mass <- sum(!is.na(data$brain.mass.g))
total_brain_volume <- length(target_data$brain_volume_mm3)
non_na_brain_volume <- sum(!is.na(target_data$brain_volume_mm3))
both_mass_and_volume <- sum(!is.na(updated_data$brain.mass.g) & !is.na(updated_data$brain_volume_mm3))

# Merge (lookup) based on Species
updated_data <- target_data %>%
  left_join(brain_data, by = "Species")  # Merging on 'Species'

#Count successful matches (non-missing BrainMass values)
matched_count <- sum(!is.na(updated_data$brain.mass.g))
total_species <- nrow(target_data)
cat("Successfully matched species:", matched_count, "out of", total_species, "\n")

# Calculate average brain density from the filtered data
# (This requires the earlier data processing steps to be run first)
data <- read.csv(file.path(getwd(), "BirdDatabase_wBrainMass.csv"))
data$brain_volume_mm3 <- as.numeric(gsub(",", "", data$brain_volume_mm3))
clean_data <- data[!is.na(data$brain.mass.g) & !is.na(data$brain_volume_mm3), ]
clean_data$brain_density_gcm3 <- (clean_data$brain.mass.g / clean_data$brain_volume_mm3) * 1000
avg_density <- mean(clean_data$brain_density_gcm3, na.rm = TRUE)

# Add calculated brain mass using average density
# Mass (g) = Density (g/cm³) × Volume (mm³) / 1000 (to convert mm³ to cm³)
updated_data$brain_mass_calculated <- ifelse(!is.na(updated_data$brain_volume_mm3),
                                             round(updated_data$brain_volume_mm3 * avg_density / 1000, 2),
                                             NA)

# Save the updated data to a new CSV file
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

# 1. Scatter plot: Brain Mass vs Density
plot(clean_data$brain.mass.g, clean_data$brain_density_gcm3,
     main = "Brain Mass vs Density",
     xlab = "Brain Mass (g)",
     ylab = "Brain Density (g/mm³)",
     pch = 16, 
     col = "blue")
# Add a trend line
abline(lm(brain_density ~ brain_mass, data = clean_data), col = "red")

# 2. Volume vs Density
plot(clean_data$brain_volume_mm3, clean_data$brain_density_gcm3,
     main = "Volume vs Density",
     xlab = "Brain Volume (mm³)",
     ylab = "Brain Density (g/cm³)",
     pch = 16, col = "green")
abline(lm(brain_density_gcm3 ~ brain_volume_mm3, data = clean_data), col = "red")

# 3. Mass vs Volume
plot(clean_data$brain.mass.g, clean_data$brain_volume_mm3,
     main = "Mass vs Volume",
     xlab = "Brain Mass (g)",
     ylab = "Brain Volume (mm³)",
     pch = 16, col = "purple")
abline(lm(brain_volume_mm3 ~ brain.mass.g, data = clean_data), col = "red")

# 4. Density Histogram (repeated for completeness)
hist(clean_data$brain_density_gcm3,
     main = "Density Distribution",
     xlab = "Brain Density (g/cm³)",
     col = "lightblue",
     breaks = 20)

# Reset plot layout
par(mfrow = c(1, 1))

# Correlation test: Brain Mass vs Density
cor_test <- cor.test(clean_data$brain.mass.g, clean_data$brain_volume_mm3, method = "pearson")
cor_test <- cor(clean_data$brain.mass.g, clean_data$brain_volume_mm3)
print("Correlation between brain mass and density:")
print(cor_test)

# Summary statistics
print("Summary Statistics:")
summary(clean_data[, c("brain.mass.g", "brain_volume_mm3", "brain_density_gcm3")])

#Sandra's model
model_lm <- lm(clean_data$brain.mass.g ~ clean_data$brain_volume_mm3, data = clean_data )
summary(model_lm)

# Scatter plot: Body Mass vs Brain Mass
plot(clean_data$Mass, clean_data$brain.mass.g,
     main = "Body Mass vs Brain Mass",
     xlab = "Body Mass (g)",
     ylab = "Brain Mass (g)",
     pch = 16,          # Solid circle points
     col = "blue")      # Point color

cor_test <- cor.test(clean_data$Mass, clean_data$brain.mass.g, method = "pearson")
cat("Correlation between body mass and brain mass:\n")
print(cor_test)

