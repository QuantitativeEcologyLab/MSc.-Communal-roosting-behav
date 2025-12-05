# install.packages("geiger")   # if you don’t have it
#install.packages("phytools")
#install.packages("pushoverr")
#install.packages("brms")
library(dplyr)
library('ape')
library('phylobase')
library('brms')
library(loo)
library("tidyverse")
library('colorspace')
library("viridis")
library('ggplot2')
library('ggtree')
library("Rcpp")
library(geiger)
library(phytools)
library("glue") #for nice string formatting in the pushover message
library(pushoverr)
#install.packages("rphylopic")
library(rphylopic)
library(png)        # for reading PNG images if needed
library(jpeg)
set_pushover_user(user = "usdw6uw28zwd2496rr2dnni3efz4ds")
set_pushover_app(token = "aa21xm4a7vi9dcfv9d46p1qypg8miu")



######### STEP 1 - DATA LOADING AND CLEANING########

#Load Data (Github)
load("Models/Consensus_Tree.Rda")
Bird_data <- read.csv("Data/Bird_data_clean.csv")
#sub <- read.csv("/home/sandracd/Downloads/sub.csv")

#Wrangle dataset
Bird_data$Species <- gsub("_+$", "", Bird_data$Species)

#Subset the variables that will go into the model
Bird_data_clean<- Bird_data %>%
  select(Species, Family, Order, Mass, Trophic_level, HWI, brain_mass_g, CRB_Final)

#Remove NA in CRB and add mass kg. #Standardize all the parameters (mass, HWI)
Bird_data_clean <- Bird_data %>%
  select(Species,Family, Order, Mass, Trophic_level, HWI, brain_mass_g, CRB_Final) %>%
  drop_na(CRB_Final)%>%
  mutate(
    mass_kg = Mass/1000,
    mass_z = as.numeric(scale(mass_kg)),
    HWI_z = as.numeric(scale(HWI))
  ) 

#Create a list with Species and their matching phylogeny to run the model
sub <- Bird_data_clean[Bird_data_clean$Species %in% phylogeny$tip.label,] 
sub$phylogeny <- sub$Species

# Number of species you want to sample - all my tree (807 total)
n_species <- 807

# Select species names (From phylogeny- column contains species names)
selected_species <- sample(unique(sub$phylogeny), size = n_species)

# Subset the data
sub_subset <- sub[sub$phylogeny %in% selected_species, ]

#Prune tree to selected spp
phy_subset <- drop.tip(
  phy    = phylogeny,
  tip    = setdiff(phylogeny$tip.label, selected_species)
)

# Create a vector of species and trait data for the calculation
crb_vec <- setNames(sub_subset$CRB_Final, sub_subset$Species)

############ CALCULATION OF PHYLOGENETIC SIGNAL ############
# sig <- phylosig(phy_subset, crb_vec, method="lambda", test=TRUE)
# sig$lambda      # ML‐estimate of Pagel’s λ

#USING fitDiscrete for binary data
library(geiger)

# Step 1: Ensure unique species in sub_subset
sub_unique <- sub_subset[!duplicated(sub_subset$Species), ]

# Step 2: Keep only species that are in the phylogeny
valid_species <- sub_unique$Species[sub_unique$Species %in% phy_subset$tip.label]
sub_unique <- sub_unique[sub_unique$Species %in% valid_species, ]

# Step 3: Create the named binary trait factor
crb_factor <- setNames(as.factor(sub_unique$CRB_Final), sub_unique$Species)

# Step 4: Fit Pagel's lambda using discrete model
fit_obs <- fitDiscrete(phy_subset, crb_factor, model = "ARD", transform = "lambda")

# Check output
fit_obs$opt$lambda

#Now use that value to Specify phylogenetic autocorrelation with Pagels lambda in the model
# say you want λ = 0.9014789
λ <- 0.9014789

# this returns a new "phylo" with all *internal* branch lengths multiplied by λ
phy_lambda <- phytools::rescale(phylogeny, model="lambda", λ)

# now get the λ-transformed variance–covariance matrix
A_lambda <- ape::vcv.phylo(phy_lambda)


# Subset the phylogenetic covariance matrix for the model
A_subset <- A_lambda[selected_species, selected_species]

# Optional: check dimensions
dim(sub_subset)
dim(A_subset)


#Randomization test against observed lambda
# --- 1. Prepare cleaned data 
# Ensure one row per species
sub_subset_unique <- sub_subset[!duplicated(sub_subset$Species), ]
matched_species <- intersect(phy_subset$tip.label, sub_subset_unique$Species)

# Prune tree and data
phy_matched <- keep.tip(phy_subset, matched_species)
sub_matched <- sub_subset_unique[sub_subset_unique$Species %in% matched_species, ]

# Create named factor
crb_factor <- setNames(as.factor(sub_matched$CRB_Final), sub_matched$Species)

# --- 2. Fit observed model
fit_obs <- fitDiscrete(phy_matched, crb_factor, model = "ARD", transform = "lambda")
obs_lambda <- fit_obs$opt$lambda

# --- 3. Run randomizations
n_iter <- 1000
rand_lambdas <- numeric(n_iter)
set.seed(42)

for (i in 1:n_iter) {
  shuffled <- sample(as.character(crb_factor))  # shuffle trait values
  names(shuffled) <- names(crb_factor)
  
  fit_rand <- tryCatch(
    fitDiscrete(phy_matched, as.factor(shuffled), model = "ARD", transform = "lambda"),
    error = function(e) return(NULL)
  )
  
  rand_lambdas[i] <- if (!is.null(fit_rand)) fit_rand$opt$lambda else NA
}

# Remove NAs (in case any model fitting failed)
rand_lambdas <- rand_lambdas[!is.na(rand_lambdas)]

# --- 4. Compute p-value and quantiles
p_val <- mean(rand_lambdas >= obs_lambda)

# Compute 95th percentile
lambda_95 <- quantile(rand_lambdas, 0.95)
print(round(lambda_95, 3))
quantile (rand_lambdas)

# --- 5. Plot
par(mar = c(5.5, 5.5, 4, 2))
hist(rand_lambdas, breaks = 30, col = "grey", main = "Random values for Pagel's λ ", xlab = "λ")
abline(v = obs_lambda, col = "red", lwd = 2)
# abline(v = lambda_95, col = "blue", lwd = 2, lty = 2)
text(obs_lambda, max(table(cut(rand_lambdas, breaks = 30))) * 0.9,
     labels = paste("Observed λ =", round(obs_lambda, 3)), col = "red", pos = 2)
# Add label
# text(lambda_95, max(table(cut(rand_lambdas, breaks = 30))) * 0.8,
#      labels = paste("95th %ile =", round(lambda_95, 3)),
#      col = "blue", pos = 4)

# --- 6. Report
cat("Observed λ:", round(obs_lambda, 3), "\n")
cat("Randomization p-value:", round(p_val, 4), "\n")







######### STEP 2 - NULL MODEL ########


start <- proc.time()["elapsed"]

#Null model with phylogeny only 
Null_model <- brm(
  CRB_Final ~ 1+ (1|gr(phylogeny, cov = A_subset)),
  data = sub_subset,
  data2 = list(A_subset = A_subset),
  family = bernoulli,
  iter = 2000,
  warmup = 1000,
  chains = 4,   
  cores = 4,
  thin= 1
)
elapsed <- round(proc.time()["elapsed"] - start)/60 # minutes
print(glue("PHylo Model of {n_species} species finished in {elapsed} min."))
pushover(message = glue("Null Model of {n_species} species finished in {elapsed} min."))
summary(Null_model)
plot(Null_model)
saveRDS(Null_model,
        file="Models/Null_model.rds")


######### STEP 3 - BASE MODEL FOR PRIORS - NO PHYLOGENY ########

start <- proc.time()["elapsed"]

#Model with no phylogeny to get the priors from
base_model_no_phylo <- brm(
  CRB_Final ~ 0 + Trophic_level + mass_kg + HWI,  # no intercept, no random effects
  data = sub_subset,
  family = bernoulli(link = "logit"),
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  save_pars = save_pars(all = TRUE),
  seed = 123  # for reproducibility
)
summary(base_model_no_phylo)
saveRDS(base_model_no_phylo,
        file="Models/base_model_no_phylo.rds")

#Use estimates for priors
priors <- c(
  set_prior("student_t(3, -0.85, 0.5)", coef = "Trophic_levelCarnivore", class = "b"),
  set_prior("student_t(3, -0.33, 0.5)", coef = "Trophic_levelHerbivore", class = "b"),
  set_prior("student_t(3, -0.14, 0.5)", coef = "Trophic_levelOmnivore", class = "b"),
  set_prior("student_t(3, -0.93, 1.0)", coef = "Trophic_levelScavenger", class = "b"),
  set_prior("student_t(3, 0.20, 0.2)", coef = "mass_kg", class = "b"),
  set_prior("student_t(3, 0.03, 0.02)", coef = "HWI", class = "b")
)



######### STEP 4 - ALL DATA MODEL########

#Model w Phylo, Priors and 0 intercept
start <- proc.time()["elapsed"]

test_model_phyl_subset_40_PRIORS_trophic <- brm(
  CRB_Final ~ 0 + Trophic_level + mass_kg + HWI + (1 | gr(phylogeny, cov = A_subset)),
  data = sub_subset,
  data2 = list(A_subset = A_subset),
  family = bernoulli(link = "logit"),
  prior = priors,
  iter = 1e6,
  warmup = 5e5,
  chains =8,
  cores = 8,
  thin = 100
)

summary(test_model_phyl_subset_40_PRIORS_trophic)

elapsed <- (proc.time()["elapsed"] - start)/60 # minutes
print(glue("With PHylo Model and Priors of {n_species} species finished in {elapsed} min."))
pushover(message = glue("With PHylo Model and Priors of {n_species} species finished in {elapsed} min."))
summary(test_model_phyl_subset_40_PRIORS_trophic)
plot(test_model_phyl_subset_40_PRIORS_trophic)

#Save model
saveRDS(test_model_phyl_subset_40_PRIORS_trophic,
        file="Models/test_model_phyl_subset_40_PRIORS_trophic.rds")




######### STEP 5 - GLOBAL MODEL ON SUBSET FOR BRAIN COMPARISON ########

#Subset of dataset with all predictor variables - no brain data yet
#Create a list with Species and their matching phylogeny to run the model for BRAIN data
#Now create another subset for species with brain data only -remove NA
Bird_data_clean_brain <- Bird_data %>%
  select(Species, Mass, Trophic_level, HWI, brain_mass_g, CRB_Final) %>%
  drop_na(CRB_Final, brain_mass_g) %>%
  mutate(
    brain_mass_g = gsub(",", "", brain_mass_g),  # Remove commas from numbers
    brain_mass_g = trimws(brain_mass_g),  # Remove extra spaces
    brain_mass_g = as.numeric(brain_mass_g)  # Convert to numeric
  ) %>%
  mutate(mass_kg = Mass / 1000)%>%
  filter(!is.na(brain_mass_g))

#subset species
sub_brain <- Bird_data_clean_brain[Bird_data_clean_brain$Species %in% phylogeny$tip.label,] 
sub_brain$phylogeny <- sub_brain$Species

# Get species in the subset
species_brain <- unique(sub_brain$phylogeny)  # or sub_brain$Species if appropriate

# Prune the matrix to only those species
A_subset_brain <- A_subset[species_brain, species_brain]

Model_2 <- brm(
  CRB_Final  ~ 0 + Trophic_level + mass_kg + HWI + (1|gr(phylogeny, cov = A_subset_brain)),
  data = sub_brain,   
  data2 = list(A_subset_brain = A_subset_brain),
  family = bernoulli(link = "logit"),
  prior = priors,
  iter = 1e6,
  warmup = 5e5,
  chains =8,
  cores = 8,
  thin = 100)

summary(Model_2)
#Save model
saveRDS(Model_2,
        file="Models/Model_2.rds")


######### STEP 6 - GLOBAL MODEL ON SUBSET WITH BRAIN DATA ########

#Model the two variables using a linear model
# Model_mass_mass<-lm(brain_mass_g ~ Mass, data =sub_brain)
Model_mass_mass <- lm(log(brain_mass_g) ~ log(Mass), data = sub_brain)
summary(Model_mass_mass)


plot(Model_mass_mass)
summary(Model_mass_mass)

#Calculate the residuals of the model and store them in the database
Residuals_body_mass <-residuals(Model_mass_mass)


# Scatterplot with fitted line
png("brain_mass_vs_body_mass.png", width = 2000, height = 1500, res = 600)
plot(sub_brain$Mass, sub_brain$brain_mass_g,
     pch = 16, cex=0.5, col = "grey40",
     xlab = "Body mass (kg)",
     ylab = "Brain mass (g)",
     main = "Brain mass vs Body mass")


# Raw scatterplot (your code)
plot(sub_brain$Mass, sub_brain$brain_mass_g,
     pch = 16, cex = 0.5, col = "grey40",
     xlab = "Body mass (kg)",
     ylab = "Brain mass (g)",
     main = "Brain mass vs Body mass")

# --- Fit allometric model (log–log) ---
Model_mass_mass <- lm(log(brain_mass_g) ~ log(Mass), data = sub_brain)

# --- Generate prediction curve on raw x-axis ---
xseq <- seq(min(sub_brain$Mass, na.rm=TRUE),
            max(sub_brain$Mass, na.rm=TRUE),
            length.out = 200)

pred_log <- predict(Model_mass_mass, 
                    newdata = data.frame(Mass = xseq))
# Back-transform (exp) to original scale
pred_raw <- exp(pred_log)

# --- Add fitted curve to the plot ---
lines(xseq, pred_raw, col = "blue", lwd = 2)
dev.off()  # close the device

# Residuals vs fitted values
plot(fitted(Model_mass_mass), Residuals_body_mass,
     pch = 19, col = "darkred",
     xlab = "Fitted brain mass (g)",
     ylab = "Residuals",
     main = "Residuals vs Fitted Values")

abline(h = 0, lwd = 2, lty = 2)

# Calculate fitted values
fitted_vals <- fitted(Model_mass_mass)

# Add residuals and fitted values to the dataset
sub_brain$Residuals <- residuals(Model_mass_mass)
sub_brain$Fitted <- fitted_vals

# Identify species above and below the regression line
species_above <- sub_brain[sub_brain$Residuals > 0, ]
species_below <- sub_brain[sub_brain$Residuals < 0, ]






#attach to df
# If residuals are a matrix, convert to vector
res_vec <- as.vector(Residuals_body_mass)

# Attach to the dataframe
sub_brain$Residuals_body_mass <- res_vec


#base mnodel for brain prior
base_model_no_phylo_brain <- brm(
  CRB_Final ~ 0 + Trophic_level + mass_kg + HWI + Residuals_body_mass,  # no intercept, no random effects
  data = sub_brain,
  family = bernoulli(link = "logit"),
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  save_pars = save_pars(all = TRUE),
  seed = 123  # for reproducibility
)

summary(base_model_no_phylo_brain)
priors_brain <- c(
  prior(student_t(3, 0.03, 0.02), class = "b", coef = "HWI"),
  prior(student_t(3, 0.28, 0.33), class = "b", coef = "mass_kg"),
  prior(student_t(3, -0.70, 0.55), class = "b", coef = "Trophic_levelCarnivore"),
  prior(student_t(3, -0.22, 0.65), class = "b", coef = "Trophic_levelHerbivore"),
  prior(student_t(3, 0.12, 0.55), class = "b", coef = "Trophic_levelOmnivore"),
  prior(student_t(3, -0.97, 2.6), class = "b", coef = "Trophic_levelScavenger"),
  prior(student_t(3, -0.69, 0.80), class = "b", coef = "Residuals_body_mass")
)


#subset model w predictors and brain mass
Model_3 <- brm(
  CRB_Final  ~ 0 + Trophic_level + mass_kg + HWI + Residuals_body_mass + (1|gr(phylogeny, cov = A_subset_brain)),
  data = sub_brain,   
  data2 = list(A_subset_brain = A_subset_brain),
  family = bernoulli,
  prior = priors_brain,
  iter = 1e6,
  warmup = 5e5,
  chains = 8,
  cores = 8,
  thin = 100)


summary(Model_3)
#Save model
saveRDS(Model_3,
        file="Models/Model_3.rds")
