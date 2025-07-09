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

#Calculate phylogenetic signal
sig <- phylosig(phy_subset, crb_vec, method="lambda", test=TRUE)
sig$lambda      # ML‐estimate of Pagel’s λ


#Now use that value to Specify phylogenetic autocorrelation with Pagels lambda in the model
# say you want λ = 0.7783954
λ <- 0.7783954

# this returns a new "phylo" with all *internal* branch lengths multiplied by λ
phy_lambda <- phytools::rescale(phylogeny, model="lambda", λ)

# now get the λ-transformed variance–covariance matrix
A_lambda <- ape::vcv.phylo(phy_lambda)


# Subset the phylogenetic covariance matrix for the model
A_subset <- A_lambda[selected_species, selected_species]

# Optional: check dimensions
dim(sub_subset)
dim(A_subset)

  # # Try using a standardized correlation matrix instead
  # A_cor <- vcv(phy_subset, corr = TRUE)
  # # And make sure your data2 list and group variable match:
  # sub_subset$phylo_name <- sub_subset$Species
  # data2_list <- list(A_cor = A_cor)

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
saveRDS(test_model_phyl_subset_40_PRIORS,
        file="Models/test_model_phyl_subset_40_PRIORS.rds")




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
Model_mass_mass<-lm(brain_mass_g ~ mass_kg, data =sub_brain)
summary(Model_mass_mass)

#Calculate the residuals of the model and store them in the database
Residuals_body_mass <-residuals(Model_mass_mass)

#attach to df
# If residuals are a matrix, convert to vector
res_vec <- as.vector(Residuals_body_mass)

# Attach to the dataframe
sub_brain$Residuals_body_mass <- res_vec

#subset model w predictors and brain mass
Model_3 <- brm(
  CRB_Final  ~ 0 + Trophic_level + mass_kg + HWI + Residuals_body_mass + (1|gr(phylogeny, cov = A_subset_brain)),
  data = sub_brain,   
  data2 = list(A_subset_brain = A_subset_brain),
  family = bernoulli,
  prior = priors,
  iter = 1e6,
  warmup = 5e5,
  chains = 8,
  cores = 8,
  thin = 100)


summary(Model_3)
#Save model
saveRDS(Model_3,
        file="Models/Model_3.rds")
