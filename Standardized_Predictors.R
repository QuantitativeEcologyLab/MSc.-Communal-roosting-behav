# install.packages("geiger")   # if you don’t have it
#install.packages("phytools")
#install.packages("pushoverr")
install.packages("brms")
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

  # Try using a standardized correlation matrix instead
  A_cor <- vcv(phy_subset, corr = TRUE)
  # And make sure your data2 list and group variable match:
  sub_subset$phylo_name <- sub_subset$Species
  data2_list <- list(A_cor = A_cor)

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

start <- proc.time()["elapsed"]
#Model no phylo for priors (updated settings)
test_model_nophyl <- brm(
  CRB_Final ~ mass_z + Trophic_level + HWI_z,
  data = sub_subset,
  #data2 = list(A_subset = A_subset),
  family = bernoulli,
  iter = 2000,
  warmup = 1000,
  chains = 4,   
  cores = 4,
  thin= 1
)
elapsed <- (proc.time()["elapsed"] - start)/60 # minutes
print(glue("No PHylo Model for Priors of {n_species} species finished in {elapsed} min."))

summary(test_model_nophyl)
plot(test_model_nophyl)

#Define priors based on nophylo models
priors <- c(
  set_prior("normal( 0   , 2.5)", class = "Intercept"), #intercept
  set_prior("normal( 0.21, 0.21)", coef = "mass_z", class = "b"),
  set_prior("normal( 0.53, 0.17)", coef = "HWI_z", class = "b"),
  set_prior("normal( 0.70, 0.20)", coef = "Trophic_levelOmnivore", class = "b"),
  set_prior("normal(-0.04, 0.77)", coef = "Trophic_levelScavenger", class = "b"),
  set_prior("normal( 0.53, 0.17)", coef = "Trophic_levelHerbivore", class = "b"),
  set_prior("exponential(1)", group = "phylogeny", class = "sd")
)

#Model w Phylo, Priors (using standardized phylogenetic matrix)
start <- proc.time()["elapsed"]
test_model_phyl_subset_40_PRIORS <- brm(
#  CRB_Final ~ mass_z + Trophic_level + HWI_z + (1 | gr(phylogeny, cov = A_cor)),
  CRB_Final ~ mass_z + HWI_z + (1 | gr(phylogeny, cov = A_cor)),
  data = sub_subset,
  data2 = data2_list, #list(A_subset = A_subset),
  family = bernoulli(link = "logit"),
#  prior = priors,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  cores = 4,
  thin = 1
)
elapsed <- (proc.time()["elapsed"] - start)/60 # minutes
print(glue("With PHylo Model and Priors of {n_species} species finished in {elapsed} min."))
pushover(message = glue("With PHylo Model and Priors of {n_species} species finished in {elapsed} min."))
summary(test_model_phyl_subset_40_PRIORS)
plot(test_model_phyl_subset_40_PRIORS)

saveRDS(test_model_phyl_subset_40_PRIORS,
        file="Models/test_model_phyl_subset_40_PRIORS.rds")

priors <- c(
  set_prior("student_t(3, -0.3, 0.2)", class = "Intercept"), #intercept
  set_prior("normal( 0.2, 0.3)", coef = "mass_z", class = "b"),
  set_prior("normal( 0.2, 0.3)", coef = "HWI_z", class = "b"),
  set_prior("normal( 0.7, 0.3)", coef = "Trophic_levelOmnivore", class = "b"),
  set_prior("normal( 0.5, 0.2)", coef = "Trophic_levelScavenger", class = "b"),
  set_prior("normal( 0, 0.8)", coef = "Trophic_levelHerbivore", class = "b"),
  set_prior("exponential(1)", group = "phylogeny", class = "sd")
)

#Model w Phylo, Default Priors
start <- proc.time()["elapsed"]
model_phyl_default_priors <- brm(
  CRB_Final ~ mass_z + Trophic_level + HWI_z + (1 | gr(phylogeny, cov = A_subset)),
  data = sub_subset,
  data2 = list(A_subset = A_subset),
  family = bernoulli,
  iter = 2000,
  warmup = 1000,
  chains = 4,
  cores = 4,
  thin = 1
)
elapsed <- (proc.time()["elapsed"] - start)/60 # minutes
print(glue("With PHylo Model and Priors of {n_species} species finished in {elapsed} min."))
pushover(message = glue("With PHylo Model and Priors of {n_species} species finished in {elapsed} min."))
summary(model_phyl_default_priors)

#Model comparison for null model vs global model
test_model_phyl_subset_40_PRIORS <- readRDS("Models/test_model_phyl_subset_40_PRIORS.rds")
Null_model <- readRDS("Models/Null_model.rds")


waic_null <- waic(Null_model)
waic_full <- waic(test_model_phyl_subset_40_PRIORS)

print(waic_null)
print(waic_full)

waic_cmp <- loo_compare(waic_null, waic_full)
print(waic_cmp)
