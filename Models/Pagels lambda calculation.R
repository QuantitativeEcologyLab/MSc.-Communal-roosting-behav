# install.packages("geiger")   # if you don’t have it
library(geiger)
library(dplyr)
library('phytools')
library('ape')
library('phylobase')
library('brms')
library('ggplot2')
library('colorspace')
library(dplyr)
library("tidyverse")
library("viridis")
library('ggtree')
library("Rcpp")
#install.packages("phytools")
library(phytools)
#install.packages("pushoverr")
library(pushoverr)
set_pushover_user(user = "usdw6uw28zwd2496rr2dnni3efz4ds")
set_pushover_app(token = "aa21xm4a7vi9dcfv9d46p1qypg8miu")

#load data
#windows
load("C:/Users/sandracd/OneDrive - UBC (1)/PhD proposal/Chapter 2/phylo_data/Consensus_Tree.Rda")
Bird_data <- read.csv("C:/Users/sandracd/OneDrive - UBC (1)/PhD proposal/Chapter 2/Brain size data and sources/Bird_data.csv")
#Linux
load("/home/sandracd/Downloads/Consensus_Tree.Rda")
Bird_data <- read.csv("/home/sandracd/Downloads/Bird_data.csv")
sub <- read.csv("/home/sandracd/Downloads/sub.csv")

#Wrangle dataset
Bird_data$Species <- gsub("_+$", "", Bird_data$Species)

#Subset the variables that will go into the model
Bird_data_clean<- Bird_data %>%
  select(Species, Family, Order, Mass, Trophic_level, HWI, brain_mass_g, CRB_Final)

#Remove NA in CRB and add mass kg
Bird_data_clean <- Bird_data %>%
  select(Species,Family, Order, Mass, Trophic_level, HWI, brain_mass_g, CRB_Final) %>%
  drop_na(CRB_Final)%>%
  mutate(mass_kg = Mass / 1000)


#Create a list with Species and their matching phylogeny to run the model
sub <- Bird_data_clean[Bird_data_clean$Species %in% phylogeny$tip.label,] 
sub$phylogeny <- sub$Species


#Calculate lambda values
#Trials with la subset spp make sure to set seed for repproducibility. I will do entire tree now
# Set seed for reproducibility
#set.seed(123)

# Number of species you want to sample - all my tree
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
phy_lambda <- rescale(phylogeny, model="lambda", λ)

# now get the λ-transformed variance–covariance matrix
A_lambda <- ape::vcv.phylo(phy_lambda)


# Subset the phylogenetic covariance matrix for the model
A_subset <- A_lambda[selected_species, selected_species]

# Optional: check dimensions
dim(sub_subset)
dim(A_subset)


#Test model Pagels lambda above from physig
test_model_phyl_subset <- brm(
  CRB_Final ~ mass_kg + Trophic_level + HWI + (1 | gr(phylogeny, cov = A_subset)),
  data = sub_subset,
  data2 = list(A_subset = A_subset),
  family = bernoulli,
  iter = 1e5,
  warmup = 5e4,
  chains = 8,   
  cores = 8,
  thin=1000
)

summary(test_model_phyl_subset)
plot(test_model_phyl_subset)




test_model_phyl_subset_100 <- brm(
  CRB_Final ~ mass_kg + Trophic_level + HWI + (1 | gr(phylogeny, cov = A_subset)),
  data = sub_subset,
  data2 = list(A_subset = A_subset),
  family = bernoulli,
  iter = 1e5,
  warmup = 5e4,
  chains = 8,   
  cores = 8,
  thin=100
)


summary(test_model_phyl_subset_100)
plot(test_model_phyl_subset_100)





#40 chains and 1000 thin#
test_model_phyl_subset_40 <- brm(
  CRB_Final ~ mass_kg + Trophic_level + HWI + (1 | gr(phylogeny, cov = A_subset)),
  data = sub_subset,
  data2 = list(A_subset = A_subset),
  family = bernoulli,
  iter = 1e5,
  warmup = 5e4,
  chains = 40,   
  cores = 40,
  thin=1000
)
pushover(message = "Model 40 chains 1000 thin finished")
summary(test_model_phyl_subset_40)
plot(test_model_phyl_subset_40)

#40 chains and 100 thin#
test_model_phyl_subset_40_2 <- brm(
  CRB_Final ~ mass_kg + Trophic_level + HWI + (1 | gr(phylogeny, cov = A_subset)),
  data = sub_subset,
  data2 = list(A_subset = A_subset),
  family = bernoulli,
  iter = 1e5,
  warmup = 5e4,
  chains = 40,   
  cores = 40,
  thin=100
)
pushover(message = "Model 40 chains 100 thin finished")

summary(test_model_phyl_subset_40_2)
plot(test_model_phyl_subset_40_2)




#model no phylo for priors#
test_model_nophyl <- brm(
  CRB_Final ~ mass_kg + Trophic_level + HWI,
  data = sub_subset,
  #data2 = list(A_subset = A_subset),
  family = bernoulli,
  iter = 1e5,
  warmup = 5e4,
  chains = 40,   
  cores = 40,
  thin=1000
)
pushover(message = "Model no phylo 1000 thin finished")

summary(test_model_nophyl)
plot(test_model_nophyl)


#Define priors based on nophylo models
priors <- c(
  set_prior("normal(-0.84, 0.23)", class = "Intercept"),
  set_prior("normal(0.19, 0.11)", coef = "mass_kg", class = "b"),
  set_prior("normal(0.53, 0.17)", coef = "Trophic_levelHerbivore", class = "b"),
  set_prior("normal(0.71, 0.20)", coef = "Trophic_levelOmnivore", class = "b"),
  set_prior("normal(-0.05, 0.76)", coef = "Trophic_levelScavenger", class = "b"),
  set_prior("normal(0.03, 0.01)", coef = "HWI", class = "b")
  )

#Test model 40 chains 1000 thin WITH PRIORS
test_model_phyl_subset_40_PRIORS <- brm(
  CRB_Final ~ mass_kg + Trophic_level + HWI + (1 | gr(phylogeny, cov = A_subset)),
  data = sub_subset,
  data2 = list(A_subset = A_subset),
  family = bernoulli,
  prior = priors,
  iter = 1e5,
  warmup = 5e4,
  chains = 40,
  cores = 40,
  thin = 1000
)
pushover(message = "Model 40 chains 1000 PRIORS thin finished")
summary(test_model_phyl_subset_40_PRIORS)
plot(test_model_phyl_subset_40_PRIORS)










# Subset the phylogenetic covariance matrix FOR BROWNIAN. WE DONT NEED THIS
A <- ape::vcv.phylo(phylogeny)

A_subset_Brownian <- A[selected_species, selected_species]

#Test model Brownian - Pagel 1
test_model_phyl_subset_Brownian <- brm(
  CRB_Final ~ Mass + Trophic_level + HWI + (1 | gr(phylogeny, cov = A_subset_Brownian)),
  data = sub_subset,
  data2 = list(A_subset_Brownian = A_subset_Brownian),
  family = bernoulli,
  iter = 1e4,
  warmup = 5e3,
  chains = 8,   
  cores = 8     
)

summary(test_model_phyl_subset_Brownian)
plot(test_model_phyl_subset_Brownian)


  
  