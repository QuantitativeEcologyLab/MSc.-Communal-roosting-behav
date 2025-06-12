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

#load data
load("C:/Users/sandracd/OneDrive - UBC (1)/PhD proposal/Chapter 2/phylo_data/Consensus_Tree.Rda")
Bird_data <- read.csv("C:/Users/sandracd/OneDrive - UBC (1)/PhD proposal/Chapter 2/Brain size data and sources/Bird_data.csv")

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

summary(test_model_phyl_subset_40)
plot(test_model_phyl_subset_40)












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


  
  