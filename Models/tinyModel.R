#This file is attempting to rerun the models with a smaller / faster data subset in order to explore different parameter settings.
#It isn't complete yet.

# install.packages("geiger")   # if you don’t have it
library(geiger)
library(loo)
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
load("Models/Consensus_Tree.Rda")
Bird_data <- read.csv("Data/Bird_data_clean.csv")
sub <- read.csv("/home/sandracd/Downloads/sub.csv") #might need to implement this later

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
set.seed(123)

# Number of species you want to sample - all my tree
n_species <- 20

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
λ <- 0.778

# this returns a new "phylo" with all *internal* branch lengths multiplied by λ
phy_lambda <- rescale(phylogeny, model="lambda", λ)

# now get the λ-transformed variance–covariance matrix
A_lambda <- ape::vcv.phylo(phy_lambda)


# Subset the phylogenetic covariance matrix for the model
A_subset <- A_lambda[selected_species, selected_species]

# Optional: check dimensions
dim(sub_subset)
dim(A_subset)


#Null model with phylogeny only 
Null_model <- brm(
  CRB_Final ~ 1+ (1|gr(phylogeny, cov = A_subset)),
  data = sub_subset,
  data2 = list(A_subset = A_subset),
  family = bernoulli,
  iter = 1000,
  warmup = 500,
  chains = 10,   
  cores = 10,
  thin=100
)

summary(Null_model)
plot(Null_model)
#saveRDS(Null_model,
#        file="Models/Null_model.rds")

priors <- c(
  set_prior("normal(0,1)", class = "b"),
  set_prior("student_t(3,0,2.5"), class = "Intercept")
)

get_prior(
  CRB_Final ~ mass_kg + Trophic_level + HWI,
  data    = sub_subset,
  family  = bernoulli()
)

#Model no phylo for priors
test_model_nophyl <- brm(
  CRB_Final ~ mass_kg + Trophic_level + HWI,
  data = sub_subset,
  #data2 = list(A_subset = A_subset),
  family = bernoulli(link="logit"),
  prior = priors,
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  thin=1,
  seed=2025,

)
pushover(message = "Model no phylo 1000 thin finished")

summary(test_model_nophyl)
plot(test_model_nophyl)


#Define priors based on nophylo models
priors <- c(
  set_prior("normal(-0.84, 1.15)", class = "Intercept"),
  set_prior("normal(0.19, 0.55)", coef = "mass_kg", class = "b"),
  set_prior("normal(0.53, 0.85)", coef = "Trophic_levelHerbivore", class = "b"),
  set_prior("normal(0.71, 1.00)", coef = "Trophic_levelOmnivore", class = "b"),
  set_prior("normal(-0.05, 3.80)", coef = "Trophic_levelScavenger", class = "b"),
  set_prior("normal(0.03, 0.05)", coef = "HWI", class = "b")
  )



#Define priors based on nophylo models
#Option A using normal priors - not using this one
priors <- c(
  set_prior("normal(-0.84, 1.15)", class = "Intercept"),
  set_prior("normal(0.19, 0.55)", coef = "mass_kg", class = "b"),
  set_prior("normal(0.53, 0.85)", coef = "Trophic_levelHerbivore", class = "b"),
  set_prior("normal(0.71, 1.00)", coef = "Trophic_levelOmnivore", class = "b"),
  set_prior("normal(-0.05, 3.80)", coef = "Trophic_levelScavenger", class = "b"),
  set_prior("normal(0.03, 0.05)", coef = "HWI", class = "b")
)

#Option B using cauchy priors. It was not converging
priors <- c(
  set_prior("cauchy(-0.84, 1.15)", class = "Intercept"),
  set_prior("cauchy(0.19, 0.55)", coef = "mass_kg", class = "b"),
  set_prior("cauchy(0.53, 0.85)", coef = "Trophic_levelHerbivore", class = "b"),
  set_prior("cauchy(0.71, 1.00)", coef = "Trophic_levelOmnivore", class = "b"),
  set_prior("cauchy(-0.05, 3.80)", coef = "Trophic_levelScavenger", class = "b"),
  set_prior("cauchy(0.03, 0.05)", coef = "HWI", class = "b")
)

#Option C using t student to promotoe convergence
priors <- c(
  set_prior("cauchy(-0.84, 1.15)", class = "Intercept"),
  set_prior("cauchy(0.19, 0.55)", coef = "mass_log", class = "b"),
  set_prior("cauchy(0.53, 0.85)", coef = "Trophic_levelHerbivore", class = "b"),
  set_prior("cauchy(0.71, 1.00)", coef = "Trophic_levelOmnivore", class = "b"),
  set_prior("cauchy(-0.05, 3.80)", coef = "Trophic_levelScavenger", class = "b"),
  set_prior("cauchy(0.03, 0.05)", coef = "HWI_z", class = "b")
)

#Model 1 with Pagels lambda included from physig calculation
#Model 40 chains 1000 thin WITH PRIORS
test_model_phyl_subset_40_PRIORS <- brm(
  CRB_Final ~ mass_kg + Trophic_level + HWI + (1 | gr(phylogeny, cov = A_subset)),
  data = sub_subset,
  data2 = list(A_subset = A_subset),
  family = bernoulli,
  prior = priors,
  iter = 1e6,
  warmup = 5e5,
  chains = 40,
  cores = 40,
  thin = 10000
)
pushover(message = "Model 40 chains 1000 PRIORS thin finished")
summary(test_model_phyl_subset_40_PRIORS)
plot(test_model_phyl_subset_40_PRIORS)

saveRDS(test_model_phyl_subset_40_PRIORS,
        file="Models/test_model_phyl_subset_40_PRIORS.rds")





#Model comparison for null model vs global model
test_model_phyl_subset_40_PRIORS <- readRDS("Models/test_model_phyl_subset_40_PRIORS.rds")
Null_model <- readRDS("Models/Null_model.rds")


#Redo from here
#fit a light model to compare becuase the chains were not saved before
light_fit <- brm(CRB_Final ~ mass_kg + Trophic_level + HWI +
                   (1 | gr(phylogeny, cov = A_subset)),
                 data = sub_subset,
                 data2 = list(A_subset = A_subset),
                 family = bernoulli,
                 prior = priors,
                 chains = 4, cores = 4, iter = 4000, warmup = 2000,
                 save_pars = save_pars(all = TRUE)   # keep latent draws!
                 )

light_fit <- add_criterion(light_fit, "loo", moment_match = TRUE)
summary(light_fit)


#Waic is working but its not ideal bc values are greater thjan 0.4 and loo is recommended
light_fit  <- add_criterion(light_fit, "waic")
Null_model <- add_criterion(Null_model, "waic")
test_model_phyl_subset_40_PRIORS <- add_criterion(test_model_phyl_subset_40_PRIORS, "waic")
summary(test_model_phyl_subset_40_PRIORS)

#compare full and null
loo::loo_compare(
  x = list(
    full_model = test_model_phyl_subset_40_PRIORS$criteria$waic,
    null_model = Null_model$criteria$waic
  )
)

#compare light and null
loo::loo_compare(
  x = list(
    full_model = light_fit$criteria$waic,
    null_model = Null_model$criteria$waic
  )
)


#loo wasnt working on my full moodel so I was suggested to rescale and refit to avoid divergent transitions..
sub_subset <- sub_subset %>%
  mutate(
    mass_log = scale(log(mass_kg)),
    HWI_z    = scale(HWI)
  )





light_fit <- brm(
  CRB_Final ~ mass_log + Trophic_level + HWI_z +
    (1 | gr(phylogeny, cov = A_subset)),
  data      = sub_subset,
  data2     = list(A_subset = A_subset),
  family    = bernoulli,
  prior     = priors,
  chains    = 4,
  iter      = 4000,
  warmup    = 2000,
  cores     = 4,
  save_pars = save_pars(all = TRUE),
  control   = list(adapt_delta = 0.99)  # helps remove divergences
)


light_fit <- add_criterion(light_fit, "loo", moment_match = TRUE)



ll <- log_lik(light_fit)
summary(as.vector(ll))  # should show no NaN or -Inf

#try with t student


priors <- c(
  set_prior("student_t(3, -0.84, 1.15)", class = "Intercept"),
  set_prior("student_t(3, 0.19, 0.55)", class = "b", coef = "mass_log"),
  set_prior("student_t(3, 0.53, 0.85)", class = "b", coef = "Trophic_levelHerbivore"),
  set_prior("student_t(3, 0.71, 1.00)", class = "b", coef = "Trophic_levelOmnivore"),
  set_prior("student_t(3, -0.05, 3.80)", class = "b", coef = "Trophic_levelScavenger"),
  set_prior("student_t(3, 0.03, 0.05)", class = "b", coef = "HWI_z"),
  set_prior("student_t(3, 0, 2.5)", class = "sd", group = "phylogeny")  # Random effect
)

#Im refitting with t student and see if divergences persist

