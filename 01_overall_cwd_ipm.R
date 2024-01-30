###
### Alison C. Ketz 05/31/2023
###
###

###########################################################
### Preliminaries
###########################################################

# rm(list = ls())

# setwd("~/Documents/ipm/cwd_ipm_sim12_chtc")

library(nimble)
library(dplyr)
library(coda)
library(splines)

#set seeds
processnum <- as.numeric(commandArgs(TRUE))
# processnum <- 100
seedy <- processnum + 1000
set.seed(seedy)

###########################################################
### Source summary function for posteriors
###########################################################

# source("support_functions.R")

###############################################################
# Load "Data" i.e. posterior summary stats from data model fit
###############################################################

source("02_load_all_data_to_run.R")

###############################################################
# Generate simulated data
###############################################################

source("04_generate_data.R")

###########################################################
### Setup consts etc for running the model
###########################################################

source("05_prelim_survival.R")

source("06_prelim_foi.R")

###########################################################
### Likelihoods
###########################################################

source("07_distributions.R")

# ###########################################################
# ### Functions for Efficient Calculations
# ###########################################################

source("08_calculations.R")

###########################################################
### Run model
###########################################################

source("09_modelcode.R")

###########################################################
### Run model
###########################################################

source("10_run_model.R")
# source("10_run_model_par.R")

###########################################################
### Post processing
###########################################################

source("11_post_process.R")
