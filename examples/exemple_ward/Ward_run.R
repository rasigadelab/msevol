################################################################################
## MSE1 RESISTRACK MODEL
##  A bacterial hospital ecosystem simulation framework
## (c)2018-2024 Jean-Philippe Rasigade, <jean-philippe.rasigade@univ-lyon1.fr>
##  University of Lyon, Centre International de Recherche en Infectiologie
################################################################################
## Author: JP Rasigade, N Lenuzza
################################################################################


setwd("C:/Users/Lenuzza/Desktop/git_msevol/Rmsevol")


##----------------------------------------------------------------------------##
## Options, packages and generic functions
##----------------------------------------------------------------------------##

## WORKING DIRECTORY
## Make sure the the working directory is Rmsevol, i.e. the root folder that
## directly contains 'bin' (with the msevolcli.exe client)

getwd()


## PACKAGES

library(data.table)
library(tidyr)
library(ggplot2)
library(visNetwork)
library(bit64)
library(Rmsevol)
library(R6)
library(gridExtra)


## WARD class
source("../examples/exemple_ward/Ward_class.R")



##----------------------------------------------------------------------------##
## RUN 1 - NO ANTIBIOTIC TREATMENT
##----------------------------------------------------------------------------##

my_params <- list(
  # graph model
  basal_growth = 0.1, 
  basal_survival = 0.99, 
  carrying_capacity = 1E6,
  antibiotic_pressure = 0.1,
  diffusion_patch_to_patch = 2.5E-11,
  # Initial state
  n_patients = 25,
  n_colonized = 0,
  # "Ward" event
  max_admission = 8,
  prop_sensitive_at_admission =0.10,
  proba_discharge = 1/6,
  screening_threshold = 100, 
  ##run param
  run_path="ward_test",
  run_n_steps = 24*10,
  run_n_writes = 1
)

my_ward <- Ward$new(params = my_params)


## Other params
compliance <- 0.9
sensitivity <- 0.9



## RUN WITH EXTENSIVE SCREENING

for(ti in 1:365){
  my_ward$run_step()
  my_ward$discharge_patient()
  my_ward$admit_patients()
  my_ward$screen_patients()
}



## SIMULATES MORE REALISTIC SCREENING

actual_sample_info <- actual_screening(my_ward$sample_info,
                                       period = 7, 
                                       compliance = 0.9, 
                                       sensitivity = 0.9)


## Saving example results

saveRDS(my_ward, 
        "../examples/exemple_ward/ward_NoATB.rds")


saveRDS(actual_sample_info,
        "../examples/exemple_ward/ward__NoATB_screening.rds")






##----------------------------------------------------------------------------##
## RUN 2 - WITH ANTIBIOTIC TREATMENT
##----------------------------------------------------------------------------##

my_params <- list(
  # graph model
  basal_growth = 0.1, 
  basal_survival = 0.99, 
  carrying_capacity = 1E6,
  antibiotic_pressure = 0.1,
  diffusion_patch_to_patch = 2.5E-11,
  # Initial state
  n_patients = 25,
  n_colonized = 0,
  # "Ward" event
  max_admission = 8,
  prop_sensitive_at_admission =0.10,
  proba_discharge = 1/6,
  screening_threshold = 100,
  proba_initiate_antibiotherapy = 0.05,
  duration_antibiotherapy = 5,
  ##run param
  run_path="ward_test",
  run_n_steps = 24*10,
  run_n_writes = 1
)

my_ward <- Ward$new(params = my_params)


## Other params
compliance <- 0.9
sensitivity <- 0.9



## RUN WITH EXTENSIVE SCREENING & monitoring of the number of treated patients 
## per day

n_treatment <-data.table(time = integer(), n_treatment = integer())

for(ti in 1:365){
  my_ward$run_step()
  my_ward$discharge_patient()
  my_ward$admit_patients()
  my_ward$stop_treatment()
  my_ward$start_treatment()
  my_ward$screen_patients()
  
  ## Monitoring of the number of patients under treatment
  n_treatment <- rbind(
    n_treatment,
    data.table(
      time = my_ward$current_time,
      n_treatment= nrow(my_ward$prm$Edge_PatchArchetype_Patch_Multiplicity[
                     source ==  my_ward$info_archid$patchArchetype_pressure, ])
      )
  )
                   
}


## SIMULATES MORE REALISTIC SCREENING

actual_sample_info <- actual_screening(my_ward$sample_info,
                                       period = 7, 
                                       compliance = 0.9, 
                                       sensitivity = 0.9)


## Saving example results

saveRDS(my_ward, 
        "../examples/exemple_ward/ward_ATB.rds")


saveRDS(actual_sample_info,
        "../examples/exemple_ward/ward__ATB_screening.rds")


saveRDS(n_treatment,
        "../examples/exemple_ward/ward__ATB_treatment.rds")

