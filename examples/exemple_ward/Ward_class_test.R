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
## MAke sure the the working directory is Rmsevol, i.e. the root folder that
## directly contains 'bin' (with the msevocli.exe client)

getwd()


## PACKAGES

library(data.table)
library(tidyr)
library(ggplot2)
library(visNetwork)
library(bit64)
library(Rmsevol)
library(R6)



## WARD class

source("C:/Users/Lenuzza/Desktop/git_msevol/examples/exemple_ward_WIP/Ward_class.R")



##----------------------------------------------------------------------------##
## Test :  initialization
##----------------------------------------------------------------------------##

init_params <- list(
  basal_growth = 0.1, 
  basal_survival = 0.99, 
  carrying_capacity = 1E6,
  antibiotic_pressure = 0.1,
  diffusion_patch_to_patch = 1E-9,
  n_patients = 10,
  n_colonized = 2
 )

my_ward <- Ward$new(params = init_params)


## check object structure
str(my_ward)


my_ward$current_time == 0
my_ward$info_archid
my_ward$patient_info
my_ward$sample_info
my_ward$params


## check prm
Rmsevol::mse1.plot_graph(my_ward$prm, layout = "fr")
Rmsevol::mse1.diagnose_model(my_ward$prm)
names(my_ward$prm)


## Test with no patient in the system at initial state

init_params <- list(
  basal_growth = 0.1, 
  basal_survival = 0.99, 
  carrying_capacity = 1E6,
  antibiotic_pressure = 0.1,
  diffusion_patch_to_patch = 1E-9,
  n_patients = 0,
  n_colonized = 0
)

my_ward <- Ward$new(params = init_params)


## check object structure
my_ward$current_time == 0
my_ward$info_archid
my_ward$patient_info
my_ward$sample_info
my_ward$params


## check prm
Rmsevol::mse1.plot_graph(my_ward$prm, layout = "fr")
Rmsevol::mse1.diagnose_model(my_ward$prm)
names(my_ward$prm)






##----------------------------------------------------------------------------##
## Test if  "run_step" is functional 
##----------------------------------------------------------------------------##

init_and_mserun_params <- list(
  basal_growth = 0.1, 
  basal_survival = 0.99, 
  carrying_capacity = 1E6,
  antibiotic_pressure = 0.1,
  diffusion_patch_to_patch = 1E-6,
  n_patients = 10,
  n_colonized = 2,
  run_path="ward_test",
  run_n_steps = 24,
  run_n_writes = 10
)

my_ward <- Ward$new(params = init_and_mserun_params)
Rmsevol::mse1.plot_graph(my_ward$prm, layout = "fr")

my_ward$run_step()

my_ward$current_time
Rmsevol::mse1.plot_graph(my_ward$prm, layout = "hr")





##----------------------------------------------------------------------------##
## Tests : treatment initiation
##----------------------------------------------------------------------------##
## Test that treatment is functionnaly attributed, ie. :
##  - patient_info is updated
##  - patients nodes actually change their archetype
##  - the modified prm is running in msevol
##  - sensitive bacteria of patient under treatment decreased and disapears
##----------------------------------------------------------------------------##

## Test 1: start of treatment at Day 0 
## All patients are initially colonised, without treatment

my_params <- list(
  # graph model
  basal_growth = 0.1, 
  basal_survival = 0.99, 
  carrying_capacity = 1E6,
  antibiotic_pressure = 0.1,
  diffusion_patch_to_patch = 0,
  # Initial state
  n_patients = 10,
  n_colonized = 10,
  # "Ward" event
  proba_initiate_antibiotherapy =0.5,
  duration_antibiotherapy = 5,
  run_path="ward_test",
  run_n_steps = 24,
  run_n_writes = 10
)

my_ward <- Ward$new(params = my_params)
str(my_ward) ## start_treatment is decleared in the object

Rmsevol::mse1.plot_graph(my_ward$prm, layout = "fr")
my_ward$patient_info

my_ward$start_treatment()
my_ward$patient_info
Rmsevol::mse1.plot_graph(my_ward$prm, layout = "fr")


## Day 1
my_ward$run_step()
my_ward$current_time
Rmsevol::mse1.plot_graph(my_ward$prm, layout = "fr",
                         edges_type = c("inclusion"))

## Day 2
my_ward$run_step()
my_ward$current_time
Rmsevol::mse1.plot_graph(my_ward$prm, layout = "fr",
                         edges_type = c("inclusion"))

## Day 3
my_ward$run_step()
my_ward$current_time
Rmsevol::mse1.plot_graph(my_ward$prm, layout = "fr",
                         edges_type = c("inclusion"))
Rmsevol::mse1.plot_graph(my_ward$prm, layout = "fr")


## Test 2 : start of treatment at Day 0 
## Only one patient is initially colonised
## Treatment is initiated in 100% of patients (prob = 1)
## Test that refill is functionnal, i.e. the sensitive cell is still 
## included in the graph at the end of a run, even if no patient are 
## colonised anymore.


my_params <- list(
  # graph model
  basal_growth = 0.1, 
  basal_survival = 0.99, 
  carrying_capacity = 1E6,
  antibiotic_pressure = 0.1,
  diffusion_patch_to_patch = 0,
  # Initial state
  n_patients = 10,
  n_colonized = 1,
  # "Ward" event
  proba_initiate_antibiotherapy =1,
  duration_antibiotherapy = 5,
  run_path="ward_test",
  run_n_steps = 24,
  run_n_writes = 10
)

my_ward <- Ward$new(params = my_params)
str(my_ward) ## start_treatment is decleared in the object

my_ward$patient_info
Rmsevol::mse1.plot_graph(my_ward$prm, layout = "fr")
my_ward$patient_info

my_ward$start_treatment()
my_ward$patient_info
Rmsevol::mse1.plot_graph(my_ward$prm, layout = "fr")


## Day 1
my_ward$run_step()
my_ward$current_time
Rmsevol::mse1.plot_graph(my_ward$prm, layout = "fr",
                         edges_type = c("inclusion"))

## Day 2
my_ward$run_step()
my_ward$current_time
Rmsevol::mse1.plot_graph(my_ward$prm, layout = "fr",
                         edges_type = c("inclusion"))

## Day 3
my_ward$run_step()
my_ward$current_time
Rmsevol::mse1.plot_graph(my_ward$prm, layout = "fr",
                         edges_type = c("inclusion"))






##----------------------------------------------------------------------------##
## Tests : treatment stop
##----------------------------------------------------------------------------##

## Test 1 : check that treatment are actually stopped both in patients & in 
## the prm graph model
## Note in this test, the functionnality of the prm model is not tested 
## (no msevol run)

my_params <- list(
  # graph model
  basal_growth = 0.1, 
  basal_survival = 0.99, 
  carrying_capacity = 1E6,
  antibiotic_pressure = 0.1,
  diffusion_patch_to_patch = 0,
  # Initial state
  n_patients = 10,
  n_colonized = 0,
  # "Ward" event
  proba_initiate_antibiotherapy =0.3,
  duration_antibiotherapy = 3,
  run_path="ward_test",
  run_n_steps = 24,
  run_n_writes = 10
)

my_ward <- Ward$new(params = my_params)
my_ward$start_treatment()
my_ward$patient_info

for(i in 1:3){
  my_ward$update_time()
  cat("\n****",  my_ward$current_time, "*****************************************\n")
  print(my_ward$patient_info)
  my_ward$stop_treatment()
  print(my_ward$patient_info)
  Rmsevol::mse1.plot_graph(my_ward$prm, layout = "fr",
                           edges_type = c("inclusion"))
}



## Test 2 : Run with msevol to test the consistency of the msevol run

my_params <- list(
  # graph model
  basal_growth = 0.1, 
  basal_survival = 0.99, 
  carrying_capacity = 1E6,
  antibiotic_pressure = 0.09,
  diffusion_patch_to_patch = 0,
  # Initial state
  n_patients = 10,
  n_colonized = 10,
  # "Ward" event
  proba_initiate_antibiotherapy =0.5,
  duration_antibiotherapy = 4,
  run_path="ward_test",
  run_n_steps = 24,
  run_n_writes = 10
)

my_ward <- Ward$new(params = my_params)
my_ward$start_treatment()
my_ward$patient_info
my_ward$current_time

for(i in 1:5){
  my_ward$run_step()
  cat("\n****",  my_ward$current_time, "*****************************************\n")
  my_ward$stop_treatment()
  print(my_ward$patient_info)
  Rmsevol::mse1.plot_graph(my_ward$prm, layout = "fr",
                           edges_type = c("inclusion"))

}


##----------------------------------------------------------------------------##
## Tests : patient_admission
##----------------------------------------------------------------------------##

my_params <- list(
  # graph model
  basal_growth = 0.1, 
  basal_survival = 0.99, 
  carrying_capacity = 1E6,
  antibiotic_pressure = 0.1,
  diffusion_patch_to_patch = 0,
  # Initial state
  n_patients = 0,
  n_colonized = 0,
    # "Ward" event
  proba_initiate_antibiotherapy =0,
  duration_antibiotherapy = 5,
  max_admission = 10,
  prop_sensitive_at_admission =0.5, 
  ##run param
  run_path="ward_test",
  run_n_steps = 24,
  run_n_writes = 10
)

my_ward <- Ward$new(params = my_params)
my_ward$patient_info
Rmsevol::mse1.plot_graph(my_ward$prm, layout = "fr")

my_ward$admit_patients()
my_ward$patient_info
Rmsevol::mse1.plot_graph(my_ward$prm, layout = "fr")
my_ward$run_step()
Rmsevol::mse1.plot_graph(my_ward$prm, layout = "fr")

my_ward$admit_patients()
my_ward$patient_info
Rmsevol::mse1.plot_graph(my_ward$prm, layout = "fr")
my_ward$run_step()
Rmsevol::mse1.plot_graph(my_ward$prm, layout = "fr")

my_ward$admit_patients()
my_ward$patient_info
Rmsevol::mse1.plot_graph(my_ward$prm, layout = "fr")
my_ward$run_step()
Rmsevol::mse1.plot_graph(my_ward$prm, layout = "fr")




##----------------------------------------------------------------------------##
## Tests : patient_discharge
##----------------------------------------------------------------------------##
source("C:/Users/Lenuzza/Desktop/git_msevol/examples/exemple_ward_WIP/Ward_class.R")
my_params <- list(
  # graph model
  basal_growth = 0.1, 
  basal_survival = 0.99, 
  carrying_capacity = 1E6,
  antibiotic_pressure = 0.1,
  diffusion_patch_to_patch = 0,
  # Initial state
  n_patients = 50,
  n_colonized = 30,
  # "Ward" event
  proba_initiate_antibiotherapy =0,
  duration_antibiotherapy = 5,
  max_admission = 10,
  prop_sensitive_at_admission =0.5,
  proba_discharge = 0.1,
  ##run param
  run_path="ward_test",
  run_n_steps = 24,
  run_n_writes = 10
)

my_ward <- Ward$new(params = my_params)
my_ward$patient_info
nrow(my_ward$patient_info[!is.na(discharge_day),])

my_ward$discharge_patient()
my_ward$patient_info
my_ward$patient_info[!is.na(discharge_day),]
nrow(my_ward$patient_info[is.na(discharge_day),])==nrow(my_ward$prm$Vertex_Patch)
all(my_ward$prm$Vertex_Patch$id %in% my_ward$patient_info[is.na(discharge_day),]$idx)

my_ward$run_step()
my_ward$discharge_patient()
my_ward$patient_info
my_ward$patient_info[!is.na(discharge_day),]
nrow(my_ward$patient_info[is.na(discharge_day),])==nrow(my_ward$prm$Vertex_Patch)
all(my_ward$prm$Vertex_Patch$id %in% my_ward$patient_info[is.na(discharge_day),]$idx)
