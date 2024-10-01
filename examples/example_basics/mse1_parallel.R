################################################################################
## MSE1 RESISTRACK MODEL                                                        
##  A bacterial hospital ecosystem simulation framework                          
## (c)2018-2024 Jean-Philippe Rasigade, <jean-philippe.rasigade@univ-lyon1.fr>
##  University of Lyon, Centre International de Recherche en Infectiologie       
################################################################################
## Complete model exemplifying how to parameter each admissible entity and edge
## Author: JP Rasigade, N Lenuzza 
################################################################################


##----------------------------------------------------------------------------##
## Options, packages and generic functions
##----------------------------------------------------------------------------##

## WORKING DIRECTORY
## Set the working directory to the root folder that directly contains 'bin' 
## (with the msevocli.exe client) and 'msvr_fun' (which contains the R 
## functions needed to use the simulator)

setwd("C:/Users/Lenuzza/Desktop/clean_msevol/msevolr") 


## PACKAGES

library(data.table)
library(tidyr)
library(ggplot2)
library(visNetwork)
library(bit64)

library(foreach)
library(doParallel)


## SOURCING FUNCTIONS

source("msvr_fun/cfg_funcs v1.0.R")
source("msvr_fun/parse_funcs.R")
source("msvr_fun/helper_funcs.R")



##----------------------------------------------------------------------------##
## PARALLEL EXECUTION - BASIC BIRTH AND DEATH
##----------------------------------------------------------------------------##


N_repetition <- 10


## Initialize cluster
n_cores <- 5
n_cores <- min(c(n_cores, detectCores()[1]-1)) ## !!
cl <- makeCluster(n_cores)
registerDoParallel(cl)

## Parallel execution of the same model

results <- foreach(i = 1:10, 
                   .packages = c("data.table", "tidyr")) %dopar% {
  
  
  ## initialize a model
  prm <- emptymodel() %>%
    chromosomeArchetype(0, fitness =0.1, survival =0.95) %>% 
    chromosome(0, 0) %>%
    cellArchetype(0) %>% cell(0, 0) %>%
    patchArchetype(0, 1E9, c(0.0, 0.0)) %>% patch(0, 0) %>%
    assign("Chromosome", 0, "Cell", 0) %>%
    assign("Cell", 0, "Patch", 0, 1E2)
  
  
  ## 
  csvres_folder_i <- paste("parallel_birth_", i, sep = "") 
  res_i <- mse_run(prm, csvres_folder_i, 1, 500, s1 = sample(1E9,1))
  
  # if needed, cleaning the csv results
  mse1.cleanup_csv(csvres_folder_i)
  
  res_i
  
  
  
}



## stop cluster
registerDoSEQ()
stopCluster(cl)


## Analysis model

d_cells <- rbindlist(lapply(results,
                            function(mylist){
                              count_from("Cell", "Patch", mylist)
                            }),
                     idcol = "repetition")


ggplot(data = d_cells,
       mapping = aes(x =step, y = multiplicity, 
                     group = interaction(Cell, Patch, repetition),
                     color = as.factor(Cell)))+
  geom_line()+
  theme_bw()+
  scale_y_log10()






##----------------------------------------------------------------------------##
## PARALLEL EXECUTION - Mutual exclusion of cells with identical features
##----------------------------------------------------------------------------##


N_repetition <- 10


## Initialize cluster
n_cores <- 5
n_cores <- min(c(n_cores, detectCores()[1]-1)) ## !!
cl <- makeCluster(n_cores)
registerDoParallel(cl)

## Parallel execution of the same model

results <- foreach(i = 1:10, 
                   .packages = c("data.table", "tidyr")) %dopar% {
                     
                     
                     ## initialize a model
                     prm <- emptymodel() %>%
                       chromosomeArchetype(0, fitness = 0.3, survival =0.95) %>% 
                       chromosome(0, 0) %>% 
                       cellArchetype(0) %>% cell(0, 0) %>%
                       cellArchetype(1) %>% cell(1,1) %>%
                       patchArchetype(0, 1E6, c(0.0, 0.0)) %>% patch(0, 0) %>%
                       assign("Chromosome", 0, "Cell", 0) %>%
                       assign("Chromosome", 0, "Cell", 1) %>%
                       assign("Cell", 0, "Patch", 0, 1E3)%>%
                       assign("Cell", 1, "Patch", 0, 1E3)
                     
                     
                     ## 
                     csvres_folder_i <- paste("parallel_exclusion_", i, sep = "") 
                     res_i <- mse_run(prm, csvres_folder_i, 1, 2000, s1 = sample(1E9,1))
                     
                     # if needed, cleaning the csv results
                     mse1.cleanup_csv(csvres_folder_i)
                     
                     res_i
 
                   }



## stop cluster
registerDoSEQ()
stopCluster(cl)


## Analysis model

d_cells <- rbindlist(lapply(results,
                            function(mylist){
                              count_from("Cell", "Patch", mylist)
                            }),
                     idcol = "repetition")


ggplot(data = d_cells,
       mapping = aes(x =step, y = multiplicity, 
                     group = interaction(Cell, Patch, repetition),
                     color = as.factor(Cell)))+
  geom_line()+
  theme_bw()+
  scale_y_log10()

