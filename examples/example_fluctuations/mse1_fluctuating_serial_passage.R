################################################################################
## MSE1 RESISTRACK MODEL                                                        
##  A bacterial hospital ecosystem simulation framework                          
## (c)2018-2024 Jean-Philippe Rasigade, <jean-philippe.rasigade@univ-lyon1.fr>
##  University of Lyon, Centre International de Recherche en Infectiologie       
################################################################################
## Example of simulations with regular "update" of parameters or configuration:
## SERIAL PASSAGE OF BACTERIA TO MONITOR THE LOSS A RESISTANT PLASMID
##
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


## SOURCING FUNCTIONS

source("msvr_fun/cfg_funcs v1.0.R")
source("msvr_fun/parse_funcs.R")
source("msvr_fun/helper_funcs.R")




#----------------------------------------------------------------------------##
## BASE EXEMPLE - WARNING !!!
##----------------------------------------------------------------------------##


## 0. global parameters
csvres_folder <- "fluctuation_serial_passage"
T_in_steps <- 24
N_periods <- 50
n_steps <- 1
n_writes <- as.integer(round(T_in_steps/n_steps))




## 1. Model initialization


prm <- emptymodel() %>%
  chromosomeArchetype(0, fitness = 0.5, survival = 1.0) %>% 
  chromosome(0, 0) %>%
  plasmidArchetype(0, loss = 1E-3, transfer = 1E-4, max_count = 1, fitness = 0.95)%>%
  plasmid(0,0)%>%
  cellArchetype(0) %>% cell(0, 0) %>%
  patchArchetype(0, 1e9, c(0.0, 0.0)) %>% patch(0, 0) %>%
  assign("Chromosome", 0, "Cell", 0) %>%
  assign("Plasmid", 0, "Cell", 0) %>%
  assign("Cell", 0, "Patch", 0, 1e6)%>%
  plasmidRange(0,0)

mse1.diagnose_model(prm)
mse1.plot_graph(prm, layout = "hr")


## 2. First iteration

res <- mse_run(prm, 
               path = csvres_folder, 
               n_steps = n_steps, n_writes = n_writes, 
               s1 = sample(1e9, 1))

mse1.cleanup_csv(csvres_folder)
rm(prm)


## 2. Iterating on the number of period

for(per_i in 2:N_periods){
  
  ## extract model from res
  new_prm <- mse1.extract_prm_from_mseres(res, 
                                          T_in_steps*(per_i-1))
  
  # Perform dilution 1:1000
  new_prm$Edge_Cell_Patch_Multiplicity[, multiplicity:=floor(multiplicity/1000)]
  
  ## run a new C++ simulation
  new_res <- mse_run(new_prm, 
                     path = csvres_folder, 
                     n_steps = n_steps, n_writes = n_writes, 
                     s1 = sample(1e9, 1))
  mse1.cleanup_csv(csvres_folder)
  
  ## update time for new_res
  new_res <- mse1.shift_step(new_res, T_in_steps*(per_i-1))
  res <-  mse1.merge_res(res, new_res)
  
  ## removed
  rm(new_prm)
  rm(new_res)
  
}


## 3. Analysis

d_cells <- count_from("Cell", "Patch", res)
ggplot(data = d_cells, 
       mapping = aes(x =step, y = multiplicity, 
                     group = Cell, color = as.factor(Cell)))+
  geom_line(linewidth = 1.2)+
  theme_bw()+
  ggtitle("Serial passage")+
  scale_y_log10()


## WARNING : PB with transitive count at duplicated time,

d_plasmids <- count_from("Plasmid", "Cell", res) %>%
  count_in("Patch", res)

ggplot(data = d_plasmids, 
       mapping = aes(x =step, y = multiplicity, 
                     group = Plasmid, color = as.factor(Plasmid)))+
  geom_line(linewidth = 1.2)+
  ggtitle("Serial passage")+
  theme_bw()+
  xlim(c(0,50))+
  geom_point(color = "blue")


## so you :
## 1) need to select which Time you want to keep and suppress the other
##    before merging res
## or, alternatively,
## 2) performed analysis individually for each period and merge the resulting
##    data instead of raw results


#----------------------------------------------------------------------------##
## BASE EXEMPLE - SOLUTION 1
##----------------------------------------------------------------------------##

## 0. global parameters
csvres_folder <- "fluctuation_serial_passage"
T_in_steps <- 24
N_periods <- 50
n_steps <- 1
n_writes <- as.integer(round(T_in_steps/n_steps))




## 1. Model initialization


prm <- emptymodel() %>%
  chromosomeArchetype(0, fitness = 0.5, survival = 1.0) %>% 
  chromosome(0, 0) %>%
  plasmidArchetype(0, loss = 1E-3, transfer = 1E-4, max_count = 1, fitness = 0.95)%>%
  plasmid(0,0)%>%
  cellArchetype(0) %>% cell(0, 0) %>%
  patchArchetype(0, 1e9, c(0.0, 0.0)) %>% patch(0, 0) %>%
  assign("Chromosome", 0, "Cell", 0) %>%
  assign("Plasmid", 0, "Cell", 0) %>%
  assign("Cell", 0, "Patch", 0, 1e6)%>%
  plasmidRange(0,0)

mse1.diagnose_model(prm)
mse1.plot_graph(prm, layout = "hr")


## 2. First iteration

res <- mse_run(prm, 
               path = csvres_folder, 
               n_steps = n_steps, n_writes = n_writes, 
               s1 = sample(1e9, 1))

mse1.cleanup_csv(csvres_folder)
rm(prm)


## 2. Iterating on the number of period

for(per_i in 2:N_periods){
  
  ## extract model from res
  new_prm <- mse1.extract_prm_from_mseres(res, 
                                          T_in_steps*(per_i-1))
  
  # Perform dilution 1:1000
  new_prm$Edge_Cell_Patch_Multiplicity[, multiplicity:=floor(multiplicity/1000)]
  
  ## run a new C++ simulation
  new_res <- mse_run(new_prm, 
                     path = csvres_folder, 
                     n_steps = n_steps, n_writes = n_writes, 
                     s1 = sample(1e9, 1))
  mse1.cleanup_csv(csvres_folder)
  
  ## update time for new_res
  new_res <- mse1.suppress_step(new_res, step_to_suppress = 0)
  new_res <- mse1.shift_step(new_res, T_in_steps*(per_i-1))
  res <-  mse1.merge_res(res, new_res)
  
  ## removed
  rm(new_prm)
  rm(new_res)
  
}


## 3. Analysis

d_cells <- count_from("Cell", "Patch", res)
ggplot(data = d_cells, 
       mapping = aes(x =step, y = multiplicity, 
                     group = Cell, color = as.factor(Cell)))+
  geom_line(linewidth = 1.2)+
  theme_bw()+
  ggtitle("Serial passage")+
  scale_y_log10()


## No more duplicated time
d_plasmids <- count_from("Plasmid", "Cell", res) %>%
  count_in("Patch", res)

ggplot(data = d_plasmids[multiplicity>0, ], 
       mapping = aes(x =step, y = multiplicity, 
                     group = Plasmid, color = as.factor(Plasmid)))+
  geom_line(linewidth = 1.2)+
  ggtitle("Serial passage : Plasmid_count")+
  #xlim(c(0,50))+
  #geom_point(color = "blue")+
  theme_bw()


