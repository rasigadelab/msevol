################################################################################
## MSE1 RESISTRACK MODEL                                                        
##  A bacterial hospital ecosystem simulation framework                          
## (c)2018-2024 Jean-Philippe Rasigade, <jean-philippe.rasigade@univ-lyon1.fr>
##  University of Lyon, Centre International de Recherche en Infectiologie       
################################################################################
## Illustration of RESISTANCE and ANTIBIOTIC STRESS
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


################################################################################
## RELATIONSHIP BETWEEN CELL NET GROWTH AND STRESS
################################################################################
## In msevol, antibiotic stress is not defined by a specific antibiotic 
## concentration in the medium, but rather by the probability of killing a 
## fully sensitive bacterium within a single time step. As a result, its effect
## on bacterial dynamics depends heavily on the inherent characteristics of 
## the bacterium, particularly its basal growth rate.
################################################################################


## Model configuration ------------------------------

time_kill_experiment_model <- function(basal_growth_probability = 0.5, 
                                       basal_survival_probability = 1.0,
                                       antibiotic_stresses =seq(0,0.5, by = 0.1),
                                       patch_capacity = 1E9,
                                       initial_popsize =1E9
){
  
  prm <- emptymodel() %>%
    chromosomeArchetype(0, 
                        fitness = basal_growth_probability,
                        survival = basal_survival_probability) %>% 
    chromosome(0, 0) %>% 
    cellArchetype(0) %>% cell(0, 0) %>%  
    assign("Chromosome", 0, "Cell", 0)
  
  for(i in 1:length(antibiotic_stresses)){
    prm <- prm %>%
      patchArchetype(i-1, patch_capacity, c(antibiotic_stresses[i], 0.0)) %>% 
      patch(i-1, i-1) %>% 
      assign("Cell", 0, "Patch", i-1, initial_popsize)
  }
    
  return(prm)
  
}


time_kill_experiment <- function(mod_args, run_args){
  
  
  
  ## model building 
  prm <- time_kill_experiment_model(mod_args$basal_growth_probability, 
                                    mod_args$basal_survival_probability,
                                    mod_args$antibiotic_stresses,
                                    mod_args$patch_capacity,
                                    mod_args$initial_popsize)
  
  mse1.plot_graph(prm, layout = "fr")
  
  
  
  ## model runing
  res <- mse_run(prm, 
                 path = run_args$path,
                 n_steps = run_args$n_steps,
                 n_writes = run_args$n_writes,
                 s1 = sample(1E9, 1))
  
   
  ## if needed, cleaning the csv results
  if(run_args$cleanup_csv){ mse1.cleanup_csv(run_args$path)}
  
  # if needed, save the R result (res) as a RDS object
  if(run_args$save_r_res){
    saveRDS(object = res, 
            file = paste(run_args$r_res_path, 
                                       "/", run_args$path, ".rds", sep = ""))
  }
  
 
  ## 
  return(res)
  
}



## Simulations A -----------------------------------------------------------------

pressures <- c(0.1, 0.2, 0.22, 0.23, 0.3, 0.4)


mse1.plot_graph(
  time_kill_experiment_model(basal_growth_probability = 0.3, 
                             basal_survival_probability = 0.99,
                             antibiotic_stresses = pressures,
                             patch_capacity = 1E9,
                             initial_popsize =1E3),
  layout = "fr")




res_A <- time_kill_experiment(mod_args = 
                               list(basal_growth_probability = 0.3, 
                                    basal_survival_probability = 0.99,
                                    antibiotic_stresses = pressures,
                                    patch_capacity = 1E9,
                                    initial_popsize =1E3),
                               run_args = 
                               list(path = "time_killed_experiment_A",
                                    n_steps = 1,
                                    n_writes = 1000,
                                    cleanup_csv = TRUE,
                                    save_r_res =TRUE,
                                    r_res_path = "example_resistance"))


cells_dynamics_A <- count_from("Cell", "Patch", res_A)
cells_dynamics_A[, stress:= pressures[Patch+1]]

plotA <- ggplot(data = cells_dynamics_A,
       mapping = aes(x =step, y = multiplicity, 
                     group = stress,
                     color = as.factor(stress)))+
  geom_line(linewidth =1.2)+
  xlab("Time")+
  ylab("Cell count")+
  theme_bw()+
  scale_y_log10()+
  ggtitle("Growth probability = 0.3")



## Simulations B -----------------------------------------------------------------

pressures_B <- c(0.04, 0.06, 0.08, 0.09, 0.1, 0.12)

res_B <- time_kill_experiment(mod_args = 
                                list(basal_growth_probability = 0.1, 
                                     basal_survival_probability = 0.99,
                                     antibiotic_stresses = pressures_B,
                                     patch_capacity = 1E9,
                                     initial_popsize =1E3),
                              run_args = 
                                list(path = "time_killed_experiment_B",
                                     n_steps = 1,
                                     n_writes = 1000,
                                     cleanup_csv = TRUE,
                                     save_r_res =TRUE,
                                     r_res_path = "example_resistance"))


cells_dynamics_B <- count_from("Cell", "Patch", res_B)
cells_dynamics_B[, stress:= pressures_B[Patch+1]]

plotB <- ggplot(data = cells_dynamics_B,
                mapping = aes(x =step, y = multiplicity, 
                              group = stress,
                              color = as.factor(stress)))+
  geom_line(linewidth =1.2)+
  xlab("Time")+
  ylab("Cell count")+
  theme_bw()+
  scale_y_log10()+
  ggtitle("Growth probability = 0.1")






################################################################################
## MODULATION BY RESISTANCE
################################################################################
## In mse1, resistance is modeled through a 'susceptibility' parameter 
## which reduces antibiotic killing under the same stress 'a'. If a stress kills 
## an average of a% of sensitive cells, it will only kill approximately a 
## proportion a*s of cells with susceptibility 's'. This implies that the 
## minimum inhibitory stress (MIS) for a cell with susceptibility 's' is
## roughly MIS/s of that of a fully sensitive cell (matching the implicit 
## formulation x2, x3, etc... for resistance defined from MIC) 
################################################################################



## Model configuration ------------------------------

resistant_time_kill_experiment_model <- 
  function(basal_growth_probability = 0.5, 
           basal_survival_probability = 1.0,
           stress_suceptibility = 1.0,
           antibiotic_stresses =seq(0,0.5, by = 0.1),
           patch_capacity = 1E9,
           initial_popsize =1E9
){
  
  prm <- emptymodel() %>%
    chromosomeArchetype(0, 
                        fitness = basal_growth_probability,
                        survival = basal_survival_probability) %>% 
    chromosome(0, 0) %>% 
    cellArchetype(0) %>% cell(0, 0) %>%  
    assign("Chromosome", 0, "Cell", 0)%>%
    geneArchetype(0, susc = c(stress_suceptibility, 1.0), fitness =1.0)%>%
    gene(0,0)%>%
    assign("Gene", 0, "Chromosome", 0)
  
  
  for(i in 1:length(antibiotic_stresses)){
    prm <- prm %>%
      patchArchetype(i-1, patch_capacity, c(antibiotic_stresses[i], 0.0)) %>% 
      patch(i-1, i-1) %>% 
      assign("Cell", 0, "Patch", i-1, initial_popsize)
  }
  
  return(prm)
  
}


resistant_time_kill_experiment <- function(mod_args, run_args){
  
  
  
  ## model building 
  prm <- resistant_time_kill_experiment_model(
    mod_args$basal_growth_probability, 
    mod_args$basal_survival_probability,
    mod_args$stress_suceptibility,
    mod_args$antibiotic_stresses,
    mod_args$patch_capacity,
    mod_args$initial_popsize)
  
  mse1.plot_graph(prm, layout = "fr")
  
  
  
  ## model runing
  res <- mse_run(prm, 
                 path = run_args$path,
                 n_steps = run_args$n_steps,
                 n_writes = run_args$n_writes,
                 s1 = sample(1E9, 1))
  
  
  ## if needed, cleaning the csv results
  if(run_args$cleanup_csv){ mse1.cleanup_csv(run_args$path)}
  
  # if needed, save the R result (res) as a RDS object
  if(run_args$save_r_res){
    saveRDS(object = res, 
            file = paste(run_args$r_res_path, 
                         "/", run_args$path, ".rds", sep = ""))
  }
  
  
  ## 
  return(res)
  
}




## Simulations C ---------------------------------------------------------------

pressures_C <- c(0.12, 0.16, 0.17,0.18, 0.20 )

mse1.plot_graph(
  resistant_time_kill_experiment_model(basal_growth_probability = 0.1, 
                                       basal_survival_probability = 0.99,
                                       stress_suceptibility = 0.5,
                                       antibiotic_stresses = pressures_C,
                                       patch_capacity = 1E9,
                                       initial_popsize =1E3),
  layout = "fr")



res_C <- resistant_time_kill_experiment(
  mod_args = 
    list(basal_growth_probability = 0.1, 
         basal_survival_probability = 0.99,
         stress_suceptibility = 0.5,
         antibiotic_stresses = pressures_C,
         patch_capacity = 1E9,
         initial_popsize =1E3),
  run_args = 
    list(path = "time_killed_experiment_C",
         n_steps = 1,
         n_writes = 1000,
         cleanup_csv = TRUE,
         save_r_res =TRUE,
         r_res_path = "example_resistance"))




cells_dynamics_C <- count_from("Cell", "Patch", res_C)
cells_dynamics_C[, stress:= pressures_C[Patch+1]]

plotC <- ggplot(data = cells_dynamics_C,
                mapping = aes(x =step, y = multiplicity, 
                              group = stress,
                              color = as.factor(stress)))+
  geom_line(linewidth =1.2)+
  xlab("Time")+
  ylab("Cell count")+
  theme_bw()+
  scale_y_log10()+
  ggtitle("Growth probability = 0.1 - Sensitivity ~ 0.5 ")




## Simulations D ---------------------------------------------------------------

pressures_D <- c(0.4, 0.6, 0.8, 0.9, 0.99)

res_D <- resistant_time_kill_experiment(
  mod_args = 
    list(basal_growth_probability = 0.1, 
         basal_survival_probability = 0.99,
         stress_suceptibility = 1/10,
         antibiotic_stresses = pressures_D,
         patch_capacity = 1E9,
         initial_popsize =1E3),
  run_args = 
    list(path = "time_killed_experiment_D",
         n_steps = 1,
         n_writes = 1000,
         cleanup_csv = TRUE,
         save_r_res =TRUE,
         r_res_path = "example_resistance"))




cells_dynamics_D <- count_from("Cell", "Patch", res_D)
cells_dynamics_D[, stress:= pressures_D[Patch+1]]

plotD <- ggplot(data = cells_dynamics_D,
                mapping = aes(x =step, y = multiplicity, 
                              group = stress,
                              color = as.factor(stress)))+
  geom_line(linewidth =1.2)+
  xlab("Time")+
  ylab("Cell count")+
  theme_bw()+
  scale_y_log10()+
  ggtitle("Growth probability = 0.1 - Sensitivity ~ 0.1 ")





## Simulations E ---------------------------------------------------------------

pressures_E <- pressures

res_E <- resistant_time_kill_experiment(
  mod_args = 
    list(basal_growth_probability = 0.1, 
         basal_survival_probability = 0.99,
         stress_suceptibility = 1/2.75,
         antibiotic_stresses = pressures_E,
         patch_capacity = 1E9,
         initial_popsize =1E3),
  run_args = 
    list(path = "time_killed_experiment_E",
         n_steps = 1,
         n_writes = 1000,
         cleanup_csv = TRUE,
         save_r_res =TRUE,
         r_res_path = "example_resistance"))




cells_dynamics_E <- count_from("Cell", "Patch", res_E)
cells_dynamics_E[, stress:= pressures_E[Patch+1]]

plotE <- ggplot(data = cells_dynamics_E,
                mapping = aes(x =step, y = multiplicity, 
                              group = stress,
                              color = as.factor(stress)))+
  geom_line(linewidth =1.2)+
  xlab("Time")+
  ylab("Cell count")+
  theme_bw()+
  scale_y_log10()+
  ggtitle("Growth probability = 0.1 - Sensitivity ~ 0.36 ")





################################################################################
## PLOTTING
################################################################################


tiff("example_resistance/time_killed_experiment_sensitive.tif", 
     units = "cm", res = 300, 
     width = 24,
     height = 10)

ggpubr::ggarrange(plotA, plotB, nrow = 1, ncol = 2)

dev.off()




tiff("example_resistance/time_killed_experiment_resistant.tif", 
     units = "cm", res = 300, 
     width = 36,
     height = 10)

ggpubr::ggarrange(plotC, plotD, plotE, nrow = 1, ncol = 3)

dev.off()
