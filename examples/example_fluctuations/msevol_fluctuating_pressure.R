################################################################################
## MSE1 RESISTRACK MODEL
##  A bacterial hospital ecosystem simulation framework
## (c)2018-2024 Jean-Philippe Rasigade, <jean-philippe.rasigade@univ-lyon1.fr>
##  University of Lyon, Centre International de Recherche en Infectiologie
################################################################################
## Example of simulations with regular "update" of parameters:
## FLUCTUATING PRESSURES
##
## Author: JP Rasigade, N Lenuzza
################################################################################


##----------------------------------------------------------------------------##
## Options, packages and generic functions
##----------------------------------------------------------------------------##

## WORKING DIRECTORY
## Make sure the the working directory is Rmsevol, i.e. the root folder that
## directly contains 'bin' (with the msevocli.exe client)

getwd()


## PACKAGES

library(data.table)
library(tidyr)
library(ggplot2)
library(visNetwork)
library(bit64)
library(Rmsevol)






#----------------------------------------------------------------------------##
## BASE EXEMPLE
##----------------------------------------------------------------------------##
## We simulate the competitive dynamics between two cell types, each
## resistant to either antibiotic A or B, competing for resources in
## a patch initially exposed to antibiotic A (for a duration T), followed
## by a period of selective pressure from antibiotic B (for the same
## duration T), and then back to antibiotic A, continuing in this
## alternating pattern.



## 0. global parameters
csvres_folder <- "fluctuation_modele"
T_in_steps <- 100
N_periods <- 10
n_steps <- 1
n_writes <- as.integer(round(T_in_steps/n_steps))


## 1. Model initialization

prm <- emptymodel() %>%
  geneArchetype(id = 0, susc = c(0.001, 1.0), fitness = 1.0) %>%
  gene(0,0)%>%
  geneArchetype(id = 1, susc = c(1.0, 0.001), fitness = 1.0) %>%
  gene(1,1)%>%
  chromosomeArchetype(0, fitness = 0.1, survival = 1.0) %>%
  chromosome(0, 0) %>% chromosome(1, 0) %>%
  cellArchetype(0) %>% cell(0, 0) %>% cell(1,0)%>%
  patchArchetype(0, 1e6, c(0.05, 0.0)) %>% patch(0, 0) %>%
  assign("Gene", 0, "Chromosome", 0) %>%
  assign("Gene", 1, "Chromosome", 1) %>%
  assign("Chromosome", 0, "Cell", 0) %>%
  assign("Chromosome", 1, "Cell", 1) %>%
  assign("Cell", 1, "Patch", 0, 5e5) %>%
  assign("Cell", 0, "Patch", 0, 5e5)

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

  ## switch pressure
  new_prm$Vertex_PatchArchetype[, pressure.0 := as.numeric(pressure.0)]
  new_prm$Vertex_PatchArchetype[, pressure.1 := as.numeric(pressure.1)]
  press_0 <- as.numeric(new_prm$Vertex_PatchArchetype[id== 0, ]$pressure.0)
  press_1 <- as.numeric(new_prm$Vertex_PatchArchetype[id== 0, ]$pressure.1)
  new_prm$Vertex_PatchArchetype[id ==0, pressure.0 := press_1]
  new_prm$Vertex_PatchArchetype[id ==0, pressure.1 := press_0]

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


## 3. Check for pressure evolution

pressure_dyn <- merge(
  res$edg$Edge_PatchArchetype_Patch_Multiplicity,
  res$vtx$Vertex_PatchArchetype,
  by.x = c("PatchArchetype", "step"),
  by.y = c("Vertex_PatchArchetype", "step")
)

ggplot(data = pressure_dyn,
       mapping = aes(x =step, group = Patch))+
  geom_line(aes(y = pressure.0), color = "red", linewidth = 1.2)+
  geom_line(aes(y = pressure.1), color = "blue", linewidth = 1.2, linetype = "dotted")+
  theme_bw()+
  ggtitle("Pressures")+
  geom_vline(xintercept = seq(0, T_in_steps*N_periods, by = T_in_steps),
             linetype = "dashed")



## 4. Cell dynamics

d_cells <- count_from("Cell", "Patch", res)
ggplot(data = d_cells,
       mapping = aes(x =step, y = multiplicity,
                     group = Cell, color = as.factor(Cell)))+
  geom_line(linewidth = 1.2)+
  theme_bw()+
  ggtitle("Fluctuating_period")+
  geom_vline(xintercept = seq(0, T_in_steps*N_periods, by = T_in_steps),
             linetype = "dashed")

