################################################################################
## MSE1 RESISTRACK MODEL
##  A bacterial hospital ecosystem simulation framework
## (c)2018-2024 Jean-Philippe Rasigade, <jean-philippe.rasigade@univ-lyon1.fr>
##  University of Lyon, Centre International de Recherche en Infectiologie
################################################################################
## Frequently errors in model configurations
## Author: JP Rasigade, N Lenuzza
################################################################################


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






##----------------------------------------------------------------------------##
## Behaviour of mse1.diagnose_model
##----------------------------------------------------------------------------##

## Correct model

correct_prm <- emptymodel() %>%
  chromosomeArchetype(id = 0, fitness = 0.5, survival = 1.0) %>%
  chromosomeArchetype(id = 1, fitness = 0.3, survival = 1.0) %>%
  chromosome(0, 0) %>% chromosome(1,1)%>%
  plasmidArchetype(0, 0.1, 0.0, 1, 1.0)%>%
  plasmid(0,0)%>%
  cellArchetype(0) %>% cell(0, 0) %>% cell(1,0) %>%
  patchArchetype(0, 1000, c(0.0, 0.0)) %>% patch(0, 0) %>%
  assign("Chromosome", 0, "Cell", 0) %>%
  assign("Chromosome", 1, "Cell", 1) %>%
  assign("Plasmid", 0, "Cell", 1) %>%
  assign("Cell", 0, "Patch", 0, 100) %>%
  assign("Cell", 1, "Patch", 0, 100) %>%
  plasmidRange(0,1)

mse1.plot_graph(correct_prm)
mse1.diagnose_model(correct_prm )



## Model 1 : two chromosome archetypes with the same id

prm <- emptymodel() %>%
  chromosomeArchetype(id = 0, fitness = 0.5, survival = 1.0) %>%
  chromosomeArchetype(id = 0, fitness = 0.3, survival = 1.0) %>%
  chromosome(0, 0) %>% chromosome(1,1)%>%
  plasmidArchetype(0, 0.1, 0.0, 1, 1.0)%>%
  plasmid(0,0)%>%
  cellArchetype(0) %>% cell(0, 0) %>% cell(1,0) %>%
  patchArchetype(0, 1000, c(0.0, 0.0)) %>% patch(0, 0) %>%
  assign("Chromosome", 0, "Cell", 0) %>%
  assign("Chromosome", 1, "Cell", 1) %>%
  assign("Plasmid", 0, "Cell", 1) %>%
  assign("Cell", 0, "Patch", 0, 100) %>%
  assign("Cell", 1, "Patch", 0, 100) %>%
  plasmidRange(0,1)


mse1.plot_graph(prm)
mse1.diagnose_model(prm)
rm(prm)






## Model 2 : Duplicated edge assignation

prm <- emptymodel() %>%
  chromosomeArchetype(id = 0, fitness = 0.5, survival = 1.0) %>%
  chromosomeArchetype(id = 1, fitness = 0.3, survival = 1.0) %>%
  chromosome(0, 0) %>% chromosome(1,1)%>%
  plasmidArchetype(0, 0.1, 0.0, 1, 1.0)%>%
  plasmid(0,0)%>%
  cellArchetype(0) %>% cell(0, 0) %>% cell(1,0) %>%
  patchArchetype(0, 1000, c(0.0, 0.0)) %>% patch(0, 0) %>%
  assign("Chromosome", 0, "Cell", 0) %>%
  assign("Chromosome", 1, "Cell", 1) %>%
  assign("Plasmid", 0, "Cell", 1) %>%
  assign("Cell", 0, "Patch", 0, 100) %>%
  assign("Cell", 0, "Patch", 0, 100) %>%
  plasmidRange(0,1)


mse1.plot_graph(prm)
mse1.plot_graph(prm, layout = "free")

mse1.diagnose_model(prm)
rm(prm)



## Model 3 : Error in cell id assignation

prm <- emptymodel() %>%
  chromosomeArchetype(id = 0, fitness = 0.5, survival = 1.0) %>%
  chromosomeArchetype(id = 1, fitness = 0.3, survival = 1.0) %>%
  chromosome(0, 0) %>% chromosome(1,1)%>%
  plasmidArchetype(0, 0.1, 0.0, 1, 1.0)%>%
  plasmid(0,0)%>%
  cellArchetype(0) %>% cell(0, 0) %>% cell(1,0) %>%
  patchArchetype(0, 1000, c(0.0, 0.0)) %>% patch(0, 0) %>%
  assign("Chromosome", 0, "Cell", 0) %>%
  assign("Chromosome", 1, "Cell", 1) %>%
  assign("Plasmid", 0, "Cell", 1) %>%
  assign("Cell", 0, "Patch", 0, 100) %>%
  assign("Cell", 2, "Patch", 0, 100) %>%
  plasmidRange(0,1)


mse1.plot_graph(prm)
mse1.diagnose_model(prm)
rm(prm)



## Model 4 : Error in archetype attribution

prm <- emptymodel() %>%
  chromosomeArchetype(id = 0, fitness = 0.5, survival = 1.0) %>%
  chromosomeArchetype(id = 1, fitness = 0.3, survival = 1.0) %>%
  chromosome(0, 0) %>% chromosome(1,1)%>%
  plasmidArchetype(0, 0.1, 0.0, 1, 1.0)%>%
  plasmid(0,1)%>%
  cellArchetype(0) %>% cell(0, 0) %>% cell(1,0) %>%
  patchArchetype(0, 1000, c(0.0, 0.0)) %>% patch(0, 0) %>%
  assign("Chromosome", 0, "Cell", 0) %>%
  assign("Chromosome", 1, "Cell", 1) %>%
  assign("Plasmid", 0, "Cell", 1) %>%
  assign("Cell", 0, "Patch", 0, 100) %>%
  assign("Cell", 1, "Patch", 0, 100) %>%
  plasmidRange(0,1)


mse1.plot_graph(prm)
mse1.diagnose_model(prm)
rm(prm)


## Model 5 : Inversion of id and type during agent creation
## Exemple of 2 plasmids
## Error detectable only in case of invalid configuration

prm <- emptymodel() %>%
  chromosomeArchetype(id = 0, fitness = 0.5, survival = 1.0) %>%
  chromosomeArchetype(id = 1, fitness = 0.3, survival = 1.0) %>%
  chromosome(0, 0) %>% chromosome(1,1)%>%
  plasmidArchetype(0, 0.1, 0.0, 1, 1.0)%>%
  plasmid(0,1)%>%plasmid(0,1)%>%
  cellArchetype(0) %>% cell(0, 0) %>% cell(1,0) %>%
  patchArchetype(0, 1000, c(0.0, 0.0)) %>% patch(0, 0) %>%
  assign("Chromosome", 0, "Cell", 0) %>%
  assign("Chromosome", 1, "Cell", 1) %>%
  assign("Plasmid", 0, "Cell", 1) %>%
  assign("Cell", 0, "Patch", 0, 100) %>%
  assign("Cell", 1, "Patch", 0, 100) %>%
  plasmidRange(0,1)


mse1.plot_graph(prm)
mse1.diagnose_model(prm)
rm(prm)




## Model 6 : Missing plasmid range for conjugative plasmid

prm <- emptymodel() %>%
  chromosomeArchetype(id = 0, fitness = 0.5, survival = 1.0) %>%
  chromosomeArchetype(id = 1, fitness = 0.3, survival = 1.0) %>%
  chromosome(0, 0) %>% chromosome(1,1)%>%
  plasmidArchetype(0, 0.0, 0.0, 1, 1.0)%>%
  plasmidArchetype(1, 0.0, 0.1, 1, 1.0)%>%
  plasmid(0,0)%>%plasmid(1,1)%>%
  cellArchetype(0) %>% cell(0, 0) %>% cell(1,0) %>%
  patchArchetype(0, 1000, c(0.0, 0.0)) %>% patch(0, 0) %>%
  assign("Chromosome", 0, "Cell", 0) %>%
  assign("Chromosome", 1, "Cell", 1) %>%
  assign("Plasmid", 0, "Cell", 1) %>%
  assign("Plasmid", 1, "Cell", 1) %>%
  assign("Cell", 0, "Patch", 0, 100) %>%
  assign("Cell", 1, "Patch", 0, 100)


mse1.plot_graph(prm)
mse1.diagnose_model(prm)
rm(prm)

