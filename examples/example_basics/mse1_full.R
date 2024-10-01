################################################################################
## MSE1 RESISTRACK MODEL                                                        
##  A bacterial hospital ecosystem simulation framework                          
## (c)2018-2024 Jean-Philippe Rasigade, <jean-philippe.rasigade@univ-lyon1.fr>
##  University of Lyon, Centre International de Recherche en Infectiologie       
################################################################################
## Complete model exemplifying how to parameter each admissible entity and edge
## allowed in a *mse1* model
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


## GlOBAL PARAMETERS, COMMON TO ALL SIMULATIONS PERFORMED IN THAT SCRIPT

cleanup_csv <- TRUE
save_r_res <- FALSE
r_res_path <- "example_basics"




##----------------------------------------------------------------------------##
## COMPLETE MODEL - STEP-BY-STEP
##----------------------------------------------------------------------------##

## MODEL CONFIGURATION -----------------------------------------

## 1. initialize with an empty model
prm <- emptymodel()


## 2. Define archetypes
prm <- prm %>% 
  geneArchetype(id = 0, susc = c(0.5, 1.0), fitness = 0.98) %>%
  geneArchetype(id = 1, susc = c(1.0, 0.01), fitness = 0.95) %>%
  plasmidArchetype(id = 0, loss = 1e-3, transfer = 1e-4, max_count = 2, fitness = 1.0)%>%
  chromosomeArchetype(id = 0, fitness = 0.1, survival = 0.99) %>%
  cellArchetype(id = 0) %>%
  patchArchetype(id = 0, capacity = 1e4, pressure = c(0.01, 0.0))  %>%
  patchArchetype(id = 1, capacity = 1e4, pressure = c(0.0, 0.05))


## 3. Define biological entities
## = any unique instance of {1 archetype + unique content}

## 3.1. Create nodes = "Agent"
prm <- prm %>% 
  gene(id = 0, type = 0) %>%
  gene(id = 1, type = 1) %>%
  plasmid(id = 0, type = 0) %>%
  chromosome(id = 0, type = 0) %>%
  cell(id = 0, type = 0)%>%
  patch(id = 0, type = 0)%>%
  patch(id =1, type = 1)


## 3.2. Add inclusion edges and multiplicity

prm <- prm %>% 
  assign("Gene", 0, "Chromosome", 0, mult = 1) %>%
  assign("Gene", 1, "Plasmid", 0, mult = 1) %>%
  assign("Chromosome", 0, "Cell", 0) %>%
  assign("Plasmid", 0, "Cell",  0) %>%
  assign("Cell", 0, "Patch", 0, mult = 10) 



## 4. if needed, Add specific agent/agent relationships

# 4.1. Plasmid range

prm <- prm %>% plasmidRange(0,0)

# 4.2. Patch-to-Patch Diffusion edges
prm <- prm %>% diffusion(0,1, diffusion = 1e-4, symmetric = TRUE)



## 5. Verification of configuration

## 5.1. Visual
mse1.plot_graph(prm)
mse1.plot_graph(prm, layout = "fr")
mse1.plot_graph(prm, layout = "free")

## 5.2. Indicative diagnostics
mse1.diagnose_model(prm)



## MODEL RUNNING -----------------------------------------------

# running parameters
csvres_folder <- "basic_full"

# run
res <- mse_run(prm, path = csvres_folder, n_steps = 1, n_writes = 2000, 
               s1 = sample(1e9, 1))

  
# if needed, clean_up csv folder
mse1.cleanup_csv(csvres_folder)


# if needed, store the R result in a RDS file for further use
if(save_r_res){
  saveRDS(object = res,
          file = paste(r_res_path, "/", csvres_folder, ".rds", sep = ""))
}






## ANALYSIS ----------------------------------------------------

## Cells diversity observed during the simulation

cells_description <- mse1.describe_cells(res)
cells_description 



## Cells in patches
d_cells <- count_from("Cell", "Patch", res)
d_cells <- merge(d_cells, cells_description[, c("Cell", "CellPrototype")],
                 by = "Cell")

ggplot(data = d_cells, 
       mapping = aes(x = step, y = multiplicity, 
                     group = interaction(CellPrototype, Patch),
                     color = as.factor(CellPrototype),
                     linetype = as.factor(Patch)
                     )
       )+
  geom_line(linewidth = 1.2)+
  xlab("Time (in step)")+
  ylab("Cell count")+
  ggtitle("Intra-patch Cell Dynamics")+
  theme_bw()



## Total cells in the ecosystem
d_cells_total <- d_cells[, .(multiplicity = sum(multiplicity)),
                         by = list(CellPrototype, step)]

ggplot(data = d_cells_total, 
       mapping = aes(x = step, y = multiplicity, 
                     group = CellPrototype,
                     color = as.factor(CellPrototype)))+
  geom_line(linewidth = 1.2)+
  xlab("Time (in step)")+
  ylab("Cell count")+
  ggtitle("Cell Dynamics (total)")+
  theme_bw()


  
## Plasmid dynamics

d_plasmids <- count_from("Plasmid", "Cell", res) %>% count_in("Patch", res)

ggplot(data = d_plasmids, 
       mapping = aes(x = step, y = multiplicity, 
                     group = interaction(Plasmid, Patch),
                     color = as.factor(Plasmid),
                     linetype = as.factor(Patch)
       )
)+
  geom_line(linewidth = 1.2)+
  xlab("Time (in step)")+
  ylab("Cell count")+
  ggtitle("Intra-patch Plasmid Dynamics")+
  theme_bw()



## Genes' dynamics

d_plas_genes <- count_from("Gene", "Plasmid", res) %>% 
  count_in("Cell", res)%>%
  count_in("Patch", res)

d_chrom_genes <- count_from("Gene", "Chromosome", res) %>% 
  count_in("Cell", res)%>%
  count_in("Patch", res)


d_genes <- rbind(d_plas_genes, 
                 d_chrom_genes)[, .(multiplicity = sum(multiplicity)),
                                by = list(Gene, Patch, step)]


ggplot(data = d_genes, 
       mapping = aes(x = step, y = multiplicity, 
                     group = interaction(Gene, Patch),
                     color = as.factor(Gene),
                     linetype = as.factor(Patch)
       )
)+
  geom_line(linewidth = 1.2)+
  xlab("Time (in step)")+
  ylab("Gene count")+
  ggtitle("Intra-patch Gene Dynamics")+
  theme_bw() 








