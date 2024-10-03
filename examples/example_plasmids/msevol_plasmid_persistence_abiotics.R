###############################################################################
## MSE1 RESISTRACK MODEL
##  A bacterial hospital ecosystem simulation framework
## (c)2018-2024 Jean-Philippe Rasigade, <jean-philippe.rasigade@univ-lyon1.fr>
##  University of Lyon, Centre International de Recherche en Infectiologie
################################################################################
## Conjugative Plasmid persistence in abiotic patch
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


library(foreach)
library(doParallel)





##----------------------------------------------------------------------------##
##  Plasmid persistence in abiotic patch
##----------------------------------------------------------------------------##
## We aim to illustrate the persistence of a conjugative plasmid in a patch
## without selective pressure, as a function of its features in terms of fitness
## cost, random segregation rate and conjugative transfer rate.
##
## This model is similar to the model of Lopatkin et al. (2019)
## DOI: 10.1038/s41467-017-01532-1 (Figure 1, with costly plasmid only), which
## demonstrated that costly plasmids can persist in the
## environment as soon as it transfer rate is above a critical threshold,
## roughly approximated by alpha*(loss+D)-D
## where :
##    * alpha is the ratio between the growth rate of plasmid-free cells and
##         the growth rate of plasmid-bearing cells
##    * loss is the random segregation rate of the plasmid
##    * D is the natural death rate of both plasmid-bearing and plasmid-free cells.
##
## In our model parameterized with probability and multiplicative penalization
## (meaning a additive cost in rate) :
##  plas_free_growth_rate = log(1+basal_growth)
##  plas_bearing_grow_rate = log(1+basal_growth*plasmid_fitness)
##  loss = -log(1-plasmid_loss)
##  D ~ -log(survival)   (slight underestimate due to the event order ...)
##
###----------------------------------------------------------------------------##



## Simulations with :
##  - growth = 0.3
##  - survival = 0.95
##  - plasmid_fitness = 0.9
##  - plasmid_loss = 1E-3
## Predicted critical transfer rate is roughly approximate by
##   log(1+0.3)/(log(1+0.3*0.9))*(-log(1-1E-3) - log(0.95)) +log(0.95)
##    ~ 0.0061086

transfer_rate <- seq(0.005, 0.01, by = 0.001)


## Initialize cluster
n_cores <- 5
n_cores <- min(c(n_cores, detectCores()[1]-1)) ## !!
cl <- makeCluster(n_cores)
registerDoParallel(cl)

## parallel simulation
results <- foreach(i = 1:length(transfer_rate),
                   .packages = c("data.table", "tidyr", "bit64", "Rmsevol")) %dopar% {

                     ##
                     lopatkin_abiotic_prm <- function(plasmid_transfer){

                       prm <- emptymodel()%>%
                         plasmidArchetype(0, loss = 1E-3,
                                          transfer = plasmid_transfer,
                                          max_count = 1,
                                          fitness = 0.90) %>%
                         plasmid(0,0)%>%
                         chromosomeArchetype(0,
                                             fitness = 0.3,
                                             survival =0.95)%>%
                         chromosome(0,0)%>%
                         plasmidRange(0,0)%>%
                         cellArchetype(0)%>% cell(0,0)%>%
                         assign("Chromosome", 0, "Cell", 0)%>%
                         assign("Plasmid", 0, "Cell", 0)%>%
                         patchArchetype(0, capacity = 1E9,
                                        pressure = c(0.0,0.0))%>%
                         patch(0,0)%>%
                         assign("Cell", 0, "Patch", 0, 1E3)

                       return(prm)


                     }

                     my_prm <- lopatkin_abiotic_prm(transfer_rate[i])


                     ## run
                     csvres_folder_i <- paste("plasmid_persistence_", i, sep = "")
                     res <- mse_run(my_prm, csvres_folder_i ,
                                    n_steps = 1, n_writes = 5000,
                                    s1 = sample(1E9, 1))
                     mse1.cleanup_csv(csvres_folder_i)

                     res


                   }
## stop cluster
registerDoSEQ()
stopCluster(cl)


## Analysis
d_cells <- rbindlist(lapply(results,
                            function(mylist){
                              cell_desc <- mse1.describe_cells(mylist)
                              cell_desc[Plasmid_0==1, label := "S1" ]
                              cell_desc[Plasmid_0==0, label := "S0" ]
                              d_cells <- count_from("Cell", "Patch", mylist)
                              d_cells <- merge(d_cells, cell_desc[, c("Cell", "label")],
                                               by = "Cell")
                              d_cells
                            }),
                     idcol = "repetition")
d_cells[, transfer_rate := transfer_rate[repetition]]


tiff(paste(dirname(getwd()),
           "/examples/example_plasmids/mse1_abiotic_critical_transfer_0.9.tiff",
           sep = ""),
     res = 300, units = "cm", width = 25, height = 15)
print(
  ggplot(d_cells,
         mapping = aes(x = step/24,
                       y = multiplicity,
                       group = interaction(label, transfer_rate),
                       color = label))+
    geom_line(linewidth = 1.2)+
    theme_bw()+
    facet_wrap("transfer_rate")
)

dev.off()



## Simulations with :
##  - growth = 0.3
##  - survival = 0.95
##  - plasmid_fitness = 0.8
##  - plasmid_loss = 1E-3
## Predicted critical transfer rate is roughly approximate by
##   log(1+0.3)/(log(1+0.3*0.8))*(-log(1-1E-3) - log(0.95)) +log(0.95)
##    ~ 0.01248772

transfer_rate <- seq(0.009, 0.02, by = 0.001)


## Initialize cluster
n_cores <- 5
n_cores <- min(c(n_cores, detectCores()[1]-1)) ## !!
cl <- makeCluster(n_cores)
registerDoParallel(cl)

## parallel simulation
results <- foreach(i = 1:length(transfer_rate),
                   .packages = c("data.table", "tidyr", "bit64", "Rmsevol")) %dopar% {

                     ##
                     lopatkin_abiotic_prm <- function(plasmid_transfer){

                       prm <- emptymodel()%>%
                         plasmidArchetype(0, loss = 1E-3,
                                          transfer = plasmid_transfer,
                                          max_count = 1,
                                          fitness = 0.80) %>%
                         plasmid(0,0)%>%
                         chromosomeArchetype(0,
                                             fitness = 0.3,
                                             survival =0.95)%>%
                         chromosome(0,0)%>%
                         plasmidRange(0,0)%>%
                         cellArchetype(0)%>% cell(0,0)%>%
                         assign("Chromosome", 0, "Cell", 0)%>%
                         assign("Plasmid", 0, "Cell", 0)%>%
                         patchArchetype(0, capacity = 1E9,
                                        pressure = c(0.0,0.0))%>%
                         patch(0,0)%>%
                         assign("Cell", 0, "Patch", 0, 1E3)

                       return(prm)


                     }

                     my_prm <- lopatkin_abiotic_prm(transfer_rate[i])


                     ## run
                     csvres_folder_i <- paste("plasmid_persistence_", i, sep = "")
                     res <- mse_run(my_prm, csvres_folder_i ,
                                    n_steps = 1, n_writes = 5000,
                                    s1 = sample(1E9, 1))
                     mse1.cleanup_csv(csvres_folder_i)

                     res


                   }
## stop cluster
registerDoSEQ()
stopCluster(cl)


## Analysis
d_cells <- rbindlist(lapply(results,
                            function(mylist){
                              cell_desc <- mse1.describe_cells(mylist)
                              cell_desc[Plasmid_0==1, label := "S1" ]
                              cell_desc[Plasmid_0==0, label := "S0" ]
                              d_cells <- count_from("Cell", "Patch", mylist)
                              d_cells <- merge(d_cells, cell_desc[, c("Cell", "label")],
                                               by = "Cell")
                              d_cells
                            }),
                     idcol = "repetition")
d_cells[, transfer_rate := transfer_rate[repetition]]



tiff(paste(dirname(getwd()),
           "/examples/example_plasmids/mse1_abiotic_critical_transfer_0.8.tiff",
           sep = ""),
     res = 300, units = "cm", width = 25, height = 15)
print(
  ggplot(d_cells,
         mapping = aes(x = step/24,
                       y = multiplicity,
                       group = interaction(label, transfer_rate),
                       color = label))+
    geom_line(linewidth = 1.2)+
    theme_bw()+
    facet_wrap("transfer_rate")
)

dev.off()




## Simulations with :
##  - growth = 0.3
##  - survival = 0.95
##  - plasmid_fitness = 0.7
##  - plasmid_loss = 1E-3
## Predicted critical transfer rate is roughly approximate by
##   log(1+0.3)/(log(1+0.3*0.7))*(-log(1-1E-3) - log(0.95)) +log(0.95)
##    ~ 0.02

transfer_rate <- seq(0.018, 0.028, by = 0.001)


## Initialize cluster
n_cores <- 5
n_cores <- min(c(n_cores, detectCores()[1]-1)) ## !!
cl <- makeCluster(n_cores)
registerDoParallel(cl)

## parallel simulation
results <- foreach(i = 1:length(transfer_rate),
                   .packages = c("data.table", "tidyr", "bit64", "Rmsevol")) %dopar% {

                     ##
                     lopatkin_abiotic_prm <- function(plasmid_transfer){

                       prm <- emptymodel()%>%
                         plasmidArchetype(0, loss = 1E-3,
                                          transfer = plasmid_transfer,
                                          max_count = 1,
                                          fitness = 0.70) %>%
                         plasmid(0,0)%>%
                         chromosomeArchetype(0,
                                             fitness = 0.3,
                                             survival =0.95)%>%
                         chromosome(0,0)%>%
                         plasmidRange(0,0)%>%
                         cellArchetype(0)%>% cell(0,0)%>%
                         assign("Chromosome", 0, "Cell", 0)%>%
                         assign("Plasmid", 0, "Cell", 0)%>%
                         patchArchetype(0, capacity = 1E9,
                                        pressure = c(0.0,0.0))%>%
                         patch(0,0)%>%
                         assign("Cell", 0, "Patch", 0, 1E3)

                       return(prm)


                     }

                     my_prm <- lopatkin_abiotic_prm(transfer_rate[i])


                     ## run
                     csvres_folder_i <- paste("plasmid_persistence_", i, sep = "")
                     res <- mse_run(my_prm, csvres_folder_i ,
                                    n_steps = 1, n_writes = 5000,
                                    s1 = sample(1E9, 1))
                     mse1.cleanup_csv(csvres_folder_i)

                     res


                   }
## stop cluster
registerDoSEQ()
stopCluster(cl)


## Analysis
d_cells <- rbindlist(lapply(results,
                            function(mylist){
                              cell_desc <- mse1.describe_cells(mylist)
                              cell_desc[Plasmid_0==1, label := "S1" ]
                              cell_desc[Plasmid_0==0, label := "S0" ]
                              d_cells <- count_from("Cell", "Patch", mylist)
                              d_cells <- merge(d_cells, cell_desc[, c("Cell", "label")],
                                               by = "Cell")
                              d_cells
                            }),
                     idcol = "repetition")
d_cells[, transfer_rate := transfer_rate[repetition]]



tiff(paste(dirname(getwd()),
           "/examples/example_plasmids/mse1_abiotic_critical_transfer_0.7.tiff",
           sep = ""),
     res = 300, units = "cm", width = 25, height = 15)
print(
  ggplot(d_cells,
         mapping = aes(x = step/24,
                       y = multiplicity,
                       group = interaction(label, transfer_rate),
                       color = label))+
    geom_line(linewidth = 1.2)+
    theme_bw()+
    facet_wrap("transfer_rate")
)

dev.off()






