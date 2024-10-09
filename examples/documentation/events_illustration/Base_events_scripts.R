################################################################################
## MSE1 RESISTRACK MODEL
##  A bacterial hospital ecosystem simulation framework
## (c)2018-2024 Jean-Philippe Rasigade, <jean-philippe.rasigade@univ-lyon1.fr>
##  University of Lyon, Centre International de Recherche en Infectiologie
################################################################################
## Basic scenarios illustrating the behavior of each events for the
## introcution of the R tutorial
##
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



## GlOBAL PARAMETERS, COMMON TO ALL SIMULATIONS PERFORMED IN THAT SCRIPT

cleanup_csv <- TRUE
path_for_fig <- paste(dirname(getwd()), "/examples/documentation/events_illustration/",
                      sep ="")







##----------------------------------------------------------------------------##
## MODEL 1 - CELL BIRTH MODULATED BY FITNESS COST INDUCED BY GENES OR PLASMIDS
##----------------------------------------------------------------------------##
## We simulate the growth of four different cells, with identical basal growth
## rate, but carrying genes and/or plasmids modulating their growth:
##  - cell with id 0 : no penalization
##  - cell with id 1 : chromosomal carriage of a costly ARG
##  - cell with id 2 : carriage of a costly plasmid without ARG
##  - cell with id 3 : carriage of a costly plasmid carrying the costly ARG.
##
## Each cell type is growing independently in their own patch, and together in
## an additional patch.
##----------------------------------------------------------------------------##


## INPUTS -----------------------

# model parameters
basal_growth_probability <- 0.1
plasmid_fitness <- 0.95
gene_fitness <- 0.8

popsize <- 1e3
capacity <- 1e6

# running parameters
csvres_folder <- "basic_fitness"
n_steps <- 1
n_writes <- 150


## MODEL CONFIGURATION ----------

prm <- emptymodel() %>%
  geneArchetype(0, c(1.0, 1.0), gene_fitness) %>% gene(0, 0) %>%
  chromosomeArchetype(0, fitness = basal_growth_probability, survival = 1.0) %>%
  chromosome(0, 0) %>% chromosome(1,0)%>%
  plasmidArchetype(0, 0.0, 0.0, 1, plasmid_fitness) %>%
  plasmid(0, 0) %>% plasmid(1, 0) %>%
  cellArchetype(0) %>% cell(0, 0) %>% cell(1,0)%>% cell(2,0)%>%cell(3,0) %>%
  patchArchetype(0, capacity, c(0.0, 0.0)) %>%
  patch(0, 0) %>%
  assign("Gene", 0, "Chromosome", 1) %>%
  assign("Gene", 0, "Plasmid", 1) %>%
  assign("Chromosome", 0, "Cell", 0) %>%
  assign("Chromosome", 1, "Cell", 1) %>%
  assign("Chromosome", 0, "Cell", 2) %>%
  assign("Chromosome", 0, "Cell", 3) %>%
  assign("Plasmid", 0, "Cell", 2) %>%
  assign("Plasmid", 1, "Cell", 3) %>%
  assign("Cell", 0, "Patch", 0, popsize)%>%
  assign("Cell", 1, "Patch", 0, popsize)%>%
  assign("Cell", 2, "Patch", 0, popsize)%>%
  assign("Cell", 3, "Patch", 0, popsize)

#mse1.plot_graph(prm)
#mse1.diagnose_model(prm)


## MODEL RUNNING ----------------

# launching simulator
res <- mse_run(prm, csvres_folder, n_steps, n_writes)

# if needed, cleaning the csv results
if(cleanup_csv){ mse1.cleanup_csv(csvres_folder)}


## RESULTS ----------------------

## computing intra-patches dynamic of cells

d_cells <- count_from("Cell", "Patch", res)

## plotting curves of independant growths
plot_cells <- ggplot(data = d_cells,
                         mapping = aes(x =step, y = multiplicity,
                                       group = interaction(Cell, Patch),
                                       color = as.factor(Cell),
                                       linetype = as.factor(Cell)))+
  geom_line(linewidth =1.2)+
  xlab("Time")+
  ylab("Cell count")+
  theme_bw(base_size = 12)+
  scale_color_manual(name = "Cell type",
                     values = RColorBrewer::brewer.pal(4, "Dark2"),
                     labels = c("WT",
                                "Chromosomal ARG",
                                "ARG-free plasmid",
                                "Plasmidic ARG"))+
  scale_linetype_manual(name = "Cell type",
                        values = c("solid", "dashed", "dotted", "dotdash"),
                        labels = c("WT",
                                   "Chromosomal ARG",
                                   "ARG-free plasmid",
                                   "Plasmidic ARG"))

png(paste( path_for_fig, "/growth.png",sep =""),
     unit = "cm", res = 300, width = 17, height = 10)
plot(plot_cells)
dev.off()




## Verification of growth rates during the exponential phases
plot(plot_cells +scale_y_log10())

step_lim <- 25
apparent_growth_rate <- rep(NA, 4)

for(i in 0:3){
  lm_mod <- lm( log(multiplicity)~step,
                  data = d_cells[step<step_lim& Cell == i, ])

  apparent_growth_rate[i+1] <-coefficients(lm_mod)["step"]
}
round(apparent_growth_rate, digit = 3)






##----------------------------------------------------------------------------##
## MODEL 2 - NATURAL DEATH AND ANTIBIOTIC KILLING
##----------------------------------------------------------------------------##


## INPUTS -----------------------

# model parameters

basal_survival_probability <- 0.95
selective_pressure <- c(0.05, 0.05)
susceptibility_0 <- c(0.0001, 1.0)
susceptibility_1 <- c(1.0, 0.2)
popsize <- 1e7
capacity <- 1e7

# running parameters
csvres_folder <- "basic_killing"
n_steps <- 1
n_writes <- 500



## MODEL CONFIGURATION ----------

prm <- emptymodel() %>%
  geneArchetype(0, susceptibility_0, 1.0) %>% gene(0, 0) %>%
  geneArchetype(1, susceptibility_1, 1.0) %>% gene(1, 1) %>%
  chromosomeArchetype(0, fitness = 0.0, survival = 0.95) %>%
  chromosome(0, 0) %>% chromosome(1,0) %>%
  chromosome(2,0) %>% chromosome(3,0)%>%
  cellArchetype(0) %>% cell(0, 0) %>% cell(1,0)%>% cell(2,0)%>%cell(3,0) %>%
  patchArchetype(0, capacity, c(0.0, 0.0)) %>% patch(0, 0) %>%
  patchArchetype(1, capacity, selective_pressure) %>% patch(1, 1) %>%
  assign("Gene", 0, "Chromosome", 1) %>%
  assign("Gene", 1, "Chromosome", 2) %>%
  assign("Gene", 0, "Chromosome", 3) %>%
  assign("Gene", 1, "Chromosome", 3) %>%
  assign("Chromosome", 0, "Cell", 0) %>%
  assign("Chromosome", 1, "Cell", 1) %>%
  assign("Chromosome", 2, "Cell", 2) %>%
  assign("Chromosome", 3, "Cell", 3) %>%
  assign("Cell", 0, "Patch", 0, popsize)%>%
  assign("Cell", 1, "Patch", 0, popsize)%>%
  assign("Cell", 2, "Patch", 0, popsize)%>%
  assign("Cell", 3, "Patch", 0, popsize)%>%
  assign("Cell", 0, "Patch", 1, popsize)%>%
  assign("Cell", 1, "Patch", 1, popsize)%>%
  assign("Cell", 2, "Patch", 1, popsize)%>%
  assign("Cell", 3, "Patch", 1, popsize)

mse1.plot_graph(prm)
mse1.diagnose_model(prm)


## MODEL RUNNING ----------------

# launching simulator
res <- mse_run(prm, csvres_folder, n_steps, n_writes)

# if needed, cleaning the csv results
if(cleanup_csv){ mse1.cleanup_csv(csvres_folder)}


## RESULTS ----------------------

## computing intra-patches dynamic of cells

d_cells <- count_from("Cell", "Patch", res)
d_cells[Patch== 0, patch := "Antibiotic free"]
d_cells[Patch== 1, patch := "Antibiotic stress (0.05, 0.05)"]

## plotting curves of independant growths
plot_cells <- ggplot(data = d_cells,
                     mapping = aes(x =step, y = multiplicity,
                                   group = interaction(Cell, patch),
                                   color = as.factor(Cell),
                                   linetype = as.factor(Cell)))+
  geom_line(linewidth =1.2)+
  xlab("Time")+
  ylab("Cell count")+
  theme_bw(base_size = 12)+
  scale_color_manual(name = "Cell type",
                     values = RColorBrewer::brewer.pal(4, "Dark2"),
                     labels = c("WT",
                                "ARG 1",
                                "ARG 2",
                                "ARG 1 + ARG 2"))+
  scale_linetype_manual(name = "Cell type",
                        values = c("solid", "dashed", "dotted", "dotdash"),
                        labels = c("WT",
                                   "ARG 1",
                                   "ARG 2",
                                   "ARG 1 + ARG 2"))+
  facet_wrap("patch")





png(paste( path_for_fig, "/death_killing.png",sep =""),
     unit = "cm", res = 300, width = 25, height = 10)
plot(plot_cells +scale_y_log10())
dev.off()




## Verification of death rates

apparent_death_rate_patch_0 <- apparent_death_rate_patch_1 <- rep(NA, 4)

for(i in 0:3){
  lm_mod <- lm( log(multiplicity)~step,
                data = d_cells[multiplicity>1E3& Cell == i & Patch ==0, ])

  apparent_death_rate_patch_0[i+1] <-coefficients(lm_mod)["step"]

  lm_mod <- lm( log(multiplicity)~step,
                data = d_cells[multiplicity>1E3& Cell == i & Patch ==1, ])

  apparent_death_rate_patch_1[i+1] <-coefficients(lm_mod)["step"]


}


round(apparent_death_rate_patch_0, digit = 3)
round(apparent_death_rate_patch_1, digit = 3)





##----------------------------------------------------------------------------##
## MODEL 3 - Plasmid transfer occuring in on specific niche
##----------------------------------------------------------------------------##

## INPUTS -----------------------

# model parameters
plasmid_transfer <- 0.01
plasmid_maxcount <- 1
capacity <- 1e6


# running parameters
csvres_folder <- "basic_plasmid_range"
n_steps <- 1
n_writes <- 1500



## MODEL CONFIGURATION ----------

prm <- emptymodel() %>%
  plasmidArchetype(0, 0.0, plasmid_transfer, plasmid_maxcount, 1.0) %>%
  plasmid(0, 0) %>%
  chromosomeArchetype(0, fitness = 0.0, survival = 1.0) %>% chromosome(0, 0) %>%
  cellArchetype(0) %>% cell(0, 0) %>% cell(1,0) %>%
  chromosomeArchetype(1, fitness = 0.0, survival = 1.0) %>% chromosome(1, 1) %>%
  cellArchetype(1) %>% cell(2, 1) %>%
  chromosomeArchetype(2, fitness = 0.0, survival = 1.0) %>% chromosome(2, 2) %>%
  cellArchetype(2) %>% cell(3, 2) %>%
  patchArchetype(0, capacity, c(0.0, 0.0)) %>% patch(0, 0) %>%
  assign("Plasmid", 0, "Cell", 0, mult = 1) %>%
  assign("Chromosome", 0, "Cell", 0) %>%
  assign("Chromosome", 0, "Cell", 1) %>%
  assign("Chromosome", 1, "Cell", 2) %>%
  assign("Chromosome", 2, "Cell", 3) %>%
  assign("Cell", 0, "Patch", 0, floor(capacity/4)) %>%
  assign("Cell", 1, "Patch", 0, floor(capacity/4)) %>%
  assign("Cell", 2, "Patch", 0, floor(capacity/5)) %>%
  assign("Cell", 3, "Patch", 0, floor(capacity/5)) %>%
  plasmidRange(0, 0)%>%
  plasmidRange(0, 1)

mse1.plot_graph(prm, layout = "fr")
mse1.diagnose_model(prm)


## MODEL RUNNING ----------------

# launching simulator
res <- mse_run(prm, csvres_folder, n_steps, n_writes)

# if needed, cleaning the csv results
if(cleanup_csv){ mse1.cleanup_csv(csvres_folder)}



## RESULTS ----------------------

## Caracterization of cells generating during the simulations
cells <- mse1.describe_cells(res)
print(cells)


## Plotting cells dynamics
d_cells <- count_from("Cell", "Patch", res)
d_cells_prototype <- merge(d_cells, cells[, c("Cell", "CellArchetype",
                                              "CellPrototype", "Plasmid_0")],
                           by = "Cell", all.x = TRUE)


plot_cells <- ggplot(data = d_cells_prototype,
                     mapping = aes(x =step, y = multiplicity,
                                   group = CellPrototype,
                                   color = as.factor(CellArchetype),
                                   linetype = as.factor(Plasmid_0)))+
  geom_line(linewidth =1.2)+
  xlab("Time")+
  ylab("Cell count")+
  theme_bw(base_size = 12)+
  scale_color_manual(name = "Cell type",
                     values =RColorBrewer::brewer.pal(4, "Dark2")[1:3],
                     labels = c(
                       "0" = "Chromosome 0",
                       "1" = "Chromosome 1",
                       "2" = "Chromosome 2"
                     ))+
  scale_linetype_manual(name = "with plasmid ?",
                        values = c(
                          "0" = "solid",
                          "1" = "dashed"
                        ),
                        labels = c("plasmid-free",
                                   "plasmid-bearing"))


png(paste( path_for_fig, "/plasmid_transfer.png",sep =""),
     unit = "cm", res = 300, width = 17, height = 10)
plot(plot_cells)
dev.off()




## Dynamics of the total number of plasmids

d_plasmids <- count_from("Plasmid", "Cell", res)
d_plasmids <- count_in(d_plasmids, "Patch", res)

n_donnors <- d_cells_prototype[Plasmid_0>0,
                               .(n_donnors = sum(multiplicity)),
                               by = step]

mean_plasmid_per_donnor <- d_cells_prototype[Plasmid_0>0,
                                             .(mean_plas = sum(multiplicity*Plasmid_0)/sum(multiplicity)),
                                             by = step]

n_recipients <- d_cells_prototype[Plasmid_0<plasmid_maxcount & CellArchetype %in% c(0,1),
                                  .(n_recipients = sum(multiplicity)),
                                  by = step]


dt <- merge(n_donnors, mean_plasmid_per_donnor, by = "step")
dt<- merge(dt, n_recipients, by = "step")
dt[, n_expected_transfer := (n_donnors/capacity)*
     (1- (1-plasmid_transfer)^mean_plas)*n_recipients]


eff_trans <- data.table(step = d_plasmids$step[-1],
                        n_transfer = diff(d_plasmids$multiplicity))

plot_transfer <- ggplot(data = eff_trans,
       mapping = aes(x = step, y = n_transfer))+
  geom_line(col = "grey")+
  theme_bw(base_size = 12)+
  geom_line(data = dt,
            mapping = aes(x = step, y = n_expected_transfer ),
            col = "red", linewidth = 1.2)+
  ggtitle("Number of transfers")+
  ylab("Number of new plasmids per time step")+
xlab("Time")

plot_transfer_cells <- ggplot(data = dt,
                              mapping = aes(x = step))+
  geom_line(aes(y = n_donnors), col = "black",linetype = "dashed",  linewidth = 1.2)+
  geom_line(aes(y = n_recipients), col = "black",
            linewidth = 1.2)+
  theme_bw(base_size = 12)+
  ggtitle("Number of donnor and recipient cells")+
  ylab("Cell count")+
  xlab("Time")


png(paste( path_for_fig, "/plasmid_transfer_2.png",sep =""),
     unit = "cm", res = 300, width = 22, height = 10)
ggpubr::ggarrange(
  plot_transfer_cells,
  plot_transfer,
  ncol = 2, nrow = 1)
dev.off()






##----------------------------------------------------------------------------##
## MODEL 4 - PLASMID LOSS
##----------------------------------------------------------------------------##

## INPUTS -----------------------

# model parameters
plasmid_loss <- 0.05
popsize <- 1e6
capacity <- 1e6

# running parameters
csvres_folder <- "basic_plasmid_loss"
n_steps <- 1
n_writes <- 500



## MODEL CONFIGURATION ----------

prm <- emptymodel() %>%
  plasmidArchetype(0, plasmid_loss, 0.0, 5, 1.0) %>% plasmid(0, 0) %>%
  chromosomeArchetype(0, fitness = 0.0, survival = 1.0) %>% chromosome(0, 0) %>%
  cellArchetype(0) %>% cell(0, 0) %>%
  patchArchetype(0, capacity, c(0.0, 0.0)) %>% patch(0, 0) %>%
  assign("Plasmid", 0, "Cell", 0, mult = 5) %>%
  assign("Chromosome", 0, "Cell", 0) %>%
  assign("Cell", 0, "Patch", 0, popsize)

mse1.plot_graph(prm)
mse1.diagnose_model(prm)


## MODEL RUNNING ----------------

# launching simulator
res <- mse_run(prm, csvres_folder, n_steps, n_writes)

# if needed, cleaning the csv results
if(cleanup_csv){ mse1.cleanup_csv(csvres_folder)}


## RESULTS ----------------------


## Caracterization of cells generating during the simulations
cells <- mse1.describe_cells(res)
print(cells)


## Plotting cells dynamics
d_cells <- count_from("Cell", "Patch", res)
d_cells <- merge(d_cells, cells[, c("Cell", "Plasmid_0")], by = "Cell")


colPal <- rainbow(6)
names(colPal) <- as.character(5:0)

plot_cells <- ggplot(data = d_cells,
                     mapping = aes(x =step, y = multiplicity,
                                   group = Plasmid_0, color = as.factor(Plasmid_0)))+
  geom_line(linewidth =1.2)+
  xlab("Time")+
  ylab("Cell count")+
  theme_bw()+
  scale_color_manual(name = "Cell type",
                     values =colPal,
                     labels = paste("Cells with ", 0:5, " plasmids", sep = ""))




png(paste( path_for_fig, "/plasmid_loss.png",sep =""),
     unit = "cm", res = 300, width = 25, height = 10)
ggpubr::ggarrange(
  plot_cells ,
  plot_cells +scale_y_log10(),
  ncol = 2, nrow = 1,
  common.legend = TRUE,
  legend = "left")
dev.off()



lm_mod <- lm(
  log(multiplicity)~step,
           data = d_cells[Cell ==0 &  step<150, ])

coefficients(lm_mod)["step"]




##----------------------------------------------------------------------------##
## MODEL 5 - PATCH-TO-PATCH DIFFUSION
##----------------------------------------------------------------------------##

## INPUTS -----------------------

# model parameters
popsize <- 1e6
capacity <- 1e6

# running parameters
csvres_folder <- "basic_diffusion"
n_steps <- 1
n_writes <- 500



## MODEL CONFIGURATION ----------

prm <- emptymodel() %>%
  chromosomeArchetype(0, fitness = 0.0, survival = 1.0) %>% chromosome(0, 0) %>%
  cellArchetype(0) %>% cell(0, 0) %>%
  patchArchetype(0, capacity, c(0.0, 0.0)) %>%
  patch(0, 0) %>%patch(1, 0) %>%patch(2, 0) %>%
  assign("Chromosome", 0, "Cell", 0) %>%
  assign("Cell", 0, "Patch", 0, popsize)%>%
  diffusion(0,1,diffusion = 0.1, symmetric = T)%>%
  diffusion(1,2,diffusion = 0.01, symmetric = F)

mse1.plot_graph(prm)
mse1.diagnose_model(prm)


## MODEL RUNNING ----------------

# launching simulator
res <- mse_run(prm, csvres_folder, n_steps, n_writes)

# if needed, cleaning the csv results
if(cleanup_csv){ mse1.cleanup_csv(csvres_folder)}


## RESULTS ----------------------

## Plotting cells dynamics
d_cells <- count_from("Cell", "Patch", res)



plot_cells <- ggplot(data = d_cells,
                     mapping = aes(x =step, y = multiplicity,
                                   group = interaction(Cell, Patch),
                                   color = as.factor(Patch),
                                   linetype = as.factor(Patch)))+
  geom_line(linewidth =1.2)+
  xlab("Time")+
  ylab("Cell count")+
  theme_bw(base_size = 12)+
  scale_color_manual(name = "Patch",
                     values = RColorBrewer::brewer.pal(4, "Dark2")[1:3],
                     labels = c("Patch 0",
                                "Patch 1",
                                "Patch 2"))+
  scale_linetype_manual(name = "Patch",
                        values = c("solid", "dashed", "dotdash"),
                        labels = c("Patch 0",
                                   "Patch 1",
                                   "Patch 2"))



png(paste( path_for_fig, "/patch_to_patch_diffusion.png",sep =""),
     unit = "cm", res = 300, width = 17, height = 10)
  plot_cells
dev.off()

