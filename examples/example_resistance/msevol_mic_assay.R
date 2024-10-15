###############################################################################
## MSE1 RESISTRACK MODEL
##  A bacterial hospital ecosystem simulation framework
## (c)2018-2024 Jean-Philippe Rasigade, <jean-philippe.rasigade@univ-lyon1.fr>
##  University of Lyon, Centre International de Recherche en Infectiologie
################################################################################
## Simulation of MIC experiment
## Author: JP Rasigade, N Lenuzza
################################################################################
## Here, we simulate a  96 well plate (8 rows, 12 columns), in which :
##    - the first column serves as a control, with no antibiotic present (patch
##      pressure = 0.0).
##    - the second column is initialized with a stress set to 1.0,
##    - in each successive column (from 3 to 12), the antibiotic concentration is
##      halved relative to the previous column.
## Each row evaluates the MIS of one strain containing a chromosomal resistant
## gene, whose susceptibility is halved relative to the previous row.
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



##----------------------------------------------------------------------------##
## Parameters
##----------------------------------------------------------------------------##

## Patches
capacity <- 1E9
popsize  <- 1E6
base_pressure <- 1
dilution_factor <- 2

## Cells lines
basal_growth <- 0.05
natural_survival <- 1.0 #no death
gene_sensitivities <- 1/2^(0:7)



##----------------------------------------------------------------------------##
## Model configuration
##----------------------------------------------------------------------------##

## Model initialisation
prm <- emptymodel()

## Create basal chromosome and cell archetypes
prm <- prm%>%
  cellArchetype(id = 0) %>%
  chromosomeArchetype(id=0,
                      fitness = basal_growth,
                      survival = 1.0)

## Create cells with given sensitivities
for(row_i in 1:8){
  prm <- prm %>%
    geneArchetype(id = row_i-1, susc = c(gene_sensitivities[row_i],
                                         1.0), fitness =1.0)%>%
    gene(row_i-1, row_i-1)%>%
    chromosome(row_i-1,0)%>%
    assign("Gene", row_i-1, "Chromosome", row_i-1)%>%
    cell(row_i-1, 0)%>%
    assign("Chromosome",row_i-1, "Cell", row_i-1)
}


## Create Patches archetypes with decreasing pressures
prm <- prm %>%
  patchArchetype(id = 0, capacity, pressure = c(0.0,0.0)) ## Abiotic control
current_stress <- 1.0
for(i in 1:11){
  prm <- prm %>%
    patchArchetype(id = i,
                   capacity = 1E9, pressure = c(current_stress ,0.0)) ## Abiotic control
  current_stress <- current_stress/2
}


## Add patches and assign cells
patch_current_id <- 0
for(patcharch_i in 0:11){
  for(cells_i in 0:7){
    prm <- prm %>%
      patch(patch_current_id, patcharch_i)%>%
      assign("Cell",cells_i, "Patch", patch_current_id, mult = popsize)
    patch_current_id <- patch_current_id+1
  }
}


## Visual checks
mse1.plot_graph(prm, vertices_type = "Patch", layout = "fr")
mse1.plot_graph(prm, vertices_type = c("Gene", "Chromosome", "Cell"), layout = "fr")

mse1.plot_graph(prm, add_archetype = FALSE, layout = "fr")

##----------------------------------------------------------------------------##
## Simulation run
##----------------------------------------------------------------------------##
## We keep only the last time (t = 3000 time step), without intermediate writes
csvres_folder <- "mic_experiment"
res <- mse_run(prm, csvres_folder , 3000, 1)
mse1.cleanup_csv(csvres_folder)


##----------------------------------------------------------------------------##
## Analysis
##----------------------------------------------------------------------------##

## counting cells in each well (= patch)

d<- count_from("Cell", "Patch", res)
d <- mse1.add_zero_steps(d)


## adding physical position of the patch in the plate
d$patch_x <- floor(d$Patch/8)
d$patch_y <- 8-d$Patch%%8


## Plot
plot(pl<- ggplot(data = d[step == max(step), ],
       mapping = aes(x = patch_x,
                     y = patch_y,
                     fill = log10(multiplicity)))+
  geom_point(shape = 21, size = 18)+
  geom_point(data = d[multiplicity ==0, ],
             aes(x = patch_x, y = patch_y),
             shape = 4, color = "black", size = 4)+
  scale_fill_gradient(
    name = "Cell count [log10]",
    low = "white",
    high = "#56B1F7",
    limits = c(0, log10(capacity))
  )+
  geom_vline(xintercept = 0.5, linetype = "dotted")+
  geom_text(data = data.frame(x = 0:11,
                              y = rep(9, 12),
                              label = c("0", paste("1:", 2^(0:10), sep = ""))),
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE)+
  geom_polygon(
    data = data.frame(
      x = c(0.5, 0.5,  11.5, 11.5),
      y = c(10.5, 9.5, 9.5, 9.5)
    ), aes(x = x, y = y),
    inherit.aes = FALSE,
    fill = "red",
    alpha = 0.5)+
  geom_text(x = 1, y = 10, label = "STRESS", color = "black", hjust = 0 )+
  geom_text(data = data.frame(x = rep(-1, 8),
                              y = 8:1,
                              label = gene_sensitivities),
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE)+
  theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )+
  ylab("Sensitivity of cells")

)

png(paste(dirname(getwd()),"/examples/example_resistance/MIS_example.png", sep = ""),
           res = 300, units = "cm", height = 15, width = 24)

plot(pl)
dev.off()

