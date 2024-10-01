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


## SOURCING FUNCTIONS

source("msvr_fun/cfg_funcs v1.0.R")
source("msvr_fun/parse_funcs.R")
source("msvr_fun/helper_funcs.R")


## GlOBAL PARAMETERS

run_for_tutorial <- TRUE
cleanup_csv <- TRUE







##----------------------------------------------------------------------------##
## COMPLETE MODEL - STEP-BY-STEP
##----------------------------------------------------------------------------##

## MODEL CONFIGURATION -----------------------------------------

## 1. initialize with an empty model
prm <- emptymodel()


## 2. Define archetypes
prm <- prm %>%
  geneArchetype(id = 0, susc = c(0.5, 1.0), fitness = 0.95)%>%
  geneArchetype(id = 1, susc = c(1.0, 0.01), fitness = 0.95)%>%
  chromosomeArchetype(id = 0, fitness = 0.1, survival = 0.95)%>%
  plasmidArchetype(id = 0, loss = 1e-3, transfer = 1e-4, max_count = 1, 
                   fitness = 1.0)%>%
  cellArchetype(id = 0)%>%
  patchArchetype(id = 0, capacity = 1e6, pressure = c(0.01, 0.0))%>%
  patchArchetype(id = 1, capacity = 1e6, pressure = c(0.0, 0.01))


## 3. Define biological entities
## = any unique instance of {1 archetype + unique content}

## 3.1. Create nodes = "Agent"
prm <- prm %>% 
  gene(id = 0, type = 0) %>%
  gene(id = 1, type = 1) %>%
  plasmid(id = 0, type = 0) %>%
  chromosome(id = 0, type = 0) %>% 
  chromosome(id = 1, type = 0) %>%
  cell(id = 0, type = 0)%>%
  cell(id = 1, type = 0)%>%
  patch(id = 0, type = 0)%>%
  patch(id =1, type = 1)


## 3.2. Add inclusion edges and multiplicity

prm <- prm %>% 
  assign("Gene", 0, "Chromosome", 1, mult = 1) %>%
  assign("Gene", 1, "Plasmid", 0, mult = 1) %>%
  assign("Chromosome", 0, "Cell", 0) %>%
  assign("Chromosome", 1, "Cell", 1) %>%
  assign("Plasmid", 0, "Cell",  1) %>%
  assign("Cell", 0, "Patch", 0, mult = 100)%>%
  assign("Cell", 1, "Patch", 0, mult = 100)



## 4. if needed, Add specific agent/agent relationships

# 4.1. Plasmid range

prm <- prm %>% plasmidRange(0,0)

# 4.2. Patch-to-Patch Diffusion edges
prm <- prm %>% diffusion(0,1, diffusion = 1e-4, symmetric = TRUE)



## 5. Verification of configuration

## 5.1. Visual
mse1.plot_graph(prm)
mse1.plot_graph(prm, layout = "fr")

## Saving files to display them in the tutorial
if(run_for_tutorial){
  visSave(mse1.plot_graph(prm, layout = "hr"), 
          file ="./documentation/toy_example_outputs/toy_exemple_initial_graph_hierachical.html")
  
  visSave(mse1.plot_graph(prm, layout = "fr"), 
          file ="./documentation/toy_example_outputs/toy_exemple_initial_graph_frlayout.html")
 
}


## Exemple of custom color and custom layout
my_col <- RColorBrewer::brewer.pal(5, "Dark2")
names(my_col) <- c("Patch","Cell", "Chromosome", "Plasmid", "Gene")
tmp <- mse1.plot_graph(prm, color_map = my_col, layout = "free")
tmp$x$nodes$x <-c(0, 100,33,67, 36.5, 20,53,0, 20, 53, 100, 67,40, 70, 53, 67)*5
tmp$x$nodes$y <-c(0,0, 0, 0, 30, 20,20, 100, 40,40, 100, 40, 80,80, 60, 60)*5
print(tmp %>% visNodes(fixed = TRUE))


## 5.2. Indicative diagnostics
mse1.diagnose_model(prm)




## MODEL RUNNING  (SINGLE RUN) --------------------------------------

# running parameters
csvres_folder <- "toy_example"

# run
res <- mse_run(prm, path = csvres_folder, 
               n_steps = 1,
               n_writes =5000, 
               s1 = sample(1e9, 1),
               s2 = sample(1e9, 1))

  
# if needed, clean_up csv folder stored in "bin"
mse1.cleanup_csv(csvres_folder)


# if needed, store the R result in a RDS file for further use
# of the result in the .rmd tutorial
if(run_for_tutorial){
  saveRDS(object = res,
          file = paste("./documentation/toy_example_outputs/",
                       csvres_folder, "_simulation.rds",
                       sep = ""))
}



## ANALYSIS ----------------------------------------------------

### CELL DYNAMICS

### Explore the diversity of cells generated during the run

cells_description <- mse1.describe_cells(res)
DT::datatable(cells_description) 

### [optional] Add labels for plotting (manual, strongly depends 
### on the use case)

cells_description[Chromosome_0== 1 & Plasmid_0 == 0,cell_label:= "Ch0"]
cells_description[Chromosome_0== 1 & Plasmid_0 > 0, 
                  cell_label:= paste0("Ch0 + ", 
                                      Plasmid_0, " plas", sep = "")
]
cells_description[Chromosome_1== 1 & Plasmid_0 == 0, cell_label:="Ch1"]
cells_description[Chromosome_1== 1 & Plasmid_0 > 0, 
                  cell_label:= paste0("Ch1 + ", 
                                      Plasmid_0, " plas", sep = "")
]

DT::datatable(cells_description)


### Within-patch cell counts (based on cell labels instead of ID)

n_cells_in_patches <- count_from("Cell", "Patch", res)
n_cells_in_patches <- merge(n_cells_in_patches,
                            cells_description[, c("Cell", "cell_label")],
                            by = "Cell")

p1 <- ggplot(data = n_cells_in_patches , 
             mapping = aes(x = step, y = multiplicity, 
                           group = interaction(cell_label, Patch),
                           color = cell_label))+
  geom_line(linewidth = 1.2)+
  xlab("Time (in step)")+
  ylab("Cell count")+
  theme_bw()+
  facet_wrap("Patch")+
  scale_colour_manual(name = "Cell type",
                      values = RColorBrewer::brewer.pal(4, "Dark2"))


plot(p1)



### Total patch cell dynamics

n_cells <- n_cells_in_patches[, 
                              .(multiplicity = sum(multiplicity)),
                              by = list(cell_label, step)]

p2 <- ggplot(data = n_cells , 
             mapping = aes(x = step, y = multiplicity, 
                           group = cell_label,
                           color = cell_label))+
  geom_line(linewidth = 1.2)+
  xlab("Time (in step)")+
  ylab("Cell count")+
  theme_bw()+
  scale_colour_manual(name = "Cell type",
                      values = RColorBrewer::brewer.pal(4, "Dark2"))

plot(p2)


if(run_for_tutorial){
  tiff("./documentation/toy_example_outputs/cell_dynamics.tiff",
       res = 300, unit = "cm", width = 25, height = 20)
  print(
    ggpubr::ggarrange(
      ggpubr::ggarrange(p2+ggtitle("Total Cell dynamics"),
                        p2+ggtitle("Total Cell dynamics (Log10)")+
                          scale_y_log10(),
                        nrow = 1, ncol = 2, common.legend = TRUE, 
                        legend = "none"),
      p1+ggtitle("Intra-patch Cell Dynamics (Log10)")+
        scale_y_log10(),
      nrow = 2, ncol = 1, common.legend = TRUE
    )
  )
  dev.off()
}



### Average between patches
## Note that the results contain only non zero multiplicity, that can bias the
## averages or any other summary statistics directy computed in data.tables 

n_cells_in_patches <- mse1.add_zero_steps(n_cells_in_patches)

mean_cells_per_patch <- n_cells_in_patches[, .(av = mean(multiplicity)), 
                                           by = list(cell_label, step)]

ggplot(data = mean_cells_per_patch, 
       mapping = aes(x = step, y = av, 
                     group = cell_label,
                     color = cell_label))+
  geom_line(linewidth = 1.2)+
  xlab("Time (in step)")+
  ylab("Cell count")+
  theme_bw()+
  scale_colour_manual(name = "Cell type",
                      values = RColorBrewer::brewer.pal(4, "Dark2"))





## PLASMID DYNAMICS 

n_plasmids_in_patches <- 
  count_from("Plasmid", "Cell",  res) %>%
  count_in("Patch", res)

ggplot(data = n_plasmids_in_patches, 
       mapping = aes(x = step, y = multiplicity,
                     group = interaction(Plasmid, Patch)))+
  geom_line(linewidth = 1.2)+
  xlab("Time (in step)")+
  ylab("Plasmid count")+
  theme_bw()+
  facet_wrap("Patch")



### GENE DYNAMICS

n_genes_in_patches <- 
  rbind( 
    ## Number of plasmidic genes
    count_from("Gene", "Plasmid",  res) %>%  
      count_in("Cell", res) %>% 
      count_in("Patch", res),
    ## Number of chromosomal genes
    count_from("Gene", "Chromosome",  res) %>%  
      count_in("Cell", res) %>% 
      count_in("Patch", res))[, ## Sum their multiplicity
                              .(multiplicity = sum(multiplicity)),
                              by = list(Gene, Patch, step)]

setcolorder(n_genes_in_patches, c("Gene", "Patch", "multiplicity", "step")) 


ggplot(data = n_genes_in_patches, 
       mapping = aes(x = step, y = multiplicity,
                     group = interaction(Gene, Patch),
                     color = as.factor(Gene)))+
  geom_line(linewidth = 1.2)+
  xlab("Time (in step)")+
  ylab("Gene count")+
  theme_bw()+
  facet_wrap("Patch")





## HOT RESTART WITH MODIFIED PRESSURE ------------------------------------------


## Extract and modify new input
new_prm <- mse1.extract_prm_from_mseres(res, step = 5000)
new_prm$Vertex_PatchArchetype[id ==0, pressure.1:= 0.01]

mse1.plot_graph(new_prm)


## Re-run the simulation with clean_up
new_res  <- mse_run(new_prm, 
                    path = "toy_example_run2", 
                    n_steps = 1,
                    n_writes =5000, 
                    s1 = sample(1e9, 1),
                    s2 = sample(1e9, 1))

mse1.cleanup_csv("toy_example_run2")

## Total results
res <- mse1.merge_res(res, 
                      mse1.shift_step(new_res, 5000))



if(run_for_tutorial){
  saveRDS(object = res,
          file = paste("./documentation/toy_example_outputs/",
                       "toy_example_simulation_with_hot_restart.rds",
                       sep = ""))
}

## Cell dynamics

full_cells_description <- mse1.describe_cells(res)
full_cells_description[Chromosome_0== 1 & Plasmid_0 == 0,cell_label:= "Ch0"]
full_cells_description[Chromosome_0== 1 & Plasmid_0 > 0, 
                  cell_label:= paste0("Ch0 + ", 
                                      Plasmid_0, " plas", sep = "")
]
full_cells_description[Chromosome_1== 1 & Plasmid_0 == 0, cell_label:="Ch1"]
full_cells_description[Chromosome_1== 1 & Plasmid_0 > 0, 
                  cell_label:= paste0("Ch1 + ", 
                                      Plasmid_0, " plas", sep = "")
]

print(full_cells_description)

full_n_cells_in_patches <- count_from("Cell", "Patch", res)
full_n_cells_in_patches <- merge(full_n_cells_in_patches,
                                 full_cells_description[, c("Cell", "cell_label")],
                            by = "Cell")

p1_full <- ggplot(data = full_n_cells_in_patches, 
             mapping = aes(x = step, y = multiplicity, 
                           group = interaction(cell_label, Patch),
                           color = cell_label))+
  geom_line(linewidth = 1.2)+
  xlab("Time (in step)")+
  ylab("Cell count")+
  theme_bw()+
  facet_wrap("Patch", ncol = 1, nrow = 2)+
  scale_colour_manual(name = "Cell type",
                      values = RColorBrewer::brewer.pal(4, "Dark2"))+
  geom_vline(xintercept = 5000, linetype = "dashed")


plot(p1_full)


if(run_for_tutorial){
  tiff("./documentation/toy_example_outputs/cell_dynamics_with_restart.tiff",
       res = 300, unit = "cm", width = 25, height = 20)
  plot(p1_full)
  dev.off()
}





## PARALLEL EXECUTION  ---------------------------------------------------------


library(foreach)
library(doParallel)


N_repetition <- 10


## Initialize cluster
n_cores <- 5
n_cores <- min(c(n_cores, detectCores()[1]-1)) ## !!
cl <- makeCluster(n_cores)
registerDoParallel(cl)


## Parallel execution of the same model

results <- foreach(i = 1:10, 
                   .packages = c("data.table", "tidyr")) %dopar% {
                     
                     csvres_folder_i <- paste("toy_example_parallel_", i, sep = "") 
                     res_i <- mse_run(prm, csvres_folder_i, 1, 5000, s1 = sample(1E9,1))
                     
                     # if needed, cleaning the csv results
                     mse1.cleanup_csv(csvres_folder_i)
                     
                     res_i
                   }



## stop cluster
registerDoSEQ()
stopCluster(cl)



if(run_for_tutorial){
  saveRDS(object = results,
          file = paste("./documentation/toy_example_outputs/",
                       "toy_example_parallel_simulations.rds",
                       sep = ""))
}





## Analyses 


n_cells_per_patches <- rbindlist(lapply(results,
                            function(mylist){
                              
                              n_cells <- count_from("Cell", "Patch", mylist)
                              
                              cell_desc <- mse1.describe_cells(mylist)
                              cell_desc[Chromosome_0== 1 & Plasmid_0 == 0,cell_label:= "Ch0"]
                              cell_desc[Chromosome_0== 1 & Plasmid_0 > 0, 
                                                     cell_label:= paste0("Ch0 + ", 
                                                                         Plasmid_0, " plas", sep = "")
                              ]
                              cell_desc[Chromosome_1== 1 & Plasmid_0 == 0, cell_label:="Ch1"]
                              cell_desc[Chromosome_1== 1 & Plasmid_0 > 0, 
                                                     cell_label:= paste0("Ch1 + ", 
                                                                         Plasmid_0, " plas", sep = "")
                              ]
                              
                              n_cells  <- merge(n_cells,
                                                cell_desc [, c("Cell", "cell_label")],
                                                               by = "Cell")
                              
                              return(n_cells)
                              
                            }),
                     idcol = "repetition")


p1_repet <- ggplot(data = n_cells_per_patches, 
       mapping = aes(x = step, y = multiplicity, 
                     group = interaction(cell_label, Patch, repetition),
                     color = cell_label))+
  geom_line()+
  xlab("Time (in step)")+
  ylab("Cell count")+
  theme_bw()+
  facet_wrap("Patch", ncol = 1, nrow = 2)+
  scale_colour_manual(name = "Cell type",
                      values = RColorBrewer::brewer.pal(4, "Dark2"))


plot(p1_repet)
if(run_for_tutorial){
  tiff("./documentation/toy_example_outputs/cell_dynamics_with_repetition.tiff",
       res = 300, unit = "cm", width = 15, height = 20)
  plot(p1_repet)
  dev.off()
}

