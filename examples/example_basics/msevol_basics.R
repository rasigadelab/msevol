################################################################################
## MSE1 RESISTRACK MODEL
##  A bacterial hospital ecosystem simulation framework
## (c)2018-2024 Jean-Philippe Rasigade, <jean-philippe.rasigade@univ-lyon1.fr>
##  University of Lyon, Centre International de Recherche en Infectiologie
################################################################################
## Basic scenarios illustrating the behavior of each event
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
save_r_res <- FALSE
r_res_path <- "example_basics"



##----------------------------------------------------------------------------##
## MODEL 1 - BASIC BIRTH
##----------------------------------------------------------------------------##
## We simulate here the logistic growth of one bacterial strain in a single
## ecological patch without selective pressure


## INPUTS -----------------------

# model parameters
basal_growth_probability <- 0.05
popsize <- 10
capacity <- 1e9

# running parameters

csvres_folder <- "basic_birth"
n_steps <- 1
n_writes <- 700



## MODEL CONFIGURATION ----------

prm <- emptymodel() %>%
  chromosomeArchetype(0, fitness = basal_growth_probability, survival = 1.0) %>%
  chromosome(0, 0) %>%
  cellArchetype(0) %>% cell(0, 0) %>%
  patchArchetype(0, capacity, c(0.0, 0.0)) %>% patch(0, 0) %>%
  assign("Chromosome", 0, "Cell", 0) %>%
  assign("Cell", 0, "Patch", 0, popsize)

mse1.plot_graph(prm, layout = "fr")



## MODEL RUNNING ----------------

# launching simulator
res <- mse_run(prm, csvres_folder, n_steps, n_writes, s1= sample(1E9,1))

# if needed, cleaning the csv results
if(cleanup_csv){ mse1.cleanup_csv(csvres_folder)}

# if needed, save the R result (res) as a RDS object
if(save_r_res){
  saveRDS(object = res, file = paste(r_res_path,
                                     "/", csvres_folder, ".rds", sep = ""))
}



## RESULTS ----------------------

## computing intra-patches dynamic of cells

d_cells <- count_from("Cell", "Patch", res)

## plotting curves
plot_mod1 <- ggplot(data = d_cells,
       mapping = aes(x =step, y = multiplicity,
                     group = Cell, color = as.factor(Cell)))+
  geom_line(linewidth =1.2)+
  xlab("Time")+
  ylab("Cell count")+
  theme_bw()

plot(plot_mod1 + ggtitle("Basic birth"))
plot(plot_mod1 + ggtitle("Basic birth (Log)")+ scale_y_continuous(trans = "log10"))


## Verification of growth rate
simulated_growth_rate <- log(1+basal_growth_probability)/1

step_lim <- 300
lm_mod <- lm( log(multiplicity)~step,
              data = d_cells[step<step_lim, ])
apparent_growth_rate <- coefficients(lm_mod)["step"]

newdata <- data.frame(step = d_cells[step<step_lim, ]$step,
                multiplicity = exp(predict(lm_mod)))

plot(plot_mod1 + ggtitle("Basic birth (Log10)")+
       scale_y_continuous(trans = "log10")+
       geom_line(data = newdata,
                 mapping = aes(x = step, y = multiplicity),
                 color = "blue",
                 linetype = "dashed",
                 inherit.aes = FALSE,
                 linewidth = 1.2)+
       geom_hline(yintercept = capacity,
                  color = "blue",
                  linetype = "dotted")+
       geom_text(data = data.frame(x = 0,  y = capacity,
                                   label = "Carrying capacity"),
                 mapping = aes(x = x, y = y, label =label),
                 color = "blue",
                 hjust = 0, vjust = 1,
                 inherit.aes = FALSE)+
       geom_text(data = cbind.data.frame(
         newdata[nrow(newdata),],
                 label = paste("Exponential growth\n with rate ",
                               round(apparent_growth_rate, digit = 3),
                               " per unit time", sep = "")
         ),
     mapping = aes(x = step, y = multiplicity, label =label),
     color = "blue",
     hjust = 0.1, vjust = 3,
     inherit.aes = FALSE)
     )





##----------------------------------------------------------------------------##
## MODEL 2 - NATURAL DEATH
##----------------------------------------------------------------------------##
## We simulate here the natural death of one single cell in a single patch
## without antibiotic pressure


## INPUTS -----------------------

# model parameters
basal_survival_probability <- 0.95
popsize <- 1e9
capacity <- 1e9

# running parameters

csvres_folder <- "basic_death"
n_steps <- 1
n_writes <- 700


## MODEL CONFIGURATION ----------

prm <- emptymodel() %>%
  chromosomeArchetype(0, fitness = 0, survival = basal_survival_probability) %>%
  chromosome(0, 0) %>%
  cellArchetype(0) %>% cell(0, 0) %>%
  patchArchetype(0, capacity, c(0.0, 0.0)) %>% patch(0, 0) %>%
  assign("Chromosome", 0, "Cell", 0) %>%
  assign("Cell", 0, "Patch", 0, popsize)

mse1.plot_graph(prm)



## MODEL RUNNING ----------------

# launching simulator
res <- mse_run(prm, csvres_folder, n_steps, n_writes, s1= sample(1E9,1))

# if needed, cleaning the csv results
if(cleanup_csv){ mse1.cleanup_csv(csvres_folder)}

# if needed, save the R result (res) as a RDS object
if(save_r_res){
  saveRDS(object = res, file = paste(r_res_path,
                                     "/", csvres_folder, ".rds", sep = ""))
}



## RESULTS ----------------------

## computing intra-patches dynamic of cells

d_cells <- count_from("Cell", "Patch", res)

## plotting curves
plot_mod1 <- ggplot(data = d_cells,
                    mapping = aes(x =step, y = multiplicity,
                                  group = Cell, color = as.factor(Cell)))+
  geom_line(linewidth =1.2)+
  xlab("Time")+
  ylab("Cell count")+
  theme_bw()

plot(plot_mod1 + ggtitle("Basic death"))
plot(plot_mod1 + ggtitle("Basic death (Log10)")+ scale_y_continuous(trans = "log10"))



## Verification of death rate
simulated_death_rate <- -log(basal_survival_probability)/1

step_lim <- 300
lm_mod <- lm(log(multiplicity)~step,
             data = d_cells[step < step_lim, ])
apparent_death_rate <- -1*coefficients(lm_mod)["step"]

newdata <- data.frame(step = d_cells[step < step_lim, ]$step,
                      multiplicity = exp(predict(lm_mod)))

plot(plot_mod1 +
       ggtitle("Basic death (Log10)")+
       scale_y_continuous(trans = "log10")+
       geom_line(data = newdata,
                 mapping = aes(x = step, y = multiplicity),
                 color = "blue",
                 linetype = "dashed",
                 inherit.aes = FALSE,
                 linewidth = 1.2)+
       geom_text(data = cbind.data.frame(
         newdata[nrow(newdata),],
         label = paste("Exponential death\n with rate ",
                       round(apparent_death_rate, digit = 3),
                       " per unit time", sep = "")
       ),
       mapping = aes(x = step, y = multiplicity, label =label),
       color = "blue",
       hjust = 0, vjust = - 1,
       inherit.aes = FALSE)
)




##----------------------------------------------------------------------------##
## MODEL 3 - NATURAL DEATH & ANTIBIOTIC KILLING
##----------------------------------------------------------------------------##


## INPUTS -----------------------

# model parameters

basal_survival_probability <- 0.95
selective_pressure <- c(0.1, 0.0)
popsize <- 1e9
capacity <- 1e9

# running parameters
csvres_folder <- "basic_killing"
n_steps <- 1
n_writes <- 700





## MODEL CONFIGURATION ----------

prm <- emptymodel() %>%
  chromosomeArchetype(0, fitness = 0, survival = basal_survival_probability) %>%
  chromosome(0, 0) %>%
  cellArchetype(0) %>% cell(0, 0) %>%
  patchArchetype(0, capacity, selective_pressure) %>% patch(0, 0) %>%
  assign("Chromosome", 0, "Cell", 0) %>%
  assign("Cell", 0, "Patch", 0, popsize)

mse1.plot_graph(prm)



## MODEL RUNNING ----------------

# launching simulator
res <- mse_run(prm, csvres_folder, n_steps, n_writes, s1= sample(1E9,1))

# if needed, cleaning the csv results
if(cleanup_csv){ mse1.cleanup_csv(csvres_folder)}

# if needed, save the R result (res) as a RDS object
if(save_r_res){
  saveRDS(object = res, file = paste(r_res_path,
                                     "/", csvres_folder, ".rds", sep = ""))
}



## RESULTS ----------------------

## computing intra-patches dynamic of cells

d_cells <- count_from("Cell", "Patch", res)

## plotting curves
plot_mod1 <- ggplot(data = d_cells,
                    mapping = aes(x =step, y = multiplicity,
                                  group = Cell, color = as.factor(Cell)))+
  geom_line(linewidth =1.2)+
  xlab("Time")+
  ylab("Cell count")+
  theme_bw()

plot(plot_mod1 + ggtitle("Basic killing"))
plot(plot_mod1 + ggtitle("Basic killing (Log10)")+ scale_y_continuous(trans = "log10"))



## Verification of death rate
simulated_death_rate <- -(log(basal_survival_probability) +
                                  log(1-selective_pressure[1])+
                                  log(1-selective_pressure[2]))/1

step_lim <- 100
lm_mod <- lm(log(multiplicity)~step,
             data = d_cells[step < step_lim, ])
apparent_death_rate <- -1*coefficients(lm_mod)["step"]

newdata <- data.frame(step = d_cells[step < step_lim, ]$step,
                      multiplicity = exp(predict(lm_mod)))

plot(plot_mod1 +
       ggtitle("Basic death (Log10)")+
       scale_y_continuous(trans = "log10")+
       geom_line(data = newdata,
                 mapping = aes(x = step, y = multiplicity),
                 color = "blue",
                 linetype = "dashed",
                 inherit.aes = FALSE,
                 linewidth = 1.2)+
       geom_text(data = cbind.data.frame(
         newdata[nrow(newdata),],
         label = paste("Exponential death\n with rate ",
                       round(apparent_death_rate, digit = 3),
                       " per unit time", sep = "")
       ),
       mapping = aes(x = step, y = multiplicity, label =label),
       color = "blue",
       hjust = 0, vjust = - 1,
       inherit.aes = FALSE)
)



##----------------------------------------------------------------------------##
## MODEL 4 - CELL RESISTANCE TO ANTIBIOTIC KILLING
##----------------------------------------------------------------------------##
## We simulate here the antibiotic killing  of resistant cells,
## (carrying a chromosomal resistance gene) vs sensitive cells, in a patch with
## one selective pressure


## INPUTS -----------------------

# model parameters

basal_survival_probability <- 0.95
selective_pressure <- c(0.1, 0.0)
antibiotic_sensitivity <- c(0.1, 1.0)
popsize <- 1e9
capacity <- 1e9

# running parameters
csvres_folder <- "basic_resistance"
n_steps <- 1
n_writes <- 700



## MODEL CONFIGURATION ----------

prm <- emptymodel() %>%
  geneArchetype(0, antibiotic_sensitivity, 1.0) %>% gene(0, 0) %>%
  chromosomeArchetype(0, fitness = 0, survival = basal_survival_probability) %>%
  chromosome(0, 0) %>%chromosome(1, 0) %>%
  cellArchetype(0) %>% cell(0, 0) %>% cell(1,0)%>%
  patchArchetype(0, capacity, selective_pressure) %>% patch(0, 0) %>%
  assign("Gene", 0, "Chromosome", 0) %>%
  assign("Chromosome", 0, "Cell", 0) %>%
  assign("Chromosome", 1, "Cell", 1) %>%
  assign("Cell", 0, "Patch", 0, floor(popsize/2))%>%
  assign("Cell", 1, "Patch", 0, floor(popsize/2))

mse1.plot_graph(prm)
mse1.diagnose_model(prm)


## MODEL RUNNING ----------------

# launching simulator
res <- mse_run(prm, csvres_folder, n_steps, n_writes, s1= sample(1E9,1))

# if needed, cleaning the csv results
if(cleanup_csv){ mse1.cleanup_csv(csvres_folder)}

# if needed, save the R result (res) as a RDS object
if(save_r_res){
  saveRDS(object = res, file = paste(r_res_path,
                                     "/", csvres_folder, ".rds", sep = ""))
}



## RESULTS ----------------------

## computing intra-patches dynamic of cells
d_cells <- count_from("Cell", "Patch", res)

## plotting curves
plot_mod1 <- ggplot(data = d_cells,
                    mapping = aes(x =step, y = multiplicity,
                                  group = Cell, color = as.factor(Cell)))+
  geom_line(linewidth =1.2)+
  xlab("Time")+
  ylab("Cell count")+
  theme_bw()+
  scale_color_manual(name = "Cell type",
                     values = c("orange", "cyan"),
                     labels = c("resistant", "sensitive"))

plot(plot_mod1 + ggtitle("Basic resistance"))
plot(plot_mod1 + ggtitle("Basic resistance (Log10)")+ scale_y_continuous(trans = "log10"))



## Computing simulated of death rates
simulated_death_rate_R <- -(log(basal_survival_probability) +
                            log(1-selective_pressure[1]*antibiotic_sensitivity[1])+
                            log(1-selective_pressure[2]*antibiotic_sensitivity[2]))/1

simulated_death_rate_S <- -(log(basal_survival_probability) +
                              log(1-selective_pressure[1])+
                              log(1-selective_pressure[2]))/1



## Verification of effective death rates

step_lim <- 200
lm_mod_R <- lm(log(multiplicity)~step,
             data = d_cells[step < step_lim & Cell == 0, ])
apparent_death_rate_R <- -1*coefficients(lm_mod_R)["step"]

newdata_R <- data.frame(step = d_cells[step < step_lim & Cell == 0, ]$step,
                      multiplicity = exp(predict(lm_mod_R)))


step_lim <- 100
lm_mod_S <- lm(log(multiplicity)~step,
               data = d_cells[step < step_lim & Cell == 1, ])
apparent_death_rate_S <- -1*coefficients(lm_mod_S)["step"]

newdata_S <- data.frame(step = d_cells[step < step_lim& Cell == 0, ]$step,
                        multiplicity = exp(predict(lm_mod_S)))


plot(plot_mod1 +
       ggtitle("Basic resistance (Log10)")+
       scale_y_continuous(trans = "log10")+
       geom_line(data = newdata_S,
                 mapping = aes(x = step, y = multiplicity),
                 color = "blue",
                 linetype = "dashed",
                 inherit.aes = FALSE,
                 linewidth = 1.2)+
       geom_text(data = cbind.data.frame(
         newdata_S[nrow(newdata_S),],
         label = paste("Exponential death\n with rate ",
                       round(apparent_death_rate_S, digit = 3),
                       " per unit time", sep = "")
       ),
       mapping = aes(x = step, y = multiplicity, label =label),
       color = "blue",
       hjust = 0, vjust = - 1,
       inherit.aes = FALSE)+
       geom_line(data = newdata_R,
                 mapping = aes(x = step, y = multiplicity),
                 color = "red",
                 linetype = "dashed",
                 inherit.aes = FALSE,
                 linewidth = 1.2)+
       geom_text(data = cbind.data.frame(
         newdata_R[nrow(newdata_R),],
         label = paste("Exponential death\n with rate ",
                       round(apparent_death_rate_R, digit = 3),
                       " per unit time", sep = "")
       ),
       mapping = aes(x = step, y = multiplicity, label =label),
       color = "red",
       hjust = 0, vjust = - 1,
       inherit.aes = FALSE)
)






##----------------------------------------------------------------------------##
## MODEL 5 - PLASMID LOSS
##----------------------------------------------------------------------------##

## INPUTS -----------------------

# model parameters
plasmid_loss <- 0.05
popsize <- 1e9
capacity <- 1e9

# running parameters
csvres_folder <- "basic_plasmid_loss"
n_steps <- 1
n_writes <- 700



## MODEL CONFIGURATION ----------

prm <- emptymodel() %>%
  plasmidArchetype(0, plasmid_loss, 0.0, 10, 1.0) %>% plasmid(0, 0) %>%
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
res <- mse_run(prm, csvres_folder, n_steps, n_writes, s1= sample(1E9,1))

# if needed, cleaning the csv results
if(cleanup_csv){ mse1.cleanup_csv(csvres_folder)}

# if needed, save the R result (res) as a RDS object
if(save_r_res){
  saveRDS(object = res, file = paste(r_res_path,
                                     "/", csvres_folder, ".rds", sep = ""))
}



## RESULTS ----------------------

## Diversification of cells thanks to plasmid loss
## Only one plasmid loose per cell per time step

mse1.plot_graph(mse1.extract_prm_from_mseres(res, 0),
                add_archetype = FALSE, edges_type = "inclusion")

mse1.plot_graph(mse1.extract_prm_from_mseres(res, 1),
                add_archetype = FALSE, edges_type = "inclusion")

mse1.plot_graph(mse1.extract_prm_from_mseres(res, 2),
                add_archetype = FALSE, edges_type = "inclusion")

mse1.plot_graph(mse1.extract_prm_from_mseres(res, 3),
                add_archetype = FALSE, edges_type = "inclusion")

mse1.plot_graph(mse1.extract_prm_from_mseres(res, 4),
                add_archetype = FALSE, edges_type = "inclusion")

mse1.plot_graph(mse1.extract_prm_from_mseres(res, 5),
                add_archetype = FALSE, edges_type = "inclusion")


## Caracterization of cells generating during the simulations
cells <- mse1.describe_cells(res)
print(cells)


## Plotting cells dynamics
d_cells <- count_from("Cell", "Patch", res)

colPal <- rainbow(6)
names(colPal) <- as.character(0:5)

plot_cells <- ggplot(data = d_cells,
                       mapping = aes(x =step, y = multiplicity,
                                     group = Cell, color = as.factor(Cell)))+
  geom_line(linewidth =1.2)+
  xlab("Time")+
  ylab("Cell count")+
  theme_bw()+
  scale_color_manual(name = "Cell type",
                     values =colPal,
                     labels = paste("Cells with ", 5:0, " plasmids", sep = ""))



plot(plot_cells + ggtitle("Basic plasmid loss : Cell dynamics"))
plot(plot_cells + ggtitle("Basic plasmid loss : Cell dynamics (log10)")+ scale_y_continuous(trans = "log10"))



## Dynamics of the total number of plasmids
d_plasmids <- count_from("Plasmid", "Cell", res)
d_plasmids[, multiplicity := as.numeric(multiplicity)] ## ! Trick to handle large number.
d_plasmids <- count_in(d_plasmids, "Patch", res)

plot_plasmid <- ggplot(data = d_plasmids,
                    mapping = aes(x =step, y = multiplicity,
                                  group = Plasmid, color = as.factor(Plasmid)))+
  geom_line(linewidth =1.2)+
  xlab("Time")+
  ylab("Plasmid count")+
  theme_bw()

plot(plot_plasmid + ggtitle("Basic plasmid loss : plasmid dynamics"))
plot(plot_plasmid + ggtitle("Basic plasmid loss : plasmid dynamics (log10)")+ scale_y_continuous(trans = "log10"))




##----------------------------------------------------------------------------##
## MODEL 6 - PLASMID TRANSFER
##----------------------------------------------------------------------------##

## INPUTS -----------------------

# model parameters
plasmid_transfer <- 0.01
plasmid_maxcount <- 5
capacity <- 1e9
popsize_plasmid <- 1e7


# running parameters
csvres_folder <- "basic_plasmid_transfer"
n_steps <- 1
n_writes <- 1500



## MODEL CONFIGURATION ----------

prm <- emptymodel() %>%
  plasmidArchetype(0, 0.0, plasmid_transfer, plasmid_maxcount, 1.0) %>%
  plasmid(0, 0) %>%
  chromosomeArchetype(0, fitness = 0.0, survival = 1.0) %>% chromosome(0, 0) %>%
  cellArchetype(0) %>% cell(0, 0) %>% cell(1,0) %>%
  patchArchetype(0, capacity, c(0.0, 0.0)) %>% patch(0, 0) %>%
  assign("Plasmid", 0, "Cell", 1, mult = 1) %>%
  assign("Chromosome", 0, "Cell", 0) %>%
  assign("Chromosome", 0, "Cell", 1) %>%
  assign("Cell", 1, "Patch", 0, popsize_plasmid) %>%
  assign("Cell", 0, "Patch", 0, capacity-popsize_plasmid) %>%
  plasmidRange(0, 0)

mse1.plot_graph(prm)
mse1.diagnose_model(prm)


## MODEL RUNNING ----------------

# launching simulator
res <- mse_run(prm, csvres_folder, n_steps, n_writes, s1= sample(1E9,1))

# if needed, cleaning the csv results
if(cleanup_csv){ mse1.cleanup_csv(csvres_folder)}

# if needed, save the R result (res) as a RDS object
if(save_r_res){
  saveRDS(object = res, file = paste(r_res_path,
                                     "/", csvres_folder, ".rds", sep = ""))
}



## RESULTS ----------------------

## Caracterization of cells generating during the simulations
cells <- mse1.describe_cells(res)
print(cells)


## Plotting cells dynamics
d_cells <- count_from("Cell", "Patch", res)
d_cells_prototype <- merge(d_cells, cells[, c("Cell", "Plasmid_0")],
                           by = "Cell", all.x = TRUE)

d_cells_prototype <- d_cells_prototype[,
                                       .(multiplicity = sum(multiplicity)),
                                       by = list(Plasmid_0,
                                                 Patch, step)]

d_cells_prototype <- mse1.add_zero_steps(d_cells_prototype)


colPal <- rainbow(6)
names(colPal) <- as.character(0:5)

plot_cells <- ggplot(data = d_cells_prototype,
                     mapping = aes(x =step, y = multiplicity,
                                   group = Plasmid_0,
                                   color = as.factor(Plasmid_0)))+
  geom_line(linewidth =1.2)+
  xlab("Time")+
  ylab("Cell count")+
  theme_bw()+
  scale_color_manual(name = "Cell type",
                     values =colPal,
                     labels = paste("Cells with ", 0:5, " plasmids", sep = ""))



plot(plot_cells + ggtitle("Basic plasmid accumulation : Cell dynamics"))
plot(plot_cells + ggtitle("Basic plasmid accumulation : Cell dynamics (log10)")+
       scale_y_continuous(trans = "log10"))



## Dynamics of the total number of plasmids
d_plasmids <- count_from("Plasmid", "Cell", res)
d_plasmids[, multiplicity := as.numeric(multiplicity)] ## ! Trick to handle large number.
d_plasmids <- count_in(d_plasmids, "Patch", res)

plot_plasmid <- ggplot(data = d_plasmids,
                       mapping = aes(x =step, y = multiplicity,
                                     group = Plasmid, color = as.factor(Plasmid)))+
  geom_line(linewidth =1.2)+
  xlab("Time")+
  ylab("Plasmid count")+
  theme_bw()

plot(plot_plasmid + ggtitle("Basic plasmid accumulation : plasmid dynamics"))
plot(plot_plasmid + ggtitle("Basic plasmid accumulation : plasmid dynamics (log10)")+ scale_y_continuous(trans = "log10"))


## Illustrating the number of transfer

n_donnors <- d_cells_prototype[Plasmid_0>0,
                               .(n_donnors = sum(multiplicity)),
                               by = step]

mean_plasmid_per_donnor <- d_cells_prototype[Plasmid_0>0,
                               .(mean_plas = sum(multiplicity*Plasmid_0)/sum(multiplicity)),
                               by = step]

n_recipients <- d_cells_prototype[Plasmid_0<plasmid_maxcount,
                                  .(n_recipients = sum(multiplicity)),
                                  by = step]


dt <- merge(n_donnors, mean_plasmid_per_donnor, by = "step")
dt<- merge(dt, n_recipients, by = "step")
dt[, n_expected_transfer := (n_donnors/capacity)*
     (1- (1-plasmid_transfer)^mean_plas)*n_recipients]


eff_trans <- data.table(step = d_plasmids$step[-1],
                         n_transfer = diff(d_plasmids$multiplicity))

ggplot(data = eff_trans,
       mapping = aes(x = step, y = n_transfer))+
  geom_line(col = "grey", linewidth = 1.2)+
  theme_bw()+
  geom_line(data = dt,
            mapping = aes(x = step, y = n_expected_transfer ),
            col = "red", linewidth = 1.2, linetype = "dotted")+
  ggtitle("Number of transfers")




##----------------------------------------------------------------------------##
## MODEL 7 - INTRA-PATCH PLASMID TRANSFER
##----------------------------------------------------------------------------##
## Similar to model 6, but with an additional patch consisting only of cells
## without plasmids. While these cells can accept the plasmid, they are unable
## to acquire it as they are not in contact with it.


## INPUTS -----------------------

# model parameters
plasmid_transfer <- 0.01
plasmid_maxcount <- 5
capacity <- 1e9
popsize_plasmid <- 1e7


# running parameters
csvres_folder <- "basic_plasmid_transfer_2"
n_steps <- 1
n_writes <- 1500



## MODEL CONFIGURATION ----------

prm <- emptymodel() %>%
  plasmidArchetype(0, 0.0, plasmid_transfer, plasmid_maxcount, 1.0) %>%
  plasmid(0, 0) %>%
  chromosomeArchetype(0, fitness = 0.0, survival = 1.0) %>% chromosome(0, 0) %>%
  cellArchetype(0) %>% cell(0, 0) %>% cell(1,0) %>%
  patchArchetype(0, capacity, c(0.0, 0.0)) %>% patch(0, 0) %>%patch(1, 0) %>%
  assign("Plasmid", 0, "Cell", 1, mult = 1) %>%
  assign("Chromosome", 0, "Cell", 0) %>%
  assign("Chromosome", 0, "Cell", 1) %>%
  assign("Cell", 1, "Patch", 0, popsize_plasmid) %>%
  assign("Cell", 0, "Patch", 0, capacity-popsize_plasmid) %>%
  assign("Cell", 0, "Patch", 1, capacity-popsize_plasmid) %>%
  plasmidRange(0, 0)

mse1.plot_graph(prm)
mse1.diagnose_model(prm)


## MODEL RUNNING ----------------

# launching simulator
res <- mse_run(prm, csvres_folder, n_steps, n_writes, s1= sample(1E9,1))

# if needed, cleaning the csv results
if(cleanup_csv){ mse1.cleanup_csv(csvres_folder)}

# if needed, save the R result (res) as a RDS object
if(save_r_res){
  saveRDS(object = res, file = paste(r_res_path,
                                     "/", csvres_folder, ".rds", sep = ""))
}



## RESULTS ----------------------

## Caracterization of cells generating during the simulations
cells <- mse1.describe_cells(res)
print(cells)


## Plotting cells dynamics
d_cells <- count_from("Cell", "Patch", res)
d_cells_prototype <- merge(d_cells, cells[, c("Cell", "Plasmid_0")],
                           by = "Cell", all.x = TRUE)

d_cells_prototype <- d_cells_prototype[,
                                       .(multiplicity = sum(multiplicity)),
                                       by = list(Plasmid_0,
                                                 Patch, step)]

d_cells_prototype <- mse1.add_zero_steps(d_cells_prototype)


colPal <- rainbow(6)
names(colPal) <- as.character(0:5)

plot_cells <- ggplot(data = d_cells_prototype,
                     mapping = aes(x =step, y = multiplicity,
                                   group = Plasmid_0,
                                   color = as.factor(Plasmid_0)))+
  geom_line(linewidth =1.2)+
  xlab("Time")+
  ylab("Cell count")+
  theme_bw()+
  scale_color_manual(name = "Cell type",
                     values =colPal,
                     labels = paste("Cells with ", 0:5, " plasmids", sep = ""))+
  facet_wrap("Patch")



plot(plot_cells + ggtitle("Basic plasmid accumulation : Cell dynamics"))
plot(plot_cells + ggtitle("Basic plasmid accumulation : Cell dynamics (log10)")+
       scale_y_continuous(trans = "log10"))





##----------------------------------------------------------------------------##
## MODEL 8 - PLASMID HOST RANGE
##----------------------------------------------------------------------------##
## Similar to model 6, but with two additional cell lines. Plasmids can be
## taken up by cells '0' (initial donors),  cells '1' (same strain as '0', but
## initially plasmid-free), by cell '2' (other strain, as indicated by the
## carriage of another chromosome, and initially plasmid-free), but not by cells
## '3' (again, other strains). Since the behavior of plasmids is independent of
## the host cell (whether cells '0' or '1'), the plasmid transfer rate remains
## the same regardless of the donor or recipient.


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
res <- mse_run(prm, csvres_folder, n_steps, n_writes, s1= sample(1E9,1))

# if needed, cleaning the csv results
if(cleanup_csv){ mse1.cleanup_csv(csvres_folder)}

# if needed, save the R result (res) as a RDS object
if(save_r_res){
  saveRDS(object = res, file = paste(r_res_path,
                                     "/", csvres_folder, ".rds", sep = ""))
}



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
  theme_bw()+
  scale_color_manual(name = "Cell type",
                     values =c(
                       "0" = "blue",
                       "1" = "green",
                       "2" = "red"
                     ),
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



plot(plot_cells + ggtitle("Basic plasmid host range"))

plot(plot_cells + ggtitle("Basic plasmid host range (log10)")+
       scale_y_continuous(trans = "log10"))



## Dynamics of the total number of plasmids
d_plasmids <- count_from("Plasmid", "Cell", res)
d_plasmids[, multiplicity := as.numeric(multiplicity)] ## ! Trick to handle large number.
d_plasmids <- count_in(d_plasmids, "Patch", res)

plot_plasmid <- ggplot(data = d_plasmids,
                       mapping = aes(x =step, y = multiplicity,
                                     group = Plasmid, color = as.factor(Plasmid)))+
  geom_line(linewidth =1.2)+
  xlab("Time")+
  ylab("Plasmid count")+
  theme_bw()

plot(plot_plasmid + ggtitle("Basic plasmid accumulation : plasmid dynamics"))
plot(plot_plasmid + ggtitle("Basic plasmid accumulation : plasmid dynamics (log10)")+ scale_y_continuous(trans = "log10"))


## Illustrating the number of transfer

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

ggplot(data = eff_trans,
       mapping = aes(x = step, y = n_transfer))+
  geom_line(col = "grey")+
  theme_bw()+
  geom_line(data = dt,
            mapping = aes(x = step, y = n_expected_transfer ),
            col = "red", linewidth = 1.2)+
  ggtitle("Number of transfers")






##----------------------------------------------------------------------------##
## MODEL 9 - PLASMID LIMITATION ~
##----------------------------------------------------------------------------##
## mse1 introduces a limitation on the number of plasmids a cell can accumulate.
## This limitation applies to the number of plasmids with a given archetype,
## regardless of their specific content. For example, two plasmids with the same
## backbone but different resistance genes cannot coexist in a cell if the max_count
## for that backbone is set to 1.
##
## This functions similarly to an incompatibility group, though under the strong
## assumption that plasmids from the same incompatibility group share identical
## characteristics in terms of host range, loss rates, and conjugation rates.
##




## INPUTS -----------------------

# model parameters
plasmid_transfer <- 0.01
capacity <- 1e6


# running parameters
csvres_folder <- "basic_plasmid_limitation"
n_steps <- 1
n_writes <- 1500



## MODEL CONFIGURATION ----------

prm <- emptymodel() %>%
  geneArchetype(0, susc = c(1.0, 1.0), fitness = 1)%>% gene(0, 0)%>%
  plasmidArchetype(0, 0.0, plasmid_transfer, max_count = 1, 1.0) %>%
  plasmid(0, 0) %>% plasmid(1,0)%>%
  plasmidArchetype(1, 0.0, plasmid_transfer, max_count = 2, 1.0) %>%
  plasmid(2, 1) %>% plasmid(3,1)%>%
  chromosomeArchetype(0, fitness = 0.0, survival = 1.0) %>% chromosome(0, 0) %>%
  cellArchetype(0) %>% cell(0, 0) %>% cell(1,0) %>%
  cell(2,0) %>%cell(3,0) %>%
  patchArchetype(0, capacity, c(0.0, 0.0)) %>% patch(0, 0) %>%
  assign("Gene", 0, "Plasmid", 0) %>%
  assign("Gene", 0, "Plasmid", 2) %>%
  assign("Plasmid", 0, "Cell", 0, mult = 1) %>%
  assign("Plasmid", 1, "Cell", 1, mult = 1) %>%
  assign("Plasmid", 2, "Cell", 2, mult = 1) %>%
  assign("Plasmid", 3, "Cell", 3, mult = 1) %>%
  assign("Chromosome", 0, "Cell", 0) %>%
  assign("Chromosome", 0, "Cell", 1) %>%
  assign("Chromosome", 0, "Cell", 2) %>%
  assign("Chromosome", 0, "Cell", 3) %>%
  assign("Cell", 0, "Patch", 0, floor(capacity/4)) %>%
  assign("Cell", 1, "Patch", 0, floor(capacity/4)) %>%
  assign("Cell", 2, "Patch", 0, floor(capacity/4)) %>%
  assign("Cell", 3, "Patch", 0, floor(capacity/4)) %>%
  plasmidRange(0, 0)%>%
  plasmidRange(1,0)


mse1.plot_graph(prm, layout = "fr")
mse1.diagnose_model(prm)



## MODEL RUNNING ----------------

# launching simulator
res <- mse_run(prm, csvres_folder, n_steps, n_writes, s1= sample(1E9,1))

# if needed, cleaning the csv results
if(cleanup_csv){ mse1.cleanup_csv(csvres_folder)}

# if needed, save the R result (res) as a RDS object
if(save_r_res){
  saveRDS(object = res, file = paste(r_res_path,
                                     "/", csvres_folder, ".rds", sep = ""))
}




## RESULTS ----------------------

## Caracterization of cells generating during the simulations
cells <- mse1.describe_cells(res)
print(cells)
DT::datatable(cells)






##----------------------------------------------------------------------------##
## MODEL 10 - PATCH-TO-PATCH BACTERIAL TRANSFER (SYMMETRICAL)
##----------------------------------------------------------------------------##


## INPUTS -----------------------

# model parameters
diffusion_probability <- 0.005
capacity <- 1e9
popsize <- 1e9


# running parameters
csvres_folder <- "basic_diffusion"
n_steps <- 1
n_writes <- 700



## MODEL CONFIGURATION ----------

prm <- emptymodel() %>%
  chromosomeArchetype(0, fitness = 0.0, survival = 1.0) %>% chromosome(0, 0) %>%
  cellArchetype(0) %>% cell(0, 0) %>%
  patchArchetype(0, capacity, c(0.0, 0.0)) %>% patch(0, 0) %>% patch(1, 0) %>%
  diffusion(0, 1, diffusion_probability, symmetric = TRUE ) %>%
  assign("Chromosome", 0, "Cell", 0) %>%
  assign("Cell", 0, "Patch", 0, popsize)

mse1.plot_graph(prm)
mse1.diagnose_model(prm)



## MODEL RUNNING ----------------

# launching simulator
res <- mse_run(prm, csvres_folder, n_steps, n_writes, s1= sample(1E9,1))

# if needed, cleaning the csv results
if(cleanup_csv){ mse1.cleanup_csv(csvres_folder)}

# if needed, save the R result (res) as a RDS object
if(save_r_res){
  saveRDS(object = res, file = paste(r_res_path,
                                     "/", csvres_folder, ".rds", sep = ""))
}



## RESULTS ----------------------

## Plotting cells dynamics
d_cells <- count_from("Cell", "Patch", res)


plot_cells <- ggplot(data = d_cells,
                     mapping = aes(x =step, y = multiplicity,
                                   group = interaction(Cell, Patch),
                                   color = as.factor(Patch)))+
  geom_line(linewidth =1.2)+
  xlab("Time")+
  ylab("Cell count")+
  theme_bw()+
  scale_color_manual(name = "Patch",
                     values =c("0" = "red", "1" ="blue"),
                     labels = c("Donnor patch", "Recipient patch"))


plot(plot_cells + ggtitle("Basic patch-to-patch cell diffusion (symmetric)"))
plot(plot_cells + ggtitle("Basic patch-to-patch cell diffusion (symmetric)(log10)")+
       scale_y_continuous(trans = "log10"))



## "Spatial evolution"
## TO DO







##----------------------------------------------------------------------------##
## MODEL 11 - FITNESS COST INDUCED BY GENES OR PLASMIDS
##----------------------------------------------------------------------------##
## We simulate four different cells, with identical basal growth rate, but
## carrying genes and/or plasmids modulating their growth:
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
capacity <- 1e9

# running parameters
csvres_folder <- "basic_fitness"
n_steps <- 1
n_writes <- 300



## MODEL CONFIGURATION ----------

prm <- emptymodel() %>%
  geneArchetype(0, c(1.0, 1.0), gene_fitness) %>% gene(0, 0) %>%
  chromosomeArchetype(0, fitness = basal_growth_probability, survival = 1.0) %>%
  chromosome(0, 0) %>% chromosome(1,0)%>%
  plasmidArchetype(0, 0.0, 0.0, 1, plasmid_fitness) %>%
  plasmid(0, 0) %>% plasmid(1, 0) %>%
  cellArchetype(0) %>% cell(0, 0) %>% cell(1,0)%>% cell(2,0)%>%cell(3,0) %>%
  patchArchetype(0, capacity, c(0.0, 0.0)) %>%
  patch(0, 0) %>% patch(1, 0) %>% patch(2, 0) %>%patch(3, 0) %>%
  patch(4, 0) %>%
  assign("Gene", 0, "Chromosome", 1) %>%
  assign("Gene", 0, "Plasmid", 1) %>%
  assign("Chromosome", 0, "Cell", 0) %>%
  assign("Chromosome", 1, "Cell", 1) %>%
  assign("Chromosome", 0, "Cell", 2) %>%
  assign("Chromosome", 0, "Cell", 3) %>%
  assign("Plasmid", 0, "Cell", 2) %>%
  assign("Plasmid", 1, "Cell", 3) %>%
  assign("Cell", 0, "Patch", 0, popsize)%>%
  assign("Cell", 1, "Patch", 1, popsize)%>%
  assign("Cell", 2, "Patch", 2, popsize)%>%
  assign("Cell", 3, "Patch", 3, popsize)%>%
  assign("Cell", 0, "Patch", 4, popsize/4)%>%
  assign("Cell", 1, "Patch", 4, popsize/4)%>%
  assign("Cell", 2, "Patch", 4, popsize/4)%>%
  assign("Cell", 3, "Patch", 4, popsize/4)

mse1.plot_graph(prm)
mse1.diagnose_model(prm)


## MODEL RUNNING ----------------

# launching simulator
res <- mse_run(prm, csvres_folder, n_steps, n_writes, s1= sample(1E9,1))

# if needed, cleaning the csv results
if(cleanup_csv){ mse1.cleanup_csv(csvres_folder)}

# if needed, save the R result (res) as a RDS object
if(save_r_res){
  saveRDS(object = res, file = paste(r_res_path,
                                     "/", csvres_folder, ".rds", sep = ""), )
}


## RESULTS ----------------------

## computing intra-patches dynamic of cells

d_cells <- count_from("Cell", "Patch", res)

## plotting curves of independant growths
plot_mod_indiv <- ggplot(data = d_cells[Patch <4, ],
                         mapping = aes(x =step, y = multiplicity,
                                       group = interaction(Cell, Patch),
                                       color = as.factor(Cell)))+
  geom_line(linewidth =1.2)+
  xlab("Time")+
  ylab("Cell count")+
  theme_bw()+
  scale_color_manual(name = "Cell type",
                     values =c("0" = "green",
                               "1" = "blue",
                               "2" = "red",
                               "3" = "purple"),
                     labels = c("WT",
                                "Chromosomal ARG",
                                "Plasmid",
                                "Plasmid + ARG"))


plot(plot_mod_indiv +
       ggtitle("Fitness cost : independant growth of cells in identical conditions"))

plot(plot_mod_indiv +
       ggtitle("Fitness cost : independant growth of cells in identical conditions")
     + scale_y_continuous(trans = "log10"))



## plotting curves of growth in a shared environment
plot_mod_all <- ggplot(data = d_cells[Patch ==4 , ],
                       mapping = aes(x =step, y = multiplicity,
                                     group = interaction(Cell, Patch),
                                     color = as.factor(Cell)))+
  geom_line(linewidth =1.2)+
  xlab("Time")+
  ylab("Cell count")+
  theme_bw()+
  scale_color_manual(name = "Cell type",
                     values =c("0" = "green",
                               "1" = "blue",
                               "2" = "red",
                               "3" = "purple"),
                     labels = c("WT",
                                "Chromosomal ARG",
                                "Plasmid",
                                "Plasmid + ARG"))


plot(plot_mod_all +
       ggtitle("Fitness cost : dependant growth of cells with shared ressources"))

plot(plot_mod_all +
       ggtitle("Fitness cost : dependant growth of cells with shared ressources (log10)")
     + scale_y_continuous(trans = "log10"))



## Verification of independant growth rates

step_lim <- 100

simulated_growth_rate_0 <- log(1+basal_growth_probability)/1
simulated_growth_rate_0

lm_mod_0 <- lm( log(multiplicity)~step,
              data = d_cells[step<step_lim& Cell ==0, ])

apparent_growth_rate_0 <- coefficients(lm_mod_0)["step"]
apparent_growth_rate_0


simulated_growth_rate_1 <- log(1+basal_growth_probability*gene_fitness)/1
simulated_growth_rate_1

lm_mod_1 <- lm( log(multiplicity)~step,
                data = d_cells[step<step_lim& Cell ==1, ])

apparent_growth_rate_1 <- coefficients(lm_mod_1)["step"]
apparent_growth_rate_1



simulated_growth_rate_2 <- log(1+basal_growth_probability*plasmid_fitness)/1
simulated_growth_rate_2

lm_mod_2 <- lm( log(multiplicity)~step,
                data = d_cells[step<step_lim& Cell ==2, ])

apparent_growth_rate_2 <- coefficients(lm_mod_2)["step"]
apparent_growth_rate_2


simulated_growth_rate_3 <- log(1+basal_growth_probability*plasmid_fitness*gene_fitness)/1
simulated_growth_rate_3

lm_mod_3 <- lm( log(multiplicity)~step,
                data = d_cells[step<step_lim& Cell ==3, ])

apparent_growth_rate_3 <- coefficients(lm_mod_3)["step"]
apparent_growth_rate_3




