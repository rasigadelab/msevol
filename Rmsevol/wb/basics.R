#################################################
# BASIC SCENARIOS
# Use as consistency checks

library(data.table)
library(tidyr)
library(bit64)
library(Rmsevol)

#################################################
# CELL BIRTH

popsize <- 1e10
capacity <- 1e12

prm <- emptymodel() %>%
  chromosomeArchetype(0, fitness = 0.05, survival = 1.0) %>% chromosome(0, 0) %>%
  cellArchetype(0) %>% cell(0, 0) %>%
  patchArchetype(0, capacity, c(0.0, 0.0)) %>% patch(0, 0) %>%
  assign("Chromosome", 0, "Cell", 0) %>%
  assign("Cell", 0, "Patch", 0, popsize)

prm

res <- mse_run(prm, "basic_birth", 1, 1000)

d <- count_from("Cell", "Patch", res)

plot(d$step, d$multiplicity, type = "l", xlab = "Time", ylab = "Cell count", log = "")



#################################################
# CELL DEATH

popsize <- 1e12
capacity <- 1e12

prm <- emptymodel() %>%
  chromosomeArchetype(0, fitness = 0, survival = 0.95) %>% chromosome(0, 0) %>%
  cellArchetype(0) %>% cell(0, 0) %>%
  patchArchetype(0, capacity, c(0.0, 0.0)) %>% patch(0, 0) %>%
  assign("Chromosome", 0, "Cell", 0) %>%
  assign("Cell", 0, "Patch", 0, popsize)

prm

res <- mse_run(prm, "basic_death", 1, 1000)

d <- count_from("Cell", "Patch", res)

plot(d$step, d$multiplicity, type = "l", xlab = "Time", ylab = "Cell count", log = "")

#################################################
# CELL KILLING

popsize <- 1e12
capacity <- 1e12

prm <- emptymodel() %>%
  chromosomeArchetype(0, fitness = 0.1, survival = 0.99) %>% chromosome(0, 0) %>%
  cellArchetype(0) %>% cell(0, 0) %>%
  patchArchetype(0, capacity, c(0.1, 0.0)) %>% patch(0, 0) %>%
  assign("Chromosome", 0, "Cell", 0) %>%
  assign("Cell", 0, "Patch", 0, popsize)

prm

res <- mse_run(prm, "basic_kill", 1, 500)

d <- count_from("Cell", "Patch", res)

plot(d$step, d$multiplicity, type = "l", xlab = "Time", ylab = "Cell count", log = "")

#################################################
# CELL RESISTANCE

popsize <- 1e10
capacity <- 1e12

# Resistant vs suceptible population

prm <- emptymodel() %>%
  geneArchetype(0, c(0.8, 0.0), 1.0) %>% gene(0, 0) %>%
  chromosomeArchetype(0, fitness = 0.1, survival = 0.99) %>% chromosome(0, 0) %>% chromosome(1, 0) %>%
  cellArchetype(0) %>% cell(0, 0) %>% cell(1, 0) %>%
  patchArchetype(0, capacity, c(0.1, 0.0)) %>% patch(0, 0) %>%
  assign("Gene", 0, "Chromosome", 1) %>%
  assign("Chromosome", 0, "Cell", 0) %>%
  assign("Chromosome", 1, "Cell", 1) %>%
  assign("Cell", 0, "Patch", 0, popsize) %>%
  assign("Cell", 1, "Patch", 0, popsize)

prm

res <- mse_run(prm, "basic_resistance", 1, 2000)

d <- count_from("Cell", "Patch", res)

plot(d[Cell == 1]$step, d[Cell == 1]$multiplicity, type = "l", xlab = "Time", ylab = "Cell count",
     log = "", col = "red", ylim = range(as.numeric(d$multiplicity)))

lines(d[Cell == 0]$step, d[Cell == 0]$multiplicity, col = "darkgreen")
legend("bottomleft", bty = "n", legend = c("S", "R"), fill = c("darkgreen", "red"))

##################################################
# CELL DIFFUSION

popsize <- 1e10
capacity <- 1e12

prm <- emptymodel() %>%
  chromosomeArchetype(0, fitness = 0.0, survival = 1.0) %>% chromosome(0, 0) %>%
  cellArchetype(0) %>% cell(0, 0) %>%
  patchArchetype(0, capacity, c(0.0, 0.0)) %>% patch(0, 0) %>% patch(1, 0) %>%
  diffusion(0, 1, 0.2) %>%
  assign("Chromosome", 0, "Cell", 0) %>%
  assign("Cell", 0, "Patch", 0, popsize)

prm

res <- mse_run(prm, "basic_diffusion", 1, 200)

d <- count_from("Cell", "Patch", res)

plot(d[Patch == 0]$step, d[Patch == 0]$multiplicity, type = "l", xlab = "Time", ylab = "Cell count", log = "", col = "red",
     ylim = range(as.numeric(d$multiplicity)))
lines(d[Patch == 1]$step, d[Patch == 1]$multiplicity, col = "darkgreen")
legend("bottomleft", bty = "n", legend = c("Donor patch", "Recipient patch"), fill = c("darkgreen", "red"))

##################################################
# PLASMID LOSS

popsize <- 1e10
capacity <- 1e12


prm <- emptymodel() %>%
  plasmidArchetype(0, 0.05, 0.0, 100, 1.0) %>% plasmid(0, 0) %>%
  chromosomeArchetype(0, fitness = 0.0, survival = 1.0) %>% chromosome(0, 0) %>%
  cellArchetype(0) %>% cell(0, 0) %>%
  patchArchetype(0, capacity, c(0.0, 0.0)) %>% patch(0, 0) %>%
  assign("Plasmid", 0, "Cell", 0, mult = 1) %>%
  assign("Chromosome", 0, "Cell", 0) %>%
  assign("Cell", 0, "Patch", 0, popsize)

prm

res <- mse_run(prm, "basic_plasmid_loss", 1, 200)

d <- count_from("Plasmid", "Cell", res) %>% count_in("Patch", res)

# d <- merge(count_from("Plasmid", "Cell", res), count_from("Cell", "Patch", res), by = c("Cell", "step"))

plot(d$step, d$multiplicity, type = "l", xlab = "Time", ylab = "Plasmid count", log = "", col = "red",
     ylim = range(as.numeric(d$multiplicity)))


##################################################
# PLASMID ACCUMULATION

popsize <- 1e10
capacity <- 1e10

prm <- emptymodel() %>%
  plasmidArchetype(0, 0.0, 0.01, 5, 1.0) %>% plasmid(0, 0) %>%
  chromosomeArchetype(0, fitness = 0.0, survival = 1.0) %>% chromosome(0, 0) %>%
  cellArchetype(0) %>% cell(0, 0) %>%
  patchArchetype(0, capacity, c(0.0, 0.0)) %>% patch(0, 0) %>%
  assign("Plasmid", 0, "Cell", 0, mult = 1) %>%
  assign("Chromosome", 0, "Cell", 0) %>%
  assign("Cell", 0, "Patch", 0, popsize) %>%
  plasmidRange(0, 0)

prm

res <- mse_run(prm, "basic_plasmid_tranfer", 1, 500)

d <- count_from("Plasmid", "Cell", res) %>% count_in("Patch", res)

# d <- merge(count_from("Plasmid", "Cell", res), count_from("Cell", "Patch", res), by = c("Cell", "step"))

plot(d$step, d$multiplicity, type = "l", xlab = "Time", ylab = "Plasmid count",
     log = "", col = "red", ylim = range(as.numeric(d$multiplicity)))


