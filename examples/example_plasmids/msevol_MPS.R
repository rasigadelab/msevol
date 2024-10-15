################################################################################
## MSE1 RESISTRACK MODEL
##  A bacterial hospital ecosystem simulation framework
## (c)2018-2024 Jean-Philippe Rasigade, <jean-philippe.rasigade@univ-lyon1.fr>
##  University of Lyon, Centre International de Recherche en Infectiologie
################################################################################
## Plasmid persistence under antibiotic pressure :
##   Determination of Minimum Persistence Threshold
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
## Persistence of costly plasmids under antibiotics resistance
##----------------------------------------------------------------------------##

single_run_mod <- function(
  killing,
  susc,
  fitness,
  loss,
  transfer,
  basal_growth,
  basal_surv,
  capacity,
  popsize,
  bin_folder,
  n_steps,
  n_writes
){

  ## Configure Model
  prm <- emptymodel()%>%
    geneArchetype(0, c(susc, 1.0), fitness) %>% gene(0, 0) %>%
    plasmidArchetype(0, loss = loss, transfer = transfer,
                     max_count = 1, fitness = 1.0) %>% plasmid(0,0) %>%
    assign("Gene", 0, "Plasmid", 0)%>%
    chromosomeArchetype(0, fitness = basal_growth, survival = basal_surv) %>%
    chromosome(0, 0) %>%
    cellArchetype(0) %>% cell(0,0)%>%
    assign("Chromosome", 0, "Cell", 0) %>%
    assign("Plasmid",0, "Cell", 0) %>%
    plasmidRange(0,0)%>%
    patchArchetype(0, capacity, c(killing, 0.0)) %>%
    patch(0, 0) %>%
    assign("Cell", 0, "Patch", 0, popsize)


  ## running simulation
  res <- mse_run(prm,
                 bin_folder,
                 n_steps,
                 n_writes, s1= sample(1E9,1))
  mse1.cleanup_csv(bin_folder)


  ## determine persistence
  celldesc <- mse1.describe_cells(res)
  celldesc[Plasmid_0 == 1, label:= "Plasmid"]
  celldesc[Plasmid_0 == 0, label:= "Sensitive"]

  d_cells <- count_from("Cell", "Patch", res)
  d_cells <- merge(d_cells, celldesc[, c("Cell", "label")], by = "Cell")



  if("Plasmid" %in% d_cells[step == max(step), label]){

    return( data.table(killing = killing,
                       persistence = TRUE,
                       level = d_cells[step == max(step)& label == "Plasmid",]$multiplicity,
                       prop = d_cells[step == max(step)& label == "Plasmid"]$multiplicity/
                         sum(d_cells[step == max(step), ]$multiplicity),
                       step = max(d_cells$step))
    )


  }else{
    return(
      data.table(killing = killing,
                 persistence = FALSE,
                 level = 0,
                 prop = 0,
                 step = max(d_cells$step))
    )
  }

}



estimate_msp <- function(args_susc,
                         args_fitness,
                         args_loss,
                         args_transfer,
                         args_basal_growth =0.1,
                         args_basal_surv =0.99,
                         args_capacity =1E9,
                         args_popsize = 1E3,
                         args_max_round = 4,
                         args_digit = 2,
                         args_bin_folder = "tmp",
                         args_n_steps = 100000,
                         args_n_writes =1
){


  ## rough estimation for killing
  ##-------------------------------------------

  current_kill <- 1
  persistence <- data.table()

  for(round in 1:args_max_round){
    for(i in 1:5){

      cat("\n[", round, ":", i, "]")
      current_kill <- current_kill/10
      persistence <- rbind(persistence,
                           single_run_mod( current_kill,
                                           args_susc,
                                           args_fitness,
                                           args_loss,
                                           args_transfer,
                                           args_basal_growth,
                                           args_basal_surv,
                                           args_capacity,
                                           args_popsize,
                                           args_bin_folder,
                                           args_n_steps,
                                           args_n_writes)
      )

    }

    ## test for another rough round
    if(!(all(persistence$persistence))){break}

  }


  ## refined estimation
  ##-------------------------------------------
  if(!(all(persistence$persistence))){

    if(any(persistence$persistence)){


    for(refined_round  in 1:args_digit){
      refined_killing <- seq(max(persistence[level ==0, ]$killing),
                             min(persistence[level >0, ]$killing),
                             length.out = 10 + ifelse(refined_round == 1, 0, 1))

      for(idx in 1:length(refined_killing)){
        cat("\n[ Refined ", round, ":", idx, "]")

        persistence <- rbind(persistence,
                             single_run_mod( refined_killing[idx],
                                             args_susc,
                                             args_fitness,
                                             args_loss,
                                             args_transfer,
                                             args_basal_growth,
                                             args_basal_surv,
                                             args_capacity,
                                             args_popsize,
                                             args_bin_folder,
                                             args_n_steps,
                                             args_n_writes)
        )
      }
    }
    }
  }

  return(persistence)
}




##----------------------------------------------------------------------------##
## Application
##----------------------------------------------------------------------------##

## Simulations (no transfer)

cost <- c(0.05, 0.10, 0.2, 0.3)
resistance = c(0.2, 0.9, 0.99)
loss = c(1E-3, 1E-5, 1E-7)
transfer = c(0)

msp_grid<- expand.grid(cost,
                       resistance,
                       loss,
                       transfer)
colnames(msp_grid) <- c("cost", "resistance", "loss", "transfer")
msp_grid <- data.table(msp_grid)
msp_grid<- msp_grid[transfer<loss, ]
msp_grid[, msp:= NA]
msp_grid[, msp_lim:= NA]
nrow(msp_grid)

for(i in 1:nrow(msp_grid)){
  cat("\n Grid :", i,  "/", nrow(msp_grid) )
  eval_msp <- estimate_msp(args_susc = 1-msp_grid$resistance[i],
                           args_fitness= 1-msp_grid$cost[i],
                           args_loss = msp_grid$loss[i],
                           args_transfer =  msp_grid$transfer[i],
                           args_basal_growth =0.1,
                           args_basal_surv =0.99,
                           args_capacity =1E9,
                           args_popsize = 1E3,
                           args_max_round = 4,
                           args_digit = 2,
                           args_bin_folder = paste("tmp", i, sep = "_"),
                           args_n_steps = 100000,
                           args_n_writes =1
                           )
  if(any(eval_msp$persistence)){
    msp_grid$msp[i] <- min(eval_msp[persistence== TRUE, ]$killing)}

  if(any(!eval_msp$persistence)){
    msp_grid$msp_lim[i] <- max(eval_msp[persistence== FALSE, ]$killing)

  }

}


saveRDS(msp_grid,
        paste(dirname(getwd()),
              "/examples/example_plasmids/msp_grid.rds",
              sep = ""))



## Plot results
msp_grid <- readRDS( paste(dirname(getwd()),
                           "/examples/example_plasmids/msp_grid.rds",
                           sep = ""))


## Plot

plot(p1 <- ggplot(data =  msp_grid,
       mapping = aes(x = cost,
                     y = msp,
                     group = interaction(resistance, loss),
                     color = as.factor(resistance),
                     linetype = as.factor(loss),
                     shape = as.factor(loss))
  )+
  geom_point()+
  geom_line()+
  theme_bw()+
  scale_y_log10()+
  labs(color = "Resistance",
       shape = "Loss",
       linetype = "Loss")
)


## Rough comparison with predicted

msp_grid[, predicted:= (loss+0.01*cost)/(resistance-cost)]
msp_grid[predicted<0, predicted := NA] ## the equilibrium cannot appends, so
                                       ## that the satbaility criteria does not
                                       ## hold anymore.

plot(
p2 <- ggplot(data =  msp_grid,
       mapping = aes(x = predicted,
                     y = msp,
                     color = as.factor(resistance),
                     shape = as.factor(cost))
)+
  geom_point()+
  theme_bw()+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", )+
  labs(color = "Resistance",
       shape = "Cost")
)



## Rough comparison with predicted

png( paste(dirname(getwd()),
               "/examples/example_plasmids/msevol_MSP.png",
               sep = ""),
     res = 300, units = "cm", width = 25,height = 10 )

ggpubr::ggarrange(p1, p2, nrow = 1, ncol = 2)

dev.off()
