################################################################################
## MSE1 RESISTRACK MODEL
##  A bacterial hospital ecosystem simulation framework
## (c)2018-2024 Jean-Philippe Rasigade, <jean-philippe.rasigade@univ-lyon1.fr>
##  University of Lyon, Centre International de Recherche en Infectiologie
################################################################################
## Author: JP Rasigade, N Lenuzza
################################################################################

##----------------------------------------------------------------------------##
## Define the class Ward
##----------------------------------------------------------------------------##


Ward <- R6Class(
  "Ward",
  
  public = list(
    
    # Attributes
    params = list(), # Parameters of the simulation
    prm = NULL,         # msevol single configuration
    patient_info = NULL, # Table of patients (current & past)
    sample_info = NULL,  # Table of samples 
    current_time = 0,  # Current time
    info_archid = NULL, ## Fixed archetype idx 

    # Constructor
    
    initialize = function(params) {
      self$params <- params
      self$current_time <- 0
      self$initialize_patient()
      self$initialize_prm()
      #self$initialize_sample()
      self$sample_info <- format_sample_info()
      self$screen_patients()
    },
    
    initialize_patient = function(){
      
      stopifnot(self$params$n_patients >= self$params$n_colonized)

      patient_info <- format_patient_info()
      
      if(self$params$n_patients>0){
        
        patient_info <- rbind(patient_info,
                              data.table(
                                idx = 0:(self$params$n_patients - 1),
                                admission_day = 0,
                                status_at_admission = "Free",
                                under_treatment = FALSE
                              ),
                              fill = TRUE)
        
        ## 
        if(self$params$n_colonized > 0){
          
          idx_to_change <- sample(patient_info$idx,
                                  self$params$n_colonized, 
                                  replace = FALSE)
          patient_info[idx %in% idx_to_change, status_at_admission:=  "Sensitive"]
        }
        
      }
      
      
      self$patient_info <- patient_info
      
      
    },
    
    initialize_prm = function(){
      
      ## Fixed archetype IDs (mandatory for refill)
      self$info_archid <- list(
        patchArchetype_nopressure = 0,
        patchArchetype_pressure = 1,
        cellArchetype_sensitive = 0,
        cell_sensitive = 0,
        chromosomeArchetype_sensitive = 0,
        chromosome_sensitive = 0)
      
      
      ## Create the "bacterial" graph model
      prm <- emptymodel() %>%
        chromosomeArchetype(self$info_archid$chromosomeArchetype_sensitive,
                            fitness = self$params$basal_growth,
                            survival = self$params$basal_survival) %>%
        chromosome(self$info_archid$chromosome_sensitive, 
                   self$info_archid$chromosomeArchetype_sensitive) %>%  
        cellArchetype(self$info_archid$cellArchetype_sensitive) %>% 
        cell(self$info_archid$cell_sensitive,
             self$info_archid$cellArchetype_sensitive) %>% 
        assign("Chromosome", self$info_archid$chromosome_sensitive, 
               "Cell", self$info_archid$cell_sensitive)
      
      
      ## Add the initial patch archetypes
      prm <- prm%>%
        patchArchetype(self$info_archid$patchArchetype_nopressure, 
                       capacity = self$params$carrying_capacity,
                       pressure = c(0, 0)) %>%
        patchArchetype(self$info_archid$patchArchetype_pressure, 
                       capacity = self$params$carrying_capacity,
                       pressure = c(self$params$antibiotic_pressure, 0))
      
      
      ## Add patients from "patients_info"
      current_patients <- self$patient_info[is.na(discharge_day), ]
      if(nrow(current_patients)){
        for(patients_i in 1:nrow(current_patients)){
          
          ## create patient
          prm <- prm %>%  
            patch(current_patients[patients_i,]$idx, 
                  ifelse(current_patients[patients_i,]$under_treatment,
                         self$info_archid$patchArchetype_pressure,
                         self$info_archid$patchArchetype_nopressure))
          
          ## connect patient to bacteria
          if(current_patients[patients_i,]$status_at_admission == "Sensitive"){
            prm <- prm %>%  
              assign("Cell",self$info_archid$cell_sensitive,
                     "Patch", current_patients[patients_i,]$idx,
                     mult = self$params$carrying_capacity) ## TO DO : moduler par growth & surv.
          }
          
          ## connect patient to existing patients
          if(patients_i >1){
            for(i in 1:(patients_i-1)){
              prm <- prm %>%  
                diffusion(current_patients[patients_i,]$idx,
                          current_patients[i,]$idx,
                          diffusion = self$params$diffusion_patch_to_patch,
                          symmetric = TRUE)
            }
          }
        }
      }
     
      
      
      
      ## update "Ward"
      self$prm <- prm
    },
    
    
    # Events
    
    update_time = function(){ ## used to test events independently of "msevol"
      self$current_time <- self$current_time + 1
    },
    
    run_step = function(){
    
        # Run msevol with clean up
        mse_res <- Rmsevol::mse_run(prm = self$prm,
                                    path = self$params$run_path,
                                    n_steps = self$params$run_n_steps,
                                    n_writes = self$params$run_n_writes,
                                    s1 = sample(1E9, 1),
                                    s2 = sample(1E9, 1))

        Rmsevol::mse1.cleanup_csv(self$params$run_path)


        ## Update prm
        self$prm <- Rmsevol::mse1.extract_prm_from_mseres(mse_res,
                                                          self$params$run_n_steps*
                                                            self$params$run_n_writes)
        ## refill prm
        self$prm <- refill_prm_function(self$prm, self$params, self$info_archid)
        
        ## Update time
        self$current_time <- self$current_time + 1
    },
    
    
    # admit_patients
    admit_patients = function(){
      
      ## Draw admission
      n_admission <- sample(0:self$params$max_admission, 1)
      
      if(n_admission){
        
        ## Draw their colonized status
        n_S <- rbinom(1, n_admission, self$params$prop_sensitive_at_admission)
        n_S_i <- c()
        if(n_S){
          n_S_i <- sample(1:n_admission, n_S)
        }
        
        
        ## actual admission
        for(admitted_patient_i in 1:n_admission){
          
          current_idx <- self$prm$Vertex_Patch$id
          last_idx <- ifelse(nrow(self$patient_info), 
                             max(self$patient_info$idx),0)
          
          ## adding patient
          self$prm <- self$prm %>%
            patch(last_idx+1, self$info_archid$patchArchetype_nopressure) 
          
          if(admitted_patient_i %in% n_S_i){
            self$prm <- self$prm %>%
              assign("Cell", self$info_archid$cell_sensitive, 
                     "Patch",last_idx+1, 
                     mult = self$prm$Vertex_PatchArchetype[
                       id ==self$info_archid$patchArchetype_nopressure, ]$capacity)
            
          }
         
          
          ## connect patients to other patients
          for(other_patient_idx in current_idx){
            self$prm <- self$prm %>%
              diffusion(last_idx+1, other_patient_idx, 
                        diffusion = self$params$diffusion_patch_to_patch, 
                        symmetric = TRUE)
          }
          
          
          ## update patient_info
          self$patient_info <- rbind(self$patient_info,
                                     data.table(idx = last_idx+1,
                                                admission_day = self$current_time, 
                                                status_at_admission = ifelse(admitted_patient_i %in% n_S_i, 
                                                                             "Sensitive", "Free"),
                                                under_treatment = FALSE),
                                     fill = TRUE
          )
          
        } ## end of admission loop
    
      } ## end if n_admission>0
      
   
    },
    
    # discharge_patient
    discharge_patient = function(){
      
      patients_to_discharge_idx <- self$patient_info[is.na(discharge_day), ]$idx
      n_possible <- length(patients_to_discharge_idx)
      n_actual <- rbinom(1,n_possible, self$params$proba_discharge )
      
      if(n_actual){
        
        patients_to_discharge_idx <- sample(patients_to_discharge_idx, n_actual)
        #cat("\n Patients to discharge = ", patients_to_discharge_idx, "\n")
        
        self$patient_info[idx %in% patients_to_discharge_idx, discharge_day := self$current_time]
        
        
        ## If needed, store the "infectious" status at discharge. 
        for(dis_i in patients_to_discharge_idx){
          
          self$patient_info[idx == dis_i, 
                            level_S_at_discharge:= ifelse(
                              dis_i %in% self$prm$Edge_Cell_Patch_Multiplicity$target,
                              self$prm$Edge_Cell_Patch_Multiplicity[
                                target == dis_i & source == self$info_archid$cell_sensitive,]$multiplicity,
                              0)
                            ]
         
          
         }
        
        ## actually delete agents
        self$prm <- delete_agents(self$prm, idxs = patients_to_discharge_idx, type = "Patch")
        
        
      }
      
    },
    
    # start_treatment
    start_treatment = function(){
      
      patients_to_start_idx <- self$patient_info[is.na(discharge_day) & 
                                                   !(under_treatment), ]$idx
      n_possible <- length(patients_to_start_idx)
      n_actual <- rbinom(1,  n_possible, self$params$proba_initiate_antibiotherapy)
      patients_to_start_idx <- patients_to_start_idx[sample(1:n_possible,n_actual)]
      
      
      ## update prm
      if(length(patients_to_start_idx)){
        self$prm$Edge_PatchArchetype_Patch_Multiplicity[
          target %in% patients_to_start_idx, 
          source:= self$info_archid$patchArchetype_pressure]
      }
      
      ## update patient_info
      self$patient_info[idx %in% patients_to_start_idx, under_treatment:=TRUE]
      self$patient_info[idx %in% patients_to_start_idx, 
                        end_current_treatment:= as.integer(self$current_time+ self$params$duration_antibiotherapy)]
      #self$patient_info[idx %in% patients_to_start_idx,
      #                  treatment_days :=c(treatment_days, self$current_time)]
      
    },
    
    # stop_treatment
    stop_treatment = function(){
      
      patients_to_stop_idx <- 
        self$patient_info[is.na(discharge_day) & end_current_treatment <= self$current_time, ]$idx
      
      if(length(patients_to_stop_idx)){
        
        self$prm$Edge_PatchArchetype_Patch_Multiplicity[target %in% patients_to_stop_idx,
                                                        source:= self$info_archid$patchArchetype_nopressure ]
        
        self$patient_info[idx %in% patients_to_stop_idx, under_treatment:= FALSE]
        self$patient_info[idx %in% patients_to_stop_idx, end_current_treatment:= NA]
        
      }
    },
    
    
    # screen_patient
    screen_patients = function(){
      
      current_sample_info <- format_sample_info()
      
      if(nrow(self$patient_info[is.na(discharge_day)])){
        
        current_sample_info <- data.table(
          time = self$current_time,
          patient = self$patient_info[is.na(discharge_day)]$idx,
          true_status = "Free")
        
        ## Note that here, we do not have to filter on the
        ##  bacteria ids since only one single strain is
        ## included in the model, with no evolution 
        current_sample_info[
          patient%in% self$prm$Edge_Cell_Patch_Multiplicity[
            multiplicity> self$params$screening_threshold, ]$target,
          true_status := "Sensitive"]
        
      }
      
      self$sample_info <- rbind(self$sample_info,
                                current_sample_info)
      
    }
    
    ) # end "Public"
) # end "R6Class"
      
    


################################################################################

format_patient_info <- function(){
  return(
    patient_info <- data.table(
      idx = integer(),
      admission_day = integer(),
      discharge_day = integer(),
      status_at_admission = character(),
      under_treatment = logical(),
      end_current_treatment = integer(),
      #treatment_days = list(),
      level_S_at_discharge = integer()
    )
  )
  
}


format_sample_info <- function(){
  return(
    patient_info <- data.table(
      time = integer(),
      patient = integer(),
      true_status = character()
    )
  )
}


refill_prm_function <- function(prm, params, info_archid){
  
  
  ## Add missing patch archetypes if needed
  
  if(!(info_archid$patchArchetype_nopressure %in% prm$Vertex_PatchArchetype$id)){
    prm <- prm %>%
       patchArchetype(info_archid$patchArchetype_nopressure,
                      capacity = params$carrying_capacity,
                      pressure = c(0,0))
  }
  
  if(!(info_archid$patchArchetype_pressure %in% prm$Vertex_PatchArchetype$id)){
    prm <- prm %>%
      patchArchetype(info_archid$patchArchetype_pressure,
                     capacity = params$carrying_capacity,
                     pressure = c(params$antibiotic_pressure,0))
  }

  
  ## Check chromosome archetype
  if(!(info_archid$chromosomeArchetype_sensitive %in% prm$Vertex_ChromosomeArchetype$id)){
    prm <- prm %>%
      chromosomeArchetype(id = info_archid$chromosomeArchetype_sensitive,
                          fitness = params$basal_growth, 
                          survival = params$basal_survival)
  }
  
  ## Check chromosome
  if(!(info_archid$chromosome_sensitive %in% prm$Vertex_Chromosome$id)){
    prm <- prm %>%
      chromosome(id = info_archid$chromosome_sensitive,
                 type = info_archid$chromosomeArchetype_sensitive)
    
  }
  
  ## Check cell archetype
  if(!(info_archid$cellArchetype_sensitive %in% prm$Vertex_CellArchetype$id)){
    prm <- prm %>%
      cellArchetype(id = info_archid$cellArchetype_sensitive)
  }
  
  ## Check cell
  if(!(info_archid$cell_sensitive %in% prm$Vertex_Cell$id)){
    prm <- prm %>%
      cell(id = info_archid$cell_sensitive,
                 type = info_archid$cellArchetype_sensitive)%>%
      assign("Chromosome", info_archid$chromosome_sensitive,
             "Cell", info_archid$cell_sensitive)
    
  }
  
  
  
  ##
  return(prm)
  
}


delete_agents <- function(prm, idxs, type = "Patch"){
  
  stopifnot(type %in% c("Patch", "Cell", "Chromosome", "Plasmid", "Gene"))
  
  ## delete agents in Vertex_'type' data.table
  prm[[paste("Vertex_", type, sep = "")]] <- 
    prm[[paste("Vertex_", type, sep = "")]][!(id %in% idxs), ]
  
  
  
  ## delete all inner and outer edges for this agents
  names_tables <- names(prm)[grep("Edge", names(prm))]
  
  for(name_table in names_tables){
    
    composantes <- strsplit(name_table, "_")[[1]]
    
    if(composantes[2] == type){
      prm[[name_table]] <- prm[[name_table]][!(source %in% idxs), ]
    }
    
    if(composantes[3] == type){
      prm[[name_table]] <- prm[[name_table]][!(target %in% idxs), ]
    }
    
  }
  
  return(prm)
  
}


################################################################################


actual_screening <- function(true_sample_data, 
                             period = 7, 
                             compliance = 0.9, 
                             sensitivity = 0.9){
  
  ## Periodic sampling only (e.g. weekly)
  actual_sample_info <-true_sample_data[time%%period == 1,] 
  
  ## Compliance
  actual_sample_info <- actual_sample_info[
    sample(1:nrow(actual_sample_info),
           rbinom(1, nrow(actual_sample_info), compliance)),]
  
  ## sensitivity
  actual_sample_info[, test_status:= true_status]
  positive_screening_idx <- which(actual_sample_info$true_status=="Sensitive")
  false_positive_idx <- sample(
    positive_screening_idx,
    rbinom(1, length(positive_screening_idx), 1-sensitivity)
  )
  actual_sample_info[false_positive_idx, test_status:="Free"]
  
  ## output
  return(actual_sample_info)

}
  







