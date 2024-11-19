################################################################################
## MSE1 RESISTRACK MODEL
##  A bacterial hospital ecosystem simulation framework
## (c)2018-2024 Jean-Philippe Rasigade, <jean-philippe.rasigade@univ-lyon1.fr>
##  University of Lyon, Centre International de Recherche en Infectiologie
################################################################################
## Author: JP Rasigade, N Lenuzza
################################################################################


setwd("C:/Users/Lenuzza/Desktop/git_msevol/Rmsevol")


##----------------------------------------------------------------------------##
## Options, packages and generic functions
##----------------------------------------------------------------------------##

## WORKING DIRECTORY
## Make sure the the working directory is Rmsevol, i.e. the root folder that
## directly contains 'bin' (with the msevolcli.exe client)

getwd()


## PACKAGES

library(data.table)
library(tidyr)
library(ggplot2)
library(visNetwork)
library(bit64)
library(Rmsevol)
library(R6)
library(gridExtra)


## WARD class
source("../examples/exemple_ward/Ward_class.R")



##----------------------------------------------------------------------------##
## LOADING DATA
##----------------------------------------------------------------------------##

## RUN 1 WITHOUT ATB
my_ward <-readRDS("../examples/exemple_ward/ward_NoATB.rds")
actual_sample_info <-readRDS("../examples/exemple_ward/ward__NoATB_screening.rds")

use_case_prefix <- "ward__NoATB"
ATB <- FALSE
period <- 7

# 
# ## RUN 2 WITH ATB
# my_ward <-readRDS("../examples/exemple_ward/ward_ATB.rds")
# actual_sample_info <-readRDS("../examples/exemple_ward/ward__ATB_screening.rds")
# n_treatment <- readRDS("../examples/exemple_ward/ward__ATB_treatment.rds")
# 
# use_case_prefix <- "ward__ATB"
# ATB <- TRUE
# period <- 7




##----------------------------------------------------------------------------##
## PATIENT DYNAMICS
##----------------------------------------------------------------------------##

## Total number of (colonized) admission
nrow(my_ward$patient_info)
nrow(my_ward$patient_info[status_at_admission == "Sensitive"])
nrow(my_ward$patient_info[status_at_admission == "Sensitive"])/nrow(my_ward$patient_info)*100


## Time dynamics of the number of patients actually in the ward
n_patients <- my_ward$sample_info[, .(n_patients = .N), by = time]

mean(n_patients$n_patients)
sd(n_patients$n_patients)
range(n_patients$n_patients)



p1 <- ggplot(n_patients, aes(x = time, y = n_patients)) +
  geom_line(color = "black") +
  geom_point(color = "black")+
  theme_bw() +
  labs(x = "Time", y = "Number of patients")+
  geom_hline(yintercept = mean(n_patients$n_patients), 
             color = "red", linetype = "dashed")

p2 <- ggplot(n_patients, aes(x =  n_patients)) +
  geom_histogram(aes(y = after_stat(density)), bins = 15, fill = "grey", color = "black") +
  stat_function(fun = dnorm, args = list(mean = mean(n_patients$n_patients), 
                                         sd = sd(n_patients$n_patients)),
                color = "red", linetype = "dashed") +
  coord_flip() + 
  theme_bw()+
  labs(x = "Density", y = "") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

png(paste("../examples/exemple_ward/", use_case_prefix, "_Figure1_npatients.png", sep = ""),
    res = 300, units = "cm", width = 21, height = 10)
grid.arrange(p1, p2, ncol = 2, widths = c(4, 1))
dev.off()


## if needed, patients_under_treatment
if(ATB){
  n_patients <- merge(n_patients, 
                      n_treatment,
                      by = "time", all = TRUE)
  n_patients[is.na(n_treatment), n_treatment:=0]
  
  mean(n_patients$n_treatment/n_patients$n_patients*100)
  sd(n_patients$n_treatment/n_patients$n_patients*100)
  range(n_patients$n_treatment/n_patients$n_patients*100)
  
  
  
  p1 <- ggplot(n_patients, aes(x = time, y = n_treatment)) +
    geom_line(color = "black") +
    geom_point(color = "black")+
    theme_bw() +
    labs(x = "Time", y = "Number of patients")+
    geom_hline(yintercept = mean(n_patients$n_treatment), 
               color = "red", linetype = "dashed")
  
  p2 <- ggplot(n_patients, aes(x =  n_treatment)) +
    geom_histogram(aes(y = after_stat(density)),  fill = "grey", color = "black") +
    stat_function(fun = dnorm, args = list(mean = mean(n_patients$n_treatment), 
                                           sd = sd(n_patients$n_treatment)),
                  color = "red", linetype = "dashed") +
    coord_flip() + 
    theme_bw()+
    labs(x = "Density", y = "") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  
  png(paste("../examples/exemple_ward/", use_case_prefix, "_Figure1bis_nATBpatients.png", sep = ""),
      res = 300, units = "cm", width = 21, height = 10)
  grid.arrange(p1, p2, ncol = 2, widths = c(4, 1))
  dev.off()
  

  
}



## Patient length-of-stay

mean(my_ward$patient_info$discharge_day - my_ward$patient_info$admission_day,
     na.rm = TRUE)
sd(my_ward$patient_info$discharge_day - my_ward$patient_info$admission_day,
   na.rm = TRUE)
range(my_ward$patient_info$discharge_day - my_ward$patient_info$admission_day,
      na.rm = TRUE)


png(paste("../examples/exemple_ward/", use_case_prefix, "_Figure2_LOS.png", sep = ""),
    res = 300, units = "cm", width = 10.5, height = 10)
ggplot(my_ward$patient_info, 
       aes(x =  discharge_day - admission_day)) +
  geom_histogram(aes(y = after_stat(density)),  binwidth = 1, 
                 fill = "grey", color = "black") +
  stat_function(fun = dexp, args = list(rate = my_ward$params$proba_discharge),
                color = "red", linetype = "dashed") +
  theme_bw()+
  labs(x = "Length of stay (in day)", y = "Density")
dev.off()




## Time dynamics of importation

importation_dynamic <- merge(
  my_ward$patient_info[status_at_admission == "Sensitive", 
                       .N, 
                       by = admission_day],
  data.table(admission_day = unique(my_ward$sample_info$time)),
  by = "admission_day",
  all.x = TRUE, all.y = TRUE)
importation_dynamic[is.na(N), N:=0]



p1 <- ggplot(importation_dynamic, aes(x = admission_day, y =N)) +
  geom_line(color = "black") +
  geom_point(color = "black")+
  theme_bw() +
  labs(x = "Time", y = "Number of importation")


p2 <- ggplot(importation_dynamic, aes(x =  N)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 1, fill = "grey", color = "black") +
  coord_flip() + 
  theme_bw()+
  labs(x = "Density", y = "") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())


png(paste("../examples/exemple_ward/", use_case_prefix, "Figure3_cases_importation.png", sep = ""),
    res = 300, units = "cm", width = 21, height = 10)
grid.arrange(p1, p2, ncol = 2, widths = c(4, 1))

dev.off()





##----------------------------------------------------------------------------##
## DYNAMICS OF COLONIZATIONS (including IMPORTATION CASES)
##----------------------------------------------------------------------------##

## TRUE PROFIL (based on extensive and 100% sensitive sampling)

n_colonized <- merge(
  data.table(time = unique(my_ward$sample_info$time)),
  my_ward$sample_info[true_status == "Sensitive", 
                      .(n_colonized = .N),
                      by = time],
  by = "time",
  all.x = TRUE, 
  all.y = TRUE)
n_colonized[is.na(n_colonized), n_colonized:= 0]


n_colonized <- merge(n_colonized, n_patients, by = "time", all = TRUE)
n_colonized[, proportion_colonized:= n_colonized/n_patients]


p1 <- ggplot(n_colonized, aes(x = time, y = n_colonized)) +
  geom_line(color = "black") +
  geom_point(color = "black")+
  theme_bw() +
  labs(x = "Time", y = "Number of current colonized patients")


p2 <- ggplot(n_colonized, aes(x =  n_colonized)) +
  geom_histogram(aes(y = ..density..), binwidth = 1, fill = "grey", color = "black") +
  coord_flip() + 
  theme_bw()+
  labs(x = "Density", y = "") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

p3 <- ggplot(n_colonized, aes(x = time, y = proportion_colonized*100)) +
  geom_line(color = "black") +
  geom_point(color = "black")+
  theme_bw() +
  labs(x = "Time", y = "Proportion of current colonized patients (%)")


p4 <- ggplot(n_colonized, aes(x =  proportion_colonized*100)) +
  geom_histogram(aes(y = ..density..), binwidth = 5, fill = "grey", color = "black") +
  coord_flip() + 
  theme_bw()+
  labs(x = "Density", y = "") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())



png(paste("../examples/exemple_ward/", use_case_prefix, "_Figure4_Colonized_TRUE.png", sep = ""),
    res = 300, units = "cm", width = 21, height = 20)
grid.arrange(p1, p2,p3, p4, nrow = 2,  ncol = 2, widths = c(4, 1))
dev.off()

mean(n_colonized$n_colonized)
sd(n_colonized$n_colonized)
range(n_colonized$n_colonized)


mean(n_colonized$proportion_colonized)
sd(n_colonized$proportion_colonized)
range(n_colonized$proportion_colonized)




## DETECTED PROFIL (assuming a 90% compliance and a 90% sensitivity)


detected_colonized <- merge(
  data.table(time = unique(actual_sample_info$time)),
  actual_sample_info[test_status == "Sensitive", 
                     .(n_colonized = .N),
                     by = time],
  by = "time",
  all.x = TRUE, 
  all.y = TRUE)
detected_colonized[is.na(n_colonized), n_colonized:= 0]


detected_colonized <- merge(detected_colonized, 
                            n_patients[time%%period ==1], by = "time", all = TRUE)
detected_colonized[, proportion_colonized:= n_colonized/n_patients]


p1 <- ggplot(detected_colonized, aes(x = time, y = n_colonized)) +
  geom_line(color = "black") +
  geom_point(color = "black")+
  theme_bw() +
  labs(x = "Time", y = "Number of current colonized patients")


p2 <- ggplot(detected_colonized, aes(x =  n_colonized)) +
  geom_histogram(aes(y = ..density..), binwidth = 1, fill = "grey", color = "black") +
  coord_flip() + 
  theme_bw()+
  labs(x = "Density", y = "") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

p3 <- ggplot(detected_colonized, aes(x = time, y = proportion_colonized*100)) +
  geom_line(color = "black") +
  geom_point(color = "black")+
  theme_bw() +
  labs(x = "Time", y = "Proportion of current colonized patients (%)")


p4 <- ggplot(detected_colonized, aes(x =  proportion_colonized*100)) +
  geom_histogram(aes(y = ..density..), binwidth = 5, fill = "grey", color = "black") +
  coord_flip() + 
  theme_bw()+
  labs(x = "Density", y = "") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())



png(paste("../examples/exemple_ward/", use_case_prefix, "_Figure5_Colonized_DETECTED.png", sep = ""),
        res = 300, units = "cm", width = 21, height = 20)
grid.arrange(p1, p2,p3, p4, nrow = 2,  ncol = 2, widths = c(4, 1))
dev.off()



mean(detected_colonized$n_colonized)
sd(detected_colonized$n_colonized)
range(detected_colonized$n_colonized)


mean(detected_colonized$proportion_colonized)
sd(detected_colonized$proportion_colonized)
range(detected_colonized$proportion_colonized)




##----------------------------------------------------------------------------##
## Dynamics of hospital-acquired contaminations
##----------------------------------------------------------------------------##

## TRUE PROFIL (based on extensive and 100% sensitive sampling)

first_positive_sample <- 
  my_ward$sample_info[true_status=="Sensitive", 
                      .(first_time = min(time)),
                      by = list(patient)]

first_positive_sample <- merge(first_positive_sample,
                               my_ward$patient_info,
                               by.x = "patient", 
                               by.y = "idx")

nrow(first_positive_sample[status_at_admission == "Free", ]) ## Number of transmission


n_hai <- merge(
  data.table(time = unique(my_ward$sample_info$time)),
  first_positive_sample[status_at_admission == "Free", 
                        .(hai = .N),
                        by = first_time],
  by.x = "time",
  by.y = "first_time",
  all.x = TRUE, 
  all.y = TRUE)

n_hai[is.na(hai), hai:= 0]



p1 <- ggplot(n_hai, aes(x = time, y = hai)) +
  geom_line(color = "black") +
  geom_point(color = "black")+
  theme_bw() +
  labs(x = "Time", y = "Number of transmissions")


p2 <- ggplot(n_hai, aes(x =  hai)) +
  geom_histogram(aes(y = ..density..), binwidth = 1, fill = "grey", color = "black") +
  coord_flip() + 
  theme_bw()+
  labs(x = "Density", y = "") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())


png(paste("../examples/exemple_ward/", use_case_prefix, "_Figure6_HAI_TRUE.png", sep = ""),
    res = 300, units = "cm", width = 21, height = 10)
grid.arrange(p1, p2 , nrow = 1,  ncol = 2, widths = c(4, 1))
dev.off()

mean(n_hai$hai)
sd(n_hai$hai)
range(n_hai$hai)





## DETECTED PROFIL 

first_detected_positive_sample <- 
  actual_sample_info[true_status=="Sensitive", 
                     .(first_time = min(time)), by = list(patient)]

first_detected_positive_sample <- merge(first_detected_positive_sample,
                                        my_ward$patient_info,
                                        by.x = "patient", 
                                        by.y = "idx")

nrow(first_detected_positive_sample[status_at_admission == "Free", ]) ## Number of transmission


n_detected_hai <- merge(
  data.table(time = unique(actual_sample_info$time)),
  first_detected_positive_sample[status_at_admission == "Free", 
                                 .(hai = .N),
                                 by = first_time],
  by.x = "time",
  by.y = "first_time",
  all.x = TRUE, 
  all.y = TRUE)

n_detected_hai[is.na(hai), hai:= 0]
sum(n_detected_hai$hai)


p1 <- ggplot(n_detected_hai, aes(x = time, y = hai)) +
  geom_line(color = "black") +
  geom_point(color = "black")+
  theme_bw() +
  labs(x = "Time", y = "Number of transmissions")

p2 <- ggplot(n_detected_hai, aes(x =  hai)) +
  geom_histogram(aes(y = ..density..), binwidth = 1, fill = "grey", color = "black") +
  coord_flip() + 
  theme_bw()+
  labs(x = "Density", y = "") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())



png(paste("../examples/exemple_ward/", use_case_prefix, "_Figure7_HAI_DETECTED.png", sep = ""),
    res = 300, units = "cm", width = 21, height = 10)
grid.arrange(p1, p2 , nrow = 1,  ncol = 2, widths = c(4, 1))
dev.off()

mean(n_detected_hai$hai)
sd(n_detected_hai$hai)
range(n_detected_hai$hai)

