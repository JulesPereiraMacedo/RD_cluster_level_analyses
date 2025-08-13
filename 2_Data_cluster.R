rm(list = ls())

########### Warning ###########
# Before starting the simulation, you have to run the workspace file first
# Save the Workspace in the "Base_file" file.
# Your "Base_file" should have : Workspace, A file with all data by scenario as "Data_file" file in github
# Then the results will be save in "Resu_file" file automatically.
# "Workspace_name" have to be the same than the one that you give in the "Workspace.R" in "Function" file




Base_file = ""
Workspace_name = "Workspace_clus_test.RData"
Workspace_file = paste(Base_file,Workspace_name,sep = "/")

load(Workspace_file)

#packages ----

library(cli)
library(Matrix)
library(MASS)
library(lme4)
library(backports)
library(dplyr)
library(parallel)
library(iterators)
library(rngtools)
library(foreach)
library(doRNG)
library(doParallel)
library(dplyr)
library(gee)
library(geepack)
library(spind)
library(arm)
library(doBy)
library(here)
library(geesmv)
library(matrixcalc)


## Data

registerDoParallel(cores = Sys.getenv('SLURM_NTASKS')) # Number of cores used for the paralellism


itt = 1000                  # Number of iteration 
Scen = rep(1:Nb_scen)
set.seed(10091958)


Total_time_d = Sys.time()

for (n in Scen) {
  dir.create(paste(Data_file,"/Scenario_",n,sep=""))
  
  Scenario_use = get(paste("Scenario",sep = "_",n))
  
  ## Simu ----
  
  
  debut = Sys.time()
  
  
  res <- foreach(i = 1:itt,
                 # .combine = rbind, #unlist
                 #.errorhandling = "remove", #s'il y a un problème enlève la ligne de
                 .packages = c("stats","arm","gee","geepack","spind","doBy","doRNG","doParallel","dplyr","here")) %dorng% fun_para_save_data_clus(itt=1,itt_para  = i,n,Scenario = Scenario_use,Pi_int,Pi_con,rho_z,OR_int,OR_con,Data_itt_File = Data_file)
  
  
  end = Sys.time()
  time_scen = end-debut
  print(paste("Time to generated data from Scenario",n,sep = " : "))
  print(time_scen)
  
}

Total_time_e = Sys.time()
Total_time = Total_time_e - Total_time_d 
print("Total time for generating Data :")
print(Total_time)
