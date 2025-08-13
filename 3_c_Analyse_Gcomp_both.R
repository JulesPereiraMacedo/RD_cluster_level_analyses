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

# File path ---- 
Resu_file = paste(Base_file,"/Res_clust_analyse",sep = "")
dir.create(Resu_file)

Resu_GC_Bin_file = paste(Resu_file,"/GC_Bin_analyse_both",sep = "")
dir.create(Resu_GC_Bin_file)

Correction_file = paste(Base_file,"/correction_gcomp",sep = "")
dir.create(Correction_file)

Scen = rep(1:Nb_scen)

## Fichier R � creer vide ----

for (n in Scen) {
  
  ## Column names for empty files ----
  
  name_CL_GC_Bin = data.frame(Risk_difference_Estimate_CL_GC_Bin = "Risk_difference_Estimate_CL_GC_Bin",
                              LL_95 = "LL_95",
                              UL_95 = "UL_95",
                              SE_itt = "SE_itt",
                              Bias_iteration_CL_GC_Bin = "Bias_iteration_CL_GC_Bin",
                              Relative_bias_iteration_CL_GC_Bin = "Relative_bias_iteration_CL_GC_Bin",
                              itt_para = "itt_para",
                              error = "error")
  # type_warning = "type_warning")
  
  
  
  ## Empty files .csv ----
  
  write.table(name_CL_GC_Bin,
              file = paste(Resu_GC_Bin_file,"/Data_output_CL_GC_both_Bin_Scenario_",n,".csv",sep = ""),
              row.names = FALSE,
              col.names = FALSE,
              sep = ";")
}




## Data

registerDoParallel(cores = Sys.getenv('SLURM_NTASKS')) # Number of cores used for the paralellism

Time = 50
itt = 1000    # Number of iteration 
# set.seed(1958)

Total_time_d = Sys.time()

for (n in Scen) {
  
  debut = Sys.time()
  
  Scenario_use = get(paste("Scenario",n,sep="_")) 
  
  res <- foreach(i = 1:itt,
                 # .combine = rbind, #unlist
                 #.errorhandling = "remove", #s'il y a un problème enlève la ligne de
                 .packages = c("stats","arm","gee","geepack","spind","doBy","doRNG","doParallel","dplyr","here","geesmv","matrixcalc")) %dorng% fun_cor_gcomp(i = i,Time,n,Base_file,Workspace_name,Data_file,Resu_file,cov_schem = "both")
  
  
  end = Sys.time()
  time_scen = end-debut
  print(paste("Time to generated data from Scenario",n,sep = " : "))
  print(time_scen)
  
}

Total_time_e = Sys.time()
Total_time = Total_time_e - Total_time_d 
print("Total time for generating Data :")
print(Total_time)
