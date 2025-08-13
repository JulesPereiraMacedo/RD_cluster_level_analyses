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
Resu_TMLE_file_indvcov = paste(Resu_file,"/TMLE_analyse_indvcov/",sep = "")
Resu_TMLE_file_bothcov = paste(Resu_file,"/TMLE_analyse_bothcov/",sep = "")


dir.create(Resu_file)
dir.create(Resu_TMLE_file_indvcov)
dir.create(Resu_TMLE_file_bothcov)



Scen = rep(1:Nb_scen)

registerDoParallel(cores = Sys.getenv('SLURM_NTASKS')) # Number of cores used for the paralellism

for (n in Scen) {
  
  ## Column names for empty files ----
  
  
  name_CL_TMLE = data.frame(RD_CL_TMLE_JAPM = "RD_CL_TMLE_JAPM",
                            SE_CL_TMLE_JAPM = "SE_CL_TMLE_JAPM",
                            LL_95_CL_TMLE_JAPM = "LL_95_CL_TMLE_JAPM",
                            UL_95_CL_TMLE_JAPM = "UL_95_CL_TMLE_JAPM",
                            Bias_JAPM = "Bias_JAPM",
                            Relative_Bias_JAPM = "Relative_Bias_JAPM",
                            Itt = "Itt")
  
  
  ## Empty files .csv ----
  
  write.table(name_CL_TMLE,
              file = paste(Resu_TMLE_file_indvcov,"/Data_output_CL_TMLE_indvcov_Scenario_",n,".csv",sep = ""),
              row.names = FALSE,
              col.names = FALSE,
              sep = ";")
  
  write.table(name_CL_TMLE,
              file = paste(Resu_TMLE_file_bothcov,"/Data_output_CL_TMLE_bothcov_Scenario_",n,".csv",sep = ""),
              row.names = FALSE,
              col.names = FALSE,
              sep = ";")
}



itt = 1000 

Total_time_d = Sys.time()

for (n in Scen) {
  
  debut = Sys.time()
  
  Scenario_use = get(paste("Scenario",n,sep="_"))
  
  res <- foreach(i = 1:itt,
                 # .combine = rbind, #unlist
                 #.errorhandling = "remove", #s'il y a un problème enlève la ligne de
                 .packages = c("stats","arm","gee","geepack","spind","doBy","doRNG","doParallel","dplyr","here","geesmv","matrixcalc")) %dorng% fun_para_analyse_cluster_TMLE(itt_para = i,Scenario = Scenario_use,n=n,Data_file=Data_file,Resu_file=Resu_file)
  
  
  end = Sys.time()
  time_scen = end-debut
  print(paste("Time to generated data from Scenario",n,sep = " : "))
  print(time_scen)
  
}

Total_time_e = Sys.time()
Total_time = Total_time_e - Total_time_d 
print("Total time for generating Data :")
print(Total_time)
