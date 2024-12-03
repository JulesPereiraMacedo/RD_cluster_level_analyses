rm(list = ls())

Base_file = ""
Workspace_name = "Workspace_clus_test.RData"
Workspace_file = paste(Base_file,Workspace_name,sep = "/")

load(Workspace_file)


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
                              # Number_of_itteration_CL_GC_Bin = "Number_of_itteration_CL_GC_Bin",
                              LL_95 = "LL_95",
                              UL_95 = "UL_95",
                              SE_itt = "SE_itt",
                              # Coverage_rate_CL_GC_Bin = "Coverage_rate_CL_GC_Bin",
                              # Abs_Bias_iteration_CL_GC_Bin = "Abs_Bias_iteration_CL_GC_Bin",
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

registerDoParallel(cores = 6) # Number of cores used for the paralellism

Time = 50
itt = 1000    # Number of iteration 

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
