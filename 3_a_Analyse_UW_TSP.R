rm(list = ls())

Base_file = ""
Workspace_name = "Workspace_clus_test.RData"
Workspace_file = paste(Base_file,Workspace_name,sep = "/")

load(Workspace_file)


# File path ---- 
Resu_file = paste(Base_file,"/Res_clust_analyse",sep = "")
dir.create(Resu_file)

Resu_UW_file = paste(Resu_file,"/UW_analyse",sep = "")
Resu_TSP_Bin_file = paste(Resu_file,"/TSP_Bin_analyse/",sep = "")
Resu_TSP_Gauss_file = paste(Resu_file,"/TSP_Gauss_analyse/",sep = "")

dir.create(Resu_UW_file)
dir.create(Resu_TSP_Bin_file)
dir.create(Resu_TSP_Gauss_file)


Scen = rep(1:Nb_scen)

## empty Excel file creation ----

for (n in Scen) {
  
  ## Column names for empty files ----
  
  name_CL_UW = data.frame(Risk_difference_Estimate_CL_UW = "Risk_difference_Estimate_CL_UW",
                          # Number_of_itteration_CL_UW = "Number_of_itteration_CL_UW",
                          LL_95 = "LL_95",
                          UL_95 = "UL_95",
                          SE_itt = "SE_itt",
                          Bias_iteration_CL_UW = "Bias_iteration_CL_UW",
                          Relative_bias_iteration_CL_UW = "Relative_bias_iteration_CL_UW",
                          itt_para = "itt_para")

  name_CL_TSP_Bin = data.frame(Risk_difference_Estimate_CL_TSP_Bin = "Risk_difference_Estimate_CL_TSP_Bin",
                               # Number_of_itteration_CL_TSP_Bin = "Number_of_itteration_CL_TSP_Bin",
                               LL_95 = "LL_95",
                               UL_95 = "UL_95",
                               SE_itt = "SE_itt",
                               Bias_iteration_CL_TSP_Bin = "Bias_iteration_CL_TSP_Bin",
                               Relative_bias_iteration_CL_TSP_Bin = "Relative_bias_iteration_CL_TSP_Bin",
                               itt_para = "itt_para")
  
  name_CL_TSP_Gauss = data.frame(Risk_difference_Estimate_CL_TSP_Gauss = "Risk_difference_Estimate_CL_TSP_Gauss",
                                 # Number_of_itteration_CL_TSP_Gauss = "Number_of_itteration_CL_TSP_Gauss",
                                 LL_95 = "LL_95",
                                 UL_95 = "UL_95",
                                 SE_itt = "SE_itt",
                                 Bias_iteration_CL_TSP_Gauss = "Bias_iteration_CL_TSP_Gauss",
                                 Relative_bias_iteration_CL_TSP_Gauss = "Relative_bias_iteration_CL_TSP_Gauss",
                                 itt_para = "itt_para")
  
  
  
  ## Empty files .csv ----
  
  write.table(name_CL_UW,
              file = paste(Resu_UW_file,"/Data_output_CL_UW_Scenario_",n,".csv",sep = ""),
              row.names = FALSE,
              col.names = FALSE,
              sep = ";")
  
  write.table(name_CL_TSP_Bin,
              file = paste(Resu_TSP_Bin_file,"/Data_output_CL_TSP_Bin_Scenario_",n,".csv",sep = ""),
              row.names = FALSE,
              col.names = FALSE,
              sep = ";")
  
  write.table(name_CL_TSP_Gauss,
              file = paste(Resu_TSP_Gauss_file,"/Data_output_CL_TSP_Gauss_Scenario_",n,".csv",sep = ""),
              row.names = FALSE,
              col.names = FALSE,
              sep = ";")
}




## Data

registerDoParallel(cores = 6) # Number of cores (= 6 here) used for the paralellism


itt = 1000    # Number of iteration 


Total_time_d = Sys.time()

for (n in Scen) {
  
  debut = Sys.time()
  
  Scenario_use = get(paste("Scenario",n,sep="_")) 
  
  res <- foreach(i = 1:itt,
                 # .combine = rbind, #unlist
                 #.errorhandling = "remove", #s'il y a un problème enlève la ligne de
                 .packages = c("stats","arm","gee","geepack","spind","doBy","doRNG","doParallel","dplyr","here","geesmv","matrixcalc")) %dorng% fun_para_analyse_cluster(itt_para=i,Scenario=Scenario_use,n=n,Data_file=Data_file,Resu_file = Resu_file)
  
  
  end = Sys.time()
  time_scen = end-debut
  print(paste("Time to generated data from Scenario",n,sep = " : "))
  print(time_scen)
  
}

Total_time_e = Sys.time()
Total_time = Total_time_e - Total_time_d 
print("Total time for generating Data :")
print(Total_time)
