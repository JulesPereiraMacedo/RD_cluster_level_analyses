rm(list = ls())

Base_file = ""
Workspace_name = "Workspace_clus_test.RData"
Workspace_file = paste(Base_file,Workspace_name,sep = "/")

load(Workspace_file)


# File path ---- 
Resu_file = paste(Base_file,"/Res_clust_analyse",sep = "")
Resu_TMLE_file_indvcov = paste(Resu_file,"/TMLE_analyse_indvcov/",sep = "")
Resu_TMLE_file_bothcov = paste(Resu_file,"/TMLE_analyse_bothcov/",sep = "")


dir.create(Resu_file)
dir.create(Resu_TMLE_file_indvcov)
dir.create(Resu_TMLE_file_bothcov)



Scen = rep(1:Nb_scen)

registerDoParallel(cores = 6) # Number of cores used for the paralellism

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
