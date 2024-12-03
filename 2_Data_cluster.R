rm(list = ls())

Base_file = ""
Workspace_name = "Workspace_clus_test.RData"
Workspace_file = paste(Base_file,Workspace_name,sep = "/")

load(Workspace_file)

#packages ----



## Data

registerDoParallel(cores = 6) # Number of cores (= 6 here) used for the paralellism


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
