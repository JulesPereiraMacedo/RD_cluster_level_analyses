rm(list = ls())

Base_file = ""
Workspace_name = "Workspace_clus_test.RData"
Workspace_file = paste(Base_file,Workspace_name,sep = "/")

load(Workspace_file)

Resu_file = paste(Base_file,"/Res_clust_analyse",sep = "")

Resu_TMLE_file_indvcov = paste(Resu_file,"/TMLE_analyse_indvcov/",sep = "")
Resu_TMLE_file_bothcov = paste(Resu_file,"/TMLE_analyse_bothcov/",sep = "")
Resu_GC_Bin_file_indvcov = paste(Resu_file,"/GC_Bin_analyse_indv",sep = "")
Resu_GC_Bin_file_bothcov = paste(Resu_file,"/GC_Bin_analyse_both",sep = "")
Resu_UW_file = paste(Resu_file,"/UW_analyse",sep = "")
Resu_TSP_Bin_file = paste(Resu_file,"/TSP_Bin_analyse/",sep = "")
Resu_TSP_Gauss_file = paste(Resu_file,"/TSP_Gauss_analyse/",sep = "")


Resume_output_file = paste(Resu_file,"/Resume_output_file",sep = "")

dir.create(Resume_output_file)

##########################################################
##########################################################
####                 Resume results                   ####
##########################################################
##########################################################


Scen = rep(1:Nb_scen)

for (n in Scen) {
  
  Scenario_use = get(paste("Scenario",sep = "_",n))
  
  True_RD = as.numeric(Scenario_use[1]) - as.numeric(Scenario_use[2])
  
  Res_UW                   <- read.csv2(paste(Resu_UW_file,"/Data_output_CL_UW_Scenario_",n,".csv",sep=""))
  Res_TSP_Bin              <- read.csv2(paste(Resu_TSP_Bin_file,"/Data_output_CL_TSP_Bin_Scenario_",n,".csv",sep=""))
  Res_TSP_Gauss            <- read.csv2(paste(Resu_TSP_Gauss_file,"/Data_output_CL_TSP_Gauss_Scenario_",n,".csv",sep=""))
  Res_Gcomp_Bin_both       <- read.csv2(paste(Resu_GC_Bin_file_bothcov,"/Data_output_CL_GC_both_Bin_Scenario_",n,".csv",sep=""))
  Res_Gcomp_Bin_indv       <- read.csv2(paste(Resu_GC_Bin_file_indvcov,"/Data_output_CL_GC_indv_Bin_Scenario_",n,".csv",sep=""))
  Res_TMLE_indv            <- read.csv2(paste(Resu_TMLE_file_indvcov,"/Data_output_CL_TMLE_indvcov_Scenario_",n,".csv",sep=""))
  Res_TMLE_both            <- read.csv2(paste(Resu_TMLE_file_bothcov,"/Data_output_CL_TMLE_bothcov_Scenario_",n,".csv",sep=""))
  
  Res_UW               <- na.omit(Res_UW)
  Res_TSP_Bin          <- na.omit(Res_TSP_Bin)
  Res_TSP_Gauss        <- na.omit(Res_TSP_Gauss)
  Res_Gcomp_Bin_both   <- na.omit(Res_Gcomp_Bin_both)
  Res_Gcomp_Bin_indv   <- na.omit(Res_Gcomp_Bin_indv)
  Res_TMLE_indv        <- na.omit(Res_TMLE_indv)
  Res_TMLE_both        <- na.omit(Res_TMLE_both)
  
  ### Resultat
  
  #### Bias ####
  Bias_Mean = c(mean(Res_UW$Bias_iteration_CL_UW),
                mean(Res_TSP_Bin$Bias_iteration_CL_TSP_Bin),
                mean(Res_TSP_Gauss$Bias_iteration_CL_TSP_Gauss),
                mean(Res_Gcomp_Bin_both$Bias_iteration_CL_GC_Bin),
                mean(Res_Gcomp_Bin_indv$Bias_iteration_CL_GC_Bin),
                mean(Res_TMLE_indv$Bias_JAPM),
                mean(Res_TMLE_both$Bias_JAPM))
  
  MC_bias = c(sqrt(1/(nrow(Res_UW)*(nrow(Res_UW)-1)) * sum((Res_UW$Risk_difference_Estimate_CL_UW - mean(Res_UW$Risk_difference_Estimate_CL_UW))**2)),
              sqrt(1/(nrow(Res_TSP_Bin)*(nrow(Res_TSP_Bin)-1)) * sum((Res_TSP_Bin$Risk_difference_Estimate_CL_TSP_Bin - mean(Res_TSP_Bin$Risk_difference_Estimate_CL_TSP_Bin))**2)),
              sqrt(1/(nrow(Res_TSP_Gauss)*(nrow(Res_TSP_Gauss)-1)) * sum((Res_TSP_Gauss$Risk_difference_Estimate_CL_TSP_Gauss - mean(Res_TSP_Gauss$Risk_difference_Estimate_CL_TSP_Gauss))**2)),
              sqrt(1/(nrow(Res_Gcomp_Bin_both)*(nrow(Res_Gcomp_Bin_both)-1)) * sum((Res_Gcomp_Bin_both$Risk_difference_Estimate_CL_GC_Bin - mean(Res_Gcomp_Bin_both$Risk_difference_Estimate_CL_GC_Bin))**2)),
              sqrt(1/(nrow(Res_Gcomp_Bin_indv)*(nrow(Res_Gcomp_Bin_indv)-1)) * sum((Res_Gcomp_Bin_indv$Risk_difference_Estimate_CL_GC_Bin - mean(Res_Gcomp_Bin_indv$Risk_difference_Estimate_CL_GC_Bin))**2)),
              sqrt(1/(nrow(Res_TMLE_indv)*(nrow(Res_TMLE_indv)-1)) * sum((Res_TMLE_indv$RD_CL_TMLE_JAPM - mean(Res_TMLE_indv$RD_CL_TMLE_JAPM))**2)),
              sqrt(1/(nrow(Res_TMLE_both)*(nrow(Res_TMLE_both)-1)) * sum((Res_TMLE_both$RD_CL_TMLE_JAPM - mean(Res_TMLE_both$RD_CL_TMLE_JAPM))**2))
  )
  
  #### Coverage Rate ####
  Coverage_rate = c(length(which( True_RD >= Res_UW$LL_95 & True_RD <= Res_UW$UL_95))/nrow(Res_UW),
                    length(which( True_RD >= Res_TSP_Bin$LL_95 & True_RD <= Res_TSP_Bin$UL_95))/nrow(Res_TSP_Bin),
                    length(which( True_RD >= Res_TSP_Gauss$LL_95 & True_RD <= Res_TSP_Gauss$UL_95))/nrow(Res_TSP_Gauss),
                    length(which( True_RD >= Res_Gcomp_Bin_both$LL_95 & True_RD <= Res_Gcomp_Bin_both$UL_95))/nrow(Res_Gcomp_Bin_both),
                    length(which( True_RD >= Res_Gcomp_Bin_indv$LL_95 & True_RD <= Res_Gcomp_Bin_indv$UL_95))/nrow(Res_Gcomp_Bin_indv),
                    length(which( True_RD >= Res_TMLE_indv$LL_95_CL_TMLE_JAPM & True_RD <= Res_TMLE_indv$UL_95_CL_TMLE_JAPM))/nrow(Res_TMLE_indv),
                    length(which( True_RD >= Res_TMLE_both$LL_95_CL_TMLE_JAPM & True_RD <= Res_TMLE_both$UL_95_CL_TMLE_JAPM))/nrow(Res_TMLE_both))
  
  # MC_CR = sqrt(Coverage_rate *(1-Coverage_rate)/c(nrow(Res_UW),nrow(Res_TSP_Bin),nrow(Res_TSP_Gauss),nrow(Res_Gcomp_Bin_both),nrow(Res_Gcomp_Bin_indv),nrow(Res_TMLE_indv),nrow(Res_TMLE_both),nrow(Res_TMLE_mixed),nrow(Res_TMLE_clus1),nrow(Res_TMLE_clus2)))
  MC_CR = sqrt(Coverage_rate *(1-Coverage_rate)/c(nrow(Res_UW),
                                                  nrow(Res_TSP_Bin),
                                                  nrow(Res_TSP_Gauss),
                                                  nrow(Res_Gcomp_Bin_both),
                                                  nrow(Res_Gcomp_Bin_indv),
                                                  nrow(Res_TMLE_indv),
                                                  nrow(Res_TMLE_both)))
  
  #### Standard Error
  
  Emp_SE = c(sqrt(
    1/(nrow(Res_UW)-1)*
      sum((Res_UW$Risk_difference_Estimate_CL_UW-mean(Res_UW$Risk_difference_Estimate_CL_UW))**2)),
    sqrt(
      1/(nrow(Res_TSP_Bin)-1)*
        sum((Res_TSP_Bin$Risk_difference_Estimate_CL_TSP_Bin-mean(Res_TSP_Bin$Risk_difference_Estimate_CL_TSP_Bin))**2)),
    sqrt(
      1/(nrow(Res_TSP_Gauss)-1)*
        sum((Res_TSP_Gauss$Risk_difference_Estimate_CL_TSP_Gauss-mean(Res_TSP_Gauss$Risk_difference_Estimate_CL_TSP_Gauss))**2)),
    sqrt(
      1/(nrow(Res_Gcomp_Bin_both)-1)*
        sum((Res_Gcomp_Bin_both$Risk_difference_Estimate_CL_GC_Bin-mean(Res_Gcomp_Bin_both$Risk_difference_Estimate_CL_GC_Bin))**2)),
    sqrt(
      1/(nrow(Res_Gcomp_Bin_indv)-1)*
        sum((Res_Gcomp_Bin_indv$Risk_difference_Estimate_CL_GC_Bin-mean(Res_Gcomp_Bin_indv$Risk_difference_Estimate_CL_GC_Bin))**2)),
    sqrt(
      1/(nrow(Res_TMLE_indv)-1)*
        sum((Res_TMLE_indv$RD_CL_TMLE_JAPM-mean(Res_TMLE_indv$RD_CL_TMLE_JAPM))**2)),
    sqrt(
      1/(nrow(Res_TMLE_both)-1)*
        sum((Res_TMLE_both$RD_CL_TMLE_JAPM-mean(Res_TMLE_both$RD_CL_TMLE_JAPM))**2))
  )
  
  ModSE = c(sqrt(mean(Res_UW$SE_itt**2)),
            sqrt(mean(Res_TSP_Bin$SE_itt**2)),
            sqrt(mean(Res_TSP_Gauss$SE_itt**2)),
            sqrt(mean(Res_Gcomp_Bin_both$SE_itt**2)),
            sqrt(mean(Res_Gcomp_Bin_indv$SE_itt**2)),
            sqrt(mean(Res_TMLE_indv$SE_CL_TMLE_JAPM**2)),
            sqrt(mean(Res_TMLE_both$SE_CL_TMLE_JAPM**2)))
  
  RE_SE = (ModSE/Emp_SE - 1)*100
  
  
  var_meth = c(var(Res_UW$SE_itt**2),
               var(Res_TSP_Bin$SE_itt**2),
               var(Res_TSP_Gauss$SE_itt**2),
               var(Res_Gcomp_Bin_both$SE_itt**2),
               var(Res_Gcomp_Bin_indv$SE_itt**2),
               var(Res_TMLE_indv$SE_CL_TMLE_JAPM**2),
               var(Res_TMLE_both$SE_CL_TMLE_JAPM**2))
  
  nrow_meth = c(nrow(Res_UW),
                nrow(Res_TSP_Bin),
                nrow(Res_TSP_Gauss),
                nrow(Res_Gcomp_Bin_both),
                nrow(Res_Gcomp_Bin_indv),
                nrow(Res_TMLE_indv),
                nrow(Res_TMLE_both))
  
  MC_RESE = 100*(ModSE/Emp_SE) * sqrt(((var_meth)/(4*nrow_meth*ModSE))+
                                        (1/(2*(nrow_meth-1))))
  
  
  # Results resume ----
  
  # MC for Monte Carlo SE of Estimate
  # RESE = Relative Error of the Standard Error
  # CR = Coverage Rate
  
  tab_res = data.frame(
    Method = c("CL_UW",
               "CL_TSP_Bin",
               "CL_TSP_Gauss",
               "CL_GC_Bin_both",
               "CL_GC_Bin_indv",
               "CL_TMLE_indiv",
               "CL_TMLE_both"),
    Mean_RD = c(mean(Res_UW$Risk_difference_Estimate_CL_UW),
                mean(Res_TSP_Bin$Risk_difference_Estimate_CL_TSP_Bin),
                mean(Res_TSP_Gauss$Risk_difference_Estimate_CL_TSP_Gauss),
                mean(Res_Gcomp_Bin_both$Risk_difference_Estimate_CL_GC_Bin),
                mean(Res_Gcomp_Bin_indv$Risk_difference_Estimate_CL_GC_Bin),
                mean(Res_TMLE_indv$RD_CL_TMLE_JAPM),
                mean(Res_TMLE_both$RD_CL_TMLE_JAPM)),
    Nb_itt = nrow_meth,
    
    Coverage_rate = Coverage_rate*100,
    MC_CR = MC_CR,
    
    Bias_Mean = Bias_Mean,
    MC_bias = MC_bias,
    Bias_relative_Mean = c(mean(Res_UW$Relative_bias_iteration_CL_UW),
                           mean(Res_TSP_Bin$Relative_bias_iteration_CL_TSP_Bin),
                           mean(Res_TSP_Gauss$Relative_bias_iteration_CL_TSP_Gauss),
                           mean(Res_Gcomp_Bin_both$Relative_bias_iteration_CL_GC_Bin),
                           mean(Res_Gcomp_Bin_indv$Relative_bias_iteration_CL_GC_Bin),
                           mean(Res_TMLE_indv$Relative_Bias_JAPM),
                           mean(Res_TMLE_both$Relative_Bias_JAPM)),
    Emp_SE = Emp_SE,
    ModSE = ModSE,
    RE_SE = RE_SE, # Relative Error of Standard Error
    MC_RESE = MC_RESE,
    
    Scenario = rep(n,9),
    k = rep(as.numeric(unlist(Scenario_use[3])),9),
    m = rep(as.numeric(unlist(Scenario_use[4])),9),
    icc = rep(as.numeric(unlist(Scenario_use[5])),9),
    nb_cov_indiv = rep(as.numeric(unlist(Scenario_use[6])),9),
    nb_cov_clus = rep(as.numeric(unlist(Scenario_use[7])),9),
    Prevalence = rep(paste(as.numeric(unlist(Scenario_use[1])),sep="/",as.numeric(unlist(Scenario_use[2]))),9),
    Convergence_rate = nrow_meth/1000*100 # Care ! 1000 -> Number of iterations of our simulation study for each scenario
                         
  )
  
  ### EXcels data resume  ----
  
  write.csv2(tab_res,here::here(Resume_output_file,paste("/Data_output_Scenario_",sep = "",n,".csv")),row.names = FALSE)
  
}
