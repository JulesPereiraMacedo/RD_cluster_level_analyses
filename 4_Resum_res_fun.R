## Fonction pour resumer les resultats

Resum_res_fun_it_spec <- function(n,Resu_file,methods,nb_itt_scen,itt_spec){
  
  if(FALSE %in% unique(methods %in% c("Unajusted",
                                      "TSP_gauss","TSP_bin",
                                      "GC_indv","GC_both",
                                      "TMLE_indv","TMLE_both",
                                      "GC_indv_indep","GC_both_indep","GC_indv_GLMM","GC_both_GLMM",
                                      "TMLE_APS","TMLE_APS_Balzer","TMLE_APS_Balzer_2","TMLE_APS_indv_Balzer"))){
    stop("ERROR IN metods USED : one method is unknown from the vector of methods, verify methods")
  }
  
  # Scenario utilisé
  
  Scenario_use = get(paste("Scenario",sep = "_",n))
  
  True_RD = as.numeric(Scenario_use[1]) - as.numeric(Scenario_use[2])
  
  n_metds = length(methods)
  
  # Différents vecteurs pour les mesures de performances qui nous intérèssent 
  
  Mean_RD = c()
  
  Bias_relative_Mean = c()
  MC_RelativeBias = c()
  
  Bias_Mean = c()
  MC_bias = c()
  
  Coverage_rate = c()
  MC_CR = c()
  CI_width = c()
  P_025_CIW = c()
  P_975_CIW = c()
  MC_CIW = c()
  
  Emp_SE = c()
  ModSE = c()
  RE_SE = c()
  
  var_meth = c()
  nrow_meth = c()
  MC_RESE = c()
  
  Methds_vec = c()
  
  
  if("Unajusted" %in% methods){
    
    Resu_UN_file = paste(Resu_file,"/UW_analyse",sep = "")
    
    Res                     <- read.csv2(paste(Resu_UN_file,"/Data_output_CL_UW_Scenario_",n,".csv",sep=""))
    Res                     <- na.omit(Res)
    Res                     <- Res[which(Res$itt_para %in% itt_spec),]
    
    Mean_RD                 = c(Mean_RD,
                                mean(Res$Risk_difference_Estimate_CL_UW))
    
    Bias_relative_Mean      = c(Bias_relative_Mean,
                                mean(Res$Relative_bias_iteration_CL_UW))
    MC_RelativeBias         = c(MC_RelativeBias,
                                sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum(((Res$Risk_difference_Estimate_CL_UW - mean(Res$Risk_difference_Estimate_CL_UW))/mean(Res$Risk_difference_Estimate_CL_UW) - mean(Res$Relative_bias_iteration_CL_UW)/100)**2)))
    
    Bias_Mean               = c(Bias_Mean,
                                mean(Res$Bias_iteration_CL_UW))
    MC_bias                 = c(MC_bias,
                                sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum((Res$Risk_difference_Estimate_CL_UW - mean(Res$Risk_difference_Estimate_CL_UW))**2)))
    
    Coverage_rate_UN        = length(which( True_RD >= Res$LL_95 & True_RD <= Res$UL_95))/nrow(Res)
    Coverage_rate           = c(Coverage_rate,
                                Coverage_rate_UN)
    
    CI_width_UN             = mean(Res$UL_95 - Res$LL_95)  
    CI_width                = c(CI_width,CI_width_UN)
    P_025_CIW               = c(P_025_CIW, as.numeric(quantile(Res$UL_95 - Res$LL_95,probs = 0.025)))
    P_975_CIW               = c(P_975_CIW, as.numeric(quantile(Res$UL_95 - Res$LL_95,probs = 0.975)))
    MC_CIW                  = c(MC_CIW,
                                sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum(((Res$UL_95 - Res$LL_95) - mean(Res$UL_95 - Res$LL_95))**2)))
                                
    
    
    MC_CR                   = c(MC_CR,
                                sqrt(Coverage_rate_UN *(1-Coverage_rate_UN)/nrow(Res)))
    
    Emp_SE_UN               = sqrt(1/(nrow(Res)-1)*sum((Res$Risk_difference_Estimate_CL_UW-mean(Res$Risk_difference_Estimate_CL_UW))**2))
    Emp_SE                  = c(Emp_SE,
                                Emp_SE_UN)
    
    ModSE_UN                = sqrt(mean(Res$SE_itt**2))
    ModSE                   = c(ModSE,
                                ModSE_UN)
    
    RE_SE_UN                = (ModSE_UN/Emp_SE_UN - 1)*100
    RE_SE                   = c(RE_SE,
                                RE_SE_UN)
    
    var_meth_UN             = var(Res$SE_itt**2)
    var_meth                = c(var_meth,
                                var_meth_UN)
    
    nrow_meth               = c(nrow_meth,
                                nrow(Res))
    
    MC_RESE                 = c(MC_RESE,
                                100*(ModSE_UN/Emp_SE_UN) * sqrt(((var_meth_UN)/(4*nrow(Res)*ModSE_UN))+
                                                                  (1/(2*(nrow(Res)-1)))))
    Methds_vec = c(Methds_vec,"CL_UN")
  }
  
  if("TSP_gauss" %in% methods){
    
    Resu_TSP_Gauss_file = paste(Resu_file,"/TSP_Gauss_analyse/",sep = "")
    
    Res                     <- read.csv2(paste(Resu_TSP_Gauss_file,"/Data_output_CL_TSP_Gauss_Scenario_",n,".csv",sep=""))
    Res                     <- na.omit(Res)
    Res                     <- Res[which(Res$itt_para %in% itt_spec),]
    
    Mean_RD                 = c(Mean_RD,
                                mean(Res$Risk_difference_Estimate_CL_TSP_Gauss))
    
    Bias_relative_Mean      = c(Bias_relative_Mean,
                                mean(Res$Relative_bias_iteration_CL_TSP_Gauss))
    MC_RelativeBias         = c(MC_RelativeBias,
                                sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum(((Res$Risk_difference_Estimate_CL_TSP_Gauss - mean(Res$Risk_difference_Estimate_CL_TSP_Gauss))/mean(Res$Risk_difference_Estimate_CL_TSP_Gauss) - mean(Res$Relative_bias_iteration_CL_TSP_Gauss)/100)**2)))
    
    
    Bias_Mean               = c(Bias_Mean,
                                mean(Res$Bias_iteration_CL_TSP_Gauss))
    MC_bias                 = c(MC_bias,
                                sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum((Res$Risk_difference_Estimate_CL_TSP_Gauss - mean(Res$Risk_difference_Estimate_CL_TSP_Gauss))**2)))
    
    Coverage_rate_TSP_G     = length(which( True_RD >= Res$LL_95 & True_RD <= Res$UL_95))/nrow(Res)
    Coverage_rate           = c(Coverage_rate,
                                Coverage_rate_TSP_G)
    
    CI_width_TSP_G          = mean(Res$UL_95 - Res$LL_95)  
    CI_width                = c(CI_width,CI_width_TSP_G)
    P_025_CIW               = c(P_025_CIW, as.numeric(quantile(Res$UL_95 - Res$LL_95,probs = 0.025)))
    P_975_CIW               = c(P_975_CIW, as.numeric(quantile(Res$UL_95 - Res$LL_95,probs = 0.975)))
    MC_CIW                  = c(MC_CIW,
                                sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum(((Res$UL_95 - Res$LL_95) - mean(Res$UL_95 - Res$LL_95))**2)))
    
    MC_CR                   = c(MC_CR,
                                sqrt(Coverage_rate_TSP_G *(1-Coverage_rate_TSP_G)/nrow(Res)))
    
    Emp_SE_TSP_G            = sqrt(1/(nrow(Res)-1)*sum((Res$Risk_difference_Estimate_CL_TSP_Gauss-mean(Res$Risk_difference_Estimate_CL_TSP_Gauss))**2))
    Emp_SE                  = c(Emp_SE,
                                Emp_SE_TSP_G)
    
    ModSE_TSP_G             = sqrt(mean(Res$SE_itt**2))
    ModSE                   = c(ModSE,
                                ModSE_TSP_G)
    
    RE_SE_TSP_G             = (ModSE_TSP_G/Emp_SE_TSP_G - 1)*100
    RE_SE                   = c(RE_SE,
                                RE_SE_TSP_G)
    
    var_meth_TSP_G          = var(Res$SE_itt**2)
    var_meth                = c(var_meth,
                                var_meth_TSP_G)
    
    nrow_meth               = c(nrow_meth,
                                nrow(Res))
    
    MC_RESE                 = c(MC_RESE,
                                100*(ModSE_TSP_G/Emp_SE_TSP_G) * sqrt(((var_meth_TSP_G)/(4*nrow(Res)*ModSE_TSP_G))+
                                                                        (1/(2*(nrow(Res)-1)))))
    # MC_RESE                 = c(MC_RESE,
    #                             100*(ModSE_TSP_G/Emp_SE_TSP_G) * sqrt( ((1/nrow(Res) * (sum(((Res$SE_itt**2) - (1/(nrow(Res)-1) * sum((Res$SE_itt**2))))**2)) )/(4*nrow(Res)*ModSE_TSP_G)) +
    #                                                                     (1/(2*(nrow(Res)-1)))))
    
    Methds_vec = c(Methds_vec,"CL_TSP_Gauss")
  }
  
  if("TSP_bin" %in% methods){
    
    Resu_TSP_Bin_file       = paste(Resu_file,"/TSP_Bin_analyse/",sep = "")
    
    Res                     <- read.csv2(paste(Resu_TSP_Bin_file,"/Data_output_CL_TSP_Bin_Scenario_",n,".csv",sep=""))
    Res                     <- na.omit(Res)
    Res                     <- Res[which(Res$itt_para %in% itt_spec),]
    
    Mean_RD                 = c(Mean_RD,
                                mean(Res$Risk_difference_Estimate_CL_TSP_Bin))
    
    Bias_relative_Mean      = c(Bias_relative_Mean,
                                mean(Res$Relative_bias_iteration_CL_TSP_Bin))
    MC_RelativeBias         = c(MC_RelativeBias,
                                sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum(((Res$Risk_difference_Estimate_CL_TSP_Bin - mean(Res$Risk_difference_Estimate_CL_TSP_Bin))/mean(Res$Risk_difference_Estimate_CL_TSP_Bin) - mean(Res$Relative_bias_iteration_CL_TSP_Bin)/100)**2)))
    
    
    Bias_Mean               = c(Bias_Mean,
                                mean(Res$Bias_iteration_CL_TSP_Bin))
    MC_bias                 = c(MC_bias,
                                sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum((Res$Risk_difference_Estimate_CL_TSP_Bin - mean(Res$Risk_difference_Estimate_CL_TSP_Bin))**2)))
    
    Coverage_rate_TSP_B     = length(which( True_RD >= Res$LL_95 & True_RD <= Res$UL_95))/nrow(Res)
    Coverage_rate           = c(Coverage_rate,
                                Coverage_rate_TSP_B)
    
    CI_width_TSP_B          = mean(Res$UL_95 - Res$LL_95)  
    CI_width                = c(CI_width,CI_width_TSP_B)
    P_025_CIW               = c(P_025_CIW, as.numeric(quantile(Res$UL_95 - Res$LL_95,probs = 0.025)))
    P_975_CIW               = c(P_975_CIW, as.numeric(quantile(Res$UL_95 - Res$LL_95,probs = 0.975)))
    MC_CIW                  = c(MC_CIW,
                                sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum(((Res$UL_95 - Res$LL_95) - mean(Res$UL_95 - Res$LL_95))**2)))
    
    MC_CR                   = c(MC_CR,
                                sqrt(Coverage_rate_TSP_B *(1-Coverage_rate_TSP_B)/nrow(Res)))
    
    Emp_SE_TSP_B            = sqrt(1/(nrow(Res)-1)*sum((Res$Risk_difference_Estimate_CL_TSP_Bin-mean(Res$Risk_difference_Estimate_CL_TSP_Bin))**2))
    Emp_SE                  = c(Emp_SE,
                                Emp_SE_TSP_B)
    
    ModSE_TSP_B             = sqrt(mean(Res$SE_itt**2))
    ModSE                   = c(ModSE,
                                ModSE_TSP_B)
    
    RE_SE_TSP_B             = (ModSE_TSP_B/Emp_SE_TSP_B - 1)*100
    RE_SE                   = c(RE_SE,
                                RE_SE_TSP_B)
    
    var_meth_TSP_B          = var(Res$SE_itt**2)
    var_meth                = c(var_meth,
                                var_meth_TSP_B)
    
    nrow_meth               = c(nrow_meth,
                                nrow(Res))
    
    MC_RESE                 = c(MC_RESE,
                                100*(ModSE_TSP_B/Emp_SE_TSP_B) * sqrt(((var_meth_TSP_B)/(4*nrow(Res)*ModSE_TSP_B))+
                                                                        (1/(2*(nrow(Res)-1)))))
    
    Methds_vec = c(Methds_vec,"CL_TSP_Bin")
  }
  
  if("GC_indv" %in% methods){
    
    Resu_GC_indv_file       = paste(Resu_file,"/GC_Bin_analyse_indv",sep = "")
    
    Res                     <- read.csv2(paste(Resu_GC_indv_file,"/Data_output_CL_GC_indv_Bin_Scenario_",n,".csv",sep=""))
    Res                     <- na.omit(Res)
    Res                     <- Res[which(Res$itt_para %in% itt_spec),]
    
    Mean_RD                 = c(Mean_RD,
                                mean(Res$Risk_difference_Estimate_CL_GC_Bin))
    
    Bias_relative_Mean      = c(Bias_relative_Mean,
                                mean(Res$Relative_bias_iteration_CL_GC_Bin))
    MC_RelativeBias         = c(MC_RelativeBias,
                                sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum(((Res$Risk_difference_Estimate_CL_GC_Bin - mean(Res$Risk_difference_Estimate_CL_GC_Bin))/mean(Res$Risk_difference_Estimate_CL_GC_Bin) - mean(Res$Relative_bias_iteration_CL_GC_Bin)/100)**2)))
    
    
    Bias_Mean               = c(Bias_Mean,
                                mean(Res$Bias_iteration_CL_GC_Bin))
    MC_bias                 = c(MC_bias,
                                sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum((Res$Risk_difference_Estimate_CL_GC_Bin - mean(Res$Risk_difference_Estimate_CL_GC_Bin))**2)))
    
    Coverage_rate_GC_indv   = length(which( True_RD >= Res$LL_95 & True_RD <= Res$UL_95))/nrow(Res)
    Coverage_rate           = c(Coverage_rate,
                                Coverage_rate_GC_indv)
    
    CI_width_GC_indv        = mean(Res$UL_95 - Res$LL_95)  
    CI_width                = c(CI_width,CI_width_GC_indv)
    P_025_CIW               = c(P_025_CIW, as.numeric(quantile(Res$UL_95 - Res$LL_95,probs = 0.025)))
    P_975_CIW               = c(P_975_CIW, as.numeric(quantile(Res$UL_95 - Res$LL_95,probs = 0.975)))
    MC_CIW                  = c(MC_CIW,
                                sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum(((Res$UL_95 - Res$LL_95) - mean(Res$UL_95 - Res$LL_95))**2)))
    
    MC_CR                   = c(MC_CR,
                                sqrt(Coverage_rate_GC_indv *(1-Coverage_rate_GC_indv)/nrow(Res)))
    
    Emp_SE_GC_indv          = sqrt(1/(nrow(Res)-1)*sum((Res$Risk_difference_Estimate_CL_GC_Bin-mean(Res$Risk_difference_Estimate_CL_GC_Bin))**2))
    Emp_SE                  = c(Emp_SE,
                                Emp_SE_GC_indv)
    
    ModSE_GC_indv           = sqrt(mean(Res$SE_itt**2))
    ModSE                   = c(ModSE,
                                ModSE_GC_indv)
    
    RE_SE_GC_indv           = (ModSE_GC_indv/Emp_SE_GC_indv - 1)*100
    RE_SE                   = c(RE_SE,
                                RE_SE_GC_indv)
    
    var_meth_GC_indv        = var(Res$SE_itt**2)
    var_meth                = c(var_meth,
                                var_meth_GC_indv)
    
    nrow_meth               = c(nrow_meth,
                                nrow(Res))
    
    MC_RESE                 = c(MC_RESE,
                                100*(ModSE_GC_indv/Emp_SE_GC_indv) * sqrt(((var_meth_GC_indv)/(4*nrow(Res)*ModSE_GC_indv))+
                                                                            (1/(2*(nrow(Res)-1)))))
    
    Methds_vec = c(Methds_vec,"CL_GC_indv")
    
  }
  
  if("GC_both" %in% methods){
    Resu_GC_both_file       = paste(Resu_file,"/GC_Bin_analyse_both",sep = "")
    
    Res                     <- read.csv2(paste(Resu_GC_both_file,"/Data_output_CL_GC_both_Bin_Scenario_",n,".csv",sep=""))
    Res                     <- na.omit(Res)
    Res                     <- Res[which(Res$itt_para %in% itt_spec),]
    
    Mean_RD                 = c(Mean_RD,
                                mean(Res$Risk_difference_Estimate_CL_GC_Bin))
    
    Bias_relative_Mean      = c(Bias_relative_Mean,
                                mean(Res$Relative_bias_iteration_CL_GC_Bin))
    MC_RelativeBias         = c(MC_RelativeBias,
                                sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum(((Res$Risk_difference_Estimate_CL_GC_Bin - mean(Res$Risk_difference_Estimate_CL_GC_Bin))/mean(Res$Risk_difference_Estimate_CL_GC_Bin) - mean(Res$Relative_bias_iteration_CL_GC_Bin)/100)**2)))
    
    Bias_Mean               = c(Bias_Mean,
                                mean(Res$Bias_iteration_CL_GC_Bin))
    MC_bias                 = c(MC_bias,
                                sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum((Res$Risk_difference_Estimate_CL_GC_Bin - mean(Res$Risk_difference_Estimate_CL_GC_Bin))**2)))
    
    Coverage_rate_GC_both   = length(which( True_RD >= Res$LL_95 & True_RD <= Res$UL_95))/nrow(Res)
    Coverage_rate           = c(Coverage_rate,
                                Coverage_rate_GC_both)
    
    CI_width_GC_both        = mean(Res$UL_95 - Res$LL_95)  
    CI_width                = c(CI_width,CI_width_GC_both)
    P_025_CIW               = c(P_025_CIW, as.numeric(quantile(Res$UL_95 - Res$LL_95,probs = 0.025)))
    P_975_CIW               = c(P_975_CIW, as.numeric(quantile(Res$UL_95 - Res$LL_95,probs = 0.975)))
    MC_CIW                  = c(MC_CIW,
                                sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum(((Res$UL_95 - Res$LL_95) - mean(Res$UL_95 - Res$LL_95))**2)))
    
    MC_CR                   = c(MC_CR,
                                sqrt(Coverage_rate_GC_both *(1-Coverage_rate_GC_both)/nrow(Res)))
    
    Emp_SE_GC_both          = sqrt(1/(nrow(Res)-1)*sum((Res$Risk_difference_Estimate_CL_GC_Bin-mean(Res$Risk_difference_Estimate_CL_GC_Bin))**2))
    Emp_SE                  = c(Emp_SE,
                                Emp_SE_GC_both)
    
    ModSE_GC_both           = sqrt(mean(Res$SE_itt**2))
    ModSE                   = c(ModSE,
                                ModSE_GC_both)
    
    RE_SE_GC_both           = (ModSE_GC_both/Emp_SE_GC_both - 1)*100
    RE_SE                   = c(RE_SE,
                                RE_SE_GC_both)
    
    var_meth_GC_both        = var(Res$SE_itt**2)
    var_meth                = c(var_meth,
                                var_meth_GC_both)
    
    nrow_meth               = c(nrow_meth,
                                nrow(Res))
    
    MC_RESE                 = c(MC_RESE,
                                100*(ModSE_GC_both/Emp_SE_GC_both) * sqrt(((var_meth_GC_both)/(4*nrow(Res)*ModSE_GC_both))+
                                                                            (1/(2*(nrow(Res)-1)))))
    
    Methds_vec = c(Methds_vec,"CL_GC_both")
  }
  
  if("GC_indv_indep" %in% methods){
    
    Resu_file_GC_indep_indv     = paste(Resu_file,"/GC_Bin_analyse_indv_indep/",sep = "")
    
    Res                         <- read.csv2(paste(Resu_file_GC_indep_indv,"Data_output_CL_GC_indv_indep_Bin_Scenario_",n,".csv",sep=""))
    Res                         <- na.omit(Res)
    Res                         <- Res[which(Res$itt_para %in% itt_spec),]
    
    Mean_RD                     = c(Mean_RD,
                                    mean(Res$Risk_difference_Estimate_CL_GC_Bin))
    
    Bias_relative_Mean          = c(Bias_relative_Mean,
                                    mean(Res$Relative_bias_iteration_CL_GC_Bin))
    MC_RelativeBias             = c(MC_RelativeBias,
                                    sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum(((Res$Risk_difference_Estimate_CL_GC_Bin - mean(Res$Risk_difference_Estimate_CL_GC_Bin))/mean(Res$Risk_difference_Estimate_CL_GC_Bin) - mean(Res$Relative_bias_iteration_CL_GC_Bin)/100)**2)))
    
    
    Bias_Mean                   = c(Bias_Mean,
                                    mean(Res$Bias_iteration_CL_GC_Bin))
    MC_bias                     = c(MC_bias,
                                    sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum((Res$Risk_difference_Estimate_CL_GC_Bin - mean(Res$Risk_difference_Estimate_CL_GC_Bin))**2)))
    
    Coverage_rate_res     = length(which( True_RD >= Res$LL_95 & True_RD <= Res$UL_95))/nrow(Res)
    Coverage_rate               = c(Coverage_rate,
                                    Coverage_rate_res)
    
    CI_width_GC_indep_indv  = mean(Res$UL_95 - Res$LL_95)  
    CI_width                = c(CI_width,CI_width_GC_indep_indv)
    P_025_CIW               = c(P_025_CIW, as.numeric(quantile(Res$UL_95 - Res$LL_95,probs = 0.025)))
    P_975_CIW               = c(P_975_CIW, as.numeric(quantile(Res$UL_95 - Res$LL_95,probs = 0.975)))
    MC_CIW                  = c(MC_CIW,
                                sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum(((Res$UL_95 - Res$LL_95) - mean(Res$UL_95 - Res$LL_95))**2)))
    
    MC_CR                       = c(MC_CR,
                                    sqrt(Coverage_rate_res *(1-Coverage_rate_res)/nrow(Res)))
    
    Emp_SE_res           = sqrt(1/(nrow(Res)-1)*sum((Res$Risk_difference_Estimate_CL_GC_Bin-mean(Res$Risk_difference_Estimate_CL_GC_Bin))**2))
    Emp_SE                      = c(Emp_SE,
                                    Emp_SE_res)
    
    ModSE_res             = sqrt(mean(Res$SE_itt**2))
    ModSE                       = c(ModSE,
                                    ModSE_res)
    
    RE_SE_res             = (ModSE_res/Emp_SE_res - 1)*100
    RE_SE                       = c(RE_SE,
                                    RE_SE_res)
    
    var_meth_res         = var(Res$SE_itt**2)
    var_meth                    = c(var_meth,
                                    var_meth_res)
    
    nrow_meth                   = c(nrow_meth,
                                    nrow(Res))
    
    MC_RESE                     = c(MC_RESE,
                                    100*(ModSE_res/Emp_SE_res) * sqrt(((var_meth_res)/(4*nrow(Res)*ModSE_res))+
                                                                        (1/(2*(nrow(Res)-1)))))
    
    Methds_vec                  = c(Methds_vec,"CL_GC_indep_indv")
  }
  
  if("GC_both_indep" %in% methods){
    
    Resu_file_GC_indep_both     = paste(Resu_file,"/GC_Bin_analyse_both_indep/",sep = "")
    
    Res                         <- read.csv2(paste(Resu_file_GC_indep_both,"Data_output_CL_GC_both_indep_Bin_Scenario_",n,".csv",sep=""))
    Res                         <- na.omit(Res)
    Res                         <- Res[which(Res$itt_para %in% itt_spec),]
    
    Mean_RD                     = c(Mean_RD,
                                    mean(Res$Risk_difference_Estimate_CL_GC_Bin))
    
    Bias_relative_Mean          = c(Bias_relative_Mean,
                                    mean(Res$Relative_bias_iteration_CL_GC_Bin))
    MC_RelativeBias             = c(MC_RelativeBias,
                                    sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum(((Res$Risk_difference_Estimate_CL_GC_Bin - mean(Res$Risk_difference_Estimate_CL_GC_Bin))/mean(Res$Risk_difference_Estimate_CL_GC_Bin) - mean(Res$Relative_bias_iteration_CL_GC_Bin)/100)**2)))
    
    
    Bias_Mean                   = c(Bias_Mean,
                                    mean(Res$Bias_iteration_CL_GC_Bin))
    MC_bias                     = c(MC_bias,
                                    sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum((Res$Risk_difference_Estimate_CL_GC_Bin - mean(Res$Risk_difference_Estimate_CL_GC_Bin))**2)))
    
    Coverage_rate_res     = length(which( True_RD >= Res$LL_95 & True_RD <= Res$UL_95))/nrow(Res)
    Coverage_rate               = c(Coverage_rate,
                                    Coverage_rate_res)
    
    CI_width_GC_indep_both  = mean(Res$UL_95 - Res$LL_95)  
    CI_width                = c(CI_width,CI_width_GC_indep_both)
    P_025_CIW               = c(P_025_CIW, as.numeric(quantile(Res$UL_95 - Res$LL_95,probs = 0.025)))
    P_975_CIW               = c(P_975_CIW, as.numeric(quantile(Res$UL_95 - Res$LL_95,probs = 0.975)))
    MC_CIW                  = c(MC_CIW,
                                sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum(((Res$UL_95 - Res$LL_95) - mean(Res$UL_95 - Res$LL_95))**2)))
    
    MC_CR                       = c(MC_CR,
                                    sqrt(Coverage_rate_res *(1-Coverage_rate_res)/nrow(Res)))
    
    Emp_SE_res           = sqrt(1/(nrow(Res)-1)*sum((Res$Risk_difference_Estimate_CL_GC_Bin-mean(Res$Risk_difference_Estimate_CL_GC_Bin))**2))
    Emp_SE                      = c(Emp_SE,
                                    Emp_SE_res)
    
    ModSE_res             = sqrt(mean(Res$SE_itt**2))
    ModSE                       = c(ModSE,
                                    ModSE_res)
    
    RE_SE_res             = (ModSE_res/Emp_SE_res - 1)*100
    RE_SE                       = c(RE_SE,
                                    RE_SE_res)
    
    var_meth_res         = var(Res$SE_itt**2)
    var_meth                    = c(var_meth,
                                    var_meth_res)
    
    nrow_meth                   = c(nrow_meth,
                                    nrow(Res))
    
    MC_RESE                     = c(MC_RESE,
                                    100*(ModSE_res/Emp_SE_res) * sqrt(((var_meth_res)/(4*nrow(Res)*ModSE_res))+
                                                                        (1/(2*(nrow(Res)-1)))))
    
    Methds_vec                  = c(Methds_vec,"CL_GC_indep_both")
    
  }
  
  if("GC_indv_GLMM" %in% methods){
    
    Resu_file_GC_GLMM_indv      = paste(Resu_file,"/GC_GLMM_Bin_analyse_indv/",sep = "")
    
    Res                         <- read.csv2(paste(Resu_file_GC_GLMM_indv,"Data_output_CL_GC_indv_Bin_glmm_Scenario_",n,".csv",sep=""))
    Res                         <- na.omit(Res)
    Res                         <- Res[which(Res$itt_para %in% itt_spec),]
    
    Mean_RD                     = c(Mean_RD,
                                    mean(Res$Risk_difference_Estimate_CL_GC_Bin))
    
    Bias_relative_Mean          = c(Bias_relative_Mean,
                                    mean(Res$Relative_bias_iteration_CL_GC_Bin))
    MC_RelativeBias             = c(MC_RelativeBias,
                                    sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum(((Res$Risk_difference_Estimate_CL_GC_Bin - mean(Res$Risk_difference_Estimate_CL_GC_Bin))/mean(Res$Risk_difference_Estimate_CL_GC_Bin) - mean(Res$Relative_bias_iteration_CL_GC_Bin)/100)**2)))
    
    
    Bias_Mean                   = c(Bias_Mean,
                                    mean(Res$Bias_iteration_CL_GC_Bin))
    MC_bias                     = c(MC_bias,
                                    sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum((Res$Risk_difference_Estimate_CL_GC_Bin - mean(Res$Risk_difference_Estimate_CL_GC_Bin))**2)))
    
    Coverage_rate_res           = length(which( True_RD >= Res$LL_95 & True_RD <= Res$UL_95))/nrow(Res)
    Coverage_rate               = c(Coverage_rate,
                                    Coverage_rate_res)
    
    CI_width_GC_GLMM_indv       = mean(Res$UL_95 - Res$LL_95)  
    CI_width                    = c(CI_width,CI_width_GC_GLMM_indv)
    P_025_CIW               = c(P_025_CIW, as.numeric(quantile(Res$UL_95 - Res$LL_95,probs = 0.025)))
    P_975_CIW               = c(P_975_CIW, as.numeric(quantile(Res$UL_95 - Res$LL_95,probs = 0.975)))
    MC_CIW                      = c(MC_CIW,
                                    sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum(((Res$UL_95 - Res$LL_95) - mean(Res$UL_95 - Res$LL_95))**2)))
    
    MC_CR                       = c(MC_CR,
                                    sqrt(Coverage_rate_res *(1-Coverage_rate_res)/nrow(Res)))
    
    Emp_SE_res                  = sqrt(1/(nrow(Res)-1)*sum((Res$Risk_difference_Estimate_CL_GC_Bin-mean(Res$Risk_difference_Estimate_CL_GC_Bin))**2))
    Emp_SE                      = c(Emp_SE,
                                    Emp_SE_res)
    
    ModSE_res                   = sqrt(mean(Res$SE_itt**2))
    ModSE                       = c(ModSE,
                                    ModSE_res)
    
    RE_SE_res                   = (ModSE_res/Emp_SE_res - 1)*100
    RE_SE                       = c(RE_SE,
                                    RE_SE_res)
    
    var_meth_res                = var(Res$SE_itt**2)
    var_meth                    = c(var_meth,
                                    var_meth_res)
    
    nrow_meth                   = c(nrow_meth,
                                    nrow(Res))
    
    MC_RESE                     = c(MC_RESE,
                                    100*(ModSE_res/Emp_SE_res) * sqrt(((var_meth_res)/(4*nrow(Res)*ModSE_res))+
                                                                        (1/(2*(nrow(Res)-1)))))
    
    Methds_vec                  = c(Methds_vec,"CL_GC_GLMM_indv")
  }
  
  if("GC_both_GLMM" %in% methods){
    
    Resu_file_GC_GLMM_both      = paste(Resu_file,"/GC_GLMM_Bin_analyse_both/",sep = "")
    
    Res                         <- read.csv2(paste(Resu_file_GC_GLMM_both,"Data_output_CL_GC_both_Bin_glmm_Scenario_",n,".csv",sep=""))
    Res                         <- na.omit(Res)
    Res                         <- Res[which(Res$itt_para %in% itt_spec),]
    
    Mean_RD                     = c(Mean_RD,
                                    mean(Res$Risk_difference_Estimate_CL_GC_Bin))
    
    Bias_relative_Mean          = c(Bias_relative_Mean,
                                    mean(Res$Relative_bias_iteration_CL_GC_Bin))
    MC_RelativeBias             = c(MC_RelativeBias,
                                    sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum(((Res$Risk_difference_Estimate_CL_GC_Bin - mean(Res$Risk_difference_Estimate_CL_GC_Bin))/mean(Res$Risk_difference_Estimate_CL_GC_Bin) - mean(Res$Relative_bias_iteration_CL_GC_Bin)/100)**2)))
    
    
    Bias_Mean                   = c(Bias_Mean,
                                    mean(Res$Bias_iteration_CL_GC_Bin))
    MC_bias                     = c(MC_bias,
                                    sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum((Res$Risk_difference_Estimate_CL_GC_Bin - mean(Res$Risk_difference_Estimate_CL_GC_Bin))**2)))
    
    Coverage_rate_res           = length(which( True_RD >= Res$LL_95 & True_RD <= Res$UL_95))/nrow(Res)
    Coverage_rate               = c(Coverage_rate,
                                    Coverage_rate_res)
    
    CI_width_GC_GLMM_both       = mean(Res$UL_95 - Res$LL_95)  
    CI_width                    = c(CI_width,CI_width_GC_GLMM_both)
    P_025_CIW               = c(P_025_CIW, as.numeric(quantile(Res$UL_95 - Res$LL_95,probs = 0.025)))
    P_975_CIW               = c(P_975_CIW, as.numeric(quantile(Res$UL_95 - Res$LL_95,probs = 0.975)))
    MC_CIW                      = c(MC_CIW,
                                    sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum(((Res$UL_95 - Res$LL_95) - mean(Res$UL_95 - Res$LL_95))**2)))
    
    MC_CR                       = c(MC_CR,
                                    sqrt(Coverage_rate_res *(1-Coverage_rate_res)/nrow(Res)))
    
    Emp_SE_res                  = sqrt(1/(nrow(Res)-1)*sum((Res$Risk_difference_Estimate_CL_GC_Bin-mean(Res$Risk_difference_Estimate_CL_GC_Bin))**2))
    Emp_SE                      = c(Emp_SE,
                                    Emp_SE_res)
    
    ModSE_res                   = sqrt(mean(Res$SE_itt**2))
    ModSE                       = c(ModSE,
                                    ModSE_res)
    
    RE_SE_res                   = (ModSE_res/Emp_SE_res - 1)*100
    RE_SE                       = c(RE_SE,
                                    RE_SE_res)
    
    var_meth_res                = var(Res$SE_itt**2)
    var_meth                    = c(var_meth,
                                    var_meth_res)
    
    nrow_meth                   = c(nrow_meth,
                                    nrow(Res))
    
    MC_RESE                     = c(MC_RESE,
                                    100*(ModSE_res/Emp_SE_res) * sqrt(((var_meth_res)/(4*nrow(Res)*ModSE_res))+
                                                                        (1/(2*(nrow(Res)-1)))))
    
    Methds_vec                  = c(Methds_vec,"CL_GC_GLMM_both")
  }
  
  if("TMLE_indv" %in% methods){
    
    Resu_TMLE_file_indv     = paste(Resu_file,"/TMLE_analyse_indvcov/",sep = "")
    
    Res                     <- read.csv2(paste(Resu_TMLE_file_indv,"/Data_output_CL_TMLE_indvcov_Scenario_",n,".csv",sep=""))
    Res                     <- na.omit(Res)
    Res                     <- Res[which(Res$Itt %in% itt_spec),]
    
    Mean_RD                 = c(Mean_RD,
                                mean(Res$RD_CL_TMLE_JAPM))
    
    Bias_relative_Mean      = c(Bias_relative_Mean,
                                mean(Res$Relative_Bias_JAPM))
    MC_RelativeBias         = c(MC_RelativeBias,
                                sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum(((Res$RD_CL_TMLE_JAPM - mean(Res$RD_CL_TMLE_JAPM))/mean(Res$RD_CL_TMLE_JAPM) - mean(Res$Relative_Bias_JAPM)/100)**2)))
    
    
    Bias_Mean               = c(Bias_Mean,
                                mean(Res$Bias_JAPM))
    MC_bias                 = c(MC_bias,
                                sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum((Res$RD_CL_TMLE_JAPM - mean(Res$RD_CL_TMLE_JAPM))**2)))
    
    Coverage_rate_TMLE_indv = length(which( True_RD >= Res$LL_95_CL_TMLE_JAPM & True_RD <= Res$UL_95_CL_TMLE_JAPM))/nrow(Res)
    Coverage_rate           = c(Coverage_rate,
                                Coverage_rate_TMLE_indv)
    
    CI_width_TMLE_indv      = mean(Res$UL_95_CL_TMLE_JAPM - Res$LL_95_CL_TMLE_JAPM)  
    CI_width                = c(CI_width,CI_width_TMLE_indv)
    P_025_CIW               = c(P_025_CIW, as.numeric(quantile(Res$UL_95_CL_TMLE_JAPM - Res$LL_95_CL_TMLE_JAPM,probs = 0.025)))
    P_975_CIW               = c(P_975_CIW, as.numeric(quantile(Res$UL_95_CL_TMLE_JAPM - Res$LL_95_CL_TMLE_JAPM,probs = 0.975)))
    MC_CIW                  = c(MC_CIW,
                                sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum(((Res$UL_95_CL_TMLE_JAPM - Res$LL_95_CL_TMLE_JAPM) - mean(Res$UL_95_CL_TMLE_JAPM - Res$LL_95_CL_TMLE_JAPM))**2)))
    
    MC_CR                   = c(MC_CR,
                                sqrt(Coverage_rate_TMLE_indv *(1-Coverage_rate_TMLE_indv)/nrow(Res)))
    
    Emp_SE_TMLE_indv        = sqrt(1/(nrow(Res)-1)*sum((Res$RD_CL_TMLE_JAPM-mean(Res$RD_CL_TMLE_JAPM))**2))
    Emp_SE                  = c(Emp_SE,
                                Emp_SE_TMLE_indv)
    
    ModSE_TMLE_indv         = sqrt(mean(Res$SE_CL_TMLE_JAPM**2))
    ModSE                   = c(ModSE,
                                ModSE_TMLE_indv)
    
    RE_SE_TMLE_indv         = (ModSE_TMLE_indv/Emp_SE_TMLE_indv - 1)*100
    RE_SE                   = c(RE_SE,
                                RE_SE_TMLE_indv)
    
    var_meth_TMLE_indv      = var(Res$SE_CL_TMLE_JAPM**2)
    var_meth                = c(var_meth,
                                var_meth_TMLE_indv)
    
    nrow_meth               = c(nrow_meth,
                                nrow(Res))
    
    MC_RESE                 = c(MC_RESE,
                                100*(ModSE_TMLE_indv/Emp_SE_TMLE_indv) * sqrt(((var_meth_TMLE_indv)/(4*nrow(Res)*ModSE_TMLE_indv))+
                                                                                (1/(2*(nrow(Res)-1)))))
    
    Methds_vec = c(Methds_vec,"CL_TMLE_indv")
  }
  
  if("TMLE_both" %in% methods){
    
    Resu_TMLE_file_both     = paste(Resu_file,"/TMLE_analyse_bothcov/",sep = "")
    
    Res                     <- read.csv2(paste(Resu_TMLE_file_both,"/Data_output_CL_TMLE_bothcov_Scenario_",n,".csv",sep=""))
    Res                     <- na.omit(Res)
    Res                     <- Res[which(Res$Itt %in% itt_spec),]
    
    Mean_RD                 = c(Mean_RD,
                                mean(Res$RD_CL_TMLE_JAPM))
    
    Bias_relative_Mean      = c(Bias_relative_Mean,
                                mean(Res$Relative_Bias_JAPM))
    MC_RelativeBias         = c(MC_RelativeBias,
                                sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum(((Res$RD_CL_TMLE_JAPM - mean(Res$RD_CL_TMLE_JAPM))/mean(Res$RD_CL_TMLE_JAPM) - mean(Res$Relative_Bias_JAPM)/100)**2)))
    
    
    Bias_Mean               = c(Bias_Mean,
                                mean(Res$Bias_JAPM))
    MC_bias                 = c(MC_bias,
                                sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum((Res$RD_CL_TMLE_JAPM - mean(Res$RD_CL_TMLE_JAPM))**2)))
    
    Coverage_rate_TMLE_both = length(which( True_RD >= Res$LL_95_CL_TMLE_JAPM & True_RD <= Res$UL_95_CL_TMLE_JAPM))/nrow(Res)
    Coverage_rate           = c(Coverage_rate,
                                Coverage_rate_TMLE_both)
    
    CI_width_TMLE_both      = mean(Res$UL_95_CL_TMLE_JAPM - Res$LL_95_CL_TMLE_JAPM)  
    CI_width                = c(CI_width,CI_width_TMLE_both)
    P_025_CIW               = c(P_025_CIW, as.numeric(quantile(Res$UL_95_CL_TMLE_JAPM - Res$LL_95_CL_TMLE_JAPM,probs = 0.025)))
    P_975_CIW               = c(P_975_CIW, as.numeric(quantile(Res$UL_95_CL_TMLE_JAPM - Res$LL_95_CL_TMLE_JAPM,probs = 0.975)))
    MC_CIW                  = c(MC_CIW,
                                sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum(((Res$UL_95_CL_TMLE_JAPM - Res$LL_95_CL_TMLE_JAPM) - mean(Res$UL_95_CL_TMLE_JAPM - Res$LL_95_CL_TMLE_JAPM))**2)))
    
    MC_CR                   = c(MC_CR,
                                sqrt(Coverage_rate_TMLE_both *(1-Coverage_rate_TMLE_both)/nrow(Res)))
    
    Emp_SE_TMLE_both        = sqrt(1/(nrow(Res)-1)*sum((Res$RD_CL_TMLE_JAPM-mean(Res$RD_CL_TMLE_JAPM))**2))
    Emp_SE                  = c(Emp_SE,
                                Emp_SE_TMLE_both)
    
    ModSE_TMLE_both         = sqrt(mean(Res$SE_CL_TMLE_JAPM**2))
    ModSE                   = c(ModSE,
                                ModSE_TMLE_both)
    
    RE_SE_TMLE_both         = (ModSE_TMLE_both/Emp_SE_TMLE_both - 1)*100
    RE_SE                   = c(RE_SE,
                                RE_SE_TMLE_both)
    
    var_meth_TMLE_both      = var(Res$SE_CL_TMLE_JAPM**2)
    var_meth                = c(var_meth,
                                var_meth_TMLE_both)
    
    nrow_meth               = c(nrow_meth,
                                nrow(Res))
    
    MC_RESE                 = c(MC_RESE,
                                100*(ModSE_TMLE_both/Emp_SE_TMLE_both) * sqrt(((var_meth_TMLE_both)/(4*nrow(Res)*ModSE_TMLE_both))+
                                                                                (1/(2*(nrow(Res)-1)))))
    
    Methds_vec = c(Methds_vec,"CL_TMLE_both")
  }
  
  if("TMLE_APS" %in% methods){
    
    Resu_TMLE_file_APS          = paste(Resu_file,"/TMLE_analyse_APS/",sep = "")
    
    Res                         <- read.csv2(paste(Resu_TMLE_file_APS,"Data_output_CL_TMLE_APS_Scenario_",n,".csv",sep=""))
    Res                         <- na.omit(Res)
    Res                         <- Res[which(Res$Itt %in% itt_spec),]
    
    Mean_RD                     = c(Mean_RD,
                                    mean(Res$RD_CL_TMLE_JAPM))
    
    Bias_relative_Mean          = c(Bias_relative_Mean,
                                    mean(Res$Relative_Bias_JAPM))
    MC_RelativeBias             = c(MC_RelativeBias,
                                    sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum(((Res$RD_CL_TMLE_JAPM - mean(Res$RD_CL_TMLE_JAPM))/mean(Res$RD_CL_TMLE_JAPM) - mean(Res$Relative_Bias_JAPM)/100)**2)))
    
    
    Bias_Mean                   = c(Bias_Mean,
                                    mean(Res$Bias_JAPM))
    MC_bias                     = c(MC_bias,
                                    sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum((Res$RD_CL_TMLE_JAPM - mean(Res$RD_CL_TMLE_JAPM))**2)))
    
    Coverage_rate_res           = length(which( True_RD >= Res$LL_95_CL_TMLE_JAPM & True_RD <= Res$UL_95_CL_TMLE_JAPM))/nrow(Res)
    Coverage_rate               = c(Coverage_rate,
                                    Coverage_rate_res)
    
    CI_width_TMLE_APS           = mean(Res$UL_95_CL_TMLE_JAPM - Res$LL_95_CL_TMLE_JAPM)  
    CI_width                    = c(CI_width,CI_width_TMLE_APS)
    P_025_CIW                   = c(P_025_CIW, as.numeric(quantile(Res$UL_95_CL_TMLE_JAPM - Res$LL_95_CL_TMLE_JAPM,probs = 0.025)))
    P_975_CIW                   = c(P_975_CIW, as.numeric(quantile(Res$UL_95_CL_TMLE_JAPM - Res$LL_95_CL_TMLE_JAPM,probs = 0.975)))
    MC_CIW                      = c(MC_CIW,
                                    sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum(((Res$UL_95_CL_TMLE_JAPM - Res$LL_95_CL_TMLE_JAPM) - mean(Res$UL_95_CL_TMLE_JAPM - Res$LL_95_CL_TMLE_JAPM))**2)))
    
    MC_CR                       = c(MC_CR,
                                    sqrt(Coverage_rate_res *(1-Coverage_rate_res)/nrow(Res)))
    
    Emp_SE_res                  = sqrt(1/(nrow(Res)-1)*sum((Res$RD_CL_TMLE_JAPM-mean(Res$RD_CL_TMLE_JAPM))**2))
    Emp_SE                      = c(Emp_SE,
                                    Emp_SE_res)
    
    ModSE_res                   = sqrt(mean(Res$SE_CL_TMLE_JAPM**2))
    ModSE                       = c(ModSE,
                                    ModSE_res)
    
    RE_SE_res                   = (ModSE_res/Emp_SE_res - 1)*100
    RE_SE                       = c(RE_SE,
                                    RE_SE_res)
    
    var_meth_res                = var(Res$SE_CL_TMLE_JAPM**2)
    var_meth                    = c(var_meth,
                                    var_meth_res)
    
    nrow_meth                   = c(nrow_meth,
                                    nrow(Res))
    
    MC_RESE                     = c(MC_RESE,
                                    100*(ModSE_res/Emp_SE_res) * sqrt(((var_meth_res)/(4*nrow(Res)*ModSE_res))+
                                                                        (1/(2*(nrow(Res)-1)))))
    
    Methds_vec                  = c(Methds_vec,"CL_TMLE_APS")
    
  }
  
  if("TMLE_APS_Balzer" %in% methods){
    
    Resu_TMLE_file_APS_Balzer   = paste(Resu_file,"/TMLE_analyse_Balzer_APS/",sep = "")
    
    Res                         <- read.csv2(paste(Resu_TMLE_file_APS_Balzer,"Data_output_CL_TMLE_Balzer_APS_Scenario_",n,".csv",sep=""))
    Res                         <- na.omit(Res)
    Res                         <- Res[which(Res$Itt %in% itt_spec),]
    
    Mean_RD                     = c(Mean_RD,
                                    mean(Res$RD_CL_TMLE_Balzer))
    
    Bias_relative_Mean          = c(Bias_relative_Mean,
                                    mean(Res$Relative_Bias_Balzer))
    MC_RelativeBias             = c(MC_RelativeBias,
                                    sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum(((Res$RD_CL_TMLE_Balzer - mean(Res$RD_CL_TMLE_Balzer))/mean(Res$RD_CL_TMLE_Balzer) - mean(Res$Relative_Bias_Balzer)/100)**2)))
    
    
    Bias_Mean                   = c(Bias_Mean,
                                    mean(Res$Bias_Balzer))
    MC_bias                     = c(MC_bias,
                                    sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum((Res$RD_CL_TMLE_Balzer - mean(Res$RD_CL_TMLE_Balzer))**2)))
    
    Coverage_rate_res           = length(which( True_RD >= Res$LL_95_CL_TMLE_Balzer & True_RD <= Res$UL_95_CL_TMLE_Balzer))/nrow(Res)
    Coverage_rate               = c(Coverage_rate,
                                    Coverage_rate_res)
    
    CI_width_TMLE_APS_Balzer    = mean(Res$UL_95_CL_TMLE_Balzer - Res$LL_95_CL_TMLE_Balzer)  
    CI_width                    = c(CI_width,CI_width_TMLE_APS_Balzer)
    P_025_CIW                   = c(P_025_CIW, as.numeric(quantile(Res$UL_95_CL_TMLE_Balzer - Res$LL_95_CL_TMLE_Balzer,probs = 0.025)))
    P_975_CIW                   = c(P_975_CIW, as.numeric(quantile(Res$UL_95_CL_TMLE_Balzer - Res$LL_95_CL_TMLE_Balzer,probs = 0.975)))
    MC_CIW                      = c(MC_CIW,
                                    sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum(((Res$UL_95_CL_TMLE_Balzer - Res$LL_95_CL_TMLE_Balzer) - mean(Res$UL_95_CL_TMLE_Balzer - Res$LL_95_CL_TMLE_Balzer))**2)))
    
    
    MC_CR                       = c(MC_CR,
                                    sqrt(Coverage_rate_res *(1-Coverage_rate_res)/nrow(Res)))
    
    Emp_SE_res                  = sqrt(1/(nrow(Res)-1)*sum((Res$RD_CL_TMLE_Balzer-mean(Res$RD_CL_TMLE_Balzer))**2))
    Emp_SE                      = c(Emp_SE,
                                    Emp_SE_res)
    
    ModSE_res                   = sqrt(mean(Res$SE_CL_TMLE_Balzer**2))
    ModSE                       = c(ModSE,
                                    ModSE_res)
    
    RE_SE_res                   = (ModSE_res/Emp_SE_res - 1)*100
    RE_SE                       = c(RE_SE,
                                    RE_SE_res)
    
    var_meth_res                = var(Res$SE_CL_TMLE_Balzer**2)
    var_meth                    = c(var_meth,
                                    var_meth_res)
    
    nrow_meth                   = c(nrow_meth,
                                    nrow(Res))
    
    MC_RESE                     = c(MC_RESE,
                                    100*(ModSE_res/Emp_SE_res) * sqrt(((var_meth_res)/(4*nrow(Res)*ModSE_res))+
                                                                        (1/(2*(nrow(Res)-1)))))
    
    Methds_vec                  = c(Methds_vec,"CL_TMLE_APS_Balzer")
  }
  
  if("TMLE_APS_Balzer_2" %in% methods){
    
    Resu_TMLE_file_APS_Balzer_2   = paste(Resu_file,"/TMLE_analyse_Balzer_APS_k_5_10/",sep = "")
    
    Res                         <- read.csv2(paste(Resu_TMLE_file_APS_Balzer,"Data_output_CL_TMLE_Balzer_APS_Scenario_",n,".csv",sep=""))
    Res                         <- na.omit(Res)
    Res                         <- Res[which(Res$Itt %in% itt_spec),]
    
    Mean_RD                     = c(Mean_RD,
                                    mean(Res$RD_CL_TMLE_Balzer))
    
    Bias_relative_Mean          = c(Bias_relative_Mean,
                                    mean(Res$Relative_Bias_Balzer))
    MC_RelativeBias             = c(MC_RelativeBias,
                                    sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum(((Res$RD_CL_TMLE_Balzer - mean(Res$RD_CL_TMLE_Balzer))/mean(Res$RD_CL_TMLE_Balzer) - mean(Res$Relative_Bias_Balzer)/100)**2)))
    
    
    Bias_Mean                   = c(Bias_Mean,
                                    mean(Res$Bias_Balzer))
    MC_bias                     = c(MC_bias,
                                    sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum((Res$RD_CL_TMLE_Balzer - mean(Res$RD_CL_TMLE_Balzer))**2)))
    
    Coverage_rate_res           = length(which( True_RD >= Res$LL_95_CL_TMLE_Balzer & True_RD <= Res$UL_95_CL_TMLE_Balzer))/nrow(Res)
    Coverage_rate               = c(Coverage_rate,
                                    Coverage_rate_res)
    
    CI_width_TMLE_APS_Balzer    = mean(Res$UL_95_CL_TMLE_Balzer - Res$LL_95_CL_TMLE_Balzer)  
    CI_width                    = c(CI_width,CI_width_TMLE_APS_Balzer)
    P_025_CIW                   = c(P_025_CIW, as.numeric(quantile(Res$UL_95_CL_TMLE_Balzer - Res$LL_95_CL_TMLE_Balzer,probs = 0.025)))
    P_975_CIW                   = c(P_975_CIW, as.numeric(quantile(Res$UL_95_CL_TMLE_Balzer - Res$LL_95_CL_TMLE_Balzer,probs = 0.975)))
    MC_CIW                      = c(MC_CIW,
                                    sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum(((Res$UL_95_CL_TMLE_Balzer - Res$LL_95_CL_TMLE_Balzer) - mean(Res$UL_95_CL_TMLE_Balzer - Res$LL_95_CL_TMLE_Balzer))**2)))
    
    
    MC_CR                       = c(MC_CR,
                                    sqrt(Coverage_rate_res *(1-Coverage_rate_res)/nrow(Res)))
    
    Emp_SE_res                  = sqrt(1/(nrow(Res)-1)*sum((Res$RD_CL_TMLE_Balzer-mean(Res$RD_CL_TMLE_Balzer))**2))
    Emp_SE                      = c(Emp_SE,
                                    Emp_SE_res)
    
    ModSE_res                   = sqrt(mean(Res$SE_CL_TMLE_Balzer**2))
    ModSE                       = c(ModSE,
                                    ModSE_res)
    
    RE_SE_res                   = (ModSE_res/Emp_SE_res - 1)*100
    RE_SE                       = c(RE_SE,
                                    RE_SE_res)
    
    var_meth_res                = var(Res$SE_CL_TMLE_Balzer**2)
    var_meth                    = c(var_meth,
                                    var_meth_res)
    
    nrow_meth                   = c(nrow_meth,
                                    nrow(Res))
    
    MC_RESE                     = c(MC_RESE,
                                    100*(ModSE_res/Emp_SE_res) * sqrt(((var_meth_res)/(4*nrow(Res)*ModSE_res))+
                                                                        (1/(2*(nrow(Res)-1)))))
    
    Methds_vec                  = c(Methds_vec,"CL_TMLE_APS_Balzer_2")
  }
  
  if("TMLE_APS_indv_Balzer" %in% methods){
    
    Resu_TMLE_file_APS_indv_Balzer   = paste(Resu_file,"/TMLE_analyse_Balzer_APS_indv/",sep = "")
    
    Res                         <- read.csv2(paste(Resu_TMLE_file_APS_indv_Balzer,"Data_output_CL_TMLE_Balzer_APS_indv_Scenario_",n,".csv",sep=""))
    Res                         <- na.omit(Res)
    Res                         <- Res[which(Res$Itt %in% itt_spec),]
    
    Mean_RD                     = c(Mean_RD,
                                    mean(Res$RD_CL_TMLE_Balzer))
    
    Bias_relative_Mean          = c(Bias_relative_Mean,
                                    mean(Res$Relative_Bias_Balzer))
    MC_RelativeBias             = c(MC_RelativeBias,
                                    sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum(((Res$RD_CL_TMLE_Balzer - mean(Res$RD_CL_TMLE_Balzer))/mean(Res$RD_CL_TMLE_Balzer) - mean(Res$Relative_Bias_Balzer)/100)**2)))
    
    
    Bias_Mean                   = c(Bias_Mean,
                                    mean(Res$Bias_Balzer))
    MC_bias                     = c(MC_bias,
                                    sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum((Res$RD_CL_TMLE_Balzer - mean(Res$RD_CL_TMLE_Balzer))**2)))
    
    Coverage_rate_res           = length(which( True_RD >= Res$LL_95_CL_TMLE_Balzer & True_RD <= Res$UL_95_CL_TMLE_Balzer))/nrow(Res)
    Coverage_rate               = c(Coverage_rate,
                                    Coverage_rate_res)
    
    CI_width_TMLE_APS_Balzer    = mean(Res$UL_95_CL_TMLE_Balzer - Res$LL_95_CL_TMLE_Balzer)  
    CI_width                    = c(CI_width,CI_width_TMLE_APS_Balzer)
    P_025_CIW                   = c(P_025_CIW, as.numeric(quantile(Res$UL_95_CL_TMLE_Balzer - Res$LL_95_CL_TMLE_Balzer,probs = 0.025)))
    P_975_CIW                   = c(P_975_CIW, as.numeric(quantile(Res$UL_95_CL_TMLE_Balzer - Res$LL_95_CL_TMLE_Balzer,probs = 0.975)))
    MC_CIW                      = c(MC_CIW,
                                    sqrt(1/(nrow(Res)*(nrow(Res)-1)) * sum(((Res$UL_95_CL_TMLE_Balzer - Res$LL_95_CL_TMLE_Balzer) - mean(Res$UL_95_CL_TMLE_Balzer - Res$LL_95_CL_TMLE_Balzer))**2)))
    
    
    MC_CR                       = c(MC_CR,
                                    sqrt(Coverage_rate_res *(1-Coverage_rate_res)/nrow(Res)))
    
    Emp_SE_res                  = sqrt(1/(nrow(Res)-1)*sum((Res$RD_CL_TMLE_Balzer-mean(Res$RD_CL_TMLE_Balzer))**2))
    Emp_SE                      = c(Emp_SE,
                                    Emp_SE_res)
    
    ModSE_res                   = sqrt(mean(Res$SE_CL_TMLE_Balzer**2))
    ModSE                       = c(ModSE,
                                    ModSE_res)
    
    RE_SE_res                   = (ModSE_res/Emp_SE_res - 1)*100
    RE_SE                       = c(RE_SE,
                                    RE_SE_res)
    
    var_meth_res                = var(Res$SE_CL_TMLE_Balzer**2)
    var_meth                    = c(var_meth,
                                    var_meth_res)
    
    nrow_meth                   = c(nrow_meth,
                                    nrow(Res))
    
    MC_RESE                     = c(MC_RESE,
                                    100*(ModSE_res/Emp_SE_res) * sqrt(((var_meth_res)/(4*nrow(Res)*ModSE_res))+
                                                                        (1/(2*(nrow(Res)-1)))))
    
    Methds_vec                  = c(Methds_vec,"CL_TMLE_APS_indv_Balzer")
  }
  
  
  tab_res = data.frame(
    Method = Methds_vec,
    
    Mean_RD = Mean_RD,
    
    Nb_itt = nrow_meth,
    
    Coverage_rate = Coverage_rate*100,
    
    MC_CR = MC_CR,
    
    CI_width = CI_width,
    
    P_025_CIW = P_025_CIW,
    P_975_CIW = P_975_CIW,
    
    MC_CIW = MC_CIW,
    
    Bias_Mean = Bias_Mean,
    
    MC_bias = MC_bias,
    
    Bias_relative_Mean = Bias_relative_Mean,
    
    MC_RelativeBias = MC_RelativeBias,
    
    Emp_SE = Emp_SE,
    
    ModSE = ModSE,
    
    RE_SE = RE_SE, # Relative Error of Standard Error
    
    MC_RESE = MC_RESE,
    
    Scenario = rep(n,n_metds),
    k = rep(as.numeric(unlist(Scenario_use[3])),n_metds),
    m = rep(as.numeric(unlist(Scenario_use[4])),n_metds),
    icc = rep(as.numeric(unlist(Scenario_use[5])),n_metds),
    nb_cov_indiv = rep(as.numeric(unlist(Scenario_use[6])),n_metds),
    nb_cov_clus = rep(as.numeric(unlist(Scenario_use[7])),n_metds),
    Prevalence = rep(paste(as.numeric(unlist(Scenario_use[1])),sep="/",as.numeric(unlist(Scenario_use[2]))),n_metds),
    Convergence_rate = nrow_meth/nb_itt_scen*100 # Care ! 1000 -> Number of iterations of our simulation study for each scenario
    
  )
  
  return(tab_res)
}


# Resum Results ----


Base_file = ""
Resu_file = paste(Base_file,"/Res_clust_analyse",sep = "")
Resume_output_conv_file = paste(Resu_file,"/Resume_output_conv_file",sep = "")
dir.create(Resume_output_conv_file)

methods = c("Unajusted",
            "TSP_gauss","TSP_bin",
            "GC_indv","GC_both",
            "TMLE_indv","TMLE_both",
            "GC_indv_indep","GC_both_indep",
            "TMLE_APS_Balzer","TMLE_APS_indv_Balzer")


Res_tot = c()
Nb_scen = 216
Scen = rep(1:Nb_scen)

for (n in Scen) {
  
  
  
  Res_n <- Resum_res_fun_it_spec(n,Resu_file,methods,nb_itt_scen = 1000,itt_spec = 1:1000)
  
  Res_tot = rbind(Res_tot,Res_n)
  
  ### EXcels data resume  ----
  
  write.csv2(Res_conv_n,here::here(Resume_output_conv_file,paste("/Data_output_Scenario_",sep = "",n,".csv")),row.names = FALSE)
  
}


write.csv2(Res_conv_tot,here::here(paste(Resu_file,"/Res_all_scen_itt_conv.csv",sep = "")),row.names = FALSE)


