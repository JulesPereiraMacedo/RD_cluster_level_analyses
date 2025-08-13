# Package functions from Balzer's package ----

TMLE.Balzer <- function(data,Scenario,itt_para,target = "indv",type_cov_adj = "both"){
  
  True_RD = as.numeric(unlist(Scenario[1])) - as.numeric(unlist(Scenario[2]))
  
  data_Balzer = data
  names(data_Balzer)[c(1,2,3)] <- c('Y','A','id')
  
  if(target == "clust"){
    N_j <- as.vector(unlist((data %>% group_by(cluster) %>% summarise(n_j = length(Outcome)))[,2]))
    data_Balzer$alpha = rep(1/N_j,N_j)
  }
  
  if(type_cov_adj=="both"){
    cov_adj =  c("IC_1","IC_2","IC_3","CC_1","CC_2")
  }
  if(type_cov_adj=="indv"){
    cov_adj =  c("IC_1","IC_2","IC_3")
  }
  
  
  Balzer_res = Stage2(goal = 'RD', target = target, data.input=data_Balzer,
                      QAdj =  cov_adj, Qform = 'glm',
                      gAdj =  cov_adj, gform = 'glm',
                      do.data.adapt = F, data.adapt.complexity = NULL,
                      cand.QAdj = NULL, cand.Qform = 'glm',
                      cand.gAdj = NULL, cand.gform = 'glm',
                      V = 5, remove.pscore = NULL, do.cv.variance = F,
                      sample.effect = F,
                      break.match = T, one.sided = F, alt.smaller = NULL,
                      verbose = F, psi = True_RD,
                      return.IC = T)
  
  
  Balzer_res$output$est.df$se
  Balzer_res$output$est.df$bias
  # Table of result
  
  Tab_res = data.frame(RD_CL_TMLE_Balzer = Balzer_res$output$est.df$est,
                       SE_CL_TMLE_Balzer = Balzer_res$output$est.df$se,
                       LL_95_CL_TMLE_JAPM = Balzer_res$output$est.df$CI.lo,
                       UL_95_CL_TMLE_JAPM = Balzer_res$output$est.df$CI.hi,
                       Bias_JAPM = Balzer_res$output$est.df$bias,
                       Relative_Bias_JAPM = Balzer_res$output$est.df$bias/Balzer_res$output$est.df$psi*100,
                       Itt = itt_para)
  return(Tab_res)
}

TMLE.Balzer.APS <- function(data,Scenario,itt_para,target = "indv"){
  
  True_RD = as.numeric(unlist(Scenario[1])) - as.numeric(unlist(Scenario[2]))
  
  data_Balzer = data
  names(data_Balzer)[c(1,2,3)] <- c('Y','A','id')
  
  if(target == "clust"){
    N_j <- as.vector(unlist((data %>% group_by(cluster) %>% summarise(n_j = length(Outcome)))[,2]))
    data_Balzer$alpha = rep(1/N_j,N_j)
  }
  
  
  
  
  Balzer_res = Stage2(goal = 'RD', target = target, data.input=data_Balzer,
                      QAdj =  NULL, Qform = NULL,
                      gAdj =  NULL, gform = NULL,
                      do.data.adapt = T, data.adapt.complexity = "high",
                      cand.QAdj = c("IC_1","IC_2","IC_3","CC_1","CC_2"), cand.Qform = c('glm',"lasso","step"),
                      cand.gAdj = c("IC_1","IC_2","IC_3","CC_1","CC_2"), cand.gform = c('glm',"lasso","step"),
                      V = 5, remove.pscore = T, do.cv.variance = T,
                      sample.effect = F,
                      break.match = T, one.sided = F, alt.smaller = NULL,
                      verbose = F, psi = True_RD,
                      return.IC = T)
  
  
  Balzer_res$output$est.df$se
  Balzer_res$output$est.df$bias
  # Table of result
  
  Tab_res = data.frame(RD_CL_TMLE_Balzer = Balzer_res$output$est.df$est,
                       SE_CL_TMLE_Balzer = Balzer_res$output$est.df$se,
                       LL_95_CL_TMLE_JAPM = Balzer_res$output$est.df$CI.lo,
                       UL_95_CL_TMLE_JAPM = Balzer_res$output$est.df$CI.hi,
                       Bias_JAPM = Balzer_res$output$est.df$bias,
                       Relative_Bias_JAPM = Balzer_res$output$est.df$bias/Balzer_res$output$est.df$psi*100,
                       Itt = itt_para)
  return(Tab_res)
}

TMLE.Balzer.APS.indv <- function(data,Scenario,itt_para,target = "indv"){
  
  True_RD = as.numeric(unlist(Scenario[1])) - as.numeric(unlist(Scenario[2]))
  
  data_Balzer = data
  names(data_Balzer)[c(1,2,3)] <- c('Y','A','id')
  
  if(target == "clust"){
    N_j <- as.vector(unlist((data %>% group_by(cluster) %>% summarise(n_j = length(Outcome)))[,2]))
    data_Balzer$alpha = rep(1/N_j,N_j)
  }
  
  
  
  
  Balzer_res = Stage2(goal = 'RD', target = target, data.input=data_Balzer,
                      QAdj =  NULL, Qform = NULL,
                      gAdj =  NULL, gform = NULL,
                      do.data.adapt = T, data.adapt.complexity = "high",
                      cand.QAdj = c("IC_1","IC_2","IC_3"), cand.Qform = c('glm',"lasso","step"),
                      cand.gAdj = c("IC_1","IC_2","IC_3"), cand.gform = c('glm',"lasso","step"),
                      V = 5, remove.pscore = T, do.cv.variance = T,
                      sample.effect = F,
                      break.match = T, one.sided = F, alt.smaller = NULL,
                      verbose = F, psi = True_RD,
                      return.IC = T)
  
  
  Balzer_res$output$est.df$se
  Balzer_res$output$est.df$bias
  # Table of result
  
  Tab_res = data.frame(RD_CL_TMLE_Balzer = Balzer_res$output$est.df$est,
                       SE_CL_TMLE_Balzer = Balzer_res$output$est.df$se,
                       LL_95_CL_TMLE_JAPM = Balzer_res$output$est.df$CI.lo,
                       UL_95_CL_TMLE_JAPM = Balzer_res$output$est.df$CI.hi,
                       Bias_JAPM = Balzer_res$output$est.df$bias,
                       Relative_Bias_JAPM = Balzer_res$output$est.df$bias/Balzer_res$output$est.df$psi*100,
                       Itt = itt_para)
  return(Tab_res)
}


fun_para_analyse_cluster_TMLE_Balzer <- function(itt_para,Scenario,n,Data_file,Resu_file){
  
  # Recuperation Data
  data = read.csv(paste(paste(Data_file,"/Scenario_",n,sep=""),sep ="",paste("/Data_itt_",sep = "",itt_para,"_scen_",n,".csv")), sep=";")
  
  ## Analyse all method
  
  res_CL_TMLE_Balzer_indv       = TMLE.Balzer(data,Scenario = Scenario,itt_para = itt_para,target = "clust",type_cov_adj = "indv")
  
  res_CL_TMLE_Balzer_both       = TMLE.Balzer(data,Scenario = Scenario,itt_para = itt_para,target = "clust",type_cov_adj = "both")
  
  # Remplissage CSV resultat pour chacune des méthodes
  
  write.table(res_CL_TMLE_Balzer_indv,file = paste(Resu_file,"/TMLE_analyse_Balzer_indvcov/Data_output_CL_TMLE_Balzer_indvcov_Scenario_",n,".csv",sep = ""),
              append = TRUE,row.names = FALSE,col.names = FALSE,sep = ";",dec = ",")
  
  write.table(res_CL_TMLE_Balzer_both,file = paste(Resu_file,"/TMLE_analyse_Balzer_bothcov/Data_output_CL_TMLE_Balzer_bothcov_Scenario_",n,".csv",sep = ""),
              append = TRUE,row.names = FALSE,col.names = FALSE,sep = ";",dec = ",")
  
}


fun_para_analyse_cluster_TMLE_Balzer_APS <- function(itt_para,Scenario,n,Data_file,Resu_file){
  
  # Recuperation Data
  data = read.csv(paste(paste(Data_file,"/Scenario_",n,sep=""),sep ="",paste("/Data_itt_",sep = "",itt_para,"_scen_",n,".csv")), sep=";")
  
  ## Analyse all method
  
  res_CL_TMLE_Balzer_APS       = TMLE.Balzer.APS(data,Scenario = Scenario,itt_para = itt_para,target = "clust")
  
  # Remplissage CSV resultat pour chacune des méthodes
  
  write.table(res_CL_TMLE_Balzer_APS,file = paste(Resu_file,"/Data_output_CL_TMLE_Balzer_APS_Scenario_",n,".csv",sep = ""),
              append = TRUE,row.names = FALSE,col.names = FALSE,sep = ";",dec = ",")
  
}


fun_para_analyse_cluster_TMLE_Balzer_APS_indv <- function(itt_para,Scenario,n,Data_file,Resu_file){
  
  # Recuperation Data
  data = read.csv(paste(paste(Data_file,"/Scenario_",n,sep=""),sep ="",paste("/Data_itt_",sep = "",itt_para,"_scen_",n,".csv")), sep=";")
  
  ## Analyse all method
  
  res_CL_TMLE_Balzer_APS       = TMLE.Balzer.APS.indv(data,Scenario = Scenario,itt_para = itt_para,target = "clust")
  
  # Remplissage CSV resultat pour chacune des méthodes
  
  write.table(res_CL_TMLE_Balzer_APS,file = paste(Resu_file,"/Data_output_CL_TMLE_Balzer_APS_indv_Scenario_",n,".csv",sep = ""),
              append = TRUE,row.names = FALSE,col.names = FALSE,sep = ";",dec = ",")
  
}

# Package functions from Balzer -----




#' Implement adaptive pre-specification
#'
#' @description function to implement adaptive prespecification as described in
#'   Balzer et al. "Adaptive pre-specification in randomized trials with and
#'   without pair-matching"
#'
#' @param goal goal of analysis ('aRR', 'OR', 'RD')
#' @param target Target of estimation and inference: cluster-level ("clust") or
#'   pooled-indv effect ("indv")
#' @param break.match Logical indicating whether to break matched pairs (if they
#'   exist)
#' @param Ldata Input data set
#' @param V
#' @param cand.QAdj Character vector of names of candidate adjustment variables
#'   for the outcome regression
#' @param cand.Qform Character vector of names of candidate adjustment
#'   algorithms for the outcome regression
#' @param cand.gAdj Character vector of names of candidate adjustment variables
#'   for the propensity score
#' @param cand.gform Character vector of names of candidate adjustment
#'   algorithms for the propensity score
#' @param remove.pscore if T, remove the variable(s) selected for adjustment in
#'   the outcome regression from candidates for the pscore... should only be
#'   used if doing adaptive prespec in RCT with few indpt units
#' @param QAdj selected adjustment variable for the outcome regression, if
#'   specifying in advannce
#' @param gAdj selected adjustment variable for the pscore regression, if
#'   specifying in advance
#' @param scale_value maximum value for outcome scaling
#' @param scale_value_min minimum value for outcome scaling
#' @param verbose
#' @param sample.effect If the target of inference the sample effect or the
#'   population effect?
#' @param data.adapt.complexity Set to 'low', 'med', or 'high'. See
#'   ?get.cand.adj for details.
#'
#' @return Working model selection for candidate TMLE
#' @export
#' 
#' @examples
do.adaptive.prespec <- function(goal, target='indv', break.match=T, Ldata, V=5,
                                cand.QAdj, cand.Qform,
                                cand.gAdj, cand.gform,
                                remove.pscore, QAdj = NULL, gAdj = NULL,
                                scale_value, scale_value_min, verbose = F,
                                sample.effect, data.adapt.complexity){
  
  # get the indpt units (will be each observation in and individual-level RCT)
  if( !break.match ){
    Ldata$indpt.unit <- Ldata$pair
  } else {
    Ldata$indpt.unit <- Ldata$id
  }
  unique.ids <- unique(Ldata$indpt.unit)
  
  # get folds
  if(length(unique.ids) > 40){
    # V-fold CV
    folds <- get.folds(V=V, Y=Ldata$Y, ids=unique.ids)
  } else {
    # leave-one-out CV
    folds <- vector("list", length(unique.ids))
    for(v in seq(length(unique.ids))){
      folds[[v]] <- unique.ids[v]
    }
  }
  
  if(is.null(remove.pscore)){
    if(length(unique.ids) > 40){
      remove.pscore <- F
    } else {
      remove.pscore <- T
    }
    message(paste0("remove.pscore not speficied; setting to ", remove.pscore, " based on # of independent units"))
  }
  
  
  #=====================================================
  # GENERATE COVARIATE SET / WORKING MODEL COMBINATIONS FOR OUTCOME REGRESSION
  # Figure out complexity level first
  if(is.null(data.adapt.complexity)){
    nvars <- length(union(cand.gAdj, cand.QAdj))
    ratio <- length(unique.ids) / nvars
    low_cutoff <- 10
    med_cutoff <- 50
    if(ratio < low_cutoff){
      message(paste0("Fewer than ", low_cutoff, " independent units per covariate. Setting data.adapt.complexity = 'low'"))
      data.adapt.complexity <- "low"
    } else if (ratio < med_cutoff){
      message(paste0("Between ", low_cutoff, " and ", med_cutoff, " independent units per covariate. Setting data.adapt.complexity = 'med'"))
      data.adapt.complexity <- "med"
    } else {
      message(paste0("More than ", med_cutoff, " independent units per covariate. Setting data.adapt.complexity = 'high'"))
      data.adapt.complexity <- "high"
    }
  }
  cand.Q <- get.cand.adj(cand.vars = cand.QAdj, cand.algos = cand.Qform,
                         data.adapt.complexity = data.adapt.complexity)
  cand.QAdj <- cand.Q$cand.varset
  cand.Qform <- cand.Q$cand.algos
  
  
  
  
  #=====================================================
  # FIGURE OUT WHICH COMBINATION FOR OUTCOME REGRESSION IS BEST
  if( is.null(QAdj) ){
    
    if(verbose) print("Examining sets for outcome regression")
    
    # do adaptive pre-specification to select from candidate approaches for Qbar
    select.Q <- suppressWarnings(CV.selector(goal=goal, target=target, break.match=break.match, Ldata=Ldata,
                                             CAND.ADJ = cand.QAdj, CAND.FORM=cand.Qform, forQ=T, verbose = verbose,
                                             scale_value=scale_value, scale_value_min=scale_value_min,
                                             folds=folds, sample.effect = sample.effect
    ))
    #if(verbose) print(select.Q)
    Q.index <- select.Q$adj.index
    QAdj <- select.Q$Adj
    Qform <- select.Q$form
    #if(verbose){print("Selected for Q:") print(QAdj); print(Qform)}
    
    # if select unadjusted estimator for QbarAW=E(Y|A,W), then stop
    if(sum('U' %in% QAdj) == 1){ 
      g.index <- -99; gAdj <- 'U'; gform <- 'glm'
      var.CV <- select.Q$var.CV
      var.CV.1 <- select.Q$var.CV.1		
      var.CV.0 <- select.Q$var.CV.0	
    }
    
    if((sum(QAdj %in% 'U')) == 0 & remove.pscore){
      if(verbose) print('remove.pscore = T; removing selected QAdj variable(s) from candidates for gAdj')
      cand.gAdj <- cand.gAdj[cand.gAdj %nin% QAdj]
      if(length(cand.gAdj) == 0){
        if(verbose) print("all covars used for outcome regression. setting cand.gAdj = U")
        cand.gAdj <- "U"
      }
    }
  }
  
  
  #=====================================================
  # GENERATE COVARIATE SET / WORKING MODEL COMBINATIONS FOR PROPENSITY SCORE
  cand.g <- get.cand.adj(cand.vars = cand.gAdj, cand.algos = cand.gform,
                         data.adapt.complexity = data.adapt.complexity)
  cand.gAdj <- cand.g$cand.varset
  cand.gform <- cand.g$cand.algos
  
  
  
  #=====================================================
  # FIGURE OUT WHICH COMBINATION FOR PROPENSITY SCORE IS BEST
  if( is.null(gAdj) ){ 		
    if(verbose) print("Examining sets for propensity score")
    
    select.G <- suppressWarnings(CV.selector(goal = goal, target = target,
                                             break.match = break.match, Ldata = Ldata,
                                             CAND.ADJ = cand.gAdj, CAND.FORM = cand.gform,
                                             forQ = F,  verbose = verbose,
                                             # input selected variables/form of the outcome regression
                                             QAdj = QAdj, Qform = Qform,
                                             scale_value = scale_value, scale_value_min = scale_value_min,
                                             folds = folds,
                                             sample.effect = sample.effect))
    
    g.index <- select.G$adj.index
    gAdj <- select.G$Adj
    gform <- select.G$form
    var.CV <- select.G$var.CV		
    var.CV.1 <- select.G$var.CV.1		
    var.CV.0 <- select.G$var.CV.0		
    #if(verbose){print(select.G);    print("Selected for g:");  print(gAdj);  print(gform)}
  }		
  output_aps <- list(Q.index = Q.index, QAdj = QAdj, Qform = Qform, 
                     g.index = g.index, gAdj = gAdj, gform = gform,
                     var.CV = var.CV, var.CV.1 = var.CV.1, var.CV.0 = var.CV.0)
  #if(verbose) print(output_aps)
  return(output_aps)
}



#-----------------------------------------------------#-----------------------------------------------------
# CV.selector: function to estimate the cross-validated risk
#		Loss function is the squared-IC; Risk is then the variance of the TMLE
#		See Balzer et al. "Adaptive pre-specification in randomized trials with and without pair-matching"
#
#	input: goal of analysis ('aRR' or 'RD), 
#   target of inference 
#		indicator to break the match (break.match)
#		dataset (Ldata)
#		candidate adjustment variables; they do not have to be at the cluster-level
#		indicator if for the conditional mean outcome (forQ)
#		selected adjustment variable for the outcome regression (QAdj)
#	output: selection for adjustment variable (corresponding to a TMLE)
#-----------------------------------------------------#-----------------------------------------------------


CV.selector <- function(goal, target, break.match, Ldata, CAND.ADJ, CAND.FORM, 
                        forQ, QAdj=NULL, Qform=NULL, verbose,
                        scale_value, scale_value_min, folds, sample.effect){
  # if(length(CAND.FORM) == 1){
  #   # if exploring only one estimation algorithm (usually GLM) then need to replicate the number forms
  #   CAND.FORM <- rep(CAND.FORM, length(CAND.ADJ))
  # }
  
  if(length(CAND.FORM) != length(CAND.ADJ)){ # After changes to get.cand.adj, this should not happen.
    stop('PROBLEM: MISMATCH LENGTHS OF ADJ VAR SETS AND MODELS')
  }
  
  # Number of candidate estimators is given by length Qform//gform
  num.tmles <- length(CAND.FORM)
  CV.risk <- var.CV <- var.CV.1 <- var.CV.0 <- rep(NA, num.tmles)
  
  for(k in 1: num.tmles){	
    
    if(forQ){
      # print("evaluating a Q option:")
      # print(CAND.ADJ[[k]])
      # print(CAND.FORM[k])
      # if selecting the adjustment approach for the outcome regression
      IC.temp <- get.IC.CV(goal = goal, target = target, break.match = break.match, Ldata = Ldata,
                           QAdj = CAND.ADJ[[k]], Qform = CAND.FORM[k],
                           gAdj = NULL, gform = 'glm', verbose = verbose,
                           scale_value = scale_value, scale_value_min = scale_value_min, 
                           folds = folds, sample.effect = sample.effect)
    } else {
      # print("evaluating g option:")
      # print(CAND.ADJ[[k]])
      # print(CAND.FORM[k])
      # if collaboratively selecting the adjustment approach for the pscore
      IC.temp <- get.IC.CV(goal = goal, target = target, break.match = break.match, Ldata = Ldata, 
                           QAdj = QAdj, Qform = Qform, verbose = verbose,
                           gAdj = CAND.ADJ[[k]], gform = CAND.FORM[k],
                           scale_value = scale_value, scale_value_min = scale_value_min, 
                           folds = folds, sample.effect = sample.effect)
    }
    # if(verbose) print(IC.temp)
    # estimating the CV risk for each candidate
    CV.risk[k]<- IC.temp$CV.risk
    # estimating the CV variance for that TMLE
    var.CV[k] <- IC.temp$var.CV
    var.CV.1[k] <- IC.temp$var.CV.1
    var.CV.0[k] <- IC.temp$var.CV.0
  }
  #print(CV.risk)
  
  # select the candidate estimator resulting in the smallest CV-risk
  adj.index <- which.min(CV.risk)
  output_best <- list(CV.risk = CV.risk[adj.index],
                      adj.index = adj.index, Adj = CAND.ADJ[[adj.index]], form = CAND.FORM[adj.index],
                      var.CV = var.CV[adj.index],
                      var.CV.1 = var.CV.1[adj.index],
                      var.CV.0 = var.CV.0[adj.index])
  return(output_best)
}




#' Calculate a cross-validated estimate of the influence curve 
#'
#' UPDATES
#' previous version only did leave-one-out (unit or pair)
#' this version generalizes to V-fold CV if V>=40
#' - can input the number of folds V (default=10)
#' - folds created stratified on binary outcomes (by default)
#' - if stratify=T and # observations in a given class is <V, 
#' then sets V=min observations in that fold
#'
#' See Balzer et al. "Adaptive pre-specification in randomized trials with and without pair-matching."
#'
#' @param goal goal of analysis ('aRR' or 'RD)
#' @param target cluster/individual effect (JN: Not implemented yet?)
#' @param break.match indicator to break the match
#' @param Ldata input dataset
#' @param QAdj adjustment variable for the outcome regression
#' @param Qform adjustment approach for outcome regression
#' @param gAdj adjustment variable for the pscore
#' @param gform adjustment approach for pscore regression
#' @param scale_value For outcomes not bounded within [0,1], the maximum value of the outcome.
#' @param scale_value_min For outcomes not bounded within [0,1], the minimum value of the outcome.
#' @param folds Number of folds in cross-validation.
#'
#' @return cross-validated estimate of the IC for pair
#' @export
#'
#' @examples
get.IC.CV <- function(goal, target, break.match, verbose,
                      Ldata, QAdj, Qform, gAdj=NULL, gform='glm', 
                      scale_value, scale_value_min, folds, sample.effect){
  
  nFolds <- length(folds)
  DY.CV <- CV.risk <- DY1.CV <- DY0.CV <- NULL
  
  # doing a cross-validated estimate
  for(i in 1:nFolds) {
    
    these <- Ldata$indpt.unit %in% folds[[i]]  ########  IMPORTANT!!!!
    valid <- Ldata[these, ]
    train <- Ldata[!these,]
    
    # run full TMLE algorithm on the training set
    train.out <- do.TMLE(goal=goal, target=target, train=train, QAdj=QAdj, Qform=Qform, 
                         gAdj=gAdj, gform=gform,
                         scale_value=scale_value, scale_value_min=scale_value_min,
                         doing.CV=T, verbose=F, sample.effect = sample.effect)	
    
    # get the relevant components of the IC for the validation set, 
    # using fits based on the training set
    valid.out <- do.TMLE.validset(goal=goal, target=target, valid=valid, train.out=train.out,
                                  scale_value=scale_value, scale_value_min=scale_value_min,
                                  sample.effect = sample.effect)	
    
    # estimating the CV risk for each candidate
    # risk = Expectation of loss with loss as IC-sq
    # risk = variance of TMLE
    if(break.match){
      DY.CV <- c(DY.CV, valid.out$DY)
      CV.risk <- c(CV.risk, mean(valid.out$DY^2))
    } else {
      DY.CV <- c(DY.CV, valid.out$DY.paired)
      CV.risk <- c(CV.risk, mean(valid.out$DY.paired^2))
    }
    DY1.CV <- c(DY1.CV, valid.out$DY1)
    DY0.CV <- c(DY0.CV, valid.out$DY0)
  }
  #if(verbose){print(valid.out)}
  
  # average across folds
  CV.risk <- mean(CV.risk)
  
  # estimating the CV variance for that TMLE
  var.CV <- stats::var(DY.CV) / length(DY.CV)
  var.CV.1 <- stats::var(DY1.CV) / length(DY1.CV)
  var.CV.0 <- stats::var(DY0.CV) / length(DY0.CV)
  
  return(list(CV.risk = CV.risk, var.CV = var.CV,
              var.CV.0 = var.CV.0, var.CV.1 = var.CV.1))
}



#-----------------------------------------------------#-----------------------------------------------------
# do.TMLE.for.valid: function to obtain a cross-validated estimate of the influence curve
#	for observations in the validation set
#		See Balzer et al. "Adaptive pre-specification in randomized trials with and without pair-matching"
#
#	input: goal of analysis ('aRR' or 'RD'),
#		validation dataset ('valid') 
#		TMLE-fits from training set (train.out)
#	output: cross-validated estimate of the IC,
#		cross-validated risk estimate (loss=IC^2)
#-----------------------------------------------------#-----------------------------------------------------

do.TMLE.validset <- function(goal, target, valid,
                             train.out, scale_value, scale_value_min,
                             sample.effect){
  
  # J <- length(unique(valid$id) )
  
  #=============================================
  # Step1 - initial estimation of E(Y|A,W)= Qbar(A,W)
  #=============================================
  valid <- do.Init.Qbar(train=valid, QAdj=train.out$QAdj, Qform=train.out$Qform, 
                        glm.out=train.out$Q.out)$train
  
  #=============================================
  # Step2: Calculate the clever covariate
  #=============================================
  valid <- get.clever.cov(train=valid, gAdj=train.out$gAdj, gform=train.out$gform, 
                          p.out=train.out$p.out)$train
  
  #=============================================
  # Step3: Targeting - 			
  #=============================================
  valid <- do.targeting(train=valid, eps=train.out$eps, goal=goal)
  
  #=============================================
  # Step5: Variance estimation using treatment-specific means from training set
  #=============================================
  get.IC.variance(goal=goal, target=target, Vdata=valid, R1=train.out$R1, R0=train.out$R0, 
                  scale_value=scale_value, scale_value_min=scale_value_min,
                  doing.CV=T, sample.effect = sample.effect)
}








#-----------------------------------------------------#-----------------------------------------------------

# adapted from .cvFolds from cvAUC package: https://CRAN.R-project.org/package=cvAUC
# by Erin LeDell 
# **** WARNING - stratify=T option is currently broken for cluster randomizzed trials! *****
get.folds <- function(V, Y, ids, stratify=F){
  
  if(stratify & length(unique(Y))==2){
    # stratify on the outcome
    classes <- tapply(1:length(Y), INDEX=Y, FUN=split, 1)
    ids.Y1 <- ids[classes$`1`]
    ids.noY1 <- ids[classes$`0`]
    if(length(ids.Y1) < V | length(ids.noY1) < V) {
      V <- min( length(ids.Y1), length(ids.noY1))
    }
    ids.Y1.split <- split(sample(length(ids.Y1)), rep(1:V, length=length(ids.Y1)))
    ids.noY1.split <- split(sample(length(ids.noY1)), rep(1:V, length=length(ids.noY1)))
    folds <- vector("list", V)
    for (v in seq(V)){
      folds[[v]] <- c(ids.Y1[ids.Y1.split[[v]]], ids.noY1[ids.noY1.split[[v]]])
    }
    
  } else {
    # dont stratify on the outcome
    ids.split <- split(sample(length(ids)), rep(1:V, length=length(ids)))
    folds <- vector("list", V)
    for (v in seq(V)){
      folds[[v]] <- ids[ids.split[[v]]]
    }
    
  }
  folds
}


#### data_generators.R ####


simulate_clustered_data <- function(treatment_arm_cluster_size = 15, treatment_arm_clusters = 50,
                                    control_arm_cluster_size = 5, control_arm_clusters = 50,
                                    txt_eff = .5, ranef_sd = 1, covar1_coef = 1, seed = NA,
                                    informative_cluster_size = F){
  if(!is.na(seed)){ set.seed(seed) }
  total_units_txt <- treatment_arm_cluster_size*treatment_arm_clusters
  total_units_ctr <- control_arm_cluster_size*control_arm_clusters
  total_units <- total_units_txt + total_units_ctr
  unit_id <- 1:(total_units)
  cluster_id_txt <- rep(1:treatment_arm_clusters, each = treatment_arm_cluster_size)
  cluster_id_ctr <- rep((max(cluster_id_txt)+1):(max(cluster_id_txt) + control_arm_clusters), each = control_arm_cluster_size)
  cluster_id <- c(cluster_id_txt, cluster_id_ctr)
  total_clusters <- max(cluster_id)
  cluster_lengths <- c(rep(treatment_arm_cluster_size, times = treatment_arm_clusters),
                       rep(control_arm_cluster_size, times = control_arm_clusters))
  cluster_cov <- sd(cluster_lengths) / mean(cluster_lengths)
  N_j <- rep(cluster_lengths, times = cluster_lengths)
  
  A <- c(rep(1, times = total_units_txt), rep(0, times = total_units_ctr))
  
  if(ranef_sd == 0){
    UE <- 0
  } else {
    UEs <- exp(rnorm(n = total_clusters, mean = 0, sd = ranef_sd))
    UE <- rep(UEs, times = cluster_lengths)
  }
  UY <- runif(n = total_units, min = 0, max = 2)
  covar0 <- rnorm(n = total_units, mean = UE, sd = .5)
  covar1 <- rbinom(n = total_units, size = 1, prob = .2)
  covar2 <- rnorm(n = total_units, mean = 0, sd = .5)
  covar3 <- rnorm(n = total_units, mean = 0, sd = .5)
  
  Y0 <- covar1_coef*covar1 + covar0^2 + covar0*covar1 + UY + UE
  Y1 <- covar1_coef*covar1 + covar0^2 + covar0*covar1 + UY + UE + txt_eff
  if(informative_cluster_size){
    Y1 <- Y1 + txt_eff*N_j/total_units
  }
  Y <- ifelse(A == 1, Y1, Y0)
  dat <- cbind.data.frame(unit_id, cluster_id, A,
                          UY, covar0, covar1, covar2, covar3, UE,
                          txt_eff, Y0, Y1, Y,
                          N_j = N_j) %>% dplyr::mutate(id = cluster_id, alpha = 1, U = 1)
  cluster_means <- dat %>% dplyr::group_by(cluster_id, A) %>% dplyr::summarise(cluster_mean0 = mean(Y0),
                                                                               cluster_mean1 = mean(Y1),
                                                                               cluster_meanY = mean(Y))
  ind_weight <- total_clusters / total_units
  ind_weighted_cluster_means <- dat %>% dplyr::group_by(cluster_id, A) %>%
    dplyr::summarise(ind_weighted_cluster_mean0 = sum(Y0)*ind_weight,
                     ind_weighted_cluster_mean1 = sum(Y1)*ind_weight,
                     ind_weighted_cluster_meanY = sum(Y)*ind_weight)
  cluster_mean_diff <- mean(cluster_means$cluster_mean1) - mean(cluster_means$cluster_mean0)
  cluster_mean_irr <- mean(cluster_means$cluster_mean1) / mean(cluster_means$cluster_mean0)
  ind_weighted_cluster_mean_diff <- mean(ind_weighted_cluster_means$ind_weighted_cluster_mean1) - mean(ind_weighted_cluster_means$ind_weighted_cluster_mean0)
  ind_weighted_cluster_mean_irr <- mean(ind_weighted_cluster_means$ind_weighted_cluster_mean1) / mean(ind_weighted_cluster_means$ind_weighted_cluster_mean0)
  
  return(cbind.data.frame(dat,
                          covar1_coef = covar1_coef,
                          true_rd_input = txt_eff,
                          true_covar_coef_input = covar1_coef,
                          ranef_sd_input = ranef_sd,
                          cluster_cov = cluster_cov,
                          sample_mean0 = mean(Y0),
                          sample_mean1 = mean(Y1),
                          unit_mean_diff = mean(Y1) - mean(Y0),
                          unit_irr = mean(Y1) / mean(Y0),
                          cluster_mean_diff = cluster_mean_diff,
                          cluster_mean_irr = cluster_mean_irr,
                          ind_weighted_cluster_mean_diff = ind_weighted_cluster_mean_diff,
                          ind_weighted_cluster_mean_irr = ind_weighted_cluster_mean_irr
  ))
}

#### stage2_functions.R ####

#' Two-stage or single-stage RCT estimation
#'
#' Takes data frame input.data and calculates arm-specific and treatment effect estimates.
#' 
#' Outcomes can be bounded continuous (in which case they are scaled to [0,1]
#' before estimation and inference, then unscaled - See Ch7 TLB) or binary. Does
#' not work for multicategory outcomes.
#' 
#' The input.data object must include treatment indicator column 'A' and outcome
#' 'Y'. If the trial pair-matched and want to keep pairs (break.match = F), the
#' column with pair IDs must be labeled as 'pair'.
#' 
#' For cluster randomized trials, should also include weights 'alpha', and cluster id 'id'. 
#' 
#' column for weights (alpha)
#'   Set = 1 for individually randomized trials.
#'   For cluster randomized trials: 
#'     value of the weights depends on the target of inference and data level 
#'     Details in Benitez et al. https://arxiv.org/abs/2110.09633v2
#'     let J=number of clusters, N_j = cluster-specific sample size, N_tot = total # participants= sum_j N_j
#'       if target='clust' with cluster-level data, alpha=1 
#'       if target='clust' with indv-level data, alpha = 1/N_j
#'       if target='indv' with cluster-level data, alpha = J/N_tot*N_j
#'       if target='indv' with indv-level data, then alpha = 1
#'     for demonstration, see sim2.R in https://github.com/LauraBalzer/Comparing_CRT_Methods
#'   weights must sum to the total # of randomized units
#' 
#' For testing, or for a prespecified (not-adaptive) estimation approach,
#' specify conditional mean outcome adjustment variables (QAdj) and method
#' (Qform), and propensity score adjustment variables (gAdj) and method (gform).
#' Note that this is NOT recommended in general--- instead use Adaptive
#' Prespecification (do.data.adapt = T).
#'
#'
#' @param goal aRR = arithmetic risk ratio;  RD=risk difference; OR= odds ratio
#'   (not recommended)
#' @param target target of inference: cluster-level ("clust") or pooled-indv
#'   effect ("indv"). With appropriate weights, this can recover cluster-level
#'   effects from individual data, or individual-level effects from clustered
#'   data (see note above re: alpha weights).
#' @param data.input Observed data with treatment column A, outcome Y. See
#'   further specifications above.
#' @param QAdj Character string of variable names to force in to the adjustment
#'   set for the outcome regression. Default is NULL. If non-null, will override
#'   'cand.QAdj'. For an unadjusted outcome regression, specify QAdj = 'U'.
#' @param Qform String pre-specifying which algorithm to use to estimate the
#'   outcome regression. Currently implemented options: 'glm', 'lasso', 'mars'
#'   (multivariate adaptive regression splines), 'mars.corP' (MARS that
#'   pre-screens in variables correlated with the outcome), 'step' (stepwise
#'   selection with main terms), and 'step.interaction' (stepwise selection with
#'   interaction terms).
#' @param gAdj Character string of variable names to force in to the adjustment
#'   set for the propensity score. Default is NULL. If non-null, will override
#'   'cand.gAdj'. For an unadjusted outcome regression, specify gAdj = 'U'.
#' @param gform String pre-specifying which algorithm to use to estimate the
#'   propensity score. Same options as 'Qform'.
#' @param do.data.adapt Should adaptive pre-specification (Balzer 2016) be used?
#'   Defaults to FALSE.
#' @param data.adapt.complexity
#' @param cand.QAdj Character vector of candidate adjustment variable(s) for
#'   conditional mean outcome
#' @param cand.Qform Character vector with candidate adjustment approache(s) for
#'   conditional mean outcome. See 'Qform' for options.
#' @param cand.gAdj Character vector of candidate adjustment variable(s) for
#'   propensity score
#' @param cand.gform Character vector with candidate adjustment approache(s) for
#'   propensity score. See 'Qform' for options.
#' @param V Number of folds for cross-validation steps. Default is 5, once the
#'   number of independent units exceeds 40. For fewer than 40 independent
#'   units, defaults to LOOCV.
#' @param remove.pscore Relevant when do.data.adapt = T. Should adjustment
#'   variables that are selected for the outcome regression be removed from the
#'   candidate set for the propensity score? See 'Collaborative TMLE'. Recommend
#'   remove.pscore = TRUE when there are few independent units; defaults to TRUE
#'   if the number of independent units is under 40.
#' @param do.cv.variance Should cross-validated variance estimate be used for
#'   inference? Default FALSE.
#' @param sample.effect Is the target the sample effect among the observed
#'   units? If FALSE, target is population effect. Default is TRUE.
#' @param break.match If pair-matched trial, should pairs be ignored during
#'   inference?
#' @param one.sided Indicator specifying whether a one-sided p-value should be
#'   returned. Defaults to FALSE.
#' @param alt.smaller Only needed if one.sided = T. Specifies whether the
#'   alternative hypothesis is that the treatment arm has a lower level of the
#'   outcome. If the outcome is 'bad' (say, incidence of an infectious disease),
#'   and the treatment is hypothesized to be 'good' and reduce incidence,
#'   alt.smaller = T.
#' @param verbose Prints more information about what's happening under the
#'   hood.
#' @param psi In (say) a simulation study, if the true value of the treatment
#'   effect is known, it can be input here and the outcome will include
#'   an indicator for rejection of the null hypothesis and coverage. Default is
#'   NA.
#' @param return.IC indicator of whether to return the influence curve
#'
#' @return Point estimate and inference
#' @export
#'
#' @examples
Stage2 <- function(goal = 'aRR', target = 'indv', data.input, 
                   QAdj = NULL, Qform = 'glm',
                   gAdj = NULL, gform = 'glm',
                   do.data.adapt = F, data.adapt.complexity = NULL,
                   cand.QAdj = NULL, cand.Qform = 'glm',
                   cand.gAdj = NULL, cand.gform = 'glm',
                   V = 5, remove.pscore = NULL, do.cv.variance = F,
                   sample.effect = T,
                   break.match = T, one.sided = F, alt.smaller = NULL,
                   verbose = F, psi = NA,
                   return.IC = F){
  
  #=====================================================
  # INPUT CHECKS
  
  # if doing a one-sided test, need to specify the alternative
  # alt.smaller = T if intervention reduces mean outcome
  # alt.smaller = F if intervention increases mean outcome
  if(one.sided & is.null(alt.smaller)){
    stop('For one-sided test, need to specify the direction of the alternative hypothesis')
  }
  if(goal %nin% c("aRR", "RD", "OR")){
    stop("Please specify goal = 'aRR' for risk/rate ratio, 'RD' for risk difference, or 'OR' for odds ratio")
  }
  if(!is.null(QAdj) & !is.null(cand.QAdj) & do.data.adapt){
    stop("You have specified both forced-in adjustment variables (QAdj) for the outcome regression, as well as a candidate set to choose from (cand.QAdj). Please specify one only and leave the other as NULL.")
  }
  if(!is.null(gAdj) & !is.null(cand.gAdj) & do.data.adapt){
    stop("You have specified both forced-in adjustment variables (gAdj) for the propensity score, as well as a candidate set to choose from (cand.gAdj). Please specify one only and leave the other as NULL.")
  }
  # Make sure data has correct structure  
  data.input <- preprocess_and_check_data(data.input, verbose = verbose)
  
  
  
  
  #=====================================================
  # TRANSFORM the outcome as in Chpt7 of TLB if not bounded in [0,1]
  if(max(data.input[,'Y']) > 1){
    scale_value <- max(data.input[,'Y'])
  } else {
    scale_value <- 1
  }
  if(min(data.input[,'Y']) < 0){
    scale_value_min <- min(data.input[,'Y'])
  } else {
    scale_value_min <- 0
  }
  data.input[,'Y'] <- (data.input[,'Y'] - scale_value_min) / (scale_value - scale_value_min)
  
  
  
  
  #=====================================================
  # ADAPTIVE PRESPECIFICATION
  # update: flexibility in CV-scheme and candidate prediction algorithms
  if(do.data.adapt){
    select <- do.adaptive.prespec(goal = goal, target = target, break.match = break.match, 
                                  Ldata = data.input, V = V,
                                  cand.QAdj = cand.QAdj, cand.Qform = cand.Qform,
                                  cand.gAdj = cand.gAdj, cand.gform = cand.gform,
                                  remove.pscore = remove.pscore,
                                  QAdj = QAdj, gAdj = gAdj,
                                  scale_value = scale_value, scale_value_min = scale_value_min,
                                  verbose = verbose, sample.effect = sample.effect, data.adapt.complexity)
    
    Q.index <- select$Q.index
    QAdj <- select$QAdj
    Qform <- select$Qform
    g.index <- select$g.index
    gAdj <- select$gAdj	
    gform <- select$gform
    #if(verbose){print(paste(QAdj, gAdj))}
  } else {
    # QAdj <- gAdj <- 'U' # THIS WAS DUMB - OVERRODE ANY USER-SPECIFIED NON-NULL QAdj or gAdj
    Q.index <- g.index <- 1
  }
  
  # RUN FULL TMLE WITH ADJUSTMENT SET 
  # runs all code for point estimation on scaled outcome
  est <- do.TMLE(goal = goal, target = target, train =  data.input,
                 QAdj = QAdj, Qform = Qform, 
                 gAdj = gAdj, gform = gform,
                 scale_value = scale_value, scale_value_min = scale_value_min,
                 doing.CV = F, verbose = verbose, sample.effect = sample.effect)  
  
  n.clust <- length(unique(data.input$id)) 
  # Get point estimates of the treatment-specific mean
  R1 <- est$R1
  R0 <- est$R0
  
  
  # Now: for the intervention effect 
  #  the point estimate on the relevant scale for getting inference
  if(goal == 'aRR') {
    psi.hat <- log(R1/R0)
  } else if (goal == 'RD'){
    psi.hat <- R1 - R0
  } else if (goal == 'OR'){
    psi.hat <- log( R1/(1-R1) * (1-R0)/R0 )
  }
  
  if(break.match){
    # if breaking the match, set df to (#clusters - 2)
    df <- n.clust - 2
    var.hat <- est$var.break
  } else{
    # if preserving the match, set df to (#pairs-1)
    df <- length(unique(data.input$pair)) -1 
    var.hat <- est$var.pair
  }
  
  ############################## GET INFERENCE for arms and   
  if(do.cv.variance){
    #print(select)
    Txt.CV <- get.inference(psi.hat = R1, se = sqrt(select$var.CV.1), # need to change SE based on 'select' object
                            df = (n.clust-2))[, c('est','CI.lo','CI.hi','se')]
    Con.CV <- get.inference(psi.hat = R0, se = sqrt(select$var.CV.0),
                            df = (n.clust-2))[, c('est','CI.lo','CI.hi','se')]
  }
  
  Txt <- get.inference(psi.hat = R1, se = sqrt(est$var.R1), df = (n.clust-2))[, c('est','CI.lo','CI.hi','se')]
  Con <- get.inference(psi.hat = R0, se = sqrt(est$var.R0), df = (n.clust-2))[, c('est','CI.lo','CI.hi','se')]
  
  inference <- get.inference(goal = goal, psi=psi, psi.hat = psi.hat,
                             se = sqrt(var.hat), df = df,
                             one.sided = one.sided, alt.smaller = alt.smaller)
  
  Q_adj_and_form_final <- make_adj_var_df(vars = QAdj, Q_or_g = "Q")
  g_adj_and_form_final <- make_adj_var_df(vars = gAdj, Q_or_g = "g")
  
  if(do.cv.variance){
    # if getting cross-validated inference
    inference.CV <- get.inference(goal=goal, psi=psi, psi.hat=psi.hat,
                                  se=sqrt(select$var.CV), df=df,
                                  one.sided=one.sided, alt.smaller = alt.smaller)
    
    est.df <- data.frame(Txt=Txt, TxtCV = Txt.CV, Con=Con, ConCV=Con.CV,  psi=psi, inference, CV=inference.CV,
                         Qform = est$Qform, #Q_adj_and_form_final,
                         gform = est$gform)#, #g_adj_and_form_final)
  } else {
    est.df <- data.frame(Txt = Txt, Con = Con, psi = psi, inference, 
                         Qform = est$Qform, #Q_adj_and_form_final,
                         gform = est$gform)#, g_adj_and_form_final)
  }
  
  if(return.IC){
    Stage2_output <- list(IC=est, est.df=est.df)
  } else{
    Stage2_output <- est.df
  }
  # Making the return a little weirder for the purposes of easy sim studies?
  return(list(output = Stage2_output,
              QAdj = QAdj, gAdj = gAdj,
              Qvars = paste0("Adjustment variable(s) for outcome regression: ",
                             paste(QAdj, sep = "", collapse = ", ")),
              Qform = paste0("Adjustment algorithm for outcome regression: ",
                             paste(est$Qform, sep = "", collapse = ", ")),
              gvars = paste0("Adjustment variable(s) for propensity score: ",
                             paste(gAdj, sep = "", collapse = ", ")),
              gform = paste0("Adjustment algorithm for propensity score: ",
                             paste(est$gform, sep = "", collapse = ", "))))
}

#### tmle_functions.R ####



#' Obtain point estimates with TMLE
#'
#' @param goal 'aRR' for arithmetic risk ratio, 'OR' for odds ratio, otherwise
#'   risk difference
#' @param target 'clust' (cluster-level effect) or 'indv' (pooled
#' individual-level effect) target of inference; see Stage2 description
#' @param train Training data
#' @param QAdj Character string of variable names for the adjustment set for the
#'   outcome regression.
#' @param Qform Prediction function/algorithm for outcome regression, see Stage2
#'   and/or do.Init.Qbar for options
#' @param gAdj Character string of variable names for the adjustment set for the
#'   propensity score.
#' @param gform Prediction function/algorithm for propensity score, see Stage2
#'   and/or do.Init.Qbar for options
#' @param Q.out Initial estimator of the conditional mean outcome
#' @param p.out initial estimator of the propensity score
#' @param scale_value For bounded continuous outcomes, the max value in the
#'   input data
#' @param scale_value_min For bounded continuous outcomes, the min value in the
#'   input data
#' @param doing.CV 
#' @param verbose Prints more information about what's happening under the
#'   hood.
#' @param sample.effect Logical indicating whether to perform inference for the
#'   sample average treatment effect (if TRUE) or the population average
#'   treatment effect (if FALSE).
#'
#' @return A list with: training data augmented with estimates, prespecified
#'   adjustment variables for the conditional mean outcome (QAdj), prespecified
#'   adjustment variables for the propensity score (gAdj), initial estimator of
#'   the conditional mean outcome (Q.out), estimator of the propensity score
#'   (p.out), estimated fluctuation coef (epsilon)
#'  
#' @export
#'
#' @examples
do.TMLE <- function(goal, target, train,
                    QAdj, Qform ='glm',
                    gAdj = NULL, gform ='glm',
                    Q.out = NULL, p.out = NULL, 
                    scale_value, scale_value_min,
                    doing.CV = F, verbose, sample.effect) {	
  
  #=====================================================
  # Step1 - initial estimation of E(Y|A,W)= Qbar(A,W)
  #=====================================================
  
  # run glm on the adjustment set
  Q <- do.Init.Qbar(train=train, QAdj=QAdj, Qform=Qform, glm.out=Q.out, verbose=verbose)
  train <- Q$train
  
  #==========================================================
  # Step2: Calculate the clever covariate
  #==========================================================	
  
  G <- get.clever.cov(train=train, gAdj=gAdj, gform=gform, p.out=p.out, verbose=verbose)
  train <- G$train
  
  #==========================================================
  # Step3: Targeting
  #==========================================================
  
  eps <- get.epsilon(train=train, goal=goal, verbose=verbose)
  
  train <- do.targeting(train=train, eps=eps, goal=goal)
  
  #==========================================================
  # Step4: Parameter  estimation
  # will unscale if appropriate
  #==========================================================
  
  if(nrow(train) > length(unique(train$id)) & target=='clust')	{
    # IF DATA ARE AT THE INDV-LEVEL, BUT GOAL IS THE CLUSTER-LEVEL EFFECT 
    # get point estimates by aggregating to the cluster level 
    #   (e.g. by taking the weighted sum)
    # then take mean of cluster-level endpoints
    if(!doing.CV & verbose) print("data row are at individual level, but cluster-level target parameter (target = 'clust') is desired")
    R1 <- mean(aggregate(data.frame(train$alpha*train$Qbar1W.star), by=list(train$id), sum)[,2])
    R0 <- mean(aggregate(data.frame(train$alpha*train$Qbar0W.star), by=list(train$id), sum)[,2])
  } else {
    # OTHERWISE, JUST TAKE THE WEIGHTED MEAN ACROSS ALL ROWS
    # future work: robustify so that dont need weights that sum to J
    R1 <- mean(train$alpha*train$Qbar1W.star)
    R0 <- mean(train$alpha*train$Qbar0W.star) 
  }
  
  # UNSCALE THE OUTCOME 
  R1 <- R1*(scale_value - scale_value_min) + scale_value_min
  R0 <- R0*(scale_value - scale_value_min) + scale_value_min
  
  #==========================================================
  # Step 5: Variance estimation
  #==========================================================
  variance.out <- get.IC.variance(goal=goal, target=target, Vdata=train, R1=R1, R0=R0,
                                  scale_value = scale_value, scale_value_min = scale_value_min, 
                                  doing.CV = doing.CV, verbose = verbose, sample.effect = sample.effect)
  
  
  RETURN<- list(train=train, 	
                QAdj=Q$QAdj, Qform=Q$Qform, Q.out=Q$glm.out,
                gAdj=G$gAdj, gform=G$gform, p.out=G$p.out, 
                eps=eps, R1=R1, R0=R0, 
                var.R1 = variance.out$var.R1, 
                var.R0 = variance.out$var.R0,
                var.pair=variance.out$var.pair, 
                var.break=variance.out$var.break)	
  RETURN
}




#' Function to do initial estimation of E\[Y|A,W\] = Qbar(A,W)
#'
#' @param train Input data.
#' @param QAdj Adjustment variable(s)
#' @param Qform outcome regression fit 
#' @param glm.out 
#' @param verbose Prints more information about what's happening under the
#'   hood.
#'
#' @return List with adjustment variable(s),	outcome regression fit data set
#'   augmented w/ initial predictions: Qbar(A,W), Qbar(1,W) and Qbar(0,W)
#'   
#' @export
#'
#' @examples
do.Init.Qbar <- function(train, QAdj, Qform = 'glm',
                         glm.out = NULL, verbose){
  
  if(is.null(QAdj)){
    QAdj <- 'U'
  }
  
  relevant_columns <- c(as.character(QAdj), 'A', 'Y')
  train.temp <- train %>% dplyr::select(dplyr::all_of(relevant_columns))
  
  X1 <- X0 <- train.temp
  X1$A <-1; X0$A <- 0	
  # needed for penalized regression
  Xl <- model.matrix(~ -1 + ., subset(train.temp, select=-Y))
  X1l <- model.matrix(~ -1 + .,  subset(X1, select=-Y))
  X0l <- model.matrix(~ -1 + .,  subset(X0, select=-Y))
  
  if( is.null(glm.out) ){
    # fit using the training data
    # run main terms regression
    glm.out<- suppressWarnings( glm( Y ~ . , family='binomial', data=train.temp, weights=train$alpha ) )	
    
    if (Qform=='step'){ # stepwise 
      glm.out <- step(glm.out, direction = "both", trace = 0, k = 2)
    } else if (Qform=='step.interaction'){
      glm.out <- step(glm.out, scope=Y~.^2, direction = "both", trace = 0, k = 2)
    } else if (Qform=='lasso'){
      familyl <- ifelse(sum(unique(train$Y))>2, 'gaussian', 'binomial')
      glm.out <- glmnet::glmnet(x=Xl,  y=train$Y, weights=train$alpha,
                                family=familyl, alpha=1, nlambda = 100)
    } else if(Qform %in% c('mars', 'mars.corP')){
      # using default settings of SL.earth in SuperLearner 
      X <- subset(train.temp, select=-Y)
      if(Qform=='mars.corP'){
        X <- X[,screen.corP(Y=train$Y, X=X, family='binomial')==T]
      } 
      
      glm.out <- earth::earth(x = X,  y=train$Y, weights=train$alpha,
                              degree = 2, penalty = 3, 
                              nk = max(21, 2 * ncol(X) + 1), pmethod = "backward", nfold = 0, 
                              ncross = 1, minspan = 0, endspan = 0,
                              glm = list(family = binomial))
    } 
    #if(verbose) print(glm.out)
  }	
  
  # get initial predictions
  if(Qform=='lasso'){
    # from LASSO 
    QbarAW <- predict(glm.out, newx=Xl, type='response', s = min(glm.out$lambda))
    Qbar1W <- predict(glm.out, newx=X1l, type='response', s = min(glm.out$lambda))
    Qbar0W <- predict(glm.out, newx=X0l, type='response', s = min(glm.out$lambda))
  } else {
    # for glms 
    QbarAW <- predict(glm.out, newdata=train.temp, type='response')
    Qbar1W <- predict(glm.out, newdata=X1, type='response')
    Qbar0W <- predict(glm.out, newdata=X0, type='response')
  }
  Qbar <- data.frame(QbarAW, Qbar1W, Qbar0W)
  colnames(Qbar) <- c('QbarAW', 'Qbar1W', 'Qbar0W')
  list(QAdj=QAdj, Qform=Qform, glm.out=glm.out, train=cbind(train, Qbar))
}




#' Function to calculate the clever covariate
#'
#' @param train data set
#' @param gAdj adjustment variable(s)
#' @param gform pscore regression form
#' @param p.out 
#' @param verbose Prints more info
#'
#' @return List with adjustment variable(s), pscore regression, data set
#'   augmented with pscore & clever covariate (H.AW, H.1W, H.0W)
#'   
#' @export
#'
#' @examples
get.clever.cov <- function(train, gAdj, gform, p.out=NULL, verbose){
  
  if(is.null(gAdj)){
    gAdj <- 'U'
  }
  
  train.temp <- train[, c(gAdj, 'A')]  
  # needed for penalized regression
  Xl <- model.matrix( ~ -1 + ., subset(train.temp, select=-A) )
  
  if( is.null(p.out) ){ # fit pscore on training set 	
    # run main terms regression
    p.out<-   suppressWarnings( glm( A~. , family='binomial', data= train.temp, weights=train$alpha) )
    
    if (gform=='step'){ # stepwise 
      p.out <- step(p.out, direction = "both", trace = 0, k = 2)
    } else if (gform=='step.interaction'){
      p.out <- step(p.out, scope=A~.^2, direction = "both", trace = 0, k = 2)
    } else if (gform=='lasso'){
      p.out <- glmnet::glmnet(x=Xl,  y=train$A, weights=train$alpha,
                              family='binomial', alpha=1, nlambda = 100)
    } else if(gform %in% c('mars', 'mars.corP')){
      # using default settings of SL.earth in SuperLearner 
      X <- subset(train.temp, select=-A)
      if(gform=='mars.corP'){
        X <- X[,screen.corP(Y=train$A, X=X, family='binomial')==T]
      } 
      p.out <- earth(x = X,  y=train$A, weights=train$alpha,
                     degree = 2, penalty = 3, 
                     nk = max(21, 2 * ncol(X) + 1), pmethod = "backward", nfold = 0, 
                     ncross = 1, minspan = 0, endspan = 0,
                     glm = list(family = binomial))
    } 
    if(verbose){ print(p.out)}
    
  }
  
  # now use p.out to get estimated pscores
  if(gform!='lasso'){
    pscore <- predict(p.out, newdata= train.temp,  type="response")
  } else {
    pscore <- predict(p.out, newx=Xl, type='response', s = min(p.out$lambda) )
  }
  
  # bound g - should not apply for a randomized trial
  pscore [pscore < 0.025] <- 0.025
  pscore [pscore > 0.975] <- 0.975
  
  A.train <- train$A
  # Clever covariate is two-dimensional; 
  H.1W <- A.train/pscore 
  H.0W <- (1-A.train)/(1-pscore )
  # via delta method
  H.AW <- H.1W - H.0W
  
  p.df <- data.frame(pscore, H.1W , H.0W , H.AW)
  colnames(p.df) <- c('pscore', 'H.1W' , 'H.0W' , 'H.AW')
  list(gAdj=gAdj, gform=gform, p.out=p.out,  train=data.frame(train, p.df) ) 
}	




#' Function to calculate the fluctuation coefficient
#'
#' @param train data set
#' @param goal 'aRR' = arithmetic RR, 'OR' = odds ratio, 'RD' risk diff
#' @param verbose 
#'
#' @return Estimated fluctuation coefficient (eps)
#' 
#' @export
#'
#' @examples
get.epsilon <- function(train, goal, verbose=F){
  
  A.train<- train$A
  Y.train<- train$Y
  
  # Skip fitting if no meaningful variance
  # in either txt or control group
  Skip.update <-  (var(Y.train[A.train==1]) < 1*10^{-4}) |
    (var(Y.train[A.train==0]) < 1*10^{-4})
  
  if(goal=='RD'){ # if going after RD, then use a 1-dim clever covariate
    if(!Skip.update){
      logitUpdate <- suppressWarnings( 
        glm(Y.train ~ -1 + offset(qlogis(train$QbarAW)) + train$H.AW, family="binomial", weights=train$alpha))
      eps <- logitUpdate$coef
    } else {
      eps <- 0
    }
    names(eps) <- 'H.AW'
  } else { # targeting the risk or odds ratio requires a two-dimensional clever covariate
    if( !Skip.update  ){
      logitUpdate <- suppressWarnings(
        glm(Y.train ~ -1 + offset(qlogis(train$QbarAW)) + train$H.0W + train$H.1W, family="binomial", weights=train$alpha))
      eps <- logitUpdate$coef
    } else {
      eps <- c(0,0)
    }
    names(eps)<- c('H.0W', 'H.1W')	
  }
  if(verbose) print(eps)
  return(eps)
}





#' Function to update initial estimators of QbarAW
#'
#' @param train data set
#' @param eps fluctuation coefficient
#' @param goal 'RD' = risk difference; otherwise RR/OR
#'
#' @return data.frame w/ targeted predictions: Qbar*(A,W), Qbar*(1,W), Qbar*(0,W)
#' 
#' @export
#'
#' @examples
do.targeting <- function(train, eps, goal){
  
  g1W <- train$pscore
  g0W <- (1 - g1W)
  
  if(goal=='RD'){
    # updated QbarAW estimates for training set. 
    QbarAW.star <- plogis( qlogis(train$QbarAW ) + eps*train$H.AW)	
    Qbar1W.star <- plogis( qlogis(train$Qbar1W ) + eps/g1W )
    Qbar0W.star <- plogis( qlogis(train$Qbar0W ) - eps/g0W )
  } else {
    # updated QbarAW estimates for training set. 
    QbarAW.star <- plogis( qlogis(train$QbarAW) + eps['H.0W']*train$H.0W + eps['H.1W']*train$H.1W)	
    Qbar0W.star <- plogis( qlogis(train$Qbar0W) + eps['H.0W']/g0W )
    Qbar1W.star <- plogis( qlogis(train$Qbar1W) + eps['H.1W']/g1W )
  }
  train <- data.frame(train, QbarAW.star, Qbar1W.star, Qbar0W.star)		
  return(train)
}


#### utility_functions.R ####

`%nin%` <- purrr::negate(`%in%`)

#' Harmonic mean
#'
#' @description Calculates the harmonic mean of a vector of numbers.
#'
#' @param a Vector of values
#'
#' @return The harmonic mean of the values: 1 / mean(1/a).
#' 
#' @export
mean.harmonic <- function(a){
  1/mean(1/a) #compute the harmonic mean
}



make_adj_var_df <- function(vars, Q_or_g = "Q"){
  if(is.null(vars)){
    stop("No adjustment variables selected. This should be 'U' in this case... something wrong? See utility_functions.R")
  }
  adj_row <- data.frame(matrix(NA, nrow = 1, ncol = length(vars)))
  adj_row[1,] <- vars
  colnames(adj_row) <- paste0(Q_or_g,"Adj",1:length(vars))
  return(adj_row)
}


#' Generate set of working models for adaptive pre-specification
#'
#' @description Function to get candidate adjustment strategies (variables +
#'   algorithms) for estimating the outcome regression and the propensity score.
#' 
#' @param cand.vars Character vector of variable names that may be used in the
#'   adjustment set. Default is NULL, in which case it returns a single
#'   unadjusted working model.
#' @param cand.algos Character vector of algorithm possibilities. See 'Qform' in
#'   Stage2 documentation for supported options. Default is NULL, in which case
#'   only GLMs are returned.
#' @param data.adapt.complexity Either 'low' (default), 'med', or 'high'. 'low'
#'   returns GLMs with a single adjustment variable, while 'med' adds all
#'   cand.algos with the full cand.vars specification. 'high' includes all
#'   non-redundant subsets of algorithms and variables, and may lead to
#'   overfitting; use carefully.
#'
#' @return A data frame where each row has a candidate model: The first column
#'   is the algorithm and the second column is a list of variables to use with
#'   it.
#' @export
#'
#' @examples
get.cand.adj <- function(cand.vars, cand.algos = NULL, data.adapt.complexity){
  unadjusted <- cbind.data.frame(cand.algos = 'glm', cand.varset = 'U')
  nvars <- length(cand.vars)
  if(nvars == 1){
    if(cand.vars == "U"){
      return(unadjusted)
    } else {
      data.adapt.complexity <- "low"
    }
  }
  
  cand.vars <- cand.vars[cand.vars %nin% 'U'] # We manually add U, so ignore if they add it as well.
  
  if(is.null(cand.algos)) {
    ###### simple Adaptive Prespec - GLMs with one adjustment variable each.
    sets <- expand.grid(cand.algos = 'glm', cand.varset = union('U', cand.vars))
  } else {                  
    ###### fancy adaptive prespec with expanded algorithms
    simple_sets <- cbind.data.frame(cand.algos = "glm", cand.varset = cand.vars)
    if(nvars > 1){
      full_set <- cbind.data.frame(cand.algos = cand.algos, cand.varset = NA) %>%
        mutate(cand.varset = list(cand.vars))
      
      # Find all of the proper subsets of length > 1
      combos <- do.call(base::c, lapply(seq_along(cand.vars), combn, x = cand.vars, simplify = FALSE))
      not_listed_yet <- intersect(which(lapply(combos, length) != 1), which(lapply(combos, length) != nvars))
      combos <- combos[not_listed_yet]
      middle_sets <- expand.grid(cand.algos = cand.algos, cand.varset = combos)
    }
    
    # Merge relevant combos (removing some unneeded duplicates) depending on desired level of complexity    
    if(data.adapt.complexity == "high"){
      message("High complexity; APS considering all unique subsets of cand.adj with each cand.form")
      sets <- rbind(simple_sets, full_set, middle_sets,
                    unadjusted)
    } else if (data.adapt.complexity == "med") {
      message("Medium complexity; APS considering GLMs with one adjustment variable, plus the full set of cand.adj with each cand.form")
      sets <- rbind(simple_sets, full_set, unadjusted)
    } else {
      if(length(cand.algos > 1)){
        message("Low complexity being used to prevent overfitting; APS considering GLMs with one adjustment variable only")
      }
      sets <- rbind(simple_sets, unadjusted)
    }
  }
  return(sets)
}




#' Check data inputs for Stage2()
#'
#' @description Performs several checks and throws warnings before running any
#'   statistical analysis.
#'
#' @param data.input Data frame input.
#' @param verbose Logical indicating if you want a lot of details to be printed
#'   to the console.
#'
#' @return A data freame with the necessary structure for Stage2().
#' @export
#'
#' @examples
preprocess_and_check_data <- function(data.input, verbose){
  if("U" %nin% colnames(data.input)){ # Will now add dummy column "U" if not in dataset already
    data.input <- data.input %>% dplyr::mutate(U = 1)
  }
  if("Y" %nin% colnames(data.input)){
    stop("No Y column in data.input")
  }
  if("id" %nin% colnames(data.input)){
    warning("data.input does not have a cluster 'id' column, assuming each row is an independent observation",
            call. = F)
    data.input <- cbind.data.frame(data.input, id = 1:nrow(data.input))
  }
  if("alpha" %nin% colnames(data.input)){
    if(verbose) print("data.input does not have an 'alpha' column, assuming each row is weighted equally.")
    data.input <- data.input %>% dplyr::mutate(alpha = 1)
  }
  if("U" %in% colnames(data.input))
    if(sum(data.input$U == 1) != nrow(data.input)){
      warning("You have a variable 'U' in your dataset that is not a column of 1s. Replacing with a column of 1s. If you have a variable named U, please rename.")
      data.input <- data.input %>% dplyr::mutate(U = 1)
    }
  return(data.input)
}






#' Aggregate influence curves to the cluster level
#'
#' Aggregates individual-level influence curve (IC) estimates to the cluster
#' level. This allows for appropriate statistical inference that respects the
#' cluster as the independent unit, whether the target parameter is an
#' individual-level or cluster-level effect.
#'
#' @param recordIC The individual-level IC estimates
#' @param id Vector of cluster membership IDs
#'
#' @return A vector of aggregated IC values with length equal to the number of
#'   unique cluster IDs.
#' @export
#'
#' @examples
aggregate_IC <- function(recordIC, id) {
  if (is.null(id)) return(recordIC)
  aggregatedIC <- as.matrix(stats::aggregate(recordIC, list(id=id), sum)[, -1, drop = FALSE])
  num.records <- nrow(recordIC)
  num.clusters <- nrow(aggregatedIC)
  aggregatedIC <- aggregatedIC * num.clusters / num.records
  return(aggregatedIC)
}



#' Calculate confidence intervals and p-values on relative or absolute scale
#'
#' @param goal String specifying the scale of the target parameter. Default is
#'   \code{RD}, risk/rate difference. Any other values assume that input values
#'   are given on the log scale, and the function will exponentiate the
#'   estimated target parameter and confidence interval bounds to output an
#'   artihmetic risk/rate ratio.
#' @param psi True value (if known) of the target parameter, for example, in a
#'   simulation study.
#' @param psi.hat Estimated value of the target parameter.
#' @param se Standard error of the estimated target parameter.
#' @param df Degrees of freedom for the Student's \emph{t} distribution as an
#'   approximation of the asymptotic normal distribution. If \code{df > 40}, the
#'   value is ignored and a normal distribution is used for inference.
#' @param sig.level Desired significance (alpha) level. Defaults to 0.05.
#' @param one.sided Logical indicating that a one-sided test is desired.
#'   Defaults to \code{FALSE}.
#' @param alt.smaller If one-sided test is desired, is the alternative
#'   hypothesis that the intervention arm will have a smaller value that the
#'   control arm? For example, if you expect a public health intervention to
#'   reduce the mean of a disease outcome, use alt.smaller = TRUE (this
#'   corresponds to a null hypothesis that the intervention did not reduce the
#'   mean disease outcome).
#'   
#' @return A one-row data frome with the estimated target parameter value
#'   (\code{est}), the (two-sided) confidence interval \code{CI.lo},
#'   \code{CI/hi}, the standard error of the estimate, the (possibly one-sided)
#'   p-value, and bias/coverage/rejection indicators (if true value or target
#'   parameter is supplied). NOTE: If \code{goal != "RD"}, the output standard
#'   error will be on the log scale.
#' @export
#'
#' @examples
get.inference <- function(goal = 'RD', psi = NA, psi.hat, se, df = 99, sig.level = 0.05, 
                          one.sided = F, alt.smaller = NULL){
  
  # test statistic
  # (on the log-transformed scale if goal is arithmetic RR or odds ratio)
  tstat <- psi.hat / se
  
  if(df > 40){
    # assume normal distribution
    cutoff <- stats::qnorm(sig.level/2, lower.tail=F)
    # one.sided hypothesis test 
    if(one.sided){
      pval <- stats::pnorm(tstat, lower.tail=alt.smaller) 
    } else{
      pval<- 2*stats::pnorm(abs(tstat), lower.tail=F) 
    }
  } else {
    # use Student's t-distribution
    # print('Using t-distribution')
    cutoff <- stats::qt(sig.level/2, df=df, lower.tail=F)
    # one.sided hypothesis test if specified
    if(one.sided){
      pval <- stats::pt(tstat, df=df, lower.tail= alt.smaller ) 
    } else{
      pval <- 2*stats::pt(abs(tstat), df=df, lower.tail=F)
    }
  }
  # 95% confidence interval 
  CI.lo <- (psi.hat - cutoff*se)
  CI.hi <- (psi.hat + cutoff*se)
  
  # If on relative (log) scale, transform back 
  if(goal != 'RD'){
    psi.hat <- exp(psi.hat)
    CI.lo <- exp(CI.lo)
    CI.hi <- exp(CI.hi)
  }  
  
  # bias
  bias <- (psi.hat - psi)
  # confidence interval coverage?
  cover <-  CI.lo <= psi & psi <= CI.hi
  # reject the null?
  reject <- as.numeric(pval < sig.level)
  return(data.frame(est = psi.hat, CI.lo, CI.hi, se = se, pval, bias, cover, reject))
}








#' Calculate the variance of the influence curve
#'
#' @description NOTE: UPDATE OF UNSCALING HAPPENS HERE
#'
#' @param goal One of 'aRR'  (arithmetic risk ratio), 'RD' (risk difference), or
#'   'OR' (odds ratio)
#' @param target Target of inference, either cluster-level ('clust') or pooled
#'   individivual-level effect ('indv')
#' @param Vdata data set
#' @param R1
#' @param R0
#' @param sample.effect Logical indicating whether to perform inference for the
#'   sample average treatment effect (if TRUE) or the population average
#'   treatment effect (if FALSE).
#' @param scale_value maximum value for outcome scaling
#' @param scale_value_min minimum value for outcome scaling
#' @param doing.CV
#'
#' @return on log scale for if goal='aRR' or 'OR' ... estimated IC & variance -
#'   preserving/breaking the match
#' @export
#'
#' @examples
get.IC.variance <- function(goal, target, Vdata, R1 = NA, R0 = NA,
                            scale_value = 1, scale_value_min = 0,
                            doing.CV = F,
                            verbose = F, sample.effect){
  
  # number of randomized units
  J <- length(unique(Vdata$id))
  
  # calculate the relevant components of the IC 
  if(sample.effect){
    # default - assume interest is in the sample effect
    DY1 <- Vdata$alpha*Vdata$H.1W*(Vdata$Y - Vdata$Qbar1W.star)
    DY0 <- Vdata$alpha*Vdata$H.0W*(Vdata$Y - Vdata$Qbar0W.star)
  } else {
    # calculate the IC for population effect (extra term for DW)
    DY1 <- Vdata$alpha*( Vdata$H.1W*(Vdata$Y - Vdata$Qbar1W.star) + Vdata$Qbar1W.star - R1 )
    DY0 <- Vdata$alpha*( Vdata$H.0W*(Vdata$Y - Vdata$Qbar0W.star) + Vdata$Qbar0W.star - R0 )	
  }
  
  # unscale 
  DY1 <- DY1*(scale_value - scale_value_min) + scale_value_min
  DY0 <- DY0*(scale_value - scale_value_min) + scale_value_min
  
  # if individual-level data, then need to aggregate the IC to the cluster-level 
  # approach for aggregation depends on the target effect
  if(length(DY1) > J) {
    if(target=='clust'){
      # Data are indv-level; target is cluster-level 
      if(!doing.CV & verbose) print('data = indv; target = clust')
      DY1 <- aggregate(DY1, by = list(Vdata$id), sum)[,-1]
      DY0 <- aggregate(DY0, by = list(Vdata$id), sum)[,-1]
    } else {
      # Data are indv-level; target is indv-level 
      if(!doing.CV & verbose) print('data = indv; target = indv')
      DY1 <- c(aggregate_IC(as.matrix(DY1), id = Vdata$id))
      DY0 <- c(aggregate_IC(as.matrix(DY0), id = Vdata$id))
    }
    # for the pair-matched IC also need to aggregate to the cluster-level
    # Vdata <- aggregate(Vdata, by=list(Vdata$id), mean)[,-1]
  } 
  
  # INFLUCENCE CURVES ARE NOW AT THE LEVEL OF THE RANDOMIZED UNIT
  if(goal=='RD'){
    # going after RD, easy IC
    DY <-  DY1 - DY0
  } else if (goal=='aRR'){ 
    # going after aRR, then get IC estimate on log scale
    #	i.e. Delta method for log(aRR) = log(R1) - log(R0)
    DY <- 1/R1*DY1 - 1/R0*DY0
  } else if(goal=='OR'){
    # Delta method for log(OR)
    DY <- 1/R1*DY1 + 1/(1-R1)*DY1 - 1/(1-R0)*DY0 - 1/R0*DY0
  }
  
  if(!doing.CV & verbose) print(paste0('Mean DY (Solve EIF): ', mean(DY) ))
  
  
  # estimated variance for txt specific means or if break the match	
  var.R1 <- stats::var(DY1) / J ## Doesn't work for LOOCV
  var.R0 <- stats::var(DY0) / J ## Doesn't work for LOOCV
  var.break <- stats::var(DY) / J
  #print(var.R1)
  
  if( 'pair' %in% colnames(Vdata) ){
    # estimated variance if preserve the match
    pairC <- aggregate(Vdata, by=list(Vdata$id), mean)[,'pair']
    pairs <- unique(pairC)
    n.pairs <- length(pairs)
    DY.paired <-  rep(NA, n.pairs)
    for(i in 1:n.pairs){		
      these<- pairC %in% pairs[i] 
      DY.paired[i]<- 0.5*sum(DY[ these] )			
    }
    
    var.pair <- stats::var(DY.paired) / n.pairs
  } else {
    DY.paired <- var.pair <- NA
  }
  
  return(list(R1=R1, R0=R0, DY1=DY1, DY0=DY0,var.R1 = var.R1, var.R0 = var.R0,
              DY=DY, var.break=var.break, 
              DY.paired=DY.paired, var.pair=var.pair))
}

## Workspace ----
########### Warning ###########
# Base_file corresponds to the file where Workspace/Data/Results will be saved
# Choose one and keep the same for all functions .R files, all simulation .R files


Base_file = ""
Workspace_name = "Workspace_Balzer_TMLE.RData"

save.image(paste(Base_file,Workspace_name,sep = "/"))