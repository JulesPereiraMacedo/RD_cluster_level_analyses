# 1) Fun_Crea_data_confusion ----

## Input ##

# pI               : Prevalence in Intervention arm
# pC               : Prevalence in Control arm
# k                : Number of cluster per arm
# icc              : Intra-cluster correlation coefficient 
# nb_cov           : Number of "individual" covariates
# p_cov            : Vector of prevalence of each "individual-level" covariates
# OR_cov           : Vector of Odd Ratio between each "individual-level" covariates and the outcome
# nb_cov_clus      : Number of "cluster" covariates
# p_cov_clus       : Vector of prevalence of each "cluster-level" covariates
# OR_cov_clus      : Vector of Odd Ratio between each "cluster-level" covariates and the outcome
# Pi_int           : Proportion of the intervention group that we will be included to add confusion
# Pi_con           : Proportion of the control group that we will be included to add confusion 

## Output ##

# Data        : A 'data.frame' which correspond to data set before selection of participant to include confusion

Crea_Data_confusion <- function(pI,pC,k,m,icc,nb_cov,p_cov,OR_cov,nb_cov_clus,p_cov_clus,OR_cov_clus,Pi_int,Pi_con){
  if (length(p_cov) < nb_cov || length(OR_cov) < nb_cov || length(p_cov_clus) < nb_cov_clus || length(OR_cov_clus) < nb_cov_clus ){stop('Vector of covariates probability or OR vector have lower size than the number of covariates')}
  
  m_int = m/Pi_int # increase the number of the mean of individual in a cluster in the intervention group because when we will apply the selection to create the confusion the number of individuals will decrease 
  m_con = m/Pi_con # increase the number of the mean of individual in a cluster in the control group because when we will apply the selection to create the confusion the number of individuals will decrease 
  
  ## Number of subject per cluster
  ## Suj_per_num[i,j] ~ Binom Negative(m,0.5)
  Suj_per_clus = matrix(0,nrow = k,ncol = 2)
  cv_int = 0.47
  cv_con = 0.47
  for (l in 1:(2*k)) {
    if(l<=k){m=m_int;cv=cv_int}else{m=m_con;cv=cv_con}
    Suj_per_clus[l] = rnbinom(n = 1,mu = m, size = m**2/((cv*m)**2-m))
  }
  
  while(TRUE %in% (Suj_per_clus <= 1)){
    for (l in 1:(2*k)) {
      if(l<=k){m=m_int;cv=cv_int}else{m=m_con;cv=cv_con}
      Suj_per_clus[l] = rnbinom(n = 1,mu = m, size = m**2/((cv*m)**2-m))
    }
  }
  
  
  ## Creation of the outcome and co variables 
  Data = c()
  Outcome = c()
  nb_indiv = sum(Suj_per_clus)
  if(nb_cov != 0){
    Covariate = matrix(0,ncol = nb_cov,nrow = nb_indiv)
    colnames_vec_Indiv_Cov = c()
    for (i in 1:nb_cov) {
      colnames_vec_Indiv_Cov[i] = paste("IC_",i,sep = "")
    }
    colnames(Covariate) <- colnames_vec_Indiv_Cov
  }
  if(nb_cov_clus != 0){
    Covariate_clus = matrix(0,nrow = nb_indiv,ncol = nb_cov_clus)
    colnames_vec_Cluster_Cov = c()
    for (i in 1:nb_cov_clus) {
      colnames_vec_Cluster_Cov[i] = paste("CC_",i,sep = "")
    }
    colnames(Covariate_clus) <- colnames_vec_Cluster_Cov
  }
  
  ## Cluster Effect
  Z_ij_int = rbinom(n = k,1,pI)
  Z_ij_con = rbinom(n = k,1,pC)
  
  Z_ij = matrix(c(Z_ij_int,Z_ij_con),nrow = k,ncol = 2)
  
  for (l in 1:(2*k)) {
    if(l <= k){p=pI}else{p=pC}
    n_ij = Suj_per_clus[l]
    
    ## Outcome : Y = (1-U)V + UZ
    ## U ~ Binom(1,sqrt(icc))
    ## V ~ Binom(1,p)  p = pI or pC (depend on which arm the individual is)
    ## Z ~ Binom(1,p)  same for individuals in one cluster
    ## Y ~ Binom(1,p)
    
    U_ij = rbinom(n_ij,1,sqrt(icc))
    V_ij = rbinom(n_ij,1,p)
    
    Y_ij = (1-U_ij)*V_ij + U_ij*Z_ij[l]
    Outcome = c(Outcome,Y_ij)
    
    
    ## Covariates i :  X_i = A_i*W_i
    ## A_ijl ~ binom(1,alpha_i)
    ## alpha_i = p_cov[i]/p  ///  p = pI or pC (depend on which arm the individual is)
    ## W_ijl = (1-R_ijl)*T_ijl + R_ijl*Y
    ## R_ijl ~ binom(1,r_i)
    ## r_i = rho_cov[i] * ( sqrt((1-p)*(1-p_cov[i])) / (sqrt(alpha_i*alpha_i+1)*(1-p)) )
    ## T_ijl ~ binom(p)
    ## Y : Outcome
    ## W_ijl ~ binom(1,p)
    
    if(nb_cov !=0){
      for (c in 1:nb_cov) {
        
        q_m = p_cov[c]
        
        # Find the good coefficient of correlation 'rho' to create the individual covariate with the correspondant OR between the covariate and the outcome
        if(OR_cov[c]!=1){
          rho <- find_rho(p,q_m,OR_cov[c]) # find 'find_rho()' in the section '2 - Correlation verification function'
          if(length(rho) == 2){
            rdm = rbinom(1,1,1/2) + 1
            rho = rho[rdm]
          }
          if(length(rho)==0){stop('No solution to the equation to find a coefficient of correlation for the individuals covariates')}
          
          c_min_max <- corr_min_max(p,pC,q_m)
          if(rho <= c_min_max[1] | rho >= c_min_max[2]){stop("Doesn't valid the hypothesis of Prentice for the correlation between two binary variables" )}
        }
        if(OR_cov[c]==1){rho = 0}
        
        
        alpha = q_m/p
        r_m = rho * ( sqrt(1-q_m) / (sqrt(alpha)*sqrt(1-p) ))
        
        R_ijl = rbinom(n_ij,1,r_m)
        T_ijl = rbinom(n_ij,1,p)
        if(l==1){
          W_ijl = (1-R_ijl)*T_ijl + R_ijl*Outcome[1:Suj_per_clus[l]]
        }else{
          W_ijl = (1-R_ijl)*T_ijl + R_ijl*Outcome[(sum(Suj_per_clus[1:(l-1)])+1):sum(Suj_per_clus[1:l])]
        }
        
        A_ijl = rbinom(n_ij,1,alpha)
        X_ijl = A_ijl*W_ijl
        if(l==1){
          Covariate[1:Suj_per_clus[l],c] = X_ijl
        }else{
          Covariate[(sum(Suj_per_clus[1:(l-1)])+1):sum(Suj_per_clus[1:l]),c] = X_ijl
        }
      }
    }
  }   
  
  if(nb_cov_clus!=0){
    for (h in 1:nb_cov_clus){
      q = p_cov_clus[h]
      
      # Aij follow a binomial distribution with parameter alpha = q/p
      # Wij follow a binomial distribution with parameter p such as W : (1-R)T + RZ
      # Rij ~ B(r) with r find with the calcul based on the OR between the outcome and the cluster level covariate
      # Tij ~ B(p)
      # Zij ~ B(p) Same as the Z variable generated for the Outcome Y
      
      if(OR_cov_clus[h]!=1){
        # Finding parameter r to simulate Rij
        r_int = find_r_clus(pI,q,OR_cov_clus[h],icc) # Function 'find_r_clus()' described in section '2 - Correlation verification function'
        r_con = find_r_clus(pC,q,OR_cov_clus[h],icc) # Function 'find_r_clus()' described in section '2 - Correlation verification function'
      }
      if(OR_cov_clus[h]==1){
        r_int = 0 # If the OR between the outcome and the cluster level covariate is 1, the r parameter should be 0
        r_con = 0
      }
      
      Rij_int = rbinom(n = k,size = 1,prob = r_int)
      Rij_con = rbinom(n = k,size = 1,prob = r_con)
      
      Rij = matrix(c(Rij_int,Rij_con),nrow = k,ncol = 2)
      
      Tij_int = rbinom(n = k,size = 1,prob = pI)
      Tij_con = rbinom(n = k,size = 1,prob = pC)
      
      T_ij = matrix(c(Tij_int,Tij_con),nrow = k,ncol = 2)
      
      Aij_int = rbinom(n = k,size = 1,prob = q/pI)
      Aij_con = rbinom(n = k,size = 1,prob = q/pC)
      
      Aij = matrix(c(Aij_int,Aij_con),nrow = k,ncol = 2)
      Wij = (1-Rij)*T_ij + Rij*Z_ij
      
      Xij = Aij*Wij
      
      verif_Xij = apply(Xij,MARGIN = 2,FUN = sum)
      while (verif_Xij[1]==0 | verif_Xij[1] == k | verif_Xij[2]==0 | verif_Xij[1] == k) {
        Rij_int = rbinom(n = k,size = 1,prob = r_int)
        Rij_con = rbinom(n = k,size = 1,prob = r_con)
        
        Rij = matrix(c(Rij_int,Rij_con),nrow = k,ncol = 2)
        
        Tij_int = rbinom(n = k,size = 1,prob = pI)
        Tij_con = rbinom(n = k,size = 1,prob = pC)
        
        T_ij = matrix(c(Tij_int,Tij_con),nrow = k,ncol = 2)
        
        Aij_int = rbinom(n = k,size = 1,prob = q/pI)
        Aij_con = rbinom(n = k,size = 1,prob = q/pC)
        
        Aij = matrix(c(Aij_int,Aij_con),nrow = k,ncol = 2)
        Wij = (1-Rij)*T_ij + Rij*Z_ij
        
        Xij = Aij*Wij
        
        verif_Xij = apply(Xij,MARGIN = 2,FUN = sum)
      }
      
      Covariate_clus[,h] = rep(Xij,Suj_per_clus)
    }
  }
  
  ## Assignation of the arm effect for each individuals
  Arm = c(rep(1,apply(Suj_per_clus, 2, sum)[1]),rep(0,apply(Suj_per_clus, 2, sum)[2]))
  
  ## Assignation of the cluster number for each individuals
  cluster = c()
  for (l in 1:(2*k)) {
    cluster = c(cluster,rep(l,Suj_per_clus[l]))
  }
  
  Data = data.frame(Outcome,Arm,cluster)
  if(nb_cov!=0){Data = data.frame(Data,Covariate)}
  if(nb_cov_clus!=0){Data = data.frame(Data,Covariate_clus)}
  
  return(Data)
}

# 1b) Crea_data_ICS ----

Crea_Data_ICS <- function(pI,pC,k_1,k_0,m,icc,nb_cov,OR_cov,nb_cov_clus,OR_cov_clus,alpha=3){
  if (length(OR_cov) < nb_cov || length(OR_cov_clus) < nb_cov_clus ){stop('Vector of covariates probability or OR vector have lower size than the number of covariates')}
  
  ##### Generate each risk for each cluster with a beta distribution with parameters B(alpha, beta) calculated in function of
  ##### the mean prevalence in each Arm.
  
  k = k_1 + k_0
  
  Risk_vec = matrix(NA,1,k)
  Arm = c(rep(1,k_1),rep(0,k_0))
  
  # alpha specified in function's argument
  
  alpha_inter <- alpha_contr <- alpha
  
  # beta_inter calculated with the prevalence of the intervention arm
  # beta_contr calculated with the prevalence of the control arm
  
  beta_inter  = (1-pI)*alpha_inter/pI
  beta_contr  = (1-pC)*alpha_contr/pC
  
  Risk_vec[which(Arm == 1)] <- rbeta(k_1,alpha_inter,beta_inter)
  Risk_vec[which(Arm == 0)] <- rbeta(k_0,alpha_contr,beta_contr)
  
  # Generate the cluster size for each cluster based on the risk of each cluster (cf: Hoffman simulation plan)
  # If the risk is under the mean of the beta variable (alpha/(alpha+beta)) we generate more individual per cluster
  # This means that the 'large' clusters will have less risk than 'small' clusters
  # The prevalence will be higher in 'small' cluster than in 'large' clusters
  
  Suj_per_clus = c()
  
  # To have the mean require 'm', we need to amplifie when we generate with the binomial distribution
  
  for (i in 1:k) {
    if(Arm[i]==1){if(Risk_vec[i] < (alpha_inter/(alpha_inter+beta_inter))){Suj_per_clus[i] = sum(rbinom(n = m*2,1,0.66))}else{Suj_per_clus[i] = sum(rbinom(n = m*2,1,0.33))}
    }else{
      if(Risk_vec[i] < (alpha_contr/(alpha_contr+beta_contr))){Suj_per_clus[i] = sum(rbinom(n = m*2,1,0.66))}else{Suj_per_clus[i] = sum(rbinom(n = m*2,1,0.33))}
    }
  }
  
  while (0 %in% Suj_per_clus) {
    for (i in 1:k) {
      if(Arm[i]==1){if(Risk_vec[i] < (alpha_inter/(alpha_inter+beta_inter))){Suj_per_clus[i] = sum(rbinom(n = m*2,1,0.66))}else{Suj_per_clus[i] = sum(rbinom(n = m*2,1,0.33))}
      }else{
        if(Risk_vec[i] < (alpha_contr/(alpha_contr+beta_contr))){Suj_per_clus[i] = sum(rbinom(n = m*2,1,0.66))}else{Suj_per_clus[i] = sum(rbinom(n = m*2,1,0.33))}
      }
    }
  }
  
  nb_indiv = sum(Suj_per_clus)
  
  # Creation de l'outcome
  
  
  
  
  if(nb_cov != 0){
    rho_vec <- matrix(NA,k,nb_cov)
    r_m_vec <- matrix(NA,k,nb_cov)
    Covariate = matrix(0,ncol = nb_cov,nrow = nb_indiv)
    colnames_vec_Indiv_Cov = c()
    for (i in 1:nb_cov) {
      colnames_vec_Indiv_Cov[i] = paste("IC_",i,sep = "")
    }
    colnames(Covariate) <- colnames_vec_Indiv_Cov
  }
  if(nb_cov_clus != 0){
    Covariate_clus = matrix(0,nrow = nb_indiv,ncol = nb_cov_clus)
    colnames_vec_Cluster_Cov = c()
    for (i in 1:nb_cov_clus) {
      colnames_vec_Cluster_Cov[i] = paste("CC_",i,sep = "")
    }
    colnames(Covariate_clus) <- colnames_vec_Cluster_Cov
  }
  
  Outcome = c()
  indiv_l = 1
  Z_j = matrix(NA,1,k)
  
  for (j in 1:k) {
    n_j = Suj_per_clus[j]
    p = Risk_vec[j]
    # p = ifelse(j<=k_1,pI,pC)
    Z_j[j] = rbinom(1,1,p)
    
    q_m = min(Risk_vec)
    # q_m = pC
    
    alpha_m = q_m/p
    
    ## Find rho for each individual covariates, cluster by cluster.
    
    if(nb_cov !=0){
      for (c in 1:nb_cov) {
        # Find the good coefficient of correlation 'rho' to create the individual covariate with the correspondant OR between the covariate and the outcome
        if(OR_cov[c]!=1){
          rho <- find_rho(pI,pC,OR_cov[c]) # find 'find_rho()' in the section '2 - Correlation verification function'
          if(length(rho)==1){rho_vec[j,c] = rho}
          if(length(rho) == 2){
            rdm = rbinom(1,1,1/2) + 1
            rho_vec[j,c] = rho[rdm]
          }
          if(length(rho)==0){stop('No solution to the equation to find a coefficient of correlation for the individuals covariates')}
          
          c_min_max <- corr_min_max(pI,pI,pC) # Here because we simulate a correlation in function of the each cluster so the validation of the prentice theorem should be verified for each cluster
          if(rho <= c_min_max[1] | rho >= c_min_max[2]){stop("Doesn't valid the hypothesis of Prentice for the correlation between two binary variables" )}
        }
        if(OR_cov[c]==1){rho_vec[j,c] = 0}
        
        r_m_vec[j,c] = rho_vec[j,c] * ( sqrt(1-pC) / (sqrt(alpha_m)*sqrt(1-pI) ))
      }
    }
    
    for (l in 1:n_j) {
      U = rbinom(1,1,sqrt(icc))
      V = rbinom(1,1,p)
      Y = (1-U)*V + U*Z_j[j]
      Outcome = c(Outcome,Y)
      
      # Covariable individuelle
      
      if(nb_cov !=0){
        for (c in 1:nb_cov) {
          
          R_ijl = rbinom(1,1,r_m_vec[j,c])
          T_ijl = rbinom(1,1,p)
          
          W_ijl = (1-R_ijl)*T_ijl + R_ijl*Y
          
          
          A_ijl = rbinom(1,1,alpha_m)
          X_ijl = A_ijl*W_ijl
          
          Covariate[indiv_l,c] = X_ijl
        }
      }
      indiv_l = indiv_l + 1
    }
  }
  
  #### Covariable au niveau cluster
  if(nb_cov_clus!=0){r_vec = matrix(NA,k,nb_cov_clus)
  Rij_vec = matrix(NA,k,nb_cov_clus)
  Tij_vec = matrix(NA,k,nb_cov_clus)
  Aij_vec = matrix(NA,k,nb_cov_clus)
  Wij_vec = matrix(NA,k,nb_cov_clus)
  Xij_vec = matrix(NA,k,nb_cov_clus)
  # Covariate_clus = matrix(NA,nb_indiv,nb_cov_clus)
  }
  
  
  if(nb_cov_clus!=0){
    for (h in 1:nb_cov_clus){
      for (j in 1:k) {
        p = Risk_vec[j]
        q = min(Risk_vec)
        # Aij_vec[j,h] = q/p
        # Aij follow a binomial distribution with parameter alpha = q/p
        # Wij follow a binomial distribution with parameter p such as W : (1-R)T + RZ
        # Rij ~ B(r) with r find with the calcul based on the OR between the outcome and the cluster level covariate
        # Tij ~ B(p)
        # Zij ~ B(p) Same as the Z variable generated for the Outcome Y
        
        if(OR_cov_clus[h]!=1){
          # Finding the parameter r
          r_vec[j,h] = find_r_clus(pI,pC,OR_cov_clus[h],icc) # Function 'find_r_clus()' described in section '2 - Correlation verification function'
        }
        if(OR_cov_clus[h]==1){
          r_vec[j,h] = 0 # If the OR between the outcome and the cluster level covariate is 1, the r parameter should be 0
        }
        Rij_vec[j,h] = rbinom(n = 1,size = 1,prob = r_vec[j,h])
        
        Tij_vec[j,h] = rbinom(n = 1,size = 1,prob = p)
        
        Aij_vec[j,h] = rbinom(n = 1,size = 1,prob = q/p)
        
        Wij_vec[j,h] = (1-Rij_vec[j,h])*Tij_vec[j,h] + Rij_vec[j,h]*Z_j[j]
        
        Xij_vec[j,h] = Aij_vec[j,h]*Wij_vec[j,h]
      }
      
      verif_Xij = sum(Xij_vec[,h])
      # while (verif_Xij==0 | verif_Xij == k) {
      #   Rij_vec[j,h] = rbinom(n = 1,size = 1,prob = r_vec[j,h])
      # 
      #   Tij_vec[j,h] = rbinom(n = 1,size = 1,prob = p)
      # 
      #   Aij_vec[j,h] = rbinom(n = 1,size = 1,prob = q/p)
      # 
      #   Wij_vec[j,h] = (1-Rij_vec[j,h])*Tij_vec[j,h] + Rij_vec[j,h]*Z_j[j]
      # 
      #   Xij_vec[j,h] = Aij_vec[j,h]*Wij_vec[j,h]
      # 
      #   verif_Xij = sum(Xij_vec[,h])
      # }
      
      Covariate_clus[,h] = rep(Xij_vec[,h],Suj_per_clus)
    }
  }
  
  ## Assignation of the arm effect for each individuals
  Arm_by_indiv = rep(Arm,Suj_per_clus)
  
  ## Assignation of the cluster number for each individuals
  cluster = rep(1:k,Suj_per_clus)
  
  Data = data.frame(Outcome,Arm = Arm_by_indiv,cluster)
  if(nb_cov!=0){Data = data.frame(Data,Covariate)}
  if(nb_cov_clus!=0){Data = data.frame(Data,Covariate_clus)}
  Data = data.frame(Data,Risk_by_cluster = rep(Risk_vec,Suj_per_clus))
  
  return(Data)
}



# 2) Confusion function ----

#Aim -> Include confusion in the data by selection in function of the covariates (Method developed by Leyrat) 

#Input

# data        : Data set (generated by the 'Crea_data_V13_confusion()' function)
# k           : Number of cluster per arm
# Pi_int      : Proportion of the intervention group that we will include to add confusion
# Pi_con      : Proportion of the control group that we will include to add confusion 
# rho_z       : Intraclass correlation coefficient for inclusion
# nb_cov      : Number of individual level covariate
# nb_cov_clus : Number of cluster level covariate
# OR_int      : Odd ratio fixed 'a priori' to create our variable of selection in the intervention group
# OR_con      : Odd ratio fixed 'a priori' to create our variable of selection in the control group

#Output

#data_fin     : Data set with confusion created by the selection of the individuals.

#Aim -> Include confusion in the data by selection in function of the covariates (Method developed by Leyrat) 

#Input

# data        : Data set (generated by the 'Crea_data_V13_confusion()' function)
# k           : Number of cluster per arm
# Pi_int      : Proportion of the intervention group that we will include to add confusion
# Pi_con      : Proportion of the control group that we will include to add confusion 
# rho_z       : Intraclass correlation coefficient for inclusion
# nb_cov      : Number of individual level covariate
# nb_cov_clus : Number of cluster level covariate
# OR_int      : Odd ratio fixed 'a priori' to create our variable of selection in the intervention group
# OR_con      : Odd ratio fixed 'a priori' to create our variable of selection in the control group

#Output

#data_fin     : Data set with confusion created by the selection of the individuals.

fun_Data_confusion <- function(pI,pC,k,m,icc,nb_cov,p_cov,OR_cov,nb_cov_clus,p_cov_clus,OR_cov_clus,Pi_int,Pi_con,rho_z,OR_int,OR_con){
  
  
  data = Crea_Data_confusion(pI,pC,k,m,icc,nb_cov,p_cov,OR_cov,nb_cov_clus,p_cov_clus,OR_cov_clus,Pi_int,Pi_con)
  
  OR_int = OR_int[1:nb_cov]
  OR_con = OR_con[1:nb_cov]
  p_cov = p_cov[1:nb_cov]
  
  # create our parameter (a,b) of the variable Tau_ij following a beta distribution Beta(a,b) in each group
  a_int = Pi_int*(1-rho_z)/rho_z
  b_int = (1-Pi_int)*(1-rho_z)/rho_z
  
  a_con = Pi_con*(1-rho_z)/rho_z
  b_con = (1-Pi_con)*(1-rho_z)/rho_z
  
  
  # Calculate Mean and Variance of the variable Tau_ij in each group
  
  E_tau_ij_int = a_int/(a_int+b_int)
  v_tau_ij_int = a_int*b_int/((a_int+b_int)**2*(a_int+b_int+1))
  
  E_tau_ij_con = a_con/(a_con+b_con)
  v_tau_ij_con = a_con*b_con/((a_con+b_con)**2*(a_con+b_con+1))
  
  
  ## Coefficient beta from de variable L_ijl
  
  Beta_pi_int = log(OR_int)
  Beta_pi_con = log(OR_con)
  
  
  # Generate the inclusion rate in each cluster of each arm
  Tau_ij_int = rbeta(k,a_int,b_int)
  Tau_ij_con = rbeta(k,a_con,b_con)
  
  while (TRUE %in% (Tau_ij_int <= 0.1) |TRUE %in% (Tau_ij_con <= 0.1) ) { # We minimize the cluster effect for each cluster to 10% (rate to be include by clus)
    Tau_ij_int = rbeta(k,a_int,b_int)
    Tau_ij_con = rbeta(k,a_con,b_con)
  }
  
  
  
  # Cluster effect corrected
  
  eff_clus_ij_int = log(Tau_ij_int/(1-Tau_ij_int)) 
  eff_clus_ij_con = log(Tau_ij_con/(1-Tau_ij_con)) 
  
  # mean(eff_clus_ij_int)
  # mean(eff_clus_ij_con)
  
  
  # Cluster effect corrected
  
  
  eff_clus_ij = c(eff_clus_ij_int,eff_clus_ij_con)
  
  # Create a vector of length the number of individuals in the data set with the value of the cluster effect for each individual corresponding to their clusters 'appartenance'
  eff_clus_ij_vec = c()
  n = length(data$Outcome) # Total number of individual in the data set
  for (i in 1:(2*k)) {
    n_ij = length(which(data$cluster==i))
    eff_clus_ij_vec = c(eff_clus_ij_vec,rep(eff_clus_ij[i],n_ij))
  }
  # length(eff_clus_ij_vec)
  
  ## L function : corresponding to the regression function to generate a probability to the individual to be include in function of their covariates and cluster effect.  
  
  L_ijl = matrix(0,nrow = n,ncol = 1)
  
  if(nb_cov>1){
    for (i in 1:n) {
      if(data[i,]$Arm==1){
        for (j in 1:nb_cov) {
          L_ijl[i] = L_ijl[i] + Beta_pi_int[j]*(data[i,3+j]-p_cov[j]) 
        }
      }else{
        for (j in 1:nb_cov) {
          L_ijl[i] = L_ijl[i] + Beta_pi_con[j]*(data[i,3+j]-p_cov[j])
        }
      }
    }
  }
  
  L_ijl = L_ijl + eff_clus_ij_vec
  
  
  # Generate the probability for each individual to be include in the trial
  
  e_z_ijl = 1/(1+exp(-L_ijl))
  
  e_z_ijl_int = e_z_ijl[which(data$Arm==1)] 
  # mean(e_z_ijl_int)
  e_z_ijl_con = e_z_ijl[which(data$Arm==0)] 
  # mean(e_z_ijl_con)
  
  e_z_ijl = c(e_z_ijl_int,e_z_ijl_con)
  
  # e_z_ijl is a probability so is in [0,1]
  # So correction will be :
  max(e_z_ijl)
  min(e_z_ijl)
  
  index_need_correction = which(e_z_ijl > 1 )
  e_z_ijl[index_need_correction] = 1
  
  
  
  # Generate the inclusion variable Z for each individuals (Z_ijl = 1 -> include / Z_ijl = 0 -> No include)
  
  Z_ijl=c()
  for (i in 1:length(e_z_ijl)) {
    Z_ijl[i] = rbinom(1,1,e_z_ijl[i])
  }
  # length(Z_ijl)
  
  Suj_per_clus = aggregate(Z_ijl ~ cluster, data = data,FUN = sum)[-1]
  
  while (TRUE %in% (Suj_per_clus <= 1)) {
    Z_ijl=c()
    for (i in 1:length(e_z_ijl)) {
      Z_ijl[i] = rbinom(1,1,e_z_ijl[i])
    }
    
    Suj_per_clus = aggregate(Z_ijl ~ cluster, data = data,FUN = sum)[-1]
  }
  
  Inclus = which(Z_ijl==1)
  
  data_fin = data[Inclus,]
  
  # while(length(unique(data_fin$cluster)) != (2*k)){ # while conditions in case if one cluster did not include at least 1 individual 
  #   # Generate the inclusion variable Z for each individuals (Z_ijl = 1 -> include / Z_ijl = 0 -> No include)
  #   num_clus = c(1:(2*k))
  #   
  #   num_clus_manquant = num_clus[-unique(data_fin$cluster)]
  #   
  #   for (i in num_clus_manquant) {
  #     index = which(data$cluster==i)
  #     for (j in index) {
  #       Z_ijl[j] = rbinom(1,1,e_z_ijl[j])
  #     }
  #   }
  #   
  #   
  #   # Z_ijl=c()
  #   # for (i in 1:length(e_z_ijl)) {
  #   #   Z_ijl[i] = rbinom(1,1,e_z_ijl[i])
  #   # }
  #   # length(Z_ijl)
  #   Inclus = which(Z_ijl==1)
  #   
  #   data_fin = data[Inclus,]
  # }
  
  
  return(data_fin)
}

# 2b) ICS Confusion ----

fun_Data_ICS_confusion <- function(pI,pC,k_1,k_0,m,icc,nb_cov,OR_cov,nb_cov_clus,OR_cov_clus,alpha = 3,Pi_int,Pi_con,rho_z,OR_int,OR_con,ICS=1){
  
  
  ## Crea data ICS (Same function as Crea_data_ICS but directly implemented here)
  if(ICS == 1){
    data = Crea_Data_ICS_V2(pI,pC,
                            k_1,k_0,m,icc,
                            nb_cov,OR_cov,
                            nb_cov_clus,OR_cov_clus,
                            alpha=3)
  }
  if(ICS == 2){
    data = Crea_Data_ICS_2(pI,pC,
                           k_1,k_0,m,icc,
                           nb_cov,OR_cov,
                           nb_cov_clus,OR_cov_clus,
                           alpha=3)
  }
  
  
  
  OR_int = OR_int[1:nb_cov]
  OR_con = OR_con[1:nb_cov]
  p_cov = min(data$Risk_by_cluster)
  k = k_1 + k_0
  
  # create our parameter (a,b) of the variable Tau_ij following a beta distribution Beta(a,b) in each group
  a_int = Pi_int*(1-rho_z)/rho_z
  b_int = (1-Pi_int)*(1-rho_z)/rho_z
  
  a_con = Pi_con*(1-rho_z)/rho_z
  b_con = (1-Pi_con)*(1-rho_z)/rho_z
  
  
  # Calculate Mean and Variance of the variable Tau_ij in each group
  
  E_tau_ij_int = a_int/(a_int+b_int)
  v_tau_ij_int = a_int*b_int/((a_int+b_int)**2*(a_int+b_int+1))
  
  E_tau_ij_con = a_con/(a_con+b_con)
  v_tau_ij_con = a_con*b_con/((a_con+b_con)**2*(a_con+b_con+1))
  
  
  ## Coefficient beta from de variable L_ijl
  
  Beta_pi_int = log(OR_int)
  Beta_pi_con = log(OR_con)
  
  
  # Generate the inclusion rate in each cluster of each arm
  Tau_ij_int = rbeta(k_1,a_int,b_int)
  Tau_ij_con = rbeta(k_0,a_con,b_con)
  
  # Cluster effect corrected
  
  eff_clus_ij_int = log(Tau_ij_int/(1-Tau_ij_int)) 
  eff_clus_ij_con = log(Tau_ij_con/(1-Tau_ij_con)) 
  
  # mean(eff_clus_ij_int)
  # mean(eff_clus_ij_con)
  
  
  # Cluster effect corrected
  
  
  eff_clus_ij = c(eff_clus_ij_int,eff_clus_ij_con)
  
  # Create a vector of length the number of individuals in the data set with the value of the cluster effect for each individual corresponding to their clusters 'appartenance'
  eff_clus_ij_vec = c()
  n = length(data$Outcome) # Total number of individual in the data set
  for (i in 1:k) {
    n_ij = length(which(data$cluster==i))
    eff_clus_ij_vec = c(eff_clus_ij_vec,rep(eff_clus_ij[i],n_ij))
  }
  # length(eff_clus_ij_vec)
  
  ## L function : corresponding to the regression function to generate a probability to the individual to be include in function of their covariates and cluster effect.  
  
  L_ijl = matrix(0,nrow = n,ncol = 1)
  
  if(nb_cov>1){
    for (i in 1:n) {
      if(data[i,]$Arm==1){
        for (j in 1:nb_cov) {
          L_ijl[i] = L_ijl[i] + Beta_pi_int[j]*(data[i,3+j]-p_cov) 
        }
      }else{
        for (j in 1:nb_cov) {
          L_ijl[i] = L_ijl[i] + Beta_pi_con[j]*(data[i,3+j]-p_cov)
        }
      }
    }
  }
  
  L_ijl = L_ijl + eff_clus_ij_vec
  
  
  # Generate the probability for each individual to be include in the trial
  
  e_z_ijl = 1/(1+exp(-L_ijl))
  
  e_z_ijl_int = e_z_ijl[which(data$Arm==1)] 
  # mean(e_z_ijl_int)
  e_z_ijl_con = e_z_ijl[which(data$Arm==0)] 
  # mean(e_z_ijl_con)
  
  e_z_ijl = c(e_z_ijl_int,e_z_ijl_con)
  
  # e_z_ijl is a probability so is in [0,1]
  # So correction will be :
  max(e_z_ijl)
  min(e_z_ijl)
  
  index_need_correction = which(e_z_ijl > 1 )
  e_z_ijl[index_need_correction] = 1
  
  
  
  # Generate the inclusion variable Z for each individuals (Z_ijl = 1 -> include / Z_ijl = 0 -> No include)
  
  Z_ijl=c()
  for (i in 1:length(e_z_ijl)) {
    Z_ijl[i] = rbinom(1,1,e_z_ijl[i])
  }
  # length(Z_ijl)
  Inclus = which(Z_ijl==1)
  
  data_fin = data[Inclus,]
  
  while(length(unique(data_fin$cluster)) != k){ # while conditions in case if one cluster did not include at least 1 individual 
    # Generate the inclusion variable Z for each individuals (Z_ijl = 1 -> include / Z_ijl = 0 -> No include)
    num_clus = c(1:k)
    
    num_clus_manquant = num_clus[-unique(data_fin$cluster)]
    
    for (i in num_clus_manquant) {
      index = which(data$cluster==i)
      for (j in index) {
        Z_ijl[j] = rbinom(1,1,e_z_ijl[j])
      }
    }
    
    
    # Z_ijl=c()
    # for (i in 1:length(e_z_ijl)) {
    #   Z_ijl[i] = rbinom(1,1,e_z_ijl[i])
    # }
    # length(Z_ijl)
    Inclus = which(Z_ijl==1)
    
    data_fin = data[Inclus,]
  }
  
  
  return(data_fin)
}




# 3) Correlation verification functions ----

# corr_min_max function
#Input :
# p1_i  : prevalence of the outcome in the intervention arm
# p1_c  : prevalence of the outcome in the intervention arm
# p2    : prevalence of the covariate

#Output :
# vector : vector of length 2, with the minimum value and maximum value of the correlation acceptable between the outcome and the covariate taking account the 
#          difference of prevalence in the two groups. Following Prentice article.

corr_min_max <- function(p1_i,p1_c,p2){
  q1_i = 1-p1_i
  q1_c = 1-p1_c
  q2 = 1-p2
  
  # Calculating for the prevalence of the outcome equals to the prevalence of the intervention arm 
  max_cor_i = min(sqrt((p1_i*q2)/(q1_i*p2)),sqrt((p2*q1_i)/(q2*p1_i)))
  min_cor_i = max(-sqrt((q1_i*q2)/(p1_i*p2)),-sqrt((p1_i*p2)/(q1_i*q2)))
  
  
  # Calculating for the prevalence of the outcome equals to the prevalence of the control arm 
  max_cor_c = min(sqrt((p1_c*q2)/(q1_c*p2)),sqrt((p2*q1_c)/(q2*p1_c)))
  min_cor_c = max(-sqrt((q1_c*q2)/(p1_c*p2)),-sqrt((p1_c*p2)/(q1_c*q2)))
  
  # Looking the max and min value of the interval for both arm conditions 
  min_cor_tot = max(min_cor_c,min_cor_i)
  max_cor_tot = min(max_cor_c,max_cor_i)
  
  return(c(min_cor_tot,max_cor_tot))
  
}

# fin_rho function
#Inpout :
# p  : prevalence of the outcome
# q  : prevalence of the individual-level covariate
# OR : odd ratio between the outcome and the covariate

#Output :
# rho : coefficient of correlation needed to create the covariate 

find_rho <- function(p,q,OR){
  ## Fin rho by using this formula of OR:
  
  # A = q*r_m+q*p-q*p*r_m
  # B = q - r_m*q - p*q + p*q*r_m
  # C = p-p*q-q*r_m+p*q*r_m
  # D = 1-q+q*r_m-p+p*q-p*q*r_m
  # 
  # OR = (A*D)/(B*C)
  
  ## OR = num / den      ==>    OR * den - num = 0
  ## with:
  # num = r_m² * num1 + r_m * num2 + num3
  # den = r_m² * den1 + r_m * den2 + den3
  # r_m = rho * coeff_rm
  
  num1 = q**2 - 2*p*(q**2) + (p**2)*(q**2)
  
  num2 = q - q**2 - 2*p*q + 3*p*(q**2) - 2*(p**2)*(q**2) + (p**2)*q
  
  num3 = p*q - p*(q**2) - (p**2)*q + (p**2)*(q**2)
  
  den1 = q**2 - 2*p*(q**2) + (p**2)*(q**2)
  
  den2 = - q**2 - p*q + 3*p*(q**2) - 2*(p**2)*(q**2) + (p**2)*q
  
  den3 = p*q - p*(q**2) - (p**2)*q + (p**2)*(q**2)
  
  coeff_rm = sqrt(1-q)/sqrt((q/p)*(1-p))
  
  # coefficient of the quadratic equation : a, b et c
  
  a = (coeff_rm**2) * (den1*OR - num1)
  b = coeff_rm * (den2*OR - num2)
  c = den3*OR - num3
  
  # delta of the quadratic equation 
  
  delta = b**2-4*a*c
  
  # 2 solutions of the quadratic equation
  
  x1 = (-b + sqrt(delta))/(2*a)
  x2 = (-b - sqrt(delta))/(2*a)
  
  rho = c(x1,x2)
  # correlation coefficient in the interval [-1;1]
  rho = rho[which(rho <= 1 & rho >= -1)]
  return(rho)
}

# find_r_clus
#Inpout :
# p     : prevalence of the outcome (corresponding to the prevalence of the outcome in each arm so different between two individual of two different arm intervention)
# q     : prevalence of the cluster-level covariate 
# OR    : Odd ratio between the cluster-level covariate and the outcome

#Output :
# r : Value of the parameter of the Rij variable to simulate the cluster-level covariate X_ij

find_r_clus <- function(p,q,OR_clus,icc){
  ## Find r by resolving this equation:
  
  # A = p**2*alpha +p*alpha*r*sqrt(icc) + p**2*alpha*r*sqrt(icc)
  # B = p*alpha - p*alpha*r*sqrt(icc) - p**2*alpha + p**2*alpha*r*sqrt(icc)
  # C = p - p*alpha*r*sqrt(icc) - p**2*alpha + p**2*alpha*r*sqrt(icc)
  # D = 1 - p - p*alpha + p*alpha*r*sqrt(icc) + p**2*alpha - p**2*alpha*r*sqrt(icc)
  # 
  # OR = (A*D)/(B*C)
  
  ## OR = num / den      ==>    OR * den - num = 0
  ## with:
  # num = r_m² * num1 + r_m * num2 + num3
  # den = r_m² * den1 + r_m * den2 + den3
  
  
  num1 = icc*(q**2) - 2*p*(q**2)*icc + (p**2)*(q**2)*icc
  
  num2 = q*sqrt(icc) - sqrt(icc)*(q**2) - 2*p*q*sqrt(icc) + 3*p*(q**2)*sqrt(icc) + (p**2)*q*sqrt(icc) - 2*(p**2)*(q**2)*sqrt(icc) 
  
  num3 = p*q - p*(q**2) - (p**2)*q + (p**2)*(q**2)
  
  den1 = icc * (q**2) - 2 * p * (q**2) * icc + (p**2)*(q**2)*icc
  
  den2 = - (q**2)*sqrt(icc) - p*q*sqrt(icc) + 3*p*(q**2)*sqrt(icc) + (p**2)*q*sqrt(icc) - 2*(p**2)*(q**2)*sqrt(icc)
  
  den3 = p*q - p*(q**2) - (p**2)*q + (p**2)*(q**2)
  
  
  # coefficient of the quadratic equation: a, b et c
  
  a = (den1*OR_clus - num1)
  b = (den2*OR_clus - num2)
  c = den3*OR_clus - num3
  
  # delta of the quadratic equation
  
  delta = b**2-4*a*c
  
  # 2 solutions of the quadratic equation
  if(delta>0){
    x1 = (-b + sqrt(delta))/(2*a)
    x2 = (-b - sqrt(delta))/(2*a)
    x = c(x1,x2)
  }
  if(delta==0){x = (-b)/(2*a)}
  if(delta<0){error(stop('No solution to the equation '))}
  
  # correlation coefficient in the interval [-1;1]
  
  x = x[which(x <= 1 & x >= 0)]
  if(length(x)==0){stop('No Solution in [0,1]')}
  return(x)
}





# 4) Function to save data ----
fun_para_save_data_clus <- function(itt,itt_para,n,Scenario,Pi_int,Pi_con,rho_z,OR_int,OR_con,Data_itt_File){
  
  
  for (i in 1:itt) {
    
    # Taking the differrent parameters set to the scenario number n
    
    pI = as.numeric(Scenario[1])                      # Prevalence in intervention arm
    pC = as.numeric(Scenario[2])                      # Prevalence in control arm
    k = as.numeric(Scenario[3])                       # Number of cluster per arm
    m = as.numeric(Scenario[4])                       # Mean cluster size
    icc = as.numeric(Scenario[5])                     # Intracluster Correlation Coefficient
    nb_cov = as.numeric(Scenario[6])                  # Number of individual level covariate
    nb_cov_clus = as.numeric(Scenario[7])             # Number of cluster level covariate
    p_cov = as.vector(unlist(Scenario[8]))            # Prevalence of the individual level covariate
    p_cov_clus = as.vector(unlist(Scenario[9]))       # Prevalence of the cluster level covariate
    OR_cov = as.vector(unlist(Scenario[10]))          # OR to measure the association between individual level covariate and the outcome
    OR_cov_clus = as.vector(unlist(Scenario[11]))     # OR to measure the association between cluster level covariate and the outcome
    
    # Depending if we have an individual level covariate or not. Because if not we do not have to include confounding, so we just have to use the regular version or the data generating function
    data = fun_Data_confusion(pI,pC,k,m,icc,nb_cov,p_cov,OR_cov,nb_cov_clus,p_cov_clus,OR_cov_clus,Pi_int,Pi_con,rho_z,OR_int,OR_con)
    
    # Saving the data set to an Excel file (.csv)
    write.csv2(data,here::here(paste(Data_itt_File,"/Scenario_",n,sep=""),paste("/Data_itt_",sep = "",itt_para,"_scen_",n,".csv")),row.names = FALSE)
  }
  
}

fun_para_save_data_clus_ICS <- function(itt,itt_para,n,Scenario,Pi_int,Pi_con,rho_z,OR_int,OR_con,alpha = 3,ICS=1,Data_itt_File){
  
  
  for (i in 1:itt) {
    
    # Taking the differrent parameters set to the scenario number n
    
    pI = as.numeric(Scenario[1])                      # Prevalence in intervention arm
    pC = as.numeric(Scenario[2])                      # Prevalence in control arm
    k = as.numeric(Scenario[3])                       # Number of cluster per arm
    m = as.numeric(Scenario[4])                       # Mean cluster size
    icc = as.numeric(Scenario[5])                     # Intracluster Correlation Coefficient
    nb_cov = as.numeric(Scenario[6])                  # Number of individual level covariate
    nb_cov_clus = as.numeric(Scenario[7])             # Number of cluster level covariate
    p_cov = as.vector(unlist(Scenario[8]))            # Prevalence of the individual level covariate
    p_cov_clus = as.vector(unlist(Scenario[9]))       # Prevalence of the cluster level covariate
    OR_cov = as.vector(unlist(Scenario[10]))          # OR to measure the association between individual level covariate and the outcome
    OR_cov_clus = as.vector(unlist(Scenario[11]))     # OR to measure the association between cluster level covariate and the outcome
    
    # Depending if we have an individual level covariate or not. Because if not we do not have to include confounding, so we just have to use the regular version or the data generating function
    # data = fun_Data_ICS_confusion(pI,pC,k_1=k,k_0=k,m=(m/0.7),icc,nb_cov,OR_cov,nb_cov_clus,OR_cov_clus,alpha = 3,Pi_int,Pi_con,rho_z,OR_int,OR_con,ICS=1)
    data = Crea_Data_ICS_V2(pI,pC,k_1=k,k_0=k,m,icc,nb_cov,OR_cov,nb_cov_clus,OR_cov_clus,alpha = 3)
    
    # Saving the data set to an Excel file (.csv)
    write.csv2(data,here::here(paste(Data_itt_File,"/Scenario_",n,sep=""),paste("Data_itt_",sep = "",itt_para,"_scen_",n,".csv")),row.names = FALSE)
  }
  
}



# 5) Fun for formula_tsp_fun ----
formula_tsp_fun <- function(nb_cov_indiv,nb_cov_clus){
  
  
  if(nb_cov_indiv > 0){
    Covariable = 'IC_1'
    if(nb_cov_indiv > 1){
      for (i in 2:nb_cov_indiv) {
        a = paste('IC_',i,sep = "")
        Covariable = paste(Covariable,a,sep = '+')
      }
    }
  }
  if(nb_cov_clus > 0){
    for (i in 1:nb_cov_clus) {
      a = paste('CC_',i,sep = "")
      Covariable = paste(Covariable,a,sep = '+')
    }
  }
  
  
  Covariable = as.factor(Covariable)
  Formula = as.formula(paste("Outcome ~ ",Covariable))
  
  return(Formula)
}

# 6) fun Two way procedure ----
CL_tsp_fun <- function(data, family = "binomial",Scenario,itt_para){
  
  ## Need scenario for the number of covariate for formula and confidence interval
  True_RD =  as.numeric(Scenario[1]) - as.numeric(Scenario[2]) 
  
  nb_cov_indiv = as.numeric(Scenario[6])
  nb_cov_clus = as.numeric(Scenario[7])
  nb_cov = nb_cov_indiv + nb_cov_clus
  
  formula_tsp = formula_tsp_fun(nb_cov_indiv,nb_cov_clus)
  
  if(family == "binomial"){
    mod_glm_without_inter = glm(formula = formula_tsp,
                                data = data,
                                family = binomial)
  } else if(family == "gaussian"){
    mod_glm_without_inter = glm(formula = formula_tsp,
                                data = data,
                                family = gaussian)
  } else if(family == "poisson"){
    mod_glm_without_inter = glm(formula = formula_tsp,
                                data = data,
                                family = poisson)
  }else{stop("ERROR IN FAMILY USED : should be a character type variable and/or choice between: binomial, gaussian and poisson")}
  
  #Number of cluster
  K = max(unique(data$cluster))
  
  #Stockage of Difference residual and Arm label for each cluster
  R_d_clus = matrix(NA,nrow = 1,ncol = K)
  Arm_mat = matrix(NA,nrow = 1,ncol = K)
  
  for (j in 1:K) {
    # Taking only the data from one cluster
    data_clus_j = data[which(data$cluster==j),]
    
    # Sum of the residuals depending on the family used in the regression ---> 'Expected' number of events in cluster j 10.5.1.2 chapter in Chapman book
    if(family == "binomial"){
      e_j = sum(invlogit(predict(mod_glm_without_inter,newdata = data_clus_j)))
    } else if(family == "gaussian"){
      e_j = sum(predict(mod_glm_without_inter,newdata = data_clus_j))
    } else if(family == "poisson"){
      e_j = sum(exp(predict(mod_glm_without_inter,newdata = data_clus_j)))
    }
    
    # d_j is the observed number of events in the cluster
    # m_j is the number of individuals in the cluster j
    d_j = sum(data_clus_j$Outcome)
    m_j = length(data_clus_j$Outcome)
    
    R_d_clus[j] = (d_j-e_j)/m_j
    Arm_mat[j] = data_clus_j[1,]$Arm
  }
  
  # Risk difference
  RD = mean(R_d_clus[which(Arm_mat==1)]) - mean(R_d_clus[which(Arm_mat==0)])
  
  # Confidence interval
  s_square = sum(c((R_d_clus[which(Arm_mat==1)] - mean(R_d_clus[which(Arm_mat==1)]))**2,
                   (R_d_clus[which(Arm_mat==0)] - mean(R_d_clus[which(Arm_mat==0)]))**2))/(K-2)
  
  SE_RD = sqrt(s_square)*sqrt(1/length(which(Arm_mat==1))+1/length(which(Arm_mat==0)))
  
  LL_95 = RD - qt(0.975,K-2-nb_cov_clus)*SE_RD   #degree of freedom following Hayes' book
  UL_95 = RD + qt(0.975,K-2-nb_cov_clus)*SE_RD
  
  # LL_95 = RD - qnorm(0.975)*SE_RD
  # UL_95 = RD + qnorm(0.975)*SE_RD
  
  RD_CI_res = tibble(RD_CL = RD,LowerLimit_95_CL = LL_95,UpperLimit_95_CL = UL_95,SE_CL = SE_RD , Bias = RD - True_RD, Relative_Bias = (RD-True_RD)/True_RD*100,Itt_para = itt_para)
  return(RD_CI_res)
}

# 7) fun un-weighted method cluster level analysis ----
CL_UW_fun <- function(data,itt_para,Scenario){
  
  True_RD = as.numeric(Scenario[1]) - as.numeric(Scenario[2]) 
  
  Prev_clus <- data %>% group_by(cluster)%>%
    summarise(Arm = mean(Arm),p_j = mean(Outcome))
  Prev_arm <- Prev_clus %>% group_by(Arm)%>%
    summarise(p_i = mean(p_j))
  
  # Number of cluster
  K = length(Prev_clus$cluster)
  k_1 = length(Prev_clus[which(Prev_clus$Arm==1),]$cluster)
  k_0 = length(Prev_clus[which(Prev_clus$Arm==0),]$cluster)
  Y_1j = Prev_clus[which(Prev_clus$Arm==1),]$p_j
  Y_0j = Prev_clus[which(Prev_clus$Arm==0),]$p_j
  
  Y_1 = Prev_arm[which(Prev_arm$Arm==1),]$p_i
  Y_0 = Prev_arm[which(Prev_arm$Arm==0),]$p_i
  
  # Risk difference CL analysis unweigthed
  RD = Prev_arm[which(Prev_arm==1),]$p_i - Prev_arm[which(Prev_arm==0),]$p_i
  
  # Confidence Interval 
  s_square = sum(c(Y_1j - Y_1,Y_0j - Y_0)**2)/(K-2)
  
  SE_CL_UW = sqrt(s_square) * sqrt(1/k_1 + 1/k_0)
  
  LL = RD - qt(0.975,K-2) * SE_CL_UW
  UL = RD + qt(0.975,K-2) * SE_CL_UW
  
  Tab_res = data.frame(RD_CL_UW = RD,
                       LL_95_CL_UW = LL,
                       UL_95_CL_UW = UL,
                       SE_CL_UW = SE_CL_UW,
                       Bias = RD - True_RD,
                       Relative_Bias = (RD-True_RD)/True_RD*100,
                       Itt = itt_para)
  return(Tab_res)
}

# 8) fun for parallelism analysis----
fun_para_analyse_cluster <- function(itt_para,Scenario,n,Data_file,Resu_file){
  
  # Recuperation Data
  data = read.csv(paste(paste(Data_file,"/Scenario_",n,sep=""),sep ="",paste("/Data_itt_",sep = "",itt_para,"_scen_",n,".csv")), sep=";")
  
  ## Analyse all method
  
  res_CL_UW = CL_UW_fun(data,itt_para = itt_para,Scenario = Scenario)
  
  res_CL_TSP_Bin = CL_tsp_fun(data = data,family = "binomial",Scenario = Scenario,itt_para = itt_para)
  
  res_CL_TSP_Gauss = CL_tsp_fun(data = data,family = "gaussian",Scenario = Scenario,itt_para = itt_para)
  
  # res_CL_Gcomp
  
  
  # Remplissage CSV resultat pour chacune des méthodes
  
  write.table(res_CL_UW,file = paste(Resu_file,"/UW_analyse/",paste("/Data_output_CL_UW_Scenario_",sep = "",n,".csv"),sep = ""),
              append = TRUE,row.names = FALSE,col.names = FALSE,sep = ";",dec = ",")
  
  write.table(res_CL_TSP_Bin,file = paste(Resu_file,"/TSP_Bin_analyse/",paste("/Data_output_CL_TSP_Bin_Scenario_",sep = "",n,".csv"),sep = ""),
              append = TRUE,row.names = FALSE,col.names = FALSE,sep = ";",dec = ",")
  
  write.table(res_CL_TSP_Gauss,file = paste(Resu_file,"/TSP_Gauss_analyse/",paste("/Data_output_CL_TSP_Gauss_Scenario_",sep = "",n,".csv"),sep = ""),
              append = TRUE,row.names = FALSE,col.names = FALSE,sep = ";",dec = ",")
  
}

fun_para_analyse_cluster_gcomp_indv <- function(itt_para,Scenario,n,Data_file,Resu_file){
  
  # Recuperation Data
  data = read.csv(paste(paste(Data_file,"/Scenario_",n,sep=""),sep ="",paste("/Data_itt_",sep = "",itt_para,"_scen_",n,".csv")), sep=";")
  
  ## Analyse all method
  
  res_CL_GC_Bin = fun_g_comp_bin(data,itt_para = itt_para,Scenario = Scenario,cov_schem = "individual")
  
  
  # Remplissage CSV resultat pour chacune des méthodes
  
  write.table(res_CL_GC_Bin,file = paste(Resu_file,"/GC_Bin_analyse_indv",paste("/Data_output_CL_GC_indv_Bin_Scenario_",sep = "",n,".csv"),sep = ""),
              append = TRUE,row.names = FALSE,col.names = FALSE,sep = ";",dec = ",")
  
}


fun_para_analyse_cluster_gcomp_both <- function(itt_para,Scenario,n,Data_file,Resu_file){
  
  # Recuperation Data
  data = read.csv(paste(paste(Data_file,"/Scenario_",n,sep=""),sep ="",paste("/Data_itt_",sep = "",itt_para,"_scen_",n,".csv")), sep=";")
  
  ## Analyse all method
  
  res_CL_GC_Bin = fun_g_comp_bin(data,itt_para = itt_para,Scenario = Scenario,cov_schem = "both")
  
  
  # Remplissage CSV resultat pour chacune des méthodes
  
  write.table(res_CL_GC_Bin,file = paste(Resu_file,"/GC_Bin_analyse_both",paste("/Data_output_CL_GC_both_Bin_Scenario_",sep = "",n,".csv"),sep = ""),
              append = TRUE,row.names = FALSE,col.names = FALSE,sep = ";",dec = ",")
  
}

fun_para_analyse_cluster_gcomp_GLMMPQL <- function(itt_para,Scenario,n,Data_file,Resu_file){
  
  # Recuperation Data
  data = read.csv(paste(paste(Data_file,"/Scenario_",n,sep=""),sep ="",paste("/Data_itt_",sep = "",itt_para,"_scen_",n,".csv")), sep=";")
  
  ## Analyse all method
  
  res_CL_GC_Bin_glmm_indv = fun_g_comp_bin_glmer(data,itt_para = itt_para,Scenario = Scenario,cov_schem = "individual")
  
  res_CL_GC_Bin_glmm_both = fun_g_comp_bin_glmer(data,itt_para = itt_para,Scenario = Scenario,cov_schem = "both")
  
  
  # Remplissage CSV resultat pour chacune des méthodes
  
  write.table(res_CL_GC_Bin_glmm_indv,file = paste(Resu_file,"/GC_GLMM_Bin_analyse_indv",paste("/Data_output_CL_GC_indv_Bin_glmm_Scenario_",sep = "",n,".csv"),sep = ""),
              append = TRUE,row.names = FALSE,col.names = FALSE,sep = ";",dec = ",")
  
  write.table(res_CL_GC_Bin_glmm_both,file = paste(Resu_file,"/GC_GLMM_Bin_analyse_both",paste("/Data_output_CL_GC_both_Bin_glmm_Scenario_",sep = "",n,".csv"),sep = ""),
              append = TRUE,row.names = FALSE,col.names = FALSE,sep = ";",dec = ",")
  
}
# 9) Fun correction for gcomp due to gee infinity laps time running ----

fun_cor_gcomp <- function(i,Time,n,Base_file,Workspace_name,Data_file,Resu_file,cov_schem = "individual"){
  
  # 0. Espace de travail
  Chemin = sprintf(paste(Base_file,"/correction_gcomp/correction_gcomp_S%d",sep = ""),n)
  dir.create(Chemin)
  setwd(Chemin)
  
  # 0.1 Nom du fichier principal
  FileName = paste("correction_gcomp_S",n,"_itt_%d.R",sep = "")
  
  
  # 1. tu crees le fichier .R qui va contenir les instructions à lancer dans ta 2eme session
  
  
  Stock_file = sprintf(paste(Base_file,"/correction_g_comp/correction_gcomp_S%d",sep = ""),n)
  
  writeLines('rm(list = ls())',
             con = sprintf(FileName,i))
  
  
  write(sprintf("Chemin  = '%s'",Stock_file),
        file  = sprintf(FileName, i),
        append = TRUE)
  
  # 1.1 tu exporteras l'id de la 2eme session pour le recuperer. Ca servira pour le killer par la suite si besoin est.
  Exp_id = paste("writeLines(as.character(Sys.getpid()), con = sprintf('%s/pid_%d.txt', getwd(),",i,"))",sep = "")
  write(Exp_id,
        file = sprintf(FileName, i),
        append = TRUE)
  
  # 1.2 les lignes de commande, pour lancer le programme qui prend du temps
  write(paste("load('",Base_file,"/",Workspace_name,"')",sep = ""), file = sprintf(FileName, i), append = TRUE)
  
  
  write("library(cli, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(Matrix, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(MASS, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(lme4, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(backports, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(dplyr, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(parallel, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(iterators, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(rngtools, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(foreach, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(doRNG, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(doParallel, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(dplyr, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(gee, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(geepack, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(spind, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(doBy, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(arm, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(here, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(geesmv, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(matrixcalc, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  
  L1 = paste("n =",n)
  write(L1,
        file = sprintf(FileName,i),
        append = TRUE)
  
  write(paste("Data_file = '",Data_file,"'",sep = ""), file = sprintf(FileName, i), append = TRUE)
  write(paste("Resu_file = '",Resu_file,"'",sep = ""), file = sprintf(FileName, i), append = TRUE)
  
  write(sprintf("Scenario_use = Scenario_%d",n),
        file = sprintf(FileName, i),
        append = TRUE)
  
  if(cov_schem == "individual"){
    write(sprintf("fun_para_analyse_cluster_gcomp_indv(itt_para=%d,Scenario=Scenario_use,n=n,Data_file=Data_file,Resu_file = Resu_file)",i),
          file = sprintf(FileName, i),
          append = TRUE)
  }
  if(cov_schem == "both"){
    write(sprintf("fun_para_analyse_cluster_gcomp_both(itt_para=%d,Scenario=Scenario_use,n=n,Data_file=Data_file,Resu_file = Resu_file)",i),
          file = sprintf(FileName, i),
          append = TRUE)
  }
  
  
  write("DONE = TRUE",file = sprintf(FileName,i),append = TRUE)
  
  # 1.3 l'export du resultat
  tace_2 = paste("save(DONE, file = sprintf('res_%d.Rdata',",i,"))",sep = "")
  write(tace_2, file = sprintf(FileName, i), append = TRUE)
  
  # 2. tu fais executer les commandes des fichiers correction_S33_itt_x.R en batch en arriere plan (wait = F)
  cmd_batch = paste("R CMD BATCH --no-restore correction_gcomp_S",n,"_itt_%d.R",sep = "")
  
  Start.T = Sys.time()
  system(sprintf(cmd_batch, i), wait = TRUE,timeout = Time)
  End.T = Sys.time()
  
  Diff.time = difftime(End.T,Start.T,units = "sec")
  if(Diff.time >= Time ){
    system(sprintf("tskill %d", scan(sprintf("pid_%d.txt", i), integer(),quiet = T)))
  }
  
}


# 10) Fun for "G-computation" method analysis ----

fun_g_comp_bin <- function(data,Scenario,itt_para,cov_schem = "individual"){
  
  ## Need scenario for the number of covariate for formula and confidence interval
  True_RD =  as.numeric(Scenario[1]) - as.numeric(Scenario[2])
  
  nb_cov_indiv = as.numeric(Scenario[6])
  nb_cov_clus = as.numeric(Scenario[7])
  nb_cov = nb_cov_indiv + nb_cov_clus
  
  if(cov_schem == "individual"){
    formula_gee <- fun_formula_for_gee(nb_cov = nb_cov_indiv,nb_cov_clus = 0)
  }
  if(cov_schem == "cluster"){
    formula_gee <- fun_formula_for_gee(nb_cov = 0,nb_cov_clus = nb_cov_clus)
  }
  if(cov_schem == "both"){
    formula_gee <- fun_formula_for_gee(nb_cov_indiv,nb_cov_clus)
  }
  
  
  cap_op2  <- capture.output(suppressMessages(mod_gee <- try(gee(formula = formula_gee,id = cluster,data = data,family = binomial,corstr = "exchangeable"))))
  
  if(inherits(mod_gee,"try-error")==FALSE){
    data_1 <- data_0 <- data
    
    data_1$Arm <- 1
    data_0$Arm <- 0
    
    # Y_1 <- invlogit(Pred_Trt_fun(mod_gee,data = data_1,nb_cov_tot = nb_cov))
    # 
    # Y_0 <- invlogit(Pred_NoTrt_fun(mod_gee,data = data_0,nb_cov_tot = nb_cov))
    if(cov_schem == "individual"){
      Y_1 <- invlogit(Pred_Trt_fun(mod_gee,data = data_1,nb_cov_tot = nb_cov_indiv))
      
      Y_0 <- invlogit(Pred_NoTrt_fun(mod_gee,data = data_0,nb_cov_tot = nb_cov_indiv))
    }
    if(cov_schem == "cluster"){
      Y_1 <- invlogit(Pred_Trt_fun(mod_gee,data = data_1,nb_cov_tot = nb_cov_clus))
      
      Y_0 <- invlogit(Pred_NoTrt_fun(mod_gee,data = data_0,nb_cov_tot = nb_cov_clus))
    }
    if(cov_schem == "both"){
      Y_1 <- invlogit(Pred_Trt_fun(mod_gee,data = data_1,nb_cov_tot = nb_cov))
      
      Y_0 <- invlogit(Pred_NoTrt_fun(mod_gee,data = data_0,nb_cov_tot = nb_cov))
    }
    
    
    # Pred_indiv <- Pred_Trt - Pred_NoTrt
    
    data_gcomp <- data.frame(data,Y_1,Y_0)
    
    Prev_clus <- data_gcomp %>% group_by(cluster)%>%
      summarise(Arm=mean(Arm),Y_1j = mean(Y_1),Y_0j = mean(Y_0))
    
    # Number of cluster
    K = length(Prev_clus$cluster)
    
    # RD by cluster
    Prev_clus$RD_ij = Prev_clus$Y_1j - Prev_clus$Y_0j 
    
    # Risk difference is the mean of all RD estimated per cluster
    RD = mean(Prev_clus$RD_ij)
    
    # Confidence Interval
    # Matrice_RobustVariance = mod_gee$robust.variance
    
    # Matrice variance covariance corrected with Fay and Graubard Method
    Mat_corr_FG_varcov = fun_var_corrected_FG(mod_gee,formula_gee,id = "cluster",family=binomial,data,corstr="exchangeable",b=0.75)
    
    Var_per_clus = c()
    
    for (i in 1:K) {
      
      if(cov_schem == "individual"){
        Var_per_clus[i] = fun_deriv_logit_tot(mod_gee,nb_cov_indiv,0,data[which(data$cluster==i),]) %*% Mat_corr_FG_varcov %*% fun_deriv_logit_tot(mod_gee,nb_cov_indiv,0,data[which(data$cluster==i),])
      }
      if(cov_schem == "cluster"){
        Var_per_clus[i] = fun_deriv_logit_tot(mod_gee,0,nb_cov_clus,data[which(data$cluster==i),]) %*% Mat_corr_FG_varcov %*% fun_deriv_logit_tot(mod_gee,0,nb_cov_clus,data[which(data$cluster==i),])
      }
      if(cov_schem == "both"){
        Var_per_clus[i] = fun_deriv_logit_tot(mod_gee,nb_cov_indiv,nb_cov_clus,data[which(data$cluster==i),]) %*% Mat_corr_FG_varcov %*% fun_deriv_logit_tot(mod_gee,nb_cov_indiv,nb_cov_clus,data[which(data$cluster==i),])
      }
      
    }
    
    Var_w = (1/K) * sum(Var_per_clus)
    
    # Var_b1 = var(Prev_clus$Y_1j-Prev_clus$Y_0j)
    # Var_b2 = sum(c((Prev_clus[1:5,]$RD_ij-mean(Prev_clus[1:5,]$RD_ij))**2,
    #                (Prev_clus[6:10,]$RD_ij-mean(Prev_clus[6:10,]$RD_ij))**2))/(k-1-2)
    # 
    # Var_b3 = sum((Prev_clus$RD_ij-mean(Prev_clus$RD_ij))**2)/(k-1)
    # var(Prev_clus$RD_ij)
    
    Var_b = sum(c((Prev_clus[Prev_clus$Arm==1,]$RD_ij-mean(Prev_clus[Prev_clus$Arm==1,]$RD_ij))**2,
                  (Prev_clus[Prev_clus$Arm==0,]$RD_ij-mean(Prev_clus[Prev_clus$Arm==0,]$RD_ij))**2))/(K-2)
    
    Var_tot = Var_w + (K +1)/K * Var_b
    
    SE_CL_GC_Bin = sqrt(Var_tot)
    
    LL = RD - qt(0.975,K-1) * SE_CL_GC_Bin
    UL = RD + qt(0.975,K-1) * SE_CL_GC_Bin
    
    # if(cov_schem == "individual"){
    #   LL = RD - qt(0.975,k-nb_cov_indiv) * SE_CL_GC_Bin
    #   UL = RD + qt(0.975,k-nb_cov_indiv) * SE_CL_GC_Bin
    # }
    # if(cov_schem == "cluster"){
    #   LL = RD - qt(0.975,k-nb_cov_clus) * SE_CL_GC_Bin
    #   UL = RD + qt(0.975,k-nb_cov_clus) * SE_CL_GC_Bin
    # }
    # if(cov_schem == "both"){
    #   LL = RD - qt(0.975,k-nb_cov) * SE_CL_GC_Bin
    #   UL = RD + qt(0.975,k-nb_cov) * SE_CL_GC_Bin
    # }
    
    
    if(mod_gee$error != 0){
      Tab_res = data.frame(RD_CL_GC_Bin = NA,
                           LL_95_CL_GC_Bin = NA,
                           UL_95_CL_GC_Bin = NA,
                           SE_CL_GC_Bin = NA,
                           Bias = NA,
                           Relative_Bias = NA,
                           Itt = itt_para,
                           Error = mod_gee$error)
    }else{
      Tab_res = data.frame(RD_CL_GC_Bin = RD,
                           LL_95_CL_GC_Bin = LL,
                           UL_95_CL_GC_Bin = UL,
                           SE_CL_GC_Bin = SE_CL_GC_Bin,
                           Bias = RD - True_RD,
                           Relative_Bias = (RD-True_RD)/True_RD*100,
                           Itt = itt_para,Error = mod_gee$error)
    }
    
  }else{
    Tab_res = data.frame(RD_CL_GC_Bin = NA,
                         LL_95_CL_GC_Bin = NA,
                         UL_95_CL_GC_Bin = NA,
                         SE_CL_GC_Bin = NA,
                         Bias = NA,
                         Relative_Bias = NA,
                         Itt = itt_para,
                         Error = 1)
  }
  
  return(Tab_res)
  
}


fun_g_comp_bin_glmer <- function(data,Scenario,itt_para,cov_schem = "individual"){
  
  ## Need scenario for the number of covariate for formula and confidence interval
  True_RD =  as.numeric(Scenario[1]) - as.numeric(Scenario[2])
  
  nb_cov_indiv = as.numeric(Scenario[6])
  nb_cov_clus = as.numeric(Scenario[7])
  nb_cov = nb_cov_indiv + nb_cov_clus
  
  if(cov_schem == "individual"){
    formula_glmer <- fun_formula_for_glmer(nb_cov = nb_cov_indiv,nb_cov_clus = 0)
  }
  if(cov_schem == "cluster"){
    formula_glmer <- fun_formula_for_glmer(nb_cov = 0,nb_cov_clus = nb_cov_clus)
  }
  if(cov_schem == "both"){
    formula_glmer <- fun_formula_for_glmer(nb_cov_indiv,nb_cov_clus)
  }
  
  
  cap_op2  <- capture.output(suppressMessages(mod_glmm <- try(mod_glmm <- glmmPQL(fixed = fun_formula_for_gee(nb_cov_indiv,nb_cov_clus),random = ~ 1|cluster,family = binomial,data = data))))
  
  if(inherits(mod_glmm,"try-error")==FALSE){
    
    fixed = mod_glmm$coefficients$fixed
    random = mod_glmm$coefficients$random$cluster
    
    
    if(cov_schem == "individual"){
      Y_1 <- invlogit(fixed["(Intercept)"]+ fixed["Arm"] + fixed["IC_1"] * data$IC_1 + fixed["IC_2"] * data$IC_2 + fixed["IC_3"] * data$IC_3 + random[data$cluster])
      Y_0 <- invlogit(fixed["(Intercept)"]+ fixed["IC_1"] * data$IC_1 + fixed["IC_2"] * data$IC_2 + fixed["IC_3"] * data$IC_3 + random[data$cluster])
    }
    if(cov_schem == "cluster"){
      Y_1 <- invlogit(fixed["(Intercept)"]+ fixed["Arm"] + fixed["CC_1"] * data$CC_1+ fixed["CC_2"] * data$CC_2 + random[data$cluster])
      Y_0 <- invlogit(fixed["(Intercept)"]+ fixed["CC_1"] * data$CC_1+ fixed["CC_2"] * data$CC_2 + random[data$cluster])
    }
    if(cov_schem == "both"){
      Y_1 <- invlogit(fixed["(Intercept)"]+ fixed["Arm"] + fixed["IC_1"] * data$IC_1 + fixed["IC_2"] * data$IC_2 + fixed["IC_3"] * data$IC_3+ fixed["CC_1"] * data$CC_1+ fixed["CC_2"] * data$CC_2 + random[data$cluster])
      Y_0 <- invlogit(fixed["(Intercept)"]+ fixed["IC_1"] * data$IC_1 + fixed["IC_2"] * data$IC_2 + fixed["IC_3"] * data$IC_3+ fixed["CC_1"] * data$CC_1+ fixed["CC_2"] * data$CC_2 + random[data$cluster])
    }
    
    
    # Pred_indiv <- Pred_Trt - Pred_NoTrt
    
    data_gcomp <- data.frame(data,Y_1,Y_0)
    
    Prev_clus <- data_gcomp %>% group_by(cluster)%>%
      summarise(Arm=mean(Arm),Y_1j = mean(Y_1),Y_0j = mean(Y_0))
    
    # Number of cluster
    K = length(Prev_clus$cluster)
    
    # RD by cluster
    Prev_clus$RD_ij = Prev_clus$Y_1j - Prev_clus$Y_0j 
    
    # Risk difference is the mean of all RD estimated per cluster
    RD = mean(Prev_clus$RD_ij)
    
    # Confidence Interval
    
    Var_per_clus = c()
    
    for (i in 1:K) {
      
      if(cov_schem == "individual"){
        Var_per_clus[i] = fun_se_delta_meth_logit_glmmPQL(mod_glmm,nb_cov_indiv,0,data[which(data$cluster==i),],adjust = cov_schem)**2
      }
      if(cov_schem == "cluster"){
        Var_per_clus[i] = fun_se_delta_meth_logit_glmmPQL(mod_glmm,0,nb_cov_clus,data[which(data$cluster==i),],adjust = cov_schem)**2
      }
      if(cov_schem == "both"){
        Var_per_clus[i] = fun_se_delta_meth_logit_glmmPQL(mod_glmm,nb_cov_indiv,nb_cov_clus,data[which(data$cluster==i),],adjust = cov_schem)**2
      }
      
    }
    
    Var_w = (1/K) * sum(Var_per_clus)
    
    Var_b = sum(c((Prev_clus[Prev_clus$Arm==1,]$RD_ij-mean(Prev_clus[Prev_clus$Arm==1,]$RD_ij))**2,
                  (Prev_clus[Prev_clus$Arm==0,]$RD_ij-mean(Prev_clus[Prev_clus$Arm==0,]$RD_ij))**2))/(K-2)
    
    Var_tot = Var_w + (K +1)/K * Var_b
    
    SE_CL_GC_Bin = sqrt(Var_tot)
    
    LL = RD - qt(0.975,K-1) * SE_CL_GC_Bin
    UL = RD + qt(0.975,K-1) * SE_CL_GC_Bin
    
   if(SE_CL_GC_Bin >1){
     Tab_res = data.frame(RD_CL_GC_Bin = NA,
                          LL_95_CL_GC_Bin = NA,
                          UL_95_CL_GC_Bin = NA,
                          SE_CL_GC_Bin = NA,
                          Bias = NA,
                          Relative_Bias = NA,
                          Itt = itt_para,
                          Error = "Standard error > 1")
   }else{
     Tab_res = data.frame(RD_CL_GC_Bin = RD,
                         LL_95_CL_GC_Bin = LL,
                         UL_95_CL_GC_Bin = UL,
                         SE_CL_GC_Bin = SE_CL_GC_Bin,
                         Bias = RD - True_RD,
                         Relative_Bias = (RD-True_RD)/True_RD*100,
                         Itt = itt_para,Error = "No error")
   }
  }else{
    Tab_res = data.frame(RD_CL_GC_Bin = NA,
                         LL_95_CL_GC_Bin = NA,
                         UL_95_CL_GC_Bin = NA,
                         SE_CL_GC_Bin = NA,
                         Bias = NA,
                         Relative_Bias = NA,
                         Itt = itt_para,
                         Error = 1)
  }
  
  return(Tab_res)
  
}

fun_se_delta_meth_logit_glmmPQL <- function(fit,nb_cov,nb_cov_clus,data,adjust){
  
  d = c()
  
  sum_fit = summary(fit)
  
  fixed = fit$coefficients$fixed
  random = fit$coefficients$random$cluster
  
  x1 <- x0 <- data
  x1$Arm <- 1
  x0$Arm <- 0
  
  x1_mod = cbind(Intercept = rep(1,nrow(x1)),x1[,-c(1,3)])
  x0_mod = cbind(Intercept = rep(1,nrow(x0)),x0[,-c(1,3)])
  
  if(adjust == "both"){
    Trt <- invlogit(fixed["(Intercept)"]+ fixed["Arm"] + fixed["IC_1"] * data$IC_1 + fixed["IC_2"] * data$IC_2 + fixed["IC_3"] * data$IC_3+ fixed["CC_1"] * data$CC_1+ fixed["CC_2"] * data$CC_2 + random[data$cluster])
    NoTrt <- invlogit(fixed["(Intercept)"]+ fixed["IC_1"] * data$IC_1 + fixed["IC_2"] * data$IC_2 + fixed["IC_3"] * data$IC_3+ fixed["CC_1"] * data$CC_1+ fixed["CC_2"] * data$CC_2 + random[data$cluster])
  }else if(adjust == "individual"){
    Trt <- invlogit(fixed["(Intercept)"]+ fixed["Arm"] + fixed["IC_1"] * data$IC_1 + fixed["IC_2"] * data$IC_2 + fixed["IC_3"] * data$IC_3 + random[data$cluster])
    NoTrt <- invlogit(fixed["(Intercept)"]+ fixed["IC_1"] * data$IC_1 + fixed["IC_2"] * data$IC_2 + fixed["IC_3"] * data$IC_3 + random[data$cluster])
  }
  
  
  fixed_beta = as.vector(sum_fit$coefficients$fixed)
  nb_var = length(fixed_beta)
  
  for (j in 1:nb_var) {
    d_j = mean(deriv_logit(Trt) * x1_mod[,j] -  deriv_logit(NoTrt) * x0_mod[,j])
    d=c(d,d_j)
  }
  
  se = sqrt(d %*% vcov(fit) %*% d)
  
  return(se)
}
# Fonction need from first paper

fun_var_corrected_FG <- function(fit,formula,id,family=gaussian,data,corstr="independence",b=0.75){
  
  #########################################################################
  # !!!! Function GEE.var.fg from the 'geesmv' package only changing 'formula' by 'fit' !!!!
  # !!!! to reduce calculation time by not doing a second gee                           !!!!
  
  # Arguments:
  # fit      specify the gee fit already performed
  # formula  specify the model of interest
  # family   "gaussian", "binomial" or "poisson"
  # data     data frame
  # corstr   Working correlation structure: "independence", "AR-M", "exchangeable", "unstructured".
  # value:   GEE returns the following elements
  #          cov.var      estimate of the variance-covariance matrix for robust variance.
  
  # !!!! Function GEE.var.fg from the 'geesmv' package only changing 'formula' by 'fit' !!!!
  # !!!! to reduce calculation time by not doing a second gee                           !!!!
  #########################################################################
  # Delete the records with missing data in predictors or outcomes;
  if (is.null(data$id)){
    index <- which(names(data)==id)
    data$id <- data[,index]}
  
  ### na.action: only na.omit is used for gee;
  init <- model.frame(formula, data)
  init$num <- 1:length(init[,1])
  if(any(is.na(init))){
    index <- na.omit(init)$num
    data <- data[index,]
    ### Get the design matrix;
    m <- model.frame(formula, data)
    mt <- attr(m, "terms") 
    data$response <- model.response(m, "numeric")
    mat <- as.data.frame(model.matrix(formula, m))
  }else{
    ### Get the design matrix;
    m <- model.frame(formula, data)
    mt <- attr(m, "terms") 
    data$response <- model.response(m, "numeric")
    mat <- as.data.frame(model.matrix(formula, m))
  }
  
  ### Fit the GEE model to get the estimate of parameters \hat{\beta};
  gee.fit <- fit
  beta_est <- gee.fit$coefficient
  alpha <- gee.fit$working.correlation[1,2]
  len <- length(beta_est)
  len_vec <- len^2
  
  ### Estimate the robust variance for \hat{\beta}
  data$id <- gee.fit$id
  cluster<-cluster.size(data$id)
  ncluster<-max(cluster$n)
  size<-cluster$m
  mat$subj <- rep(unique(data$id), cluster$n)
  if(is.character(corstr)){
    var <- switch(corstr,
                  "independence"=cormax.ind(ncluster),
                  "exchangeable"=cormax.exch(ncluster, alpha),
                  "AR-M"=cormax.ar1(ncluster, alpha),
                  "unstructured"=summary(gee.fit)$working.correlation,)
  }else{
    print(corstr)
    stop("'working correlation structure' not recognized")
  }   
  if(is.character(family)){
    family <- switch(family,
                     "gaussian"="gaussian",
                     "binomial"="binomial",
                     "poisson"="poisson")
  }else{ 
    if(is.function(family)){
      family <- family()[[1]]
    }else{
      print(family)
      stop("'family' not recognized")
    }    
  }
  
  cov.beta<-unstr<-matrix(0,nrow=len,ncol=len)
  step11<-matrix(0, nrow=len, ncol=len)
  for (i in 1:size){
    y<-as.matrix(data$response[data$id==unique(data$id)[i]])
    covariate<-as.matrix(subset(mat[,-length(mat[1,])], mat$subj==unique(data$id)[i]))
    ncluster=cluster$n[i]
    var1=var[1:ncluster,1:ncluster] 
    if (family=="gaussian"){ 
      Vi=gee.fit$scale*var1
      xx<-t(covariate)%*%solve(Vi)%*%covariate
      step11<-step11+xx  
    }else if (family=="poisson") {
      D<-mat.prod(covariate, exp(covariate%*%beta_est))
      Vi <- diag(sqrt(c(exp(covariate%*%beta_est))),ncluster)%*%var1%*%diag(sqrt(c(exp(covariate%*%beta_est))),ncluster)
      xx<-t(D)%*%solve(Vi)%*%D
      step11<-step11+xx
    }else if (family=="binomial"){
      D<-mat.prod(covariate, exp(covariate%*%beta_est)/((1+exp(covariate%*%beta_est))^2))
      Vi <- diag(sqrt(c(exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est))^2)),ncluster)%*%var1%*%diag(sqrt(c(exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est))^2)),ncluster)
      xx<-t(D)%*%solve(Vi)%*%D
      step11<-step11+xx 
    }
  }
  step12<-matrix(0,nrow=len,ncol=len)
  p<-matrix(0,nrow=len_vec,ncol=size)
  
  for (i in 1:size){
    y<-as.matrix(data$response[data$id==unique(data$id)[i]])
    covariate<-as.matrix(subset(mat[,-length(mat[1,])], mat$subj==unique(data$id)[i]))
    ncluster=cluster$n[i]
    var1=var[1:ncluster,1:ncluster]
    if (family=="gaussian"){ 
      ## set up the scale parameter;
      Vi=gee.fit$scale*var1
      xx<-t(covariate)%*%solve(Vi)%*%covariate
      Qi <- xx%*%solve(step11)
      Ai<-diag((1-pmin(b,diag(Qi)))^(-0.5))
      xy<-Ai%*%t(covariate)%*%solve(Vi)%*%(y-covariate%*%beta_est)
      step12<-step12+xy%*%t(xy)
      p[,i]<-vec(xy%*%t(xy))
    }else if (family=="poisson") {
      ## set up the scale parameter;
      D<-mat.prod(covariate, exp(covariate%*%beta_est))
      Vi <- diag(sqrt(c(exp(covariate%*%beta_est))),ncluster)%*%var1%*%diag(sqrt(c(exp(covariate%*%beta_est))),ncluster)
      xx<-t(D)%*%solve(Vi)%*%D
      Qi <- xx%*%solve(step11)
      Ai<-diag((1-pmin(b,diag(Qi)))^(-0.5))
      xy<-Ai%*%t(D)%*%solve(Vi)%*%(y-exp(covariate%*%beta_est))
      step12<-step12+xy%*%t(xy)
      p[,i]<-vec(xy%*%t(xy))
    }else if (family=="binomial"){
      ## set up the scale parameter;
      D<-mat.prod(covariate, exp(covariate%*%beta_est)/((1+exp(covariate%*%beta_est))^2))
      Vi <- diag(sqrt(c(exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est))^2)),ncluster)%*%var1%*%diag(sqrt(c(exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est))^2)),ncluster)
      xx<-t(D)%*%solve(Vi)%*%D
      Qi <- xx%*%solve(step11)
      Ai<-diag((1-pmin(b,diag(Qi)))^(-0.5))
      xy<-Ai%*%t(D)%*%solve(Vi)%*%(y-exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est)))
      step12<-step12+xy%*%t(xy)
      p[,i]<-vec(xy%*%t(xy)) 
    }    
  }
  cov.beta<-solve(step11)%*%(step12)%*%solve(step11)
  return(cov.beta)
}


Pred_Trt_fun <- function(fit,data,nb_cov_tot){
  coef = fit$coefficients
  Pred_Trt = rep(coef[1]+ coef[2],nrow(data)) 
  if(nb_cov_tot>0){
    for (i in 1:nb_cov_tot) {
      Pred_Trt = Pred_Trt + data[,3+i]*coef[2+i]
    }
  }
  
  Pred_Trt <- as.vector(Pred_Trt)
  return(Pred_Trt)
}

Pred_NoTrt_fun <- function(fit,data,nb_cov_tot){
  coef = fit$coefficients
  Pred_NoTrt = rep(coef[1],nrow(data))
  if(nb_cov_tot>0){
    for (i in 1:nb_cov_tot) {
      Pred_NoTrt = Pred_NoTrt + data[,3+i]*coef[2+i]
    }
  }
  
  Pred_NoTrt <- as.vector(Pred_NoTrt)
  return(Pred_NoTrt)
}

fun_formula_for_glm <- function(nb_cov,nb_cov_clus){
  
  Covariable = ''
  
  if(nb_cov > 0){
    for (i in 1:nb_cov) {
      a = paste('IC_',i,sep = "")
      Covariable = paste(Covariable,a,sep = '+')
    }
  }
  if(nb_cov_clus > 0){
    for (i in 1:nb_cov_clus) {
      a = paste('CC_',i,sep = "")
      Covariable = paste(Covariable,a,sep = '+')
    }
  }
  
  
  Covariable = as.factor(Covariable)
  Formula = as.formula(paste("Outcome ~ Arm ",Covariable))
  
  return(Formula)
}

fun_formula_for_gee <- function(nb_cov,nb_cov_clus){
  
  Covariable = ''
  
  if(nb_cov > 0){
    for (i in 1:nb_cov) {
      a = paste('IC_',i,sep = "")
      Covariable = paste(Covariable,a,sep = '+')
    }
  }
  if(nb_cov_clus > 0){
    for (i in 1:nb_cov_clus) {
      a = paste('CC_',i,sep = "")
      Covariable = paste(Covariable,a,sep = '+')
    }
  }
  
  
  Covariable = as.factor(Covariable)
  Formula = as.formula(paste("Outcome ~ Arm ",Covariable))
  
  return(Formula)
}

fun_formula_for_glmer <- function(nb_cov,nb_cov_clus){
  
  Covariable = ''
  
  if(nb_cov > 0){
    for (i in 1:nb_cov) {
      a = paste('IC_',i,sep = "")
      Covariable = paste(Covariable,a,sep = '+')
    }
  }
  if(nb_cov_clus > 0){
    for (i in 1:nb_cov_clus) {
      a = paste('CC_',i,sep = "")
      Covariable = paste(Covariable,a,sep = '+')
    }
  }
  
  
  Covariable = as.factor(Covariable)
  Formula = as.formula(paste("Outcome ~ Arm ",Covariable," + (1|cluster)"))
  
  return(Formula)
}


# Input

# fit : model of the 'geese' function 
# nb_cov : number of individual level covariates of the data set
# nb_cov_tot : number of cluster level covariates of the data set
# data = data set generated by the function "Crea_data_V12" or "fun_Data_confusion_V7()" (see below section "7-Confusion function")

# Output 

# d: The resulting matrix calculating for the delta method corresponding with the link function logit or log.


fun_deriv_logit_tot <- function(fit,nb_cov,nb_cov_clus,data){
  d = c()
  nb_cov_tot = nb_cov + nb_cov_clus
  x = get_data(data,nb_cov,nb_cov_clus)
  x1 <- x0 <- x
  x1[,2]<-1
  x0[,2]<-0
  
  betas <- as.vector(fit$coefficients)
  nb_var = length(betas)
  Trt = betas[1] + betas[2]
  NoTrt = betas[1]
  if(nb_cov_tot>0){
    for (i in 1:nb_cov_tot) {
      Trt = Trt + betas[2+i]*x[,2+i]            # it is the estimation of each individual considering in the intervention group
      NoTrt = NoTrt + betas[2+i]*x[,2+i]        # it is the estimation of each individual considering in the control group
    }
    for (j in 1:nb_var) {
      d_j = mean(deriv_logit(Trt) * x1[,j] -  deriv_logit(NoTrt) * x0[,j]) # Difference between intervention and control estimation for each individual as the function used in the delta method
      d[j] = d_j
    }
  }
  else{
    for (j in 1:nb_var) {
      d_j = mean(deriv_logit(Trt) * x1[,j] -  deriv_logit(NoTrt) * x0[,j])
      d[j] = d_j
    }
  }
  return(d)
}

deriv_logit <- function(x){
  exp(-x)/(1+exp(-x))^2
}

get_data <- function(data,nb_cov,nb_cov_clus){
  nb_cov_tot = nb_cov + nb_cov_clus
  intercept = rep(1,nrow(data))
  x = data.frame(intercept,data$Arm)
  if(nb_cov_tot>0){
    for (i in 1:nb_cov_tot) {
      x = data.frame(x,data[,3+i])
    }
  }
  return(as.matrix(x))
}

### 11) TMLE function ----

#### TMLE Indiv cov ----
TMLE.JAPM.Indvcov <- function(data,Scenario,itt_para){
  
  # Step 0: Need scenario for the number of covariate for formula and confidence interval ----
  
  True_RD =  as.numeric(Scenario[1]) - as.numeric(Scenario[2])
  
  nb_cov_indiv = as.numeric(Scenario[6])
  nb_cov_clus = as.numeric(Scenario[7])
  nb_cov = nb_cov_indiv + nb_cov_clus
  K <- length(unique(data$cluster))
  N <- nrow(data)
  N_j <- as.vector(unlist((data %>% group_by(cluster) %>% summarise(n_j = length(Outcome)))[,2]))
  
  
  # Step 1: Initial Outcome, fit the model, predict with status set to "treatment" and "control" ----
  
  # Formula used in GEE considering only individual-level covariates
  formula_glm <- fun_formula_for_glm(nb_cov = nb_cov_indiv,nb_cov_clus = 0)
  
  
  # Fit the model with GEE approach using a exchangeable correlation matrix
  # Trying to catch a error, (convergence problem, suspect results etc)
  
  weight_data = rep(1/N_j,N_j)
  data$weight_data = rep(1/N_j,N_j)
  # cap_op2  <- capture.output(suppressMessages(mod_glm <- try(glm(formula = formula_glm,
  #                                                                data = data,
  #                                                                family = binomial,weights = data$weight_data))))
  
  mod_glm <- glm(formula = formula_glm,data = data,family = binomial,weights = weight_data)
  
  ## If "Error_gee" is TRUE, we have automatically a convergence problem
  ## If "Error_gee" is FALSE but its value is different of 0, GEE estimates the parameters but they are "suspect" so we admit a convergence problem.
  ## In this last case, it can be a non positive estimation of the correlation matrix
  
  # Prediction with status set to "treatment" -> data_1 and "control" -> data_0
  
  data_1 <- data_0 <- data
  data_1$Arm <- 1
  data_0$Arm <- 0
  
  # Prediction of the outcome if all individuals where in treatment group
  data$Y_1 <- predict(mod_glm,newdata = data_1,type = "response")
  
  # Prediction of the outcome if all individuals where in control group
  data$Y_0 <- predict(mod_glm,newdata = data_0,type = "response")
  
  # Prediction of the outcome if all individuals in there respective group
  data$Y_Arm <- invlogit(mod_glm$linear.predictors)
  
  
  # Step 2: Probability of treatment, estimation of "inverse probability of receiving treatment" and "Negative inverse probability of not receiving treatment"  ----
  
  # Fit the model with a glm using binomial family
  g_fit <- glm(Arm ~ IC_1 + IC_2 + IC_3 ,family = binomial,data = data,weights = weight_data )
  
  
  ## Individual-level
  
  pscore = as.vector(predict(g_fit,newdata = data,type = "response"))
  
  pscore = ifelse(pscore>0.975,0.975,pscore)
  pscore = ifelse(pscore<0.025,0.025,pscore)
  
  data$pscore = pscore
  
  # We maximize and minimize the propensity score.
  
  # Step : Aggregating data to the cluster level ----
  
  Mod_data_clus <- aggregate(data,by = list(data$cluster),mean)[,-1]
  
  # Inverse probability of revceiving treatment
  
  Mod_data_clus$H_1 <- Mod_data_clus$Arm/Mod_data_clus$pscore
  
  # Negative inverse probability of not revceiving treatment
  
  Mod_data_clus$H_0 <- (1-Mod_data_clus$Arm)/(1-Mod_data_clus$pscore)
  
  # Clever covariate
  
  Mod_data_clus$H_Arm <- Mod_data_clus$Arm/Mod_data_clus$pscore - (1-Mod_data_clus$Arm)/(1-Mod_data_clus$pscore)
  
  # Step 3: Fluctuation parameter ----
  
  eps_fit <- glm(Outcome ~ -1 + offset(qlogis(Y_Arm)) + H_1 + H_0,family = binomial,data = Mod_data_clus)
  eps <- coef(eps_fit)
  
  
  # Step 4: Update Initial Outcome ----
  
  Mod_data_clus$Y_star_1 = plogis(qlogis(Mod_data_clus$Y_1) + eps[1]/(Mod_data_clus$pscore))
  Mod_data_clus$Y_star_0 = plogis(qlogis(Mod_data_clus$Y_0) + eps[2]/(1-Mod_data_clus$pscore))
  Mod_data_clus$Y_star_A = Mod_data_clus$Arm*Mod_data_clus$Y_star_1 + (1-Mod_data_clus$Arm)*Mod_data_clus$Y_star_0
  
  # Step 5: Compute ATE & Inference JAPM ----
  
  
  Psi_1 = mean(Mod_data_clus$Y_star_1)
  
  Psi_0 = mean(Mod_data_clus$Y_star_0)
  
  # RD_JAPM = mean(agg_Y1$x) - mean(agg_Y0$x)
  
  RD_JAPM = Psi_1 - Psi_0
  
  
  D1_ij = (Mod_data_clus$H_1*(Mod_data_clus$Outcome-Mod_data_clus$Y_star_1) + Mod_data_clus$Y_star_1 - Psi_1)
  D0_ij = (Mod_data_clus$H_0*(Mod_data_clus$Outcome-Mod_data_clus$Y_star_0) + Mod_data_clus$Y_star_0 - Psi_0)
  
  # sqrt(var(D1_ij)/K)
  # sqrt(var(D0_ij)/K)
  # sqrt(var(D1_ij - D0_ij)/K)
  
  # Method dans le visual guide de Hoffman
  
  SE_CL_TMLE_JAPM = as.numeric(sqrt(var(D1_ij - D0_ij)/K))
  
  LL_JAPM = RD_JAPM - qt(0.975,K-2) * SE_CL_TMLE_JAPM
  UL_JAPM = RD_JAPM + qt(0.975,K-2) * SE_CL_TMLE_JAPM
  
  # Table of result
  
  Tab_res = data.frame(RD_CL_TMLE_JAPM = RD_JAPM,
                       SE_CL_TMLE_JAPM = SE_CL_TMLE_JAPM,
                       LL_95_CL_TMLE_JAPM = LL_JAPM,
                       UL_95_CL_TMLE_JAPM = UL_JAPM,
                       Bias_JAPM = RD_JAPM - True_RD,
                       Relative_Bias_JAPM = (RD_JAPM-True_RD)/True_RD*100,
                       Itt = itt_para)
  
  
  return(Tab_res)
}


#### TMLE both cov ----
TMLE.JAPM.bothcov <- function(data,Scenario,itt_para){
  
  # Step 0: Need scenario for the number of covariate for formula and confidence interval ----
  
  True_RD =  as.numeric(Scenario[1]) - as.numeric(Scenario[2])
  
  nb_cov_indiv = as.numeric(Scenario[6])
  nb_cov_clus = as.numeric(Scenario[7])
  nb_cov = nb_cov_indiv + nb_cov_clus
  K <- length(unique(data$cluster))
  N <- nrow(data)
  N_j <- as.vector(unlist((data %>% group_by(cluster) %>% summarise(n_j = length(Outcome)))[,2]))
  
  
  # Step 1: Initial Outcome, fit the model, predict with status set to "treatment" and "control" ----
  
  # Formula used in GEE considering only individual-level covariates
  formula_glm <- fun_formula_for_glm(nb_cov_indiv,nb_cov_clus)
  
  
  # Fit the model with GEE approach using a exchangeable correlation matrix
  # Trying to catch a error, (convergence problem, suspect results etc)
  
  weight_data = rep(1/N_j,N_j)
  data$weight_data = rep(1/N_j,N_j)
  
  cap_op2  <- capture.output(suppressMessages(mod_glm <- try(glm(formula = formula_glm,
                                                                 data = data,
                                                                 family = binomial,
                                                                 weights = weight_data))))
  
  
  ## If "Error_gee" is TRUE, we have automatically a convergence problem
  ## If "Error_gee" is FALSE but its value is different of 0, GEE estimates the parameters but they are "suspect" so we admit a convergence problem.
  ## In this last case, it can be a non positive estimation of the correlation matrix
  
  # Prediction with status set to "treatment" -> data_1 and "control" -> data_0
  
  data_1 <- data_0 <- data
  data_1$Arm <- 1
  data_0$Arm <- 0
  
  # Prediction of the outcome if all individuals where in treatment group
  data$Y_1 <- predict(mod_glm,newdata = data_1,type = "response")
  
  # Prediction of the outcome if all individuals where in control group
  data$Y_0 <- predict(mod_glm,newdata = data_0,type = "response")
  
  # Prediction of the outcome if all individuals in there respective group
  data$Y_Arm <- invlogit(mod_glm$linear.predictors)
  
  
  # Step 2: Probability of treatment, estimation of "inverse probability of receiving treatment" and "Negative inverse probability of not receiving treatment"  ----
  
  # Fit the model with a glm using binomial family
  
  g_fit <- glm(Arm ~ IC_1 + IC_2 + IC_3 + CC_1 + CC_2 ,family = binomial,data = data,weights = weight_data)
  
  
  ## Individual-level
  
  pscore = as.vector(predict(g_fit,newdata = data,type = "response"))
  
  pscore = ifelse(pscore>0.975,0.975,pscore)
  pscore = ifelse(pscore<0.025,0.025,pscore)
  
  data$pscore = pscore
  
  # We maximize and minimize the propensity score.
  
  # Step : Aggregating data to the cluster level ----
  
  Mod_data_clus <- aggregate(data,by = list(data$cluster),mean)[,-1]
  
  # Inverse probability of revceiving treatment
  
  Mod_data_clus$H_1 <- Mod_data_clus$Arm/Mod_data_clus$pscore
  
  # Negative inverse probability of not revceiving treatment
  
  Mod_data_clus$H_0 <- (1-Mod_data_clus$Arm)/(1-Mod_data_clus$pscore)
  
  # Clever covariate
  
  Mod_data_clus$H_Arm <- Mod_data_clus$Arm/Mod_data_clus$pscore - (1-Mod_data_clus$Arm)/(1-Mod_data_clus$pscore)
  
  # Step 3: Fluctuation parameter ----
  
  eps_fit <- glm(Outcome ~ -1 + offset(qlogis(Y_Arm)) + H_1 + H_0,family = binomial,data = Mod_data_clus)
  eps <- coef(eps_fit)
  
  
  # Step 4: Update Initial Outcome ----
  
  Mod_data_clus$Y_star_1 = plogis(qlogis(Mod_data_clus$Y_1) + eps[1]/(Mod_data_clus$pscore))
  Mod_data_clus$Y_star_0 = plogis(qlogis(Mod_data_clus$Y_0) + eps[2]/(1-Mod_data_clus$pscore))
  Mod_data_clus$Y_star_A = Mod_data_clus$Arm*Mod_data_clus$Y_star_1 + (1-Mod_data_clus$Arm)*Mod_data_clus$Y_star_0
  
  # Step 5: Compute ATE & Inference JAPM ----
  
  
  Psi_1 = mean(Mod_data_clus$Y_star_1)
  
  Psi_0 = mean(Mod_data_clus$Y_star_0)
  
  # RD_JAPM = mean(agg_Y1$x) - mean(agg_Y0$x)
  
  RD_JAPM = Psi_1 - Psi_0
  
  
  D1_ij = (Mod_data_clus$H_1*(Mod_data_clus$Outcome-Mod_data_clus$Y_star_1) + Mod_data_clus$Y_star_1 - Psi_1)
  D0_ij = (Mod_data_clus$H_0*(Mod_data_clus$Outcome-Mod_data_clus$Y_star_0) + Mod_data_clus$Y_star_0 - Psi_0)
  
  sqrt(var(D1_ij)/K)
  sqrt(var(D0_ij)/K)
  sqrt(var(D1_ij - D0_ij)/K)
  
  # Method dans le visual guide de Hoffman
  
  SE_CL_TMLE_JAPM = as.numeric(sqrt(var(D1_ij - D0_ij)/K))
  
  LL_JAPM = RD_JAPM - qt(0.975,K-2) * SE_CL_TMLE_JAPM
  UL_JAPM = RD_JAPM + qt(0.975,K-2) * SE_CL_TMLE_JAPM
  
  # Table of result
  
  Tab_res = data.frame(RD_CL_TMLE_JAPM = RD_JAPM,
                       SE_CL_TMLE_JAPM = SE_CL_TMLE_JAPM,
                       LL_95_CL_TMLE_JAPM = LL_JAPM,
                       UL_95_CL_TMLE_JAPM = UL_JAPM,
                       Bias_JAPM = RD_JAPM - True_RD,
                       Relative_Bias_JAPM = (RD_JAPM-True_RD)/True_RD*100,
                       Itt = itt_para)
  
  
  return(Tab_res)
}


TMLE.JAPM.ML <- function(data,Scenario,itt_para,ML.method = "step"){
  
  # Step 0: Need scenario for the number of covariate for formula and confidence interval ----
  
  True_RD =  as.numeric(Scenario[1]) - as.numeric(Scenario[2])
  
  nb_cov_indiv = as.numeric(Scenario[6])
  nb_cov_clus = as.numeric(Scenario[7])
  nb_cov = nb_cov_indiv + nb_cov_clus
  K <- length(unique(data$cluster))
  N <- nrow(data)
  N_j <- as.vector(unlist((data %>% group_by(cluster) %>% summarise(n_j = length(Outcome)))[,2]))
  
  
  # Step 1: Initial Outcome, fit the model, predict with status set to "treatment" and "control" ----
  
  # Formula used in GEE considering only individual-level covariates
  # formula_glm <- fun_formula_for_glm(nb_cov_indiv,nb_cov_clus)
  # 
  
  weight_data = rep(1/N_j,N_j)
  data$weight_data = rep(1/N_j,N_j)
  
  
  data_Step1 = data[,c("Outcome","Arm","IC_1","IC_2","IC_3","CC_1","CC_2")]
  
  data_1 <- data_0 <- data_Step1
  data_1$Arm <- 1
  data_0$Arm <- 0
  
  cap_op2  <- capture.output(suppressMessages(mod_glm <- try(glm(formula = Outcome~.,
                                                                 data = data_Step1,
                                                                 family = binomial,
                                                                 weights = data$weight_data))))
  if(ML.method == 'step'){
    mod_glm <- step(mod_glm,scope = list(lower = ~Arm), direction = "both", trace = 0, k = 2)
  }else if(ML.method=='lasso'){
    
    Adj_Step1 <- model.matrix(~-1 +. , subset(data_Step1, select=-Outcome) )
    
    data_1 <- model.matrix(~-1 +. , subset(data_1, select=-Outcome) )
    data_0 <- model.matrix(~-1 +. , subset(data_0, select=-Outcome) )
    
    mod_glm <- glmnet(x=Adj_Step1,  y=data_Step1$Outcome, weights=data$weight_data,
                      family=binomial(), alpha=1, nlambda = 100)
  }else if(ML.method %in% c('mars', 'mars.corP')){
    # using default settings of SL.earth in SuperLearner 
    X <- subset(data_Step1, select=-Outcome)
    if(ML.method=='mars.corP'){
      Vec_True_var = which(screen.corP(Y=data_Step1$Outcome, X=X, family='binomial')==T)
      if(1 %in% Vec_True_var){X <- X[Vec_True_var]}else{X <- X[c(1,Vec_True_var)]}
    } 
    
    mod_glm <- earth::earth(x = X,  y=data_Step1$Outcome, weights=data$weight_data,
                            degree = 2, penalty = 3, 
                            nk = max(21, 2 * ncol(X) + 1), pmethod = "backward", nfold = 0, 
                            ncross = 1, minspan = 0, endspan = 0,
                            glm = list(family = binomial))
  }
  
  # Prediction with status set to "treatment" -> data_1 and "control" -> data_0
  
  if(ML.method=='lasso'){
    # from LASSO 
    data$Y_Arm <- as.vector(predict(mod_glm, newx=Adj_Step1, type='response', s = min(mod_glm$lambda)))
    data$Y_1 <- as.vector(predict(mod_glm, newx=data_1, type='response', s = min(mod_glm$lambda)))
    data$Y_0 <- as.vector(predict(mod_glm, newx=data_0, type='response', s = min(mod_glm$lambda)))
  } else {
    # for glms 
    data$Y_Arm <- as.vector(predict(mod_glm,newdata = data,type = "response"))
    data$Y_1 <- as.vector(predict(mod_glm,newdata = data_1,type = "response"))
    data$Y_0 <- as.vector(predict(mod_glm,newdata = data_0,type = "response"))
  }
  
  # Step 2: Probability of treatment, estimation of "inverse probability of receiving treatment" and "Negative inverse probability of not receiving treatment"  ----
  
  # Fit the model with a glm using binomial family
  
  # g_fit <- glm(Arm ~ IC_1 + IC_2 + IC_3 + CC_1 + CC_2 ,family = binomial,data = data,weights = weight_data)
  
  data_Step2 = data[,c("Arm","IC_1","IC_2","IC_3","CC_1","CC_2")]
  
  
  cap_op2  <- capture.output(suppressMessages(g_fit <- try(glm(formula = Arm~.,
                                                               data = data_Step2,
                                                               family = binomial,
                                                               weights = data$weight_data))))
  if(ML.method == 'step'){
    g_fit <- step(g_fit, direction = "both", trace = 0, k = 2)
  }else if(ML.method=='lasso'){
    
    Adj_Step2 <- model.matrix(~-1 +. , subset(data_Step2, select=-Arm) )
    
    g_fit <- glmnet(x=Adj_Step2,  y=data_Step2$Arm, weights=data$weight_data,
                    family=binomial(), alpha=1, nlambda = 100)
  }else if(ML.method %in% c('mars', 'mars.corP')){
    # using default settings of SL.earth in SuperLearner 
    X <- subset(data_Step2, select=-Arm)
    if(ML.method=='mars.corP'){
      X <- X[screen.corP(Y=data_Step2$Arm, X=X, family='binomial')==T]
    } 
    
    g_fit <- earth::earth(x = X,  y=data_Step2$Arm, weights=data$weight_data,
                          degree = 2, penalty = 3, 
                          nk = max(21, 2 * ncol(X) + 1), pmethod = "backward", nfold = 0, 
                          ncross = 1, minspan = 0, endspan = 0,
                          glm = list(family = binomial))
  }
  
  
  
  ## Individual-level
  
  if(ML.method=='lasso'){
    # from LASSO 
    pscore <- as.vector(predict(g_fit, newx=Adj_Step2, type='response', s = min(g_fit$lambda)))
  } else {
    # for glms 
    pscore = as.vector(predict(g_fit,newdata = data,type = "response"))
    # pscore = as.vector(g_fit$fitted.values)
  }
  
  
  # pscore = as.vector(predict(g_fit,newdata = data,type = "response"))
  
  pscore = ifelse(pscore>0.975,0.975,pscore)
  pscore = ifelse(pscore<0.025,0.025,pscore)
  
  data$pscore = pscore
  
  # We maximize and minimize the propensity score.
  
  # Step : Aggregating data to the cluster level ----
  
  Mod_data_clus <- aggregate(data,by = list(data$cluster),mean)[,-1]
  
  # Inverse probability of revceiving treatment
  
  Mod_data_clus$H_1 <- Mod_data_clus$Arm/Mod_data_clus$pscore
  
  # Negative inverse probability of not revceiving treatment
  
  Mod_data_clus$H_0 <- (1-Mod_data_clus$Arm)/(1-Mod_data_clus$pscore)
  
  # Clever covariate
  
  Mod_data_clus$H_Arm <- Mod_data_clus$Arm/Mod_data_clus$pscore - (1-Mod_data_clus$Arm)/(1-Mod_data_clus$pscore)
  
  # Step 3: Fluctuation parameter ----
  
  if(Inf %in% qlogis(Mod_data_clus$Y_Arm)){
    Tab_res = data.frame(RD_CL_TMLE_JAPM = NA,
                         SE_CL_TMLE_JAPM = NA,
                         LL_95_CL_TMLE_JAPM = NA,
                         UL_95_CL_TMLE_JAPM = NA,
                         Bias_JAPM = NA,
                         Relative_Bias_JAPM = NA,
                         Itt = itt_para,
                         ML.method = ML.method)
  }else{
    eps_fit <- glm(Outcome ~ -1 + offset(qlogis(Y_Arm)) + H_1 + H_0,family = binomial,data = Mod_data_clus)
    eps <- coef(eps_fit)
    
    
    # Step 4: Update Initial Outcome ----
    
    Mod_data_clus$Y_star_1 = plogis(qlogis(Mod_data_clus$Y_1) + eps[1]/(Mod_data_clus$pscore))
    Mod_data_clus$Y_star_0 = plogis(qlogis(Mod_data_clus$Y_0) + eps[2]/(1-Mod_data_clus$pscore))
    Mod_data_clus$Y_star_A = Mod_data_clus$Arm*Mod_data_clus$Y_star_1 + (1-Mod_data_clus$Arm)*Mod_data_clus$Y_star_0
    
    # Step 5: Compute ATE & Inference JAPM ----
    
    
    Psi_1 = mean(Mod_data_clus$Y_star_1)
    
    Psi_0 = mean(Mod_data_clus$Y_star_0)
    
    # RD_JAPM = mean(agg_Y1$x) - mean(agg_Y0$x)
    
    RD_JAPM = Psi_1 - Psi_0
    
    
    D1_ij = (Mod_data_clus$H_1*(Mod_data_clus$Outcome-Mod_data_clus$Y_star_1) + Mod_data_clus$Y_star_1 - Psi_1)
    D0_ij = (Mod_data_clus$H_0*(Mod_data_clus$Outcome-Mod_data_clus$Y_star_0) + Mod_data_clus$Y_star_0 - Psi_0)
    
    sqrt(var(D1_ij)/K)
    sqrt(var(D0_ij)/K)
    sqrt(var(D1_ij - D0_ij)/K)
    
    # Method dans le visual guide de Hoffman
    
    SE_CL_TMLE_JAPM = as.numeric(sqrt(var(D1_ij - D0_ij)/K))
    
    LL_JAPM = RD_JAPM - qt(0.975,K-2) * SE_CL_TMLE_JAPM
    UL_JAPM = RD_JAPM + qt(0.975,K-2) * SE_CL_TMLE_JAPM
    
    # Table of result
    
    Tab_res = data.frame(RD_CL_TMLE_JAPM = RD_JAPM,
                         SE_CL_TMLE_JAPM = SE_CL_TMLE_JAPM,
                         LL_95_CL_TMLE_JAPM = LL_JAPM,
                         UL_95_CL_TMLE_JAPM = UL_JAPM,
                         Bias_JAPM = RD_JAPM - True_RD,
                         Relative_Bias_JAPM = (RD_JAPM-True_RD)/True_RD*100,
                         Itt = itt_para,
                         ML.method = ML.method)
  }
  
  return(Tab_res)
}

R_model <- function(data,model_reg){
  
  if(model_reg=="glm_both"){
    
    res <- capture.output(suppressMessages(model <- glm(formula = Outcome~Arm + IC_1 +IC_2+ IC_3+ CC_1 + CC_2,
                                                        family = binomial(),
                                                        data = data,
                                                        weights = data$weight_data)))
  }else if(model_reg=="glm_indv"){
    
    res <- capture.output(suppressMessages(model <- glm(formula = Outcome~Arm + IC_1 +IC_2+ IC_3,
                                                        family = binomial(),
                                                        data = data,
                                                        weights = data$weight_data)))
  }else if(model_reg=="step"){
    
    res <- capture.output(suppressMessages(model <- glm(formula = Outcome ~ Arm + IC_1 +IC_2+ IC_3+ CC_1 + CC_2,
                                                        family = binomial(),
                                                        data = data,
                                                        weights = data$weight_data)))
    capture.output(suppressMessages(model <- step(model,scope = list(lower=as.formula(Outcome ~ Arm )),direction = "both")))
  }else if(model_reg=="lasso"){
    
    Adj_Step1_data <- model.matrix(~-1 +. , subset(data, select=-c(Outcome,cluster,weight_data)))
    res <- capture.output(suppressMessages(model <- glmnet(x=Adj_Step1_data,
                                                           y=data$Outcome,
                                                           weights=data$weight_data,
                                                           family=binomial(),
                                                           alpha=1, nlambda = 100)))
  }else if(model_reg=="mars"){
    
    X <- subset(data, select=-c(Outcome,cluster,weight_data))
    
    model <- earth::earth(x = X,  y=data$Outcome, weights=data$weight_data,
                          degree = 2, penalty = 3, 
                          nk = max(21, 2 * ncol(X) + 1), pmethod = "backward", nfold = 0, 
                          ncross = 1, minspan = 0, endspan = 0,
                          glm = list(family = binomial))
  }
  
  # model prédit
  return(model)
}

Pred_step1 <- function(data,model_reg,fit_reg){
  
  Adj_Step1 <- model.matrix(~-1 +. , subset(data, select=-c(Outcome,cluster,weight_data)))
  
  if(model_reg == "lasso"){
    pred_reg_s1 <- as.vector(predict(fit_reg, newx=Adj_Step1, type='response', s = min(fit_reg$lambda)))
  }
  if(model_reg %in% c("glm_indv","glm_both","step")){
    pred_reg_s1 <- predict.glm(fit_reg,data,type = "response")
  }
  if(model_reg == "mars"){
    pred_reg_s1 <- predict(fit_reg,data,type = "response")
  }
  return(pred_reg_s1)
}

P_model <- function(data,model_pscore){
  
  if(model_pscore=="glm_both"){
    
    res <- capture.output(suppressMessages(model <- glm(formula = Arm ~ IC_1 +IC_2+ IC_3+ CC_1 + CC_2,
                                                        family = binomial(),
                                                        data = data,
                                                        weights = data$weight_data)))
  }else if(model_pscore=="glm_indv"){
    
    res <- capture.output(suppressMessages(model <- glm(formula = Arm ~ IC_1 +IC_2+ IC_3,
                                                        family = binomial(),
                                                        data = data,
                                                        weights = data$weight_data)))
  }else if(model_pscore=="step"){
    
    res <- capture.output(suppressMessages(model <- glm(formula = Arm ~ IC_1 +IC_2+ IC_3+ CC_1 + CC_2,
                                                        family = binomial(),
                                                        data = data,
                                                        weights = data$weight_data)))
    
    capture.output(suppressMessages(model <-step(model,direction = "both")))
    
  }else if(model_pscore=="lasso"){
    
    Adj_Step2 <- model.matrix(~-1 +. , subset(data, select=-c(Outcome,Arm,cluster,weight_data)))
    res <- capture.output(suppressMessages(model <- glmnet(x=Adj_Step2,
                                                           y=data$Arm,
                                                           weights=data$weight_data,
                                                           family=binomial(),
                                                           alpha=1, nlambda = 100)))
  }else if(model_pscore=="mars"){
    
    X <- subset(data, select=-c(Outcome,Arm,cluster,weight_data))
    
    model <- earth::earth(x = X,  y=data$Arm, weights=data$weight_data,
                          degree = 2, penalty = 3, 
                          nk = max(21, 2 * ncol(X) + 1), pmethod = "backward", nfold = 0, 
                          ncross = 1, minspan = 0, endspan = 0,
                          glm = list(family = binomial))
  }
  
  # model prédit
  return(model)
}

Pred_step2 <- function(data,model_pscore,fit_pscore){
  if(model_pscore=='lasso'){
    Adj_Step2 <- model.matrix(~-1 +. , subset(data, select=-c(Outcome,Arm,cluster,weight_data)) )
    # from LASSO 
    pscore <- as.vector(predict(fit_pscore, newx=Adj_Step2, type='response', s = min(fit_pscore$lambda)))
  } else {
    # for glms 
    pscore = as.vector(predict(fit_pscore,newdata = data,type = "response"))
  }
  return(pscore)
}

TMLE.APS.Step1.CV <- function(data,model_reg){
  
  K = length(unique(data$cluster))
  IC_CV_S1 = c()
  
  for (j in 1:K) {
    
    train.data <- data[-which(data$cluster==j),]
    train.data.0 <- train.data.1 <- data[-which(data$cluster==j),]
    train.data.0$Arm <- 0
    train.data.1$Arm <- 1
    
    valid.data <- data[which(data$cluster==j),]
    valid.data.0 <- valid.data.1 <- data[which(data$cluster==j),]
    valid.data.0$Arm <- 0
    valid.data.1$Arm <- 1
    
    # Train data ----
    fit_reg <- R_model(data = train.data,model_reg = model_reg)
    
    # model_step <- step(model,scope = list(lower=as.formula(Outcome ~ Arm)),direction = "both")
    # predict.glm(model_step,train.data,type = "response")
    
    
    ## Step 1 TMLE ----
    
    train.data$Y_A <- Pred_step1(train.data,model_reg,fit_reg)
    train.data$Y_0 <- Pred_step1(train.data.0,model_reg,fit_reg)
    train.data$Y_1 <- Pred_step1(train.data.1,model_reg,fit_reg)
    
    ## Step 2 on fait l'hypothese d'un score de propension de 0.5 ----
    train.data$pscore <- 0.5
    
    ## Step : Aggregating data to the cluster level ----
    
    train.data.clus <- aggregate(train.data,by = list(train.data$cluster),mean)[,-1]
    
    # Inverse probability of revceiving treatment
    
    train.data.clus$H_1 <- train.data.clus$Arm/train.data.clus$pscore
    train.data.clus$H_0 <- (1-train.data.clus$Arm)/(1-train.data.clus$pscore)
    train.data.clus$H_Arm <- train.data.clus$Arm/train.data.clus$pscore - (1-train.data.clus$Arm)/(1-train.data.clus$pscore)
    
    ## Step 3: Fluctuation parameter ----
    
    res <- capture.output(suppressMessages(eps_fit <- glm(Outcome ~ -1 + offset(qlogis(Y_A)) + H_1 + H_0,
                                                          family = binomial,
                                                          data = train.data.clus)))
    eps <- coef(eps_fit)
    
    
    # Valid data estimation ----
    
    ## Step 1 ----
    valid.data$Y_A <- Pred_step1(valid.data,model_reg,fit_reg)
    valid.data$Y_0 <- Pred_step1(valid.data.0,model_reg,fit_reg)
    valid.data$Y_1 <- Pred_step1(valid.data.1,model_reg,fit_reg)
    
    
    ## Step 2 on fait l'hypothese d'un score de propension de 0.5 ----
    valid.data$pscore <- 0.5
    
    ## Step : Aggregating data to the cluster level ----
    
    valid.data.clus <- aggregate(valid.data,by = list(valid.data$cluster),mean)[,-1]
    
    
    valid.data.clus$H_1 <- valid.data.clus$Arm/valid.data.clus$pscore
    valid.data.clus$H_0 <- (1-valid.data.clus$Arm)/(1-valid.data.clus$pscore)
    valid.data.clus$H_Arm <- valid.data.clus$Arm/valid.data.clus$pscore - (1-valid.data.clus$Arm)/(1-valid.data.clus$pscore)
    
    
    ## Step 4: Update Initial Outcome ----
    
    valid.data.clus$Y_star_1 = plogis(qlogis(valid.data.clus$Y_1) + eps[1]/(valid.data.clus$pscore))
    valid.data.clus$Y_star_0 = plogis(qlogis(valid.data.clus$Y_0) + eps[2]/(1-valid.data.clus$pscore))
    valid.data.clus$Y_star_A = valid.data.clus$Arm*valid.data.clus$Y_star_1 + (1-valid.data.clus$Arm)*valid.data.clus$Y_star_0
    
    ## Step 5: Compute ATE & Inference JAPM ----
    
    
    Psi_1 = mean(valid.data.clus$Y_star_1)
    Psi_0 = mean(valid.data.clus$Y_star_0)
    Psi_A = mean(valid.data.clus$Y_star_A)
    
    # RD_JAPM = mean(agg_Y1$x) - mean(agg_Y0$x)
    
    RD_JAPM = Psi_1 - Psi_0
    
    DA_ij = (valid.data.clus$H_Arm*(valid.data.clus$Outcome-valid.data.clus$Y_star_A) + valid.data.clus$Y_star_A - Psi_A)
    # D1_ij = (valid.data.clus$H_1*(valid.data.clus$Outcome-valid.data.clus$Y_star_1) + valid.data.clus$Y_star_1 - Psi_1)
    # D0_ij = (valid.data.clus$H_0*(valid.data.clus$Outcome-valid.data.clus$Y_star_0) + valid.data.clus$Y_star_A - Psi_0)
    
    
    # Method dans le visual guide de Hoffman
    
    IC_CV_S1[j] = as.numeric(DA_ij)
  }
  return(IC_CV_S1)
}

TMLE.APS.Step2.CV <- function(data,model_reg,model_pscore){
  K = length(unique(data$cluster))
  IC_CV_S2 = c()
  
  for (j in 1:K) {
    
    train.data <- data[-which(data$cluster==j),]
    train.data.0 <- train.data.1 <- data[-which(data$cluster==j),]
    train.data.0$Arm <- 0
    train.data.1$Arm <- 1
    
    valid.data <- data[which(data$cluster==j),]
    valid.data.0 <- valid.data.1 <- data[which(data$cluster==j),]
    valid.data.0$Arm <- 0
    valid.data.1$Arm <- 1
    
    # Train data ----
    fit_reg <- R_model(data = train.data,model_reg = model_reg)
    
    
    ## Step 1 TMLE ----
    
    train.data$Y_A <- Pred_step1(train.data,model_reg,fit_reg)
    train.data$Y_0 <- Pred_step1(train.data.0,model_reg,fit_reg)
    train.data$Y_1 <- Pred_step1(train.data.1,model_reg,fit_reg)
    
    ## Step 2 on fait l'hypothese d'un score de propension de 0.5 ----
    
    fit_pscore <- P_model(train.data,model_pscore)
    
    train.data$pscore <- Pred_step2(train.data,model_pscore,fit_pscore)
    
    ## Step : Aggregating data to the cluster level ----
    
    train.data.clus <- aggregate(train.data,by = list(train.data$cluster),mean)[,-1]
    
    # Inverse probability of revceiving treatment
    
    train.data.clus$H_1 <- train.data.clus$Arm/train.data.clus$pscore
    train.data.clus$H_0 <- (1-train.data.clus$Arm)/(1-train.data.clus$pscore)
    train.data.clus$H_Arm <- train.data.clus$Arm/train.data.clus$pscore - (1-train.data.clus$Arm)/(1-train.data.clus$pscore)
    
    ## Step 3: Fluctuation parameter ----
    
    res <- capture.output(suppressMessages(eps_fit <- glm(Outcome ~ -1 + offset(qlogis(Y_A)) + H_1 + H_0,
                                                          family = binomial,
                                                          data = train.data.clus)))
    eps <- coef(eps_fit)
    
    
    # Valid data estimation ----
    
    ## Step 1 ----
    valid.data$Y_A <- Pred_step1(valid.data,model_reg,fit_reg)
    valid.data$Y_0 <- Pred_step1(valid.data.0,model_reg,fit_reg)
    valid.data$Y_1 <- Pred_step1(valid.data.1,model_reg,fit_reg)
    
    
    ## Step 2 on fait l'hypothese d'un score de propension de 0.5 ----
    valid.data$pscore <- Pred_step2(valid.data,model_pscore,fit_pscore)
    
    valid.data$pscore = ifelse(valid.data$pscore>0.975,0.975,valid.data$pscore)
    valid.data$pscore = ifelse(valid.data$pscore<0.025,0.025,valid.data$pscore)
    
    ## Step : Aggregating data to the cluster level ----
    
    valid.data.clus <- aggregate(valid.data,by = list(valid.data$cluster),mean)[,-1]
    
    
    valid.data.clus$H_1 <- valid.data.clus$Arm/valid.data.clus$pscore
    valid.data.clus$H_0 <- (1-valid.data.clus$Arm)/(1-valid.data.clus$pscore)
    valid.data.clus$H_Arm <- valid.data.clus$Arm/valid.data.clus$pscore - (1-valid.data.clus$Arm)/(1-valid.data.clus$pscore)
    
    
    ## Step 4: Update Initial Outcome ----
    
    valid.data.clus$Y_star_1 = plogis(qlogis(valid.data.clus$Y_1) + eps[1]/(valid.data.clus$pscore))
    valid.data.clus$Y_star_0 = plogis(qlogis(valid.data.clus$Y_0) + eps[2]/(1-valid.data.clus$pscore))
    valid.data.clus$Y_star_A = valid.data.clus$Arm*valid.data.clus$Y_star_1 + (1-valid.data.clus$Arm)*valid.data.clus$Y_star_0
    
    ## Step 5: Compute ATE & Inference JAPM ----
    
    
    Psi_1 = mean(valid.data.clus$Y_star_1)
    Psi_0 = mean(valid.data.clus$Y_star_0)
    Psi_A = mean(valid.data.clus$Y_star_A)
    
    # RD_JAPM = mean(agg_Y1$x) - mean(agg_Y0$x)
    
    RD_JAPM = Psi_1 - Psi_0
    
    DA_ij = (valid.data.clus$H_Arm*(valid.data.clus$Outcome-valid.data.clus$Y_star_A) + valid.data.clus$Y_star_A - Psi_A)
    # D1_ij = (valid.data.clus$H_1*(valid.data.clus$Outcome-valid.data.clus$Y_star_1) + valid.data.clus$Y_star_1 - Psi_1)
    # D0_ij = (valid.data.clus$H_0*(valid.data.clus$Outcome-valid.data.clus$Y_star_0) + valid.data.clus$Y_star_A - Psi_0)
    
    
    # Method dans le visual guide de Hoffman
    
    IC_CV_S2[j] = as.numeric(DA_ij)
  }
  return(IC_CV_S2)
}

fun_selec_MLmeth <- function(data){
  var_s1 <- c(
    lasso_S1    <- mean(TMLE.APS.Step1.CV(data,"lasso")^2),
    # mars_S1     <- mean(TMLE.APS.Step1.CV(data,"mars")^2),
    glm_indv_S1 <- mean(TMLE.APS.Step1.CV(data,"glm_indv")^2),
    glm_both_S1 <- mean(TMLE.APS.Step1.CV(data,"glm_both")^2),
    step_S1     <- mean(TMLE.APS.Step1.CV(data,"step")^2)
  )
  
  method_S1_index <- which(var_s1==min(var_s1))
  if(method_S1_index==1){model_reg="lasso"}
  # if(method_S1_index==2){model_reg="mars"}
  if(method_S1_index==2){model_reg="glm_indv"}
  if(method_S1_index==3){model_reg="glm_both"}
  if(method_S1_index==4){model_reg="step"}
  
  
  var_s2 <- c(
    lasso_S2    <- mean(TMLE.APS.Step2.CV(data,model_reg,"lasso")^2),
    # mars_S2     <- mean(TMLE.APS.Step2.CV(data,model_reg,"mars")^2),
    glm_indv_S2 <- mean(TMLE.APS.Step2.CV(data,model_reg,"glm_indv")^2),
    glm_both_S2 <- mean(TMLE.APS.Step2.CV(data,model_reg,"glm_both")^2),
    step_S2     <- mean(TMLE.APS.Step2.CV(data,model_reg,"step")^2)
  )
  
  method_S2_index <- which(var_s2==min(var_s2))
  if(method_S2_index==1){model_pscore="lasso"}
  # if(method_S2_index==2){model_pscore="mars"}
  if(method_S2_index==2){model_pscore="glm_indv"}
  if(method_S2_index==3){model_pscore="glm_both"}
  if(method_S2_index==4){model_pscore="step"}
  
  return(c(model_reg,model_pscore))
}

# fun_selec_MLmeth(data)

TMLE.APS <- function(data,Scenario,itt_para){
  
  # Step 0: Need scenario for the true risk difference ----
  
  True_RD =  as.numeric(Scenario[1]) - as.numeric(Scenario[2])
  
  K <- length(unique(data$cluster))
  N <- nrow(data)
  N_j <- as.vector(unlist((data %>% group_by(cluster) %>% summarise(n_j = length(Outcome)))[,2]))
  weight_data = rep(1/N_j,N_j)
  data$weight_data = rep(1/N_j,N_j)
  
  # Models selection with the cross validation method
  
  method_selec  = fun_selec_MLmeth(data)
  method_reg    = method_selec[1]
  method_pscore = method_selec[2]
  
  
  # Step 1: Initial Outcome, fit the model, predict with status set to "treatment" and "control" using ML method selected for regression model ----
  
  
  fit_reg <- R_model(data,method_reg)
  
  # Prediction with status set to "treatment" -> data_1 and "control" -> data_0
  
  data_1 <- data_0 <- data
  data_1$Arm <- 1
  data_0$Arm <- 0
  
  # Prediction of the outcome if all individuals in there respective group
  data$Y_Arm <- Pred_step1(data,method_reg,fit_reg)
  
  # Prediction of the outcome if all individuals where in treatment group
  data$Y_1 <- Pred_step1(data_1,method_reg,fit_reg)
  
  # Prediction of the outcome if all individuals where in control group
  data$Y_0 <- Pred_step1(data_0,method_reg,fit_reg)
  
  
  # Step 2: Probability of treatment, estimation of "inverse probability of receiving treatment" and "Negative inverse probability of not receiving treatment"  ----
  
  # Fit the model with the ML method selected using cross validation
  g_fit <- P_model(data,method_pscore)
  
  
  ## Individual-level
  
  pscore = Pred_step2(data,method_pscore,g_fit)
  
  pscore = ifelse(pscore>0.975,0.975,pscore)
  pscore = ifelse(pscore<0.025,0.025,pscore)
  
  data$pscore = pscore
  
  # We maximize and minimize the propensity score.
  
  # Step : Aggregating data to the cluster level ----
  
  Mod_data_clus <- aggregate(data,by = list(data$cluster),mean)[,-1]
  
  # Inverse probability of revceiving treatment
  
  Mod_data_clus$H_1 <- Mod_data_clus$Arm/Mod_data_clus$pscore
  
  # Negative inverse probability of not revceiving treatment
  
  Mod_data_clus$H_0 <- (1-Mod_data_clus$Arm)/(1-Mod_data_clus$pscore)
  
  # Clever covariate
  
  Mod_data_clus$H_Arm <- Mod_data_clus$Arm/Mod_data_clus$pscore - (1-Mod_data_clus$Arm)/(1-Mod_data_clus$pscore)
  
  # Step 3: Fluctuation parameter ----
  
  eps_fit <- glm(Outcome ~ -1 + offset(qlogis(Y_Arm)) + H_1 + H_0,family = binomial,data = Mod_data_clus)
  eps <- coef(eps_fit)
  
  
  # Step 4: Update Initial Outcome ----
  
  Mod_data_clus$Y_star_1 = plogis(qlogis(Mod_data_clus$Y_1) + eps[1]/(Mod_data_clus$pscore))
  Mod_data_clus$Y_star_0 = plogis(qlogis(Mod_data_clus$Y_0) + eps[2]/(1-Mod_data_clus$pscore))
  Mod_data_clus$Y_star_A = Mod_data_clus$Arm*Mod_data_clus$Y_star_1 + (1-Mod_data_clus$Arm)*Mod_data_clus$Y_star_0
  
  # Step 5: Compute ATE & Inference JAPM ----
  
  
  Psi_1 = mean(Mod_data_clus$Y_star_1)
  
  Psi_0 = mean(Mod_data_clus$Y_star_0)
  
  # RD_JAPM = mean(agg_Y1$x) - mean(agg_Y0$x)
  
  RD_JAPM = Psi_1 - Psi_0
  
  
  D1_ij = (Mod_data_clus$H_1*(Mod_data_clus$Outcome-Mod_data_clus$Y_star_1) + Mod_data_clus$Y_star_1 - Psi_1)
  D0_ij = (Mod_data_clus$H_0*(Mod_data_clus$Outcome-Mod_data_clus$Y_star_0) + Mod_data_clus$Y_star_0 - Psi_0)
  
  # sqrt(var(D1_ij)/K)
  # sqrt(var(D0_ij)/K)
  # sqrt(var(D1_ij - D0_ij)/K)
  
  # Method dans le visual guide de Hoffman
  
  SE_CL_TMLE_JAPM = as.numeric(sqrt(var(D1_ij - D0_ij)/K))
  
  LL_JAPM = RD_JAPM - qt(0.975,K-2) * SE_CL_TMLE_JAPM
  UL_JAPM = RD_JAPM + qt(0.975,K-2) * SE_CL_TMLE_JAPM
  
  # Table of result
  
  Tab_res = data.frame(RD_CL_TMLE_JAPM = RD_JAPM,
                       SE_CL_TMLE_JAPM = SE_CL_TMLE_JAPM,
                       LL_95_CL_TMLE_JAPM = LL_JAPM,
                       UL_95_CL_TMLE_JAPM = UL_JAPM,
                       Bias_JAPM = RD_JAPM - True_RD,
                       Relative_Bias_JAPM = (RD_JAPM-True_RD)/True_RD*100,
                       Itt = itt_para)
  
  
  return(Tab_res)
  
  
  
}

fun_para_analyse_cluster_TMLE_APS <- function(itt_para,Scenario,n,Data_file,Resu_file){
  
  # Recuperation Data
  data = read.csv(paste(paste(Data_file,"/Scenario_",n,sep=""),sep ="",paste("/Data_itt_",sep = "",itt_para,"_scen_",n,".csv")), sep=";")
  
  ## Analyse all method
  
  TMLE_res_APS = TMLE.APS(data,Scenario,itt_para)
  
  # Remplissage CSV resultat pour chacune des méthodes
  
  write.table(TMLE_res_APS,file = paste(Resu_file,"/TMLE_analyse_APS",paste("/Data_output_CL_TMLE_APS_Scenario_",sep = "",n,".csv"),sep = ""),
              append = TRUE,row.names = FALSE,col.names = FALSE,sep = ";",dec = ",")
  
  
}




fun_para_analyse_cluster_TMLE <- function(itt_para,Scenario,n,Data_file,Resu_file){
  
  # Recuperation Data
  data = read.csv(paste(paste(Data_file,"/Scenario_",n,sep=""),sep ="",paste("/Data_itt_",sep = "",itt_para,"_scen_",n,".csv")), sep=";")
  
  ## Analyse all method
  
  TMLE_res_glm                      = TMLE.JAPM.ML(data,Scenario,itt_para,ML.method = "glm")
  
  TMLE_res_lasso                    = TMLE.JAPM.ML(data,Scenario,itt_para,ML.method = "lasso")
  
  TMLE_res_step                     = TMLE.JAPM.ML(data,Scenario,itt_para,ML.method = "step")
  
  TMLE_res_mars                     = TMLE.JAPM.ML(data,Scenario,itt_para,ML.method = "mars")
  
  TMLE_res_mars.corP                = TMLE.JAPM.ML(data,Scenario,itt_para,ML.method = "mars.corP")
  
  # res_CL_TMLE_indvcov       = TMLE.JAPM.Indvcov(data,Scenario = Scenario,itt_para = itt_para)
  # 
  # res_CL_TMLE_bothcov       = TMLE.JAPM.bothcov(data,Scenario = Scenario,itt_para = itt_para)
  
  # Remplissage CSV resultat pour chacune des méthodes
  
  write.table(TMLE_res_glm,file = paste(Resu_file,"/TMLE_analyse_glm",paste("/Data_output_CL_TMLE_glm_Scenario_",sep = "",n,".csv"),sep = ""),
              append = TRUE,row.names = FALSE,col.names = FALSE,sep = ";",dec = ",")
  
  write.table(TMLE_res_lasso,file = paste(Resu_file,"/TMLE_analyse_lasso",paste("/Data_output_CL_TMLE_lasso_Scenario_",sep = "",n,".csv"),sep = ""),
              append = TRUE,row.names = FALSE,col.names = FALSE,sep = ";",dec = ",")
  
  write.table(TMLE_res_step,file = paste(Resu_file,"/TMLE_analyse_step",paste("/Data_output_CL_TMLE_step_Scenario_",sep = "",n,".csv"),sep = ""),
              append = TRUE,row.names = FALSE,col.names = FALSE,sep = ";",dec = ",")
  
  write.table(TMLE_res_mars,file = paste(Resu_file,"/TMLE_analyse_mars",paste("/Data_output_CL_TMLE_mars_Scenario_",sep = "",n,".csv"),sep = ""),
              append = TRUE,row.names = FALSE,col.names = FALSE,sep = ";",dec = ",")
  
  write.table(TMLE_res_mars.corP,file = paste(Resu_file,"/TMLE_analyse_marscorP",paste("/Data_output_CL_TMLE_marscorP_Scenario_",sep = "",n,".csv"),sep = ""),
              append = TRUE,row.names = FALSE,col.names = FALSE,sep = ";",dec = ",")
  
}


fun_para_analyse_cluster_TMLE_marsML <- function(itt_para,Scenario,n,Data_file,Resu_file){
  
  # Recuperation Data
  data = read.csv(paste(paste(Data_file,"/Scenario_",n,sep=""),sep ="",paste("/Data_itt_",sep = "",itt_para,"_scen_",n,".csv")), sep=";")
  
  ## Analyse all method
  
  #TMLE_res_glm                      = TMLE.JAPM.ML(data,Scenario,itt_para,ML.method = "glm")
  
  #TMLE_res_lasso                    = TMLE.JAPM.ML(data,Scenario,itt_para,ML.method = "lasso")
  
  #TMLE_res_step                     = TMLE.JAPM.ML(data,Scenario,itt_para,ML.method = "step")
  
  TMLE_res_mars                     = TMLE.JAPM.ML(data,Scenario,itt_para,ML.method = "mars")
  
  TMLE_res_mars.corP                = TMLE.JAPM.ML(data,Scenario,itt_para,ML.method = "mars.corP")
  
  # res_CL_TMLE_indvcov       = TMLE.JAPM.Indvcov(data,Scenario = Scenario,itt_para = itt_para)
  # 
  # res_CL_TMLE_bothcov       = TMLE.JAPM.bothcov(data,Scenario = Scenario,itt_para = itt_para)
  
  # Remplissage CSV resultat pour chacune des méthodes
  
  #write.table(TMLE_res_glm,file = paste(Resu_file,"/TMLE_analyse_glm",paste("/Data_output_CL_TMLE_glm_Scenario_",sep = "",n,".csv"),sep = ""),
  #            append = TRUE,row.names = FALSE,col.names = FALSE,sep = ";",dec = ",")
  
  #write.table(TMLE_res_lasso,file = paste(Resu_file,"/TMLE_analyse_lasso",paste("/Data_output_CL_TMLE_lasso_Scenario_",sep = "",n,".csv"),sep = ""),
  #            append = TRUE,row.names = FALSE,col.names = FALSE,sep = ";",dec = ",")
  
  #write.table(TMLE_res_step,file = paste(Resu_file,"/TMLE_analyse_step",paste("/Data_output_CL_TMLE_step_Scenario_",sep = "",n,".csv"),sep = ""),
  #            append = TRUE,row.names = FALSE,col.names = FALSE,sep = ";",dec = ",")
  
  write.table(TMLE_res_mars,file = paste(Resu_file,"/TMLE_analyse_mars",paste("/Data_output_CL_TMLE_mars_Scenario_",sep = "",n,".csv"),sep = ""),
              append = TRUE,row.names = FALSE,col.names = FALSE,sep = ";",dec = ",")
  
  write.table(TMLE_res_mars.corP,file = paste(Resu_file,"/TMLE_analyse_marscorP",paste("/Data_output_CL_TMLE_marscorP_Scenario_",sep = "",n,".csv"),sep = ""),
              append = TRUE,row.names = FALSE,col.names = FALSE,sep = ";",dec = ",")
  
}






fun_g_comp_bin_independent <- function(data,Scenario,itt_para,cov_schem = "individual"){
  
  ## Need scenario for the number of covariate for formula and confidence interval
  True_RD =  as.numeric(Scenario[1]) - as.numeric(Scenario[2])
  
  nb_cov_indiv = as.numeric(Scenario[6])
  nb_cov_clus = as.numeric(Scenario[7])
  nb_cov = nb_cov_indiv + nb_cov_clus
  
  if(cov_schem == "individual"){
    formula_gee <- fun_formula_for_gee(nb_cov = nb_cov_indiv,nb_cov_clus = 0)
  }
  if(cov_schem == "cluster"){
    formula_gee <- fun_formula_for_gee(nb_cov = 0,nb_cov_clus = nb_cov_clus)
  }
  if(cov_schem == "both"){
    formula_gee <- fun_formula_for_gee(nb_cov_indiv,nb_cov_clus)
  }
  
  
  cap_op2  <- capture.output(suppressMessages(mod_gee <- try(gee(formula = formula_gee,id = cluster,data = data,family = binomial,corstr = "independence"))))
  
  if(inherits(mod_gee,"try-error")==FALSE){
    data_1 <- data_0 <- data
    
    data_1$Arm <- 1
    data_0$Arm <- 0
    
    # Y_1 <- invlogit(Pred_Trt_fun(mod_gee,data = data_1,nb_cov_tot = nb_cov))
    # 
    # Y_0 <- invlogit(Pred_NoTrt_fun(mod_gee,data = data_0,nb_cov_tot = nb_cov))
    if(cov_schem == "individual"){
      Y_1 <- invlogit(Pred_Trt_fun(mod_gee,data = data_1,nb_cov_tot = nb_cov_indiv))
      
      Y_0 <- invlogit(Pred_NoTrt_fun(mod_gee,data = data_0,nb_cov_tot = nb_cov_indiv))
    }
    if(cov_schem == "cluster"){
      Y_1 <- invlogit(Pred_Trt_fun(mod_gee,data = data_1,nb_cov_tot = nb_cov_clus))
      
      Y_0 <- invlogit(Pred_NoTrt_fun(mod_gee,data = data_0,nb_cov_tot = nb_cov_clus))
    }
    if(cov_schem == "both"){
      Y_1 <- invlogit(Pred_Trt_fun(mod_gee,data = data_1,nb_cov_tot = nb_cov))
      
      Y_0 <- invlogit(Pred_NoTrt_fun(mod_gee,data = data_0,nb_cov_tot = nb_cov))
    }
    
    
    # Pred_indiv <- Pred_Trt - Pred_NoTrt
    
    data_gcomp <- data.frame(data,Y_1,Y_0)
    
    Prev_clus <- data_gcomp %>% group_by(cluster)%>%
      summarise(Arm=mean(Arm),Y_1j = mean(Y_1),Y_0j = mean(Y_0))
    
    # Number of cluster
    K = length(Prev_clus$cluster)
    
    # RD by cluster
    Prev_clus$RD_ij = Prev_clus$Y_1j - Prev_clus$Y_0j 
    
    # Risk difference is the mean of all RD estimated per cluster
    RD = mean(Prev_clus$RD_ij)
    
    # Confidence Interval
    # Matrice_RobustVariance = mod_gee$robust.variance
    
    # Matrice variance covariance corrected with Fay and Graubard Method
    Mat_corr_FG_varcov = fun_var_corrected_FG(mod_gee,formula_gee,id = "cluster",family=binomial,data,corstr="exchangeable",b=0.75)
    
    Var_per_clus = c()
    
    for (i in 1:K) {
      
      if(cov_schem == "individual"){
        Var_per_clus[i] = fun_deriv_logit_tot(mod_gee,nb_cov_indiv,0,data[which(data$cluster==i),]) %*% Mat_corr_FG_varcov %*% fun_deriv_logit_tot(mod_gee,nb_cov_indiv,0,data[which(data$cluster==i),])
      }
      if(cov_schem == "cluster"){
        Var_per_clus[i] = fun_deriv_logit_tot(mod_gee,0,nb_cov_clus,data[which(data$cluster==i),]) %*% Mat_corr_FG_varcov %*% fun_deriv_logit_tot(mod_gee,0,nb_cov_clus,data[which(data$cluster==i),])
      }
      if(cov_schem == "both"){
        Var_per_clus[i] = fun_deriv_logit_tot(mod_gee,nb_cov_indiv,nb_cov_clus,data[which(data$cluster==i),]) %*% Mat_corr_FG_varcov %*% fun_deriv_logit_tot(mod_gee,nb_cov_indiv,nb_cov_clus,data[which(data$cluster==i),])
      }
      
    }
    
    Var_w = (1/K) * sum(Var_per_clus)
    
    # Var_b1 = var(Prev_clus$Y_1j-Prev_clus$Y_0j)
    # Var_b2 = sum(c((Prev_clus[1:5,]$RD_ij-mean(Prev_clus[1:5,]$RD_ij))**2,
    #                (Prev_clus[6:10,]$RD_ij-mean(Prev_clus[6:10,]$RD_ij))**2))/(k-1-2)
    # 
    # Var_b3 = sum((Prev_clus$RD_ij-mean(Prev_clus$RD_ij))**2)/(k-1)
    # var(Prev_clus$RD_ij)
    
    Var_b = sum(c((Prev_clus[Prev_clus$Arm==1,]$RD_ij-mean(Prev_clus[Prev_clus$Arm==1,]$RD_ij))**2,
                  (Prev_clus[Prev_clus$Arm==0,]$RD_ij-mean(Prev_clus[Prev_clus$Arm==0,]$RD_ij))**2))/(K-2)
    
    Var_tot = Var_w + (K +1)/K * Var_b
    
    SE_CL_GC_Bin = sqrt(Var_tot)
    
    LL = RD - qt(0.975,K-1) * SE_CL_GC_Bin
    UL = RD + qt(0.975,K-1) * SE_CL_GC_Bin
    
    # if(cov_schem == "individual"){
    #   LL = RD - qt(0.975,k-nb_cov_indiv) * SE_CL_GC_Bin
    #   UL = RD + qt(0.975,k-nb_cov_indiv) * SE_CL_GC_Bin
    # }
    # if(cov_schem == "cluster"){
    #   LL = RD - qt(0.975,k-nb_cov_clus) * SE_CL_GC_Bin
    #   UL = RD + qt(0.975,k-nb_cov_clus) * SE_CL_GC_Bin
    # }
    # if(cov_schem == "both"){
    #   LL = RD - qt(0.975,k-nb_cov) * SE_CL_GC_Bin
    #   UL = RD + qt(0.975,k-nb_cov) * SE_CL_GC_Bin
    # }
    
    
    if(mod_gee$error != 0){
      Tab_res = data.frame(RD_CL_GC_Bin = NA,
                           LL_95_CL_GC_Bin = NA,
                           UL_95_CL_GC_Bin = NA,
                           SE_CL_GC_Bin = NA,
                           Bias = NA,
                           Relative_Bias = NA,
                           Itt = itt_para,
                           Error = mod_gee$error)
    }else{
      Tab_res = data.frame(RD_CL_GC_Bin = RD,
                           LL_95_CL_GC_Bin = LL,
                           UL_95_CL_GC_Bin = UL,
                           SE_CL_GC_Bin = SE_CL_GC_Bin,
                           Bias = RD - True_RD,
                           Relative_Bias = (RD-True_RD)/True_RD*100,
                           Itt = itt_para,Error = mod_gee$error)
    }
    
  }else{
    Tab_res = data.frame(RD_CL_GC_Bin = NA,
                         LL_95_CL_GC_Bin = NA,
                         UL_95_CL_GC_Bin = NA,
                         SE_CL_GC_Bin = NA,
                         Bias = NA,
                         Relative_Bias = NA,
                         Itt = itt_para,
                         Error = 1)
  }
  
  return(Tab_res)
  
}

fun_para_analyse_cluster_gcomp_indv_independence <- function(itt_para,Scenario,n,Data_file,Resu_file){
  
  # Recuperation Data
  data = read.csv(paste(paste(Data_file,"/Scenario_",n,sep=""),sep ="",paste("/Data_itt_",sep = "",itt_para,"_scen_",n,".csv")), sep=";")
  
  ## Analyse all method
  
  # res_CL_GC_Bin_indep_both = fun_g_comp_bin_independent(data,itt_para = itt_para,Scenario = Scenario,cov_schem = "both")
  
  res_CL_GC_Bin_indep_indv = fun_g_comp_bin_independent(data,itt_para = itt_para,Scenario = Scenario,cov_schem = "individual")
  
  
  # Remplissage CSV resultat pour chacune des méthodes
  
  # write.table(res_CL_GC_Bin_indep_both,file = paste(Resu_file,"/GC_Bin_analyse_both_indep",paste("/Data_output_CL_GC_both_indep_Bin_Scenario_",sep = "",n,".csv"),sep = ""),
  #             append = TRUE,row.names = FALSE,col.names = FALSE,sep = ";",dec = ",")
  
  write.table(res_CL_GC_Bin_indep_indv,file = paste(Resu_file,"/GC_Bin_analyse_indv_indep",paste("/Data_output_CL_GC_indv_indep_Bin_Scenario_",sep = "",n,".csv"),sep = ""),
              append = TRUE,row.names = FALSE,col.names = FALSE,sep = ";",dec = ",")
  
}

fun_para_analyse_cluster_gcomp_both_independence <- function(itt_para,Scenario,n,Data_file,Resu_file){
  
  # Recuperation Data
  data = read.csv(paste(paste(Data_file,"/Scenario_",n,sep=""),sep ="",paste("/Data_itt_",sep = "",itt_para,"_scen_",n,".csv")), sep=";")
  
  ## Analyse all method
  
  res_CL_GC_Bin_indep_both = fun_g_comp_bin_independent(data,itt_para = itt_para,Scenario = Scenario,cov_schem = "both")
  
  # res_CL_GC_Bin_indep_indv = fun_g_comp_bin_independent(data,itt_para = itt_para,Scenario = Scenario,cov_schem = "individual")
  
  
  # Remplissage CSV resultat pour chacune des méthodes
  
  write.table(res_CL_GC_Bin_indep_both,file = paste(Resu_file,"/GC_Bin_analyse_both_indep",paste("/Data_output_CL_GC_both_indep_Bin_Scenario_",sep = "",n,".csv"),sep = ""),
              append = TRUE,row.names = FALSE,col.names = FALSE,sep = ";",dec = ",")
  
  # write.table(res_CL_GC_Bin_indep_indv,file = paste(Resu_file,"/GC_Bin_analyse_indv_indep",paste("/Data_output_CL_GC_indv_indep_Bin_Scenario_",sep = "",n,".csv"),sep = ""),
  #             append = TRUE,row.names = FALSE,col.names = FALSE,sep = ";",dec = ",")
  
}

fun_cor_gcomp_indep <- function(i,Time,n,Base_file,Workspace_name,Data_file,Resu_file,cov_schem = "individual"){
  
  # 0. Espace de travail
  Chemin = sprintf(paste(Base_file,"/correction_gcomp/correction_gcomp_S%d",sep = ""),n)
  dir.create(Chemin)
  setwd(Chemin)
  
  # 0.1 Nom du fichier principal
  FileName = paste("correction_gcomp_S",n,"_itt_%d.R",sep = "")
  
  
  # 1. tu crees le fichier .R qui va contenir les instructions à lancer dans ta 2eme session
  
  
  Stock_file = sprintf(paste(Base_file,"/correction_g_comp/correction_gcomp_S%d",sep = ""),n)
  
  writeLines('rm(list = ls())',
             con = sprintf(FileName,i))
  
  
  write(sprintf("Chemin  = '%s'",Stock_file),
        file  = sprintf(FileName, i),
        append = TRUE)
  
  # 1.1 tu exporteras l'id de la 2eme session pour le recuperer. Ca servira pour le killer par la suite si besoin est.
  Exp_id = paste("writeLines(as.character(Sys.getpid()), con = sprintf('%s/pid_%d.txt', getwd(),",i,"))",sep = "")
  write(Exp_id,
        file = sprintf(FileName, i),
        append = TRUE)
  
  # 1.2 les lignes de commande, pour lancer le programme qui prend du temps
  write(paste("load('",Base_file,"/",Workspace_name,"')",sep = ""), file = sprintf(FileName, i), append = TRUE)
  
  
  write("library(cli, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(Matrix, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(MASS, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(lme4, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(backports, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(dplyr, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(parallel, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(iterators, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(rngtools, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(foreach, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(doRNG, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(doParallel, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(dplyr, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(gee, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(geepack, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(spind, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(doBy, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(arm, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(here, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(geesmv, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  write("library(matrixcalc, lib.loc='/home/jpereira/z.libs_Pierre/common')",file = sprintf(FileName, i), append = TRUE)
  
  L1 = paste("n =",n)
  write(L1,
        file = sprintf(FileName,i),
        append = TRUE)
  
  write(paste("Data_file = '",Data_file,"'",sep = ""), file = sprintf(FileName, i), append = TRUE)
  write(paste("Resu_file = '",Resu_file,"'",sep = ""), file = sprintf(FileName, i), append = TRUE)
  
  write(sprintf("Scenario_use = Scenario_%d",n),
        file = sprintf(FileName, i),
        append = TRUE)
  
  if(cov_schem == "individual"){
    write(sprintf("fun_para_analyse_cluster_gcomp_indv_independence(itt_para=%d,Scenario=Scenario_use,n=n,Data_file=Data_file,Resu_file = Resu_file)",i),
          file = sprintf(FileName, i),
          append = TRUE)
  }
  if(cov_schem == "both"){
    write(sprintf("fun_para_analyse_cluster_gcomp_both_independence(itt_para=%d,Scenario=Scenario_use,n=n,Data_file=Data_file,Resu_file = Resu_file)",i),
          file = sprintf(FileName, i),
          append = TRUE)
  }
  
  
  write("DONE = TRUE",file = sprintf(FileName,i),append = TRUE)
  
  # 1.3 l'export du resultat
  tace_2 = paste("save(DONE, file = sprintf('res_%d.Rdata',",i,"))",sep = "")
  write(tace_2, file = sprintf(FileName, i), append = TRUE)
  
  # 2. tu fais executer les commandes des fichiers correction_S33_itt_x.R en batch en arriere plan (wait = F)
  cmd_batch = paste("R CMD BATCH --no-restore correction_gcomp_S",n,"_itt_%d.R",sep = "")
  
  Start.T = Sys.time()
  system(sprintf(cmd_batch, i), wait = TRUE,timeout = Time)
  End.T = Sys.time()
  
  Diff.time = difftime(End.T,Start.T,units = "sec")
  if(Diff.time >= Time ){
    system(sprintf("tskill %d", scan(sprintf("pid_%d.txt", i), integer(),quiet = T)))
  }
  
}
###################################################################################################################################################a
###################################################################################################################################################a
###################################################################################################################################################a
###################################################################################################################################################a
###################################################################################################################################################a

# Generation scenario ----

########### Warning ###########
# Base_file corresponds to the file where Workspace/Data/Results will be saved
# Choose one and keep the same for all functions .R files, all simulation .R files

Base_file = ""
# Base_file = "/home/jpereira/Simulation_Article" # Exemple

# Create a file where the Data will be stock for all scenarii
# We recommend to create this file in the Base_file

Data_file = paste(Base_file,"/Data_itt",sep = "")


# We set a seed to generate 405 seed for each scenario
# Because of the convergence problem of the gee function we had to simulated each scenario independently
# And problem could happen at every time for a scenario so, we chose to have a seed for each scenario
# so if the program stop a one point we could restart the program starting from the scenario whom crashed
# It does not crashed a second time if we start from this scenario, we have no explanation to this phenomenon

Pi_int = 0.7                # Prevalence of inclusion in the intervention arm when we add confounding
Pi_con = 0.7                # Prevalence of inclusion in the control arm when we add confounding
rho_z = 0.2                 # Inclusion ICC
OR_int = c(2.5,2.5,2.5)     # Beta_m_1 parameter in our article, corresponding the parameter quantifying the associaation between the inclusion and the individual level covariate in the intervention arm
OR_con = c(1,1,1)           # Beta_m_0 parameter in our article, corresponding the parameter quantifying the associaation between the inclusion and the individual level covariate in the control arm


# Parameters set in our simulation study, Pair of prevalence, ICC, k number of clusters per arm,
# cov_para number of individual level covariate,
# cov_clus_para number of cluster level covariate,
# list_OR = Association measured between the individual level covariate and the outcome
# OR_cov_clus = Association measured between the cluster level covariate and the outcome

list_prev = list(c(0.5,0.4),c(0.2,0.15),c(0.1,0.05),c(0.5,0.5),c(0.2,0.2),c(0.1,0.1))

Nb_prev_paire = length(list_prev)
Mat_prev = matrix(0,nrow = 2,ncol = Nb_prev_paire)
for (i in 1:Nb_prev_paire) {
  Mat_prev[,i] = as.double(unlist(list_prev[i]))
}
Paire_prevalance = Mat_prev

icc_para = c(0.2,0.05,0.01,0.001)

Nb_icc = length(icc_para)

k_para = c(5,10,20,30,40)

Nb_k = length(k_para)

cov_para = c(3)

Nb_cov = length(cov_para)

cov_clus_para = c(2)

Nb_cov_clus = length(cov_clus_para)

list_OR = list(c(2.5,1.5,1))


num = 1
for (a in 1:Nb_prev_paire) {
  for (b in 1:Nb_icc) {
    for (c in 1:Nb_k) {
      if(k_para[c]==5){
        m_para = c(25,280)
        Nb_m = length(m_para)
      }
      if(k_para[c]==10){
        m_para = c(15,200)
        Nb_m = length(m_para)
      }
      if(k_para[c]==20){
        m_para = c(10,150)
        Nb_m = length(m_para)
      }
      if(k_para[c]==30){
        m_para = c(10,150)
        Nb_m = length(m_para)
      }
      if(k_para[c]==40){
        m_para = c(50)
        Nb_m = length(m_para)
      }
      for (d in 1:Nb_m) {
        for (e in 1:Nb_cov) {
          for (f in 1:Nb_cov_clus) {
            pI = Mat_prev[1,a]
            pC = Mat_prev[2,a]
            k = k_para[c]
            m = m_para[d]
            icc = icc_para[b]
            nb_cov = cov_para[e]
            nb_cov_clus = cov_clus_para[f]
            p_cov = rep(pC,max(cov_para))
            p_cov_clus = rep(pC,max(cov_clus_para))
            OR_cov = as.vector(unlist(list_OR))
            
            # Vecteur OR pour les cluster-level covariables 
            if(nb_cov_clus == 0){OR_cov_clus = c(0)}
            if(nb_cov_clus == 1){
              if(icc == 0.05 || icc == 0.01){OR_cov_clus = c(1.2)}else{OR_cov_clus = c(1.1)}
            }
            if(nb_cov_clus == 2){
              if(icc == 0.05 || icc == 0.01){OR_cov_clus = c(1.2,1)}else{OR_cov_clus = c(1.1,1)}
            }
            
            
            assign(paste("Scenario",num,sep = "_"),list(pI,
                                                        pC,
                                                        k,
                                                        m,
                                                        icc,
                                                        nb_cov,
                                                        nb_cov_clus,
                                                        p_cov,
                                                        p_cov_clus,
                                                        OR_cov,
                                                        OR_cov_clus))
            num = num+1
          }
        }
      }
    }
  }
}


Nb_scen = num - 1



### SAVE THE WORKSPACE TO CREATE THE 'Workspace_Fun' ----
# We save it and then we will download it in each R files 
# !care it will be save in the Base_file path

Workspace_name = "Workspace_clus_test.RData"

save.image(paste(Base_file,Workspace_name,sep = "/"))
