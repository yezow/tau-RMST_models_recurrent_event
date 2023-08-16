##################
#Obtain initial estimated paras
##################
MI_impute_0 <- function(){
  #step1:obtain initial estimated paras
  complete_data<-longi_data_long_censored[-row_MI,]
  complete_data$I_tau <- ifelse(complete_data$censored_measurement==tau,1,0)
  
  fm <- geem2(I_tau ~ tr+age+male+fev1pp+smoker+center_2+center_3+center_4+center_5+center_6, id=subject, data =complete_data, family = FunList,
              corstr="userdefined",corr.mat=mycorr.mat,scale.fix = TRUE,waves=complete_data$waves )
  
  complete_data<-complete_data[-which(complete_data$censored_measurement==tau),]
  complete_data$meas_tau<-complete_data$censored_measurement/tau
  
  respbeta <- geem2(meas_tau ~ tr+age+male+fev1pp+smoker+center_2+center_3+center_4+center_5+center_6, id=subject, data = complete_data, family = FunList,
                    corstr="userdefined",corr.mat=mycorr.mat,waves=complete_data$waves) #Z1+Z2 for simulation
  coef_0 <- c(respbeta$beta,fm$beta,(1/respbeta$phi-1)) ##****add one more estimate
  return(coef_0)
}



#######################
#ES estimation (get converged para estimates for MI)
#######################
ES_converge_example <- function(){
  initial_para <- MI_impute_0() #get initial para estimates
  mu_est<- inv.logit(initial_para[1]+Z_mu%*%initial_para[2:(ncol(Z_mu)+1)])
  pi_est<- inv.logit(initial_para[(ncol(Z_mu)+2)]+Z_pi%*%initial_para[(ncol(Z_mu)+3):(ncol(Z_mu)+ncol(Z_pi)+2)])
  alpha_b <- mu_est*initial_para[ncol(Z_mu)+ncol(Z_pi)+3]
  beta_b <- (1-mu_est)*initial_para[ncol(Z_mu)+ncol(Z_pi)+3]
  
  diff_coef_converge<-3
  while (diff_coef_converge>1e-4 ) { 
    ####calculate expectation of b_i(t)
    E_b_i<- as.numeric(lapply(row_MI,function(x) pi_est[x]/(pi_est[x]+(1-pi_est[x])*(1-pbeta(longi_data_long_censored$censored_measurement[x]/tau, alpha_b[x],beta_b[x])))))
    w<-as.numeric(lapply(1:nrow(longi_data_long_censored), function(x) ifelse(longi_data_long_censored$censored_measurement[x]>=tau,1,0)))
    w[row_MI]<-E_b_i
    
    ####calculate expectation of y_i(t) 
    E_y_i <- numeric()
    for ( j in row_MI){
      alpha_i = alpha_b[j]
      beta_i = beta_b[j]
      integrand <- function(x) {x*x^(alpha_i-1)*(1-x)^(beta_i-1)/beta(alpha_i,beta_i)}
      integrand_1 <- function(x) {x^(alpha_i-1)*(1-x)^(beta_i-1)/beta(alpha_i,beta_i)}
      E_y_i <-append(E_y_i,as.numeric(integrate(integrand, lower = longi_data_long_censored$censored_measurement[j]/tau, upper = 1)[1])/
                       as.numeric(integrate(integrand_1, lower = longi_data_long_censored$censored_measurement[j]/tau, upper = 1)[1]))
    }
    v<-longi_data_long_censored$censored_measurement/tau
    v[row_MI]<-E_y_i
    
    #######estimate current model parameters
    data_ES <- longi_data_long_censored
    data_ES$w <- w
    
    #######Please edit models as needed
    fm<-geem2(w ~ tr+age+male+fev1pp+smoker+center_2+center_3+center_4+center_5+center_6, id=subject, data = data_ES,
              family=FunList,corstr="userdefined",corr.mat=mycorr.mat,scale.fix = TRUE)

    respbeta<-geem2(v ~ tr+age+male+fev1pp+smoker+center_2+center_3+center_4+center_5+center_6, id=subject,  data = data_ES,
                    family=FunList, corstr="userdefined",corr.mat=mycorr.mat,weights = 1-w) #geeglm can handle weights but doesn't give correct correlation parameters using binomial family
    
    coef_converge<- c(respbeta$beta,fm$beta,(1/respbeta$phi-1),respbeta$alpha,fm$alpha)#alpha is correlation parameter
    diff_coef_converge <- max(abs(coef_converge[1:(2+ncol(Z_mu)+ncol(Z_pi))]-initial_para[1:(2+ncol(Z_mu)+ncol(Z_pi))]))
    diff_coef_coverge_pi <- max(abs(coef_converge[(2+ncol(Z_mu)):(2+ncol(Z_mu)+ncol(Z_pi))]-initial_para[(2+ncol(Z_mu)):(2+ncol(Z_mu)+ncol(Z_pi))]))
    diff_coef_coverge_mu <- max(abs(coef_converge[1:(1+ncol(Z_mu))]-initial_para[1:(1+ncol(Z_mu))]))
    print(diff_coef_converge)
    initial_para <- coef_converge
    mu_est<- inv.logit(initial_para[1]+Z_mu%*%initial_para[2:(ncol(Z_mu)+1)])
    pi_est<- inv.logit(initial_para[(ncol(Z_mu)+2)]+Z_pi%*%initial_para[(ncol(Z_mu)+3):(ncol(Z_mu)+ncol(Z_pi)+2)])
    alpha_b <- mu_est*initial_para[length(initial_para)]
    beta_b <- (1-mu_est)*initial_para[length(initial_para)]
  }
  #update w and E(Y|Y<X) using the last EM para estimates
  E_b_i<- as.numeric(lapply(row_MI,function(x) pi_est[x]/(pi_est[x]+(1-pi_est[x])*(1-pbeta(longi_data_long_censored$censored_measurement[x]/tau, alpha_b[x],beta_b[x])))))
  w<-as.numeric(lapply(1:nrow(longi_data_long_censored), function(x) ifelse(longi_data_long_censored$censored_measurement[x]>=tau,1,0)))
  w[row_MI]<-E_b_i
  
  E_y_i <- numeric()
  for ( j in row_MI){
    alpha_i = alpha_b[j]
    beta_i = beta_b[j]
    integrand <- function(x) {x*x^(alpha_i-1)*(1-x)^(beta_i-1)/beta(alpha_i,beta_i)}
    integrand_1 <- function(x) {x^(alpha_i-1)*(1-x)^(beta_i-1)/beta(alpha_i,beta_i)}
    E_y_i <-append(E_y_i,as.numeric(integrate(integrand, lower = longi_data_long_censored$censored_measurement[j]/tau, upper = 1)[1])/
                     as.numeric(integrate(integrand_1, lower = longi_data_long_censored$censored_measurement[j]/tau, upper = 1)[1]))
  }
  v<-longi_data_long_censored$censored_measurement/tau
  v[row_MI]<-E_y_i
  
  return(list(coef_converge,w,v))
}


ES_var_example <- function(coef_ES,w,v){
  ####Calculate the derivative of mu w.r.t alpha and pi w.r.t beta for each individual
  e_alpha_Z<-exp(coef_ES[1]+Z_mu%*%coef_ES[2:(ncol(Z_mu)+1)])
  e_beta_Z<-exp(coef_ES[ncol(Z_mu)+2]+Z_pi%*%coef_ES[(ncol(Z_mu)+3):(ncol(Z_mu)+ncol(Z_pi)+2)])
  Z_mu_1<-cbind(1,Z_mu)
  Z_pi_1<-cbind(1,Z_pi)
  deriv_mu_alpha_i<-list()
  deriv_pi_beta_i<-list()
  for(i in 1:nrow(Z_mu_1)) {                
    deriv_mu_alpha_i[[i]] <- Z_mu_1[i,]*e_alpha_Z[i]/(1+e_alpha_Z[i])^2       
    deriv_pi_beta_i[[i]] <- Z_pi_1[i,]*e_beta_Z[i]/(1+e_beta_Z[i])^2     
  }
  
  ####Calculate each component of M, which is quasi-score function
  M <- matrix(0,ncol=ncol(Z_mu)+ncol(Z_pi)+2,nrow = ncol(Z_mu)+ncol(Z_pi)+2)
  B1 <- matrix(0,ncol=ncol(Z_mu)+ncol(Z_pi)+2,nrow = ncol(Z_mu)+ncol(Z_pi)+2)
  B2 <- matrix(0,ncol=ncol(Z_mu)+ncol(Z_pi)+2,nrow = ncol(Z_mu)+ncol(Z_pi)+2)
  for ( i in unique(longi_data_long_censored$subject)){#nrow(Z_mu_1) 
    ##p*n_i derivative matrix for alpha
    obs_i<-which(longi_data_long_censored$subject==i) #which row is for indiviudal i
    deriv_mu_alpha_matrix<-matrix(unlist(deriv_mu_alpha_i[obs_i]),nrow=dim(Z_mu_1)[2],ncol=length(obs_i))
    ##Estimated variance matrix of Y_i
    estimated_mu<-inv.logit(coef_ES[1]+Z_mu[obs_i,]%*%coef_ES[2:(ncol(Z_mu)+1)])
    estimated_nu <- coef_ES[ncol(Z_mu)+ncol(Z_pi)+3]
    A_mu <- ifelse(length(obs_i)>1,list(diag(as.numeric((estimated_nu+1)^(-1)*estimated_mu*(1-estimated_mu)))),list((estimated_nu+1)^(-1)*estimated_mu*(1-estimated_mu)))[[1]]
    corr_para<-coef_ES[(ncol(Z_mu)+ncol(Z_pi)+4):(ncol(Z_mu)+ncol(Z_pi)+4+3)]
    R_mu <- matrix(corr_para[b],b,b)
    diag(R_mu)<- 1
    for(p in 1:b){for (q in 1:b){if (abs(p-q)==1){R_mu[p,q]=corr_para[2]} else if (abs(p-q)==2){R_mu[p,q]=corr_para[3]}}}
    R_mu<-R_mu[1:length(obs_i),1:length(obs_i)]
    var_Y<-sqrt(A_mu)%*%R_mu%*%sqrt(A_mu)
    
    ##Estimated quasi score function of Y
    w_i <- w[longi_data_long_censored$subject==i]
    v_i <- v[longi_data_long_censored$subject==i]
    quasi_Y<-deriv_mu_alpha_matrix%*%solve(var_Y)%*%
      ifelse(length(obs_i)>1,list(diag(1-w_i)),list(1-w_i))[[1]]%*%
      (v_i-estimated_mu)
    
    ##p*n_i derivative matrix for beta
    deriv_pi_beta_matrix<-matrix(unlist(deriv_pi_beta_i[obs_i]),nrow=dim(Z_mu_1)[2],ncol=length(obs_i))
    
    ##Estimated variance matrix of B_i
    estimated_pi<-inv.logit(coef_ES[ncol(Z_mu)+2]+Z_pi[obs_i,]%*%coef_ES[(ncol(Z_mu)+3):(ncol(Z_mu)+ncol(Z_pi)+2)])
    A_pi <-ifelse(length(obs_i)>1,list(diag(as.numeric(estimated_pi*(1-estimated_pi)))),list((estimated_pi*(1-estimated_pi))))[[1]] 
    corr_para<-coef_ES[(ncol(Z_mu)+ncol(Z_pi)+8):(ncol(Z_mu)+ncol(Z_pi)+8+3)]
    R_pi <- matrix(corr_para[b],b,b)
    diag(R_pi)<- 1
    for(p in 1:b){for (q in 1:b){if (abs(p-q)==1){R_pi[p,q]=corr_para[2]} else if (abs(p-q)==2){R_pi[p,q]=corr_para[3]}}}
    R_pi<-R_pi[1:length(obs_i),1:length(obs_i)]
    var_B<-sqrt(A_pi)%*%R_pi%*%sqrt(A_pi)
    
    ##Estimated quasi score function of Y
    quasi_B<-deriv_pi_beta_matrix%*%solve(var_B)%*%
      (w_i-estimated_pi)
    
    M_i<- rbind(quasi_Y,quasi_B)%*%t(rbind(quasi_Y,quasi_B))
    M<- M+M_i
    
    ##B_1i in hessian matrix
    component_1 <- -1*deriv_mu_alpha_matrix%*%solve(var_Y)%*%
      ifelse(length(obs_i)>1,list(diag(1-w_i)),list(1-w_i))[[1]]%*%t(deriv_mu_alpha_matrix)
    component_2 <- -1*deriv_pi_beta_matrix%*%solve(var_B)%*%t(deriv_pi_beta_matrix)
    B_1i<-matrix(0,ncol=ncol(Z_mu)+ncol(Z_pi)+2,nrow = ncol(Z_mu)+ncol(Z_pi)+2)
    B_1i[1:(ncol(Z_mu)+1),1:(ncol(Z_mu)+1)] <- component_1
    B_1i[(ncol(Z_mu)+2):(ncol(Z_mu)+ncol(Z_pi)+2),(ncol(Z_mu)+2):(ncol(Z_mu)+ncol(Z_pi)+2)]<-component_2
    
    ##B_2i in hessian matrix
    component_3 <- -1*deriv_mu_alpha_matrix%*%solve(var_Y)%*%
      ifelse(length(obs_i)>1,list(diag(as.numeric(v_i-estimated_mu))),list(v_i-estimated_mu))[[1]]
    component_4 <- deriv_pi_beta_matrix%*%solve(var_B)
    component_5 <- sqrt(ifelse(length(obs_i)>1,list(diag(w_i*(1-w_i))),list(w_i*(1-w_i)))[[1]])%*%R_pi%*%sqrt(ifelse(length(obs_i)>1,list(diag(w_i*(1-w_i))),list(w_i*(1-w_i)))[[1]])
    B_2i <- rbind(component_3,component_4) %*% component_5 %*% t(rbind(component_3,component_4))
    
    B1<- B1+B_1i
    B2<- B2+B_2i
  }
  var_ES<-solve(B1+B2)%*%M%*%solve(B1+B2)
  return(var_ES)
}

####################################
#MI Step1: Risk set
####################################
MI_impute_risk_set <- function(coef_MI){

  LINP1_MI <- coef_MI[ncol(Z_mu)+2]  + Z_pi%*%coef_MI[(ncol(Z_mu)+3):(ncol(Z_mu)+ncol(Z_pi)+2)]
  PI_MI <- inv.logit(LINP1_MI)
  LINP_MI <- coef_MI[1]+Z_mu%*%coef_MI[2:(1+ncol(Z_mu))]
  MU_MI <- inv.logit(LINP_MI)
  
  Risk_set <- list()
  epsilon_set<-numeric()
  m <- 0
  
  for (i in subj_imputed){ #from subj id
    risk_row <- numeric()
    epsilon <- 0.05 #You can choose another reasonable number
    set <- which(loop_num==i)
    m <- m+1
    row_index_i<-which(longi_data_long_censored$subject==i)[last_window[m]]
    imputation_window <- (last_window[m]-length(rows_imputed_subj[[m]])+1):last_window[m]
    subj_imputed_prior <- longi_data_long_censored$censored_measurement[longi_data_long_censored$subject==i][-imputation_window] #history obs for imputed subject i
    subj_imputed_prior_mean <- mean(subj_imputed_prior) #mean of obs in prior windows for imputed subject i
    
    while (length(risk_row) < 15 & epsilon <=0.5){
      for (j in loop_num[-set] ){ #j is subj id excluded i
        row_index_j<-which(longi_data_long_censored$subject==j)[last_window[m]]
        abs_diff <- max(abs(PI_MI[row_index_j]-PI_MI[row_index_i]),abs(MU_MI[row_index_j]-MU_MI[row_index_i])) #PI_MU[loop_num=1] also works for example because the order of PI_MU is 1:1111->1:1111...
        last_window_obs_time <- longi_data_long_censored$censored_measurement[row_index_j]#observed time of candidate in the last window of the subject required imputation
        
        #Criteria of selecting subjects in risk set
        if (abs_diff <= epsilon&!is.na(last_window_obs_time)&
            last_window_obs_time>censor_last_window[m]){
          imputation_window_candidate <- ((last_window[m]-length(rows_imputed_subj[[m]])+1):length(longi_data_long_censored$censored_measurement[longi_data_long_censored$subject==j]))
          subj_risk_set_prior <- longi_data_long_censored$censored_measurement[longi_data_long_censored$subject==j][-imputation_window_candidate] #obs history for subject j in the risk set
          subj_risk_set_prior_mean <- mean(subj_risk_set_prior) #mean of obs in prior windows for subject j in risk set
          if (length(subj_risk_set_prior)>0 & all(tau %in% subj_risk_set_prior,tau %in% subj_imputed_prior)){
            risk_row <- append(risk_row,which(loop_num==j))
          } else if (length(subj_risk_set_prior)>0 & all(!(tau %in% subj_risk_set_prior),!(tau %in% subj_imputed_prior)) &
                     abs(subj_risk_set_prior_mean-subj_imputed_prior_mean)< tau*0.6 ){ #You can specify your own number instead of 0.6
            risk_row <- append(risk_row,which(loop_num==j))
          } else if (length(subj_risk_set_prior)==0){
            risk_row <- append(risk_row,which(loop_num==j))
          }
        }
      }
      set <- append(which(loop_num==i),risk_row)
      epsilon <- epsilon+0.005
    }
    
    
    ###Second restriction: the last event time should be failure or attain tau, otherwise increase sample size
    event_list<-numeric()
    delta_list<-numeric()
    for (item in loop_num[risk_row]){ #still subj id
      event_list <- append(event_list,longi_data_long_censored$censored_measurement[longi_data_long_censored$subject==item][last_window[m]])
      delta_list <- append(delta_list,longi_data_long_censored$Delta[longi_data_long_censored$subject==item][last_window[m]])
    }#get the last observed time in the risk set
    last_event_time_in_R <- max(event_list)
    last_censor_in_R <- delta_list[which(event_list==last_event_time_in_R)] #get the last censoring status in risk set
    while(any(last_censor_in_R==0) &last_event_time_in_R<tau){
      epsilon <- epsilon+0.001
      for (j in loop_num[-set] ){
        row_index_j<-which(longi_data_long_censored$subject==j)[last_window[m]]
        abs_diff <- max(abs(PI_MI[row_index_j]-PI_MI[row_index_i]),abs(MU_MI[row_index_j]-MU_MI[row_index_i])) 
        last_window_obs_time <- longi_data_long_censored$censored_measurement[row_index_j]#observed time of candidate in the last window of the subject required imputation
        
        #criteria of selecting subjects in risk set
        if (abs_diff <= epsilon&!is.na(last_window_obs_time)
            &last_window_obs_time>censor_last_window[m]){
          imputation_window_candidate <- ((last_window[m]-length(rows_imputed_subj[[m]])+1):length(longi_data_long_censored$censored_measurement[longi_data_long_censored$subject==j]))
          subj_risk_set_prior <- longi_data_long_censored$censored_measurement[longi_data_long_censored$subject==j][-imputation_window_candidate]#obs history for subject j in the risk set
          subj_risk_set_prior_mean <- mean(subj_risk_set_prior) #mean of obs in prior windows for subject j in risk set
          if (length(subj_risk_set_prior)>0 & all(tau %in% subj_risk_set_prior,tau %in% subj_imputed_prior)){
            risk_row <- append(risk_row,which(loop_num==j))
          } else if (length(subj_risk_set_prior)>0 & all(!(tau %in% subj_risk_set_prior),!(tau %in% subj_imputed_prior))){
            risk_row <- append(risk_row,which(loop_num==j))
          } else if (length(subj_risk_set_prior)==0){
            risk_row <- append(risk_row,which(loop_num==j))
          }
        }
      }
      
      #check the second restriction again
      event_list<-numeric()
      delta_list<-numeric()
      for (item in loop_num[risk_row]){
        event_list <- append(event_list,longi_data_long_censored$censored_measurement[longi_data_long_censored$subject==item][last_window[m]])
        delta_list <- append(delta_list,longi_data_long_censored$Delta[longi_data_long_censored$subject==item][last_window[m]])
      }#get the last observed time in the risk set
      last_event_time_in_R <- max(event_list)
      last_censor_in_R <- delta_list[which(event_list==last_event_time_in_R)] #get the last censoring status in risk set
      set <- append(which(loop_num==i),risk_row)
      print("add more subjs in risk set")
    }
    
    epsilon_set[m] <- epsilon
    Risk_set[[m]] <- risk_row
  }
  return(Risk_set)
}


##########################################################################
#MI Step 2: IPTI imputation (impute values for the last window)
##########################################################################
MI_impute_I_Y<- function(Risk_set){
  I_i_impute<-numeric() 
  Y_j_impute <- numeric()
  for (i in 1:length(subj_imputed)){
    risk_set <- Risk_set[[i]]
    risk_row <- numeric()
    for (j in risk_set){
      risk_row <- append(risk_row,which(longi_data_long_censored$subject==loop_num[j])[last_window[i]]) #risk_set[j]
    } #risk_row is the rows used in KM curve
    
    data_impute <- longi_data_long_censored[risk_row,]
    KM <- survfit(Surv(censored_measurement, Delta) ~ 1,  type="kaplan-meier", conf.type="log", data=data_impute)
    
    #inverse transform
    v=0
    while (v<=censor_last_window[i]){
      u <- runif(1,0,1)
      v <- KM$time[KM$surv<=u][1]
      if(v>=tau){
        I_i_impute[i]<-1
      }else{
        I_i_impute[i]<-0
        Y_j_impute<-append(Y_j_impute,v)
      }
    }
    #print(i)
  }
  return(list(I_i_impute,Y_j_impute))
}

########################################################################
# impute and get the complete dataset using the result from MI_impute_I_Y
########################################################################
impute_func <- function(result_last_win){
  data_complete <- longi_data_long_censored
  data_complete$censored_measurement[last_rows_imputed][result_last_win[[1]]==1] <- tau #impute for the last window
  data_complete$censored_measurement[last_rows_imputed][result_last_win[[1]]==0] <- result_last_win[[2]] #impute for the last window
  #impute for other windows for subjects
  for(i in 1:length(subj_imputed)){
    impute_num <- length(rows_imputed_subj[[i]]) #the number of windows need imputation
    while(impute_num>1){
      last_impute_value <- data_complete$censored_measurement[rows_imputed_subj[[i]][impute_num]]
      data_complete$censored_measurement[rows_imputed_subj[[i]][impute_num-1]]<-ifelse(last_impute_value>=(tau-a),tau,a+last_impute_value)
      impute_num <- impute_num-1
    }
  }
  return(data_complete)
}

####################################################################
#MI Step3: Final step: Get final estimated parameters and its var matrix
####################################################################

MI_impute <- function(coef_ES){
  #repeat above steps until parameter estimators converge
  coef_MI <- coef_ES 
  
  #create a new complete survival time used to fit model
  Risk_set_1 <- MI_impute_risk_set(coef_MI)
  coef_MI_para_10 <- list()
  var_MI_para_10 <- list()
  
  #repeat bernoulli and uniform steps and get para estimators using 10 complete datasets
  M_imputed_results <-foreach (i = 1:10, .combine=c) %dopar% {
    impute_result_1 <- MI_impute_I_Y(Risk_set_1)
    data_no_NA <- impute_func(impute_result_1)
    
    #re-estimate parameters
    data_no_NA$I_tau <- ifelse(data_no_NA$censored_measurement==tau,1,0)
    fm <- geem2(I_tau ~ tr+age+male+fev1pp+smoker+center_2+center_3+center_4+center_5+center_6,
                id=subject, data =data_no_NA, family = FunList,
                corstr="userdefined",corr.mat=mycorr.mat,scale.fix = TRUE,waves=data_no_NA$waves )
  
    data_no_NA_1<-data_no_NA[-which(data_no_NA$censored_measurement==tau),]
    data_no_NA_1$meas_tau<-data_no_NA_1$censored_measurement/tau
    respbeta <- geem2(meas_tau ~tr+age+male+fev1pp+smoker+center_2+center_3+center_4+center_5+center_6, 
                      id=subject, data = data_no_NA_1, family = FunList,
                      corstr="userdefined",corr.mat=mycorr.mat,waves=data_no_NA_1$waves)
    coef_MI_para_10 <- c(respbeta$beta,fm$beta)
    var_MI_para_10<-matrix(0,nrow=(2+ncol(Z_mu)+ncol(Z_pi)),ncol=(2+ncol(Z_mu)+ncol(Z_pi)))#7  #create a combined vcov matrix for calculation of var of restricted mean
    var_MI_para_10[1:(1+ncol(Z_mu)),1:(1+ncol(Z_mu))] <- as.matrix(respbeta$var)
    var_MI_para_10[(2+ncol(Z_mu)):(2+ncol(Z_mu)+ncol(Z_pi)),(2+ncol(Z_mu)):(2+ncol(Z_mu)+ncol(Z_pi))]<-as.matrix(fm$var)
    list(coef_MI_para_10,var_MI_para_10,data_no_NA$censored_measurement)
  }
  coef_MI_para_10<-M_imputed_results[c(1,4,7,10,13,16,19,22,25,28)]
  var_MI_para_10<-M_imputed_results[c(2,5,8,11,14,17,20,23,26,29)]
  MI_para_mean<-MI_combine(coef_MI_para_10,var_MI_para_10,10)[[1]]
  vcov_mean <- MI_combine(coef_MI_para_10,var_MI_para_10,10)[[2]]
  imputed_data<-M_imputed_results[c(3,6,9,12,15,18,21,24,27,30)]
  return(list(MI_para_mean,vcov_mean,imputed_data))
}


######################################################################
#Other used functions in MI method: combine multiple imputation results
######################################################################
MI_combine <- function(coef_list,vcov_list,m){
  para_data<-data.frame(coef_list)
  para_bar<-rowMeans(para_data)
  vcov_bar <- vcov_list[[1]]
  B <- (coef_list[[1]]-para_bar)%*%t(coef_list[[1]]-para_bar) #calculate the variance between the MI data sets
  for (i in 2:m){
    vcov_bar <- vcov_bar+vcov_list[[i]]
    B <- B+(coef_list[[i]]-para_bar)%*%t(coef_list[[i]]-para_bar)
  }
  vcov_bar <- vcov_bar/m
  B <-B/(m-1)
  
  #calculate the total variance
  comb.vcov <- vcov_bar+(1+1/m)*B
  return(list(para_bar,comb.vcov))
}

#################################################
#calculate the variance of restricted mean
#################################################
Var_restricted_mean_func <- function(coef_vcov,coef,Z_mu,Z_pi){
  beta_vcov <- coef_vcov[(ncol(Z_mu)+2):(ncol(Z_pi)+ncol(Z_mu)+2),(ncol(Z_mu)+2):(ncol(Z_pi)+ncol(Z_mu)+2)]
  alpha_vcov<- coef_vcov[1:(ncol(Z_mu)+1),1:(ncol(Z_mu)+1)]
  var_restricted_mean <- numeric()
  
  for (i in 1:nrow(Z_pi)){
    #calculate variance of pi
    logit_mu_est <- coef[1]+Z_mu[i,]%*%coef[2:(ncol(Z_mu)+1)]
    logit_pi_est <- coef[(ncol(Z_mu)+2)]+Z_pi[i,]%*%coef[(ncol(Z_mu)+3):(ncol(Z_mu)+ncol(Z_pi)+2)]
    X_pi <- matrix(c(1,Z_pi[i,]), nrow = 1+ncol(Z_pi))
    var_logit_pi <- t(X_pi) %*% beta_vcov %*% X_pi
    var_pi <- var_logit_pi*(exp(-logit_pi_est)^2)/(1+exp(-logit_pi_est))^4
    
    #calculate the first term in variance of restricted mean
    var_term_1<-tau^2*(1-inv.logit(logit_mu_est))^2*var_pi
    
    #calculate variance of logit of pi
    X_mu <- matrix(c(1,Z_mu[i,]),nrow=1+ncol(Z_mu))
    var_logit_mu <- t(X_mu) %*% alpha_vcov %*% X_mu
    var_mu <- var_logit_mu*(exp(-logit_mu_est)^2)/(1+exp(-logit_mu_est))^4
    
    #calculate the derivative of transformation function
    g_deriv <- tau*(1-1/(1+exp(-logit_pi_est)))
    
    #calculate the second term in variance of restricted mean
    var_term_2<-g_deriv^2*var_mu
    
    var_restricted_mean[i] <- var_term_1+var_term_2
  }
  return(var_restricted_mean)
}

##########################################
#calculate the covariance of restricted mean
##########################################
Cov_restricted_mean_func <- function(coef_vcov,coef,Z_i_mu,Z_i_pi,Z_j_mu,Z_j_pi){
  beta_vcov <- coef_vcov[(ncol(Z_mu)+2):(ncol(Z_pi)+ncol(Z_mu)+2),(ncol(Z_mu)+2):(ncol(Z_pi)+ncol(Z_mu)+2)]
  alpha_vcov<- coef_vcov[1:(ncol(Z_mu)+1),1:(ncol(Z_mu)+1)]
  var_restricted_mean <- numeric()
  
  #calculate variance of pi and mu
  logit_mu_est_i <- coef[1]+Z_i_mu%*%coef[2:(ncol(Z_mu)+1)]
  logit_pi_est_i <- coef[(ncol(Z_mu)+2)]+Z_i_pi%*%coef[(ncol(Z_mu)+3):(ncol(Z_mu)+ncol(Z_pi)+2)]
  logit_mu_est_j <- coef[1]+Z_j_mu%*%coef[2:(ncol(Z_pi)+1)]
  logit_pi_est_j <- coef[(ncol(Z_mu)+2)]+Z_j_pi%*%coef[(ncol(Z_mu)+3):(ncol(Z_mu)+ncol(Z_pi)+2)]
  X_pi_i <- matrix(c(1,Z_i_pi), nrow = 1+ncol(Z_pi))
  X_pi_j <- matrix(c(1,Z_j_pi), nrow = 1+ncol(Z_pi))
  var_logit_pi <- t(X_pi_i) %*% beta_vcov %*% X_pi_j
  var_pi <- var_logit_pi*exp(-logit_pi_est_i)*exp(-logit_pi_est_j)/(1+exp(-logit_pi_est_i))^2/(1+exp(-logit_pi_est_j))^2
  
  X_mu_i <- matrix(c(1,Z_i_mu), nrow = 1+ncol(Z_mu))
  X_mu_j <- matrix(c(1,Z_j_mu), nrow = 1+ncol(Z_mu))
  var_logit_mu <- t(X_mu_i) %*% alpha_vcov %*% X_mu_j
  var_mu <- var_logit_mu*exp(-logit_mu_est_i)*exp(-logit_mu_est_j)/(1+exp(-logit_mu_est_i))^2/(1+exp(-logit_mu_est_j))^2
  
  #calculate the first term in variance of restricted mean
  var_term_1<-tau^2*(1-inv.logit(logit_mu_est_i))*(1-inv.logit(logit_mu_est_j))*var_pi
  
  #calculate the second term in variance of restricted mean
  var_term_2<-tau^2*(1-inv.logit(logit_pi_est_i))*(1-inv.logit(logit_pi_est_j))*var_mu
  
  cov_restricted_mean <- var_term_1+var_term_2
  return(cov_restricted_mean)
}


################################################
#fold change and its variance computation
################################################
fold_change <- function(coef){
  beta<-coef[2]
  exp_sum <- exp(sum(coef))
  fold<-(exp(beta)+exp_sum)/(1+exp_sum)
  return(fold)
}

fold_change_var<-function(coef,coef_vcov){
  exp_sum <- exp(sum(coef))
  beta<-coef[2]
  gradient_beta <- (exp(beta)+exp_sum)/(1+exp_sum)-(exp(beta)+exp_sum)*exp_sum/(1+exp_sum)^2
  gradient_intercept<-exp_sum/(1+exp_sum)-(exp(beta)+exp_sum)*exp_sum/(1+exp_sum)^2
  gradient <- matrix(c(gradient_intercept,gradient_beta),nrow=1)
  var <- gradient%*%coef_vcov%*%t(gradient)
  return(var)
}
