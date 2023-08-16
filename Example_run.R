library(fastDummies)
library(dplyr)
library(tidyr)
library(boot)
library(MASS)
library(purrr)
library(doBy)
library(bindata)
library(coda)
library(numDeriv)
library(survival)
library(gdata)
library(MuMIn)
library(matrixStats)
library(base)
library(geepack)
library(foreach)
library(doParallel)
library(mmmgee)
library(expm)
library(latex2exp)

registerDoParallel(6)
source("/Users/yizhuowang/Downloads/Yizhuo\ Wang/thesis/Paper2\ code/Github_data_code/Github.R")

#####################################################
# User-defined link and variance functions used in geem2 
#####################################################
LinkFun <- function(arg){log(arg/(1-arg))}
InvLink <- function(arg){exp(arg)/(1+exp(arg))}
InvLinkDeriv <- function(arg){exp(arg)/(1+exp(arg))^2}
VarFun <- function(arg){arg*(1-arg)}
FunList <- list(LinkFun, VarFun, InvLink, InvLinkDeriv)

#######################################
# example analysis
#######################################
data_1<- read.csv("/Users/yizhuowang/Downloads/Yizhuo\ Wang/thesis/Paper2\ code/Github_data_code/Example_data.csv")

##################################################################
#Specify follow-up windows
##################################################################
s=380
a=60
tau=180 #Same with paper1
t=seq(from=0, to=s-tau, by=a)
b=length(t)

##################################################################
#get T(t), time to the first event observed for each window starting at t
##################################################################
n <- nrow(data_1) #Number of subjects
X <- data_1$Days_In_Study #follow-up days
delta_X <- data_1$diedle380 #if died before 380 days
R_star <- array(cbind(data_1$Days_To_Onset1,data_1$Days_To_Onset2,data_1$Days_To_Onset3,data_1$Days_To_Onset4,data_1$Days_To_Onset5,data_1$Days_To_Onset6,data_1$Days_To_Onset7,data_1$Days_To_Onset8,data_1$Days_To_Onset9,data_1$Days_To_Onset10,data_1$Days_To_Onset11),c(n,11))

T_t=array(NA,c(n,b))
delta=array(NA,c(n,b))
for(j in 1:b)
{
  R_star_j=R_star-t[j] #recurrent events time from t[j]
  R_star_j[R_star_j<=0]=NA #recurrent events observed before time t[j]
  X_j=X-t[j] #total follow-up time from t[j]
  X_j[X_j<=0]=NA
  R_star_j_prime=ifelse(is.na(X_j), NA, apply(cbind(X_j,R_star_j), 1, 
                 function(x) {ifelse(sum(!is.na(x))>0, min(x, na.rm = T), NA)})) #time to next recurrent event
  T_t[,j]=R_star_j_prime
  delta[,j]=ifelse(R_star_j_prime==X_j, delta_X, 1) #time to next event time = remained follow up time then delta = diedled status, otherwise =1
}
min_T_tau <- apply(T_t, c(1,2), function(x) min(x, tau))

##################################################################
#Format to long data 
##################################################################
data_wide <- data.frame(cbind(1:n, data_1$tr, data_1$age, data_1$male, data_1$smoker, data_1$fev1pp, data_1$center,min_T_tau))
colnames(data_wide) <- c("ID","tr","age","male","smoker","fev1pp","center",paste0("min_T_tau_",1:b))
#Transform from wide to long
data_long <- gather(data_wide, tj, min_T_tau, paste0("min_T_tau_",1:b), factor_key=TRUE)
data_long <- data_long[!is.na(data_long$min_T_tau),] #remove missing
#Get delta
delta_long=array(delta)
data_long$delta=delta_long[!is.na(delta_long)]
#data_long$delta=ifelse(data_long$min_T_tau==tau,0,data_long$delta) #used in pseudo observation
data_long$delta_MI <- ifelse(data_long$min_T_tau==tau,1,data_long$delta) #if min_T_tau=tau,then this window doesn't need imputation

##################################################################
#define Toeplitz correlation structure
##################################################################
mycorr.mat <- matrix(4,b,b)
diag(mycorr.mat)<- 1
for(p in 1:b){for (q in 1:b){if (abs(p-q)==1){mycorr.mat[p,q]=2} else if (abs(p-q)==2){mycorr.mat[p,q]=3}}}

############################################################################################
#ES & MI
############################################################################################
set.seed(123)
#number of subjs
nsubj <- 1111

#variables
data_long$center<-factor(data_long$center)
longi_data_long_censored <- dummy_cols(data_long, select_columns = "center") #create dummy center variabls
longi_data_long_censored <- dummy_cols(longi_data_long_censored,select_columns = "tj")
longi_data_long_censored$waves[longi_data_long_censored$tj=="min_T_tau_1"]<-1 
longi_data_long_censored$waves[longi_data_long_censored$tj=="min_T_tau_2"]<-2
longi_data_long_censored$waves[longi_data_long_censored$tj=="min_T_tau_3"]<-3
longi_data_long_censored$waves[longi_data_long_censored$tj=="min_T_tau_4"]<-4

longi_data_long_censored<-longi_data_long_censored[order(longi_data_long_censored$ID),]

Z_mu <-as.matrix(longi_data_long_censored[,c("tr","age","male","fev1pp","smoker","center_2","center_3","center_4","center_5","center_6")]) 
Z_pi<-Z_mu

longi_data_long_censored<-longi_data_long_censored %>% #rename so that can use func
  rename(
    censored_measurement = min_T_tau,
    subject = ID,
    Delta = delta_MI
)

######################################################
#Get variables needed in MI functions
######################################################
row_MI <- which(longi_data_long_censored$Delta==0) #which rows need imputation
loop_num <- unique(longi_data_long_censored$subject) #all subj id in the dataset
subj_imputed <- unique(longi_data_long_censored$subject[row_MI]) #determine the subj id who needs to be imputed
m=0
censor_last_window<-numeric()
last_window<-numeric()
last_rows_imputed <- numeric()
rows_imputed_subj<-list()
for (i in subj_imputed){
  #get the censoring time at the last window for each subj needed to be imputed
  m = m+1
  T_last_window <- longi_data_long_censored$censored_measurement[longi_data_long_censored$subject==i]
  censor_last_window[m] <- tail(T_last_window,n=1)
  last_window[m] <- tail(which(T_last_window==censor_last_window[m]),n=1)
  rows_imputed_subj[[m]] <- which(longi_data_long_censored$subject==i & longi_data_long_censored$Delta==0) #determine which rows need imputation for subject i 
  last_rows_imputed <- append(last_rows_imputed,tail(rows_imputed_subj[[m]],n=1)) #determine last rows that need imputation
}

######################################################
#Fit model using ES algorithm
######################################################
ES_result <- ES_converge_example()
coef_ES <- ES_result[[1]] #Besides model coefficients (1:22), also get scale parameter and corr parameters estimates
var_ES <- ES_var_example(coef_ES,ES_result[[2]],ES_result[[3]])

#######Get fold change for specific covariate
fold_change_trt<-fold_change(coef_ES[c(1,2)]) 
var_fold_change_trt <- fold_change_var(coef_ES[c(1,2)],var_ES[c(1,2),c(1,2)])
fold_change_CI_lower<-fold_change_trt-1.96*sqrt(var_fold_change_trt)
fold_change_CI_upper<-fold_change_trt+1.96*sqrt(var_fold_change_trt)
fold_change_trt
fold_change_CI_lower
fold_change_CI_upper


######################################################
#Fit model using MI algorithm
######################################################
set.seed(123)
MI_result <- MI_impute(coef_ES) #get MI results
mu_est_MI <- inv.logit(MI_result[[1]][1]+Z_mu%*%c(MI_result[[1]][2:(ncol(Z_mu)+1)]))
pi_est_MI <- inv.logit(MI_result[[1]][ncol(Z_mu)+2]+Z_pi%*%c(MI_result[[1]][(ncol(Z_mu)+3):(ncol(Z_mu)+ncol(Z_pi)+2)]))
Mean_tau_T_est_MI <- mu_est_MI*tau*(1-pi_est_MI)+tau*pi_est_MI
CI_lower <- MI_result[[1]]-1.964739*sqrt(diag(MI_result[[2]]))
CI_upper <- MI_result[[1]]+1.964739*sqrt(diag(MI_result[[2]]))
var_restricted_mean_MI_est<-Var_restricted_mean_func(MI_result[[2]],MI_result[[1]],Z_mu,Z_pi)
coef_MI<-MI_result[[1]]
var_MI<-MI_result[[2]]

