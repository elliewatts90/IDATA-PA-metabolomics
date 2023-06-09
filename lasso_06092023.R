# Script to create LASSO models splitting by training and testing set and evaluate model performance - created JUNE 2023 by ELW

setwd("C:/Users/wattsel/OneDrive - National Institutes of Health/iDATA/Data/")

# # Function to load packages 
install_load <- function(packages){
  for (p in packages) {
    if (p %in% rownames(installed.packages())) {
      library(p, character.only=TRUE)
    } else {
      install.packages(p)
      library(p,character.only = TRUE)
    }
  }
}

cbind.fill <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

# function to split dataframe into training and testing set
splitdf <- function(dataframe, seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  index <- 1:nrow(dataframe)
  trainindex <- sample(index, trunc(length(index)/5))
  trainset <- dataframe[-trainindex,]
  testset <- dataframe[trainindex,]
  list(trainset=trainset,testset=testset)
} 

lasso_function <- function(x.1, y.1, x.2, y.2, adj_n=0 , nfold=10) {
  
  # Set remove penalty if adjustments > 0
  if (adj_n > 0){ 
    pfac <- rep(1, ncol(x.1))
    n_col <- ncol(x.1)
    i <-n_col + 1 - adj_n
    pfac[c(i:n_col)] <- 0
  } else {
    pfac <-rep (1, ncol(x.1))
  }
  
  npf       <- ceiling(nrow(x.1)/nfold)
  foldid    <- rep(1:nfold,each=npf)[1:nrow(x.1)]
  
  fit1 <- glmnet(x = x.1, y = y.1 , family='gaussian')
  lambdanew = seq(max(fit1$lambda),min(fit1$lambda),length.out=10000)
  
  #' cross validation (determination of lambda and the optimal number of lipid species to keep)
  cv.fit1 <- cv.glmnet(x = x.1, y = y.1, family = "gaussian", alpha = 1, foldid = foldid, penalty.factor = pfac, lambda= lambdanew) 
  coefficients <- coef(cv.fit1) %>% as.matrix() %>% as.data.frame() %>% filter(s1 != 0) %>% dplyr::rename(beta = s1) %>%
    tibble::rownames_to_column("Variable")
  
  #' fit to new data 
  prediction_model <- predict(cv.fit1, newx = x.2, s = cv.fit1$lambda.min) 
  
  # R-squared
  #calculate SST and SSE
  sst <- sum((y.2 - mean(y.2))^2)
  sse <- sum((prediction_model - y.2)^2)
  
  #find R-Squared
  rsq <- 1 - sse/sst
  rsq<-cbind("r-sq", rsq) %>% as.data.frame() %>% dplyr::rename("Variable" = V1, beta="rsq")
  combined<-rbind(coefficients, rsq)
  
  return(combined)
}

`%ni%` <- Negate(`%in%`)

package_list<-c("glmnet","dplyr", "data.table", "fastDummies", "logr", "tibble", "stringr")
install_load(package_list)

## Load data 
# Loop through serum, 24-hr urine and 1st morning void urine to estimate partial correlations
 dta_list <- c("serum","urine",'fmv')

# ACTIVITY LIST
activity_var <- c(
  # TEE, PAEE, PAL - DLW
  "dlw_tee_avg_ln_sd", "dlw_paee_avg_ln_sd", "dlw_pal_avg_ln_sd", 
  
  # Total active time (questionnaire and devices)
  "ap_methrs_ln_sd", "ag_soj3x_methrs_ln_sd", "a24_total_methrs_ln_sd",
  
  # Total active time by intensity 
  "a24_sed_hrs_ln_sd", "ag_soj3x_sed_hrs_ln_sd", "ap_sed_hrs_ln_sd",
  "a24_light_hrs_ln_sd", "ag_soj3x_light_hrs_ln_sd", "ap_light_hrs_ln_sd",
  "a24_mod_hrs_ln_sd", "ag_soj3x_mod_hrs_ln_sd",
  "a24_vig_hrs_ln_sd", "ag_soj3x_vig_hrs_ln_sd",
  "a24_mvpa_hrs_ln_sd", "ag_soj3x_mvpa_hrs_ln_sd", "ap_mvpa_hrs_ln_sd" )


for (s in dta_list ){
  
  full_dat<-fread(paste0("C:/Users/wattsel/OneDrive - National Institutes of Health/iDATA/Data/PA_avg_met_avg_",s,".csv"))
  
for (i in activity_var) {

  # drop participants with missing variable 
  dat <- full_dat %>% filter(!is.na(!!rlang::sym(i))) 
  
  # drop participants with missing body comp data
  dat <- dat %>% filter(!is.na(pct_fat_avg) )
  
  dat$raceG<-as.numeric(dat$raceG)
  dat$sex<-as.numeric(dat$sex)
  dat$smok_curr<-as.numeric(dat$smok_curr)
 
  set.seed(123)
  
  split <- splitdf(dat, seed=123)
  
  # N observations in each data frame
  lapply(split,nrow)
  
  #view the first few columns in each data frame
  lapply(split,head)
  
  # save the training and testing sets as data frames
  training = split$trainset
  testing = split$testset
  
### 1 ) select metabolites only 
x.1 <- training %>%
  dplyr::select((matches("[A-Z]", ignore.case = FALSE)| starts_with("_") )) %>% 
  dplyr::select(-c(AGE, ageG, raceG, COTININE )) %>% data.matrix

x.2 <- testing %>%
  dplyr::select((matches("[A-Z]", ignore.case = FALSE)| starts_with("_") )) %>% 
  dplyr::select(-c(AGE, ageG, raceG, COTININE )) %>% data.matrix

adjustment_var <- c("AGE","pct_fat_avg","raceG", "sex", "smok_curr")
adjustments_training <- training %>% select(match(adjustment_var, names(training)))
adjustments_testing <- testing %>% select(match(adjustment_var, names(testing)))

x.1_adj <- cbind(x.1, adjustments_training) %>% data.matrix()
x.2_adj <- cbind(x.2, adjustments_testing) %>% data.matrix()

#create target variables
y.1 <- training %>% select(paste0(i)) %>% data.matrix()
y.2 <- testing %>% select(paste0(i)) %>% data.matrix()

# no adjustments 
lasso_model_noadj <- lasso_function(x.1 = x.1, y.1 = y.1, x.2 = x.2, y.2 = y.2)

# adjusting for "adjustment_var". NB- update adjustment number if changes 
lasso_model_adj <- lasso_function(x.1 = x.1_adj, y.1=y.1, x.2=x.2_adj, y.2=y.2, adj_n=5 )

combined<-cbind.fill(lasso_model_noadj, lasso_model_adj)

fwrite(combined, paste0("C:/Users/wattsel/OneDrive - National Institutes of Health/iDATA/Results/LASSO_results_",s,"_",i,".csv"))

}
}          
           



