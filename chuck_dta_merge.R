# SCRIPT CREATED TO merge in iDATA data from the various data sources - last edited June 9th 2023
# This creates separate serum, 24hr urine and FMV datasets 

# EW MERGE IN CHUCK'S generated data
### IDATA MERGE DLW AND METABOLOMICS DATA 
rm(list=ls())

## Set working directory
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

######## INSTALL PACKAGES
package_list<-c("dplyr", "data.table", "haven", "lubridate",  "fastDummies", "mice")
install_load(package_list)

# # fix plyr/dplyr issue
# if("dplyr" %in% (.packages())){
#   detach("package:dplyr", unload=TRUE) 
#   detach("package:plyr", unload=TRUE) 
# } 
# library(plyr)
# library(dplyr)

### FUNCTION TO LOG STANDARDIZE VARIABLES FOR PA SOME ACTIVITY TYPES ARE = 0, add small number to enable log transformation
log_standardize <- function(variable) {
  # Filter out missing values
  non_missing_values <- variable[!is.na(variable)]
  non_missing_values[non_missing_values == 0] <- non_missing_values[non_missing_values == 0] + 0.00001
  # Log transform the non-missing values
  log_variable <-log(non_missing_values)
  # Standardize the log-transformed variable
  standardized_variable <- (log_variable - mean(log_variable)) / sd(log_variable)
  # Create a vector with missing values preserved
  result <- rep(NA, length(variable))
  result[!is.na(variable)] <- standardized_variable
  return(result)
}


activity_var <- c(
  # TEE, PAEE, PAL - DLW
  "dlw_tee_avg", "dlw_paee_avg", "dlw_pal_avg", 
  
  # Total active time (questionnaire and devices)
  "ap_methrs", "ag_soj3x_methrs", "a24_total_methrs",
  
  # Total active time by intensity 
  "a24_sed_hrs", "ag_soj3x_sed_hrs", "ap_sed_hrs",
  "a24_light_hrs", "ag_soj3x_light_hrs", "ap_light_hrs",
  "a24_mod_hrs", "ag_soj3x_mod_hrs",
  "a24_vig_hrs", "ag_soj3x_vig_hrs",
  "a24_mvpa_hrs", "ag_soj3x_mvpa_hrs", "ap_mvpa_hrs" )

averages_dat <- read_sas("./Chuck PA data/Data final/avgs_all.sas7bdat") %>%
  mutate(ageG = case_when(AGE < 55 ~ 0, 
                          AGE >=55 & AGE< 60 ~ 1,
                          AGE >=60 & AGE< 65 ~ 2,
                          AGE >=65 & AGE< 70 ~ 3,
                          AGE >=70 & AGE< 75 ~ 4,
                          AGE >=75 ~ 5 )) %>%
  mutate(raceG = case_when(race > 1 ~ 1, TRUE ~ 0)) %>%
  # for two people their RMR is the same as their total energy expenditure, drop for now
  mutate(dlw_tee = if_else(dlw_paee_mif < 0, NA_real_, dlw_tee)) %>%
  mutate(dlw_tee2nd = if_else(dlw_paee_mif2nd < 0, NA_real_, dlw_tee2nd)) %>%
  mutate(dlw_pal_mif = if_else(dlw_paee_mif < 0, NA_real_, dlw_pal_mif)) %>%
  mutate(dlw_pal_mif2nd = if_else(dlw_paee_mif2nd < 0, NA_real_, dlw_pal_mif2nd)) %>%
  mutate(dlw_paee_mif = if_else(dlw_paee_mif < 0, NA_real_, dlw_paee_mif)) %>%
  mutate(dlw_paee_mif2nd = if_else(dlw_paee_mif2nd < 0, NA_real_, dlw_paee_mif2nd)) %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  rowwise() %>% 
  dplyr::mutate(weight_avg = mean(c(cv1_avg_wgtf, cv2_avg_wgtf , cv3_avg_wgtf), na.rm=TRUE)) %>%
  dplyr::mutate(dlw_wt_avg = mean(c(cv_wt_dlw1, cv_wt_dlw2), na.rm=TRUE)) %>%
  dplyr::mutate(dlw_tee_avg = mean (c(dlw_tee, dlw_tee2nd), na.rm=TRUE)) %>%
  dplyr::mutate(dlw_pal_avg = mean(c(dlw_pal_mif, dlw_pal_mif2nd), na.rm=TRUE)) %>%
  dplyr::mutate(dlw_ree_est_avg = mean(c(cv_ree_mif, cv_ree_mif2nd), na.rm=TRUE)) %>%
  dplyr::mutate(dlw_paee_avg = mean(c(dlw_paee_mif, dlw_paee_mif2nd), na.rm=TRUE)) %>%
  dplyr::mutate(dlw_pai_avg = mean(c(dlw_pai_mif, dlw_pai_mif2nd), na.rm=TRUE)) %>%
  ungroup 

# log transform and standardize physical activity variables
for (i in activity_var) {
  x <-log_standardize(averages_dat[[i]])
  new_col <- paste0(i, "_ln_sd")
  averages_dat[[new_col]] <- x
}


# Pull out dlw body composition data (not in Chuck's data)
dlw_comp <-fread("./DLW Energy Expenditure/Primary/dlw.jul22.d070122.csv") %>%
  dplyr::select(c(iid, ffmc, ffmc2nd, pct_fat, pct_fat2nd)) %>%
  na_if("M") %>%
  mutate_at(c("ffmc", "ffmc2nd", "pct_fat", "pct_fat2nd"), as.numeric) %>%
  rowwise() %>% 
  dplyr::mutate(ffmc_avg = mean(c(ffmc, ffmc2nd), na.rm=TRUE)) %>%
  dplyr::mutate(pct_fat_avg = mean(c(pct_fat, pct_fat2nd), na.rm=TRUE)) %>%
  ungroup() 


### Other clinical assessments - this pulls out standing height and average bmi over the 3 clinical assessments 
clin_vis <- fread("C:/Users/wattsel/OneDrive - National Institutes of Health/iDATA/Data/Participant/Primary/participant.may17.d051117.csv") %>% 
  dplyr::select(c(iid, stand_hgt, avg_bmi)) %>%
  mutate(avg_bmi = ifelse(avg_bmi =="N", NA, avg_bmi)) %>% mutate_at(vars(2:3),as.numeric)

# Merge dlw and averages
PA_COMBINED<-merge(averages_dat, dlw_comp, by="iid", all.x=TRUE )  %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  merge(clin_vis, by="iid", all.x=TRUE)

# impute pct fat and FFM - revisit whether we need to do the pooling etc. due to the imputation but kept simple for now. 
impute_bodycomp <- PA_COMBINED %>% dplyr::select(ffmc_avg, pct_fat_avg, AGE, race, sex , avg_bmi, stand_hgt ) %>% mice(seed=123) %>% complete() %>% select(c(pct_fat_avg, ffmc_avg)) %>%
  rename(pct_fat_avg_imp = pct_fat_avg,
         ffmc_avg_imp = ffmc_avg )

# # fit a linear model to predict pct_fat - use MICE FOR NOW
# fat_mod <- lm(pct_fat_avg ~ AGE + race + sex + avg_bmi, data = PA_COMBINED)
# summary(fat_mod)
# PA_COMBINED$pct_fat_avg[is.na(PA_COMBINED$pct_fat_avg)] <- predict(fat_mod, newdata = PA_COMBINED)[is.na(PA_COMBINED$pct_fat_avg)]
# test<-PA_COMBINED %>% select(pct_fat_avg, avg_bmi, AGE ,race , sex)

# Merge into main dataset 
PA_COMBINED <- cbind(PA_COMBINED , impute_bodycomp )
  
##########################  
## METABOLOMIC data     ##
##########################
##########################
## 1) SERUM 
##########################
serum_metabolites<-fread("./IDATA_WAS_serum_urine_update/IDATA_WAS_NCIA030221_serum_rslts.csv") 

# EXTRACT COTININE AS THIS IS A MARKER OF SMOKING THAT IS MISSING FOR > 50% OF PARTICIPANTS BUT WANTED FOR ADJUSTMENTS LATER 
participant_dat<-serum_metabolites %>% dplyr::select(c(1:8, "COTININE")) %>% colnames() %>% unlist

# Count number of metabolites missing for > 50% participants
met_missing<- serum_metabolites %>% filter(QC_FLAG!=1 & RSLT_TYPE=="Raw Ion")  %>% dplyr::select(c(9:1472))
non_miss_met<- (colMeans(is.na(met_missing))*100) %>% as.data.frame() %>% dplyr::rename(v1=".") %>%
  # Extract list of metabolites with > 50% missing - at the moment this is for across the 2 visits
  filter(v1 < 50) %>%
  # sum(counts_tab$v1 > 50)
  # # 259 missing for > 50% of participants
  # sum(counts_tab$v1 > 90)
  # # 148 missing for > 90% participants
  # EXTRACT LIST METABOLITES WITH LOW AMOUNTS OF MISSING DATA
  rownames() %>% unlist()
colnames_keep <- append(participant_dat, non_miss_met)

# remove metabolite columns based on if not in list non-miss-met list
serum_metabolites <- serum_metabolites %>% filter(RSLT_TYPE == "Log" & QC_FLAG!=1) 
serum_metabolites <- serum_metabolites[, colnames_keep, with = FALSE]

serum_v1 <- serum_metabolites %>% filter(VISITTYPE == "Visit 1") %>% dplyr::select(c(SUBJECT_ID, COTININE:X_26119)) %>%
  # Average over duplicate samples
  group_by(SUBJECT_ID) %>% 
  dplyr::summarise(across(where(is.numeric), mean))

serum_v2 <- serum_metabolites %>% filter(VISITTYPE == "Visit 2") %>% dplyr::select(c(SUBJECT_ID, COTININE:X_26119)) %>%
  # Average over duplicate samples
  group_by(SUBJECT_ID) %>% 
  dplyr::summarise(across(where(is.numeric), mean))

# metabolites average over both visits
serum_metabolites_avg <-serum_metabolites %>% dplyr::select(c(SUBJECT_ID,COTININE:X_26119)) %>%
  group_by(SUBJECT_ID) %>% 
  dplyr::summarise(across(where(is.numeric), mean)) 

### MERGE serum with PA data 
combined_dta_met_serum<- merge(PA_COMBINED, serum_metabolites_avg, by.x="iid", by.y="SUBJECT_ID" ) %>% 
  mutate(smok_curr = case_when(COTININE == min(COTININE) ~ 0, TRUE ~ 1)) 
combined_dta_met_serum$smok_curr<-as.numeric(combined_dta_met_serum$smok_curr)

fwrite(combined_dta_met_serum,"PA_avg_met_avg_serum.csv")

#############################
## 2) URINE (FMV AND 24-HR) #
#############################

urine<-fread("./IDATA_WAS_serum_urine_update/IDATA_WAS_NCIA030321_urine_rslts.csv")
urine_sample<-fread("./IDATA_WAS_serum_urine_update/IDATA_WAS_NCIA030321_urine_samp_meta.csv")  

twenty_four_hr_urine_ids <- urine_sample %>% filter(MATERIAL_CAT=="24 hour urine") %>% dplyr::select(CLIENT_SAMPLE_ID) %>% unlist
fmv_urine_ids <- urine_sample %>% filter(MATERIAL_CAT=="FMV urine") %>% dplyr::select(CLIENT_SAMPLE_ID) %>% unlist
participant_urine<-urine %>% dplyr::select(c(1:8, "COTININE")) %>% colnames() %>% unlist

#########################
### A)  24-hr samples ###
#########################
twenty_four_hr_urine <- urine %>%  filter(RSLT_TYPE == "Osmo" & QC_FLAG!=1) %>% subset(CLIENT_SAMPLE_ID %in% twenty_four_hr_urine_ids)
# Missing metabolites
missings_24_hr <- urine %>% filter(QC_FLAG!=1 & RSLT_TYPE=="Raw Ion") %>%  subset(CLIENT_SAMPLE_ID %in% twenty_four_hr_urine_ids )  %>% dplyr::select(c(9:1581))
non_miss_met_24hr<- (colMeans(is.na(missings_24_hr))*100) %>% as.data.frame() %>% rename(v1=".") %>%
  # Extract list of metabolites with > 50% missing
  filter(v1 < 50) %>%
  rownames() %>% unlist()
colnames_keep <- append(participant_urine, non_miss_met_24hr)

# remove metabolite columns based on if not in list non-miss-met list
twenty_four_hr_urine <- twenty_four_hr_urine[, colnames_keep, with = FALSE]

urine_v1 <- twenty_four_hr_urine %>% filter(VISITTYPE == "Visit 1") %>% dplyr::select(c(SUBJECT_ID, COTININE:`_999926193`)) %>%
  # Average over duplicate samples
  group_by(SUBJECT_ID) %>% 
  summarise(across(where(is.numeric), mean))

urine_v2 <- twenty_four_hr_urine %>% filter(VISITTYPE == "Visit 2") %>% dplyr::select(c(SUBJECT_ID, COTININE:`_999926193`)) %>%
  # Average over duplicate samples
  group_by(SUBJECT_ID) %>% 
  summarise(across(where(is.numeric), mean))

# metabolites average over both visits
twenty_four_urine_avg <-twenty_four_hr_urine %>% dplyr::select(c(SUBJECT_ID,COTININE:`_999926193`)) %>%
  group_by(SUBJECT_ID) %>% 
  summarise(across(where(is.numeric), mean))

### MERGE 24 hr urine with dlw and baseline data 
combined_dta_met_urine<- merge(PA_COMBINED, twenty_four_urine_avg, by.x="iid", by.y="SUBJECT_ID" ) %>% 
  mutate(smok_curr = case_when(COTININE == min(COTININE) ~ 0, TRUE ~ 1)) 
combined_dta_met_urine$smok_curr<-as.numeric(combined_dta_met_urine$smok_curr)

fwrite(combined_dta_met_urine,"PA_avg_met_avg_urine.csv")


########################
### B)   FMV        ####
########################
fmv_urine_ids <- urine_sample %>% filter(MATERIAL_CAT=="FMV urine") %>% dplyr::select(CLIENT_SAMPLE_ID) %>% unlist

# separate out 24hr urine samples
fmv_urine <- urine %>%  filter(RSLT_TYPE == "Osmo" & QC_FLAG!=1) %>% subset(CLIENT_SAMPLE_ID %in% fmv_urine_ids)
# Missing metabolites
missings_fmv <- urine %>% filter(QC_FLAG!=1 & RSLT_TYPE=="Raw Ion") %>%  subset(CLIENT_SAMPLE_ID %in% fmv_urine_ids )  %>% dplyr::select(c(9:1581))
non_miss_met_fmvr<- (colMeans(is.na(missings_fmv))*100) %>% as.data.frame() %>% rename(v1=".") %>%
  # Extract list of metabolites with > 50% missing
  filter(v1 < 50) %>%
  rownames() %>% unlist()
colnames_keep <- append(participant_urine, non_miss_met_fmvr)

# remove metabolite columns based on if not in list non-miss-met list
fmv_urine <- fmv_urine[, colnames_keep, with = FALSE]

urine_v1 <- fmv_urine %>% filter(VISITTYPE == "Visit 1") %>% dplyr::select(c(SUBJECT_ID, COTININE:`_999926193`)) %>%
  # Average over duplicate samples
  group_by(SUBJECT_ID) %>% 
  summarise(across(where(is.numeric), mean))

urine_v2 <- fmv_urine %>% filter(VISITTYPE == "Visit 2") %>% dplyr::select(c(SUBJECT_ID, COTININE:`_999926193`)) %>%
  # Average over duplicate samples
  group_by(SUBJECT_ID) %>% 
  summarise(across(where(is.numeric), mean))

# metabolites average over both visits
fmv_avg <-fmv_urine %>% dplyr::select(c(SUBJECT_ID,COTININE:`_999926193`)) %>%
  group_by(SUBJECT_ID) %>% 
  summarise(across(where(is.numeric), mean))

### MERGE FMV urine with dlw and baseline data 
combined_dta_met_fmv<- merge(PA_COMBINED, fmv_avg, by.x="iid", by.y="SUBJECT_ID" ) %>%
  mutate(smok_curr = case_when(COTININE == min(COTININE) ~ 0, TRUE ~ 1)) 

combined_dta_met_fmv$smok_curr<-as.numeric(combined_dta_met_fmv$smok_curr)

fwrite(combined_dta_met_fmv,"PA_avg_met_avg_fmv.csv")


