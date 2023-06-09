rm(list=ls())

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


corr_plots <- function(dat, x_name, y_name, sample, savefolder=NULL){
  # overall correlation 
  r_sq<-round(cor(z$estimate.2.x, z$estimate.2.y, use = "complete.obs", method = "pearson"), digits = 2)
  # Create figure
  emf(file = paste("C:/Users/wattsel/OneDrive - National Institutes of Health/iDATA/Results/correlations/plots/",savefolder,"/Corr_",x_name,"_",y_name,"_",sample,".emf", sep=""), 
      width = 9, height = 9)
  print(ggplot(z, aes(x=estimate.2.x, y=estimate.2.y)) +
          geom_point(color="#9999CC", size=5, alpha=0.6)+
          annotate("text", x=0.2, y=-0.2, label= paste("Overall correlation=", r_sq, sep=""),size = 7) + 
          geom_abline(slope=1, intercept=0, color="grey", size=1) +
          xlab(paste(x_name," correlation coefficient", sep="")) + 
          scale_x_continuous(limits=c(-0.3,0.3)) +
          scale_y_continuous(limits=c(-0.3,0.3), name = paste(y_name," correlation coefficient", sep="")) +
          theme_minimal() +
          theme(text = element_text(size = 25, family="TT Arial"), legend.text=element_text(size=25), axis.title=element_text(face="bold")) 
        
  )       
  dev.off()
  
}


package_list<-c("dplyr", "data.table", "ggplot2", "stringr",  "tidyverse", "devEMF")
install_load(package_list)

dta_list <- c("serum","urine",'fmv')

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

full_dta <- fread("C:/Users/wattsel/OneDrive - National Institutes of Health/iDATA/Results/correlations/corr_all_ln_sd.csv")

serum <- full_dta %>% filter(str_detect(sample, paste("serum"))) %>% dplyr::select(c(metabolite, estimate.2, p.value.2, act_var)) 
urine <- full_dta %>% filter(str_detect(sample, paste("urine"))) %>% dplyr::select(c(metabolite, estimate.2, p.value.2, act_var)) 
fmv <- full_dta %>% filter(str_detect(sample, paste("fmv"))) %>% dplyr::select(c(metabolite, estimate.2, p.value.2, act_var)) 

# WITHIN SAMPLE TYPES
for (x in dta_list){
  for (i in 1:(length(activity_var) - 1)){
    for (d in (i + 1):length(activity_var)) {
      
    temp1<-full_dta %>% filter(str_detect(sample, paste0(x))) %>% filter(str_detect(act_var, paste(activity_var[[i]]))) %>% dplyr::select(metabolite, estimate.2)
    temp2<-full_dta %>% filter(str_detect(sample, paste0(x))) %>% filter(str_detect(act_var, paste(activity_var[[d]]))) %>% dplyr::select(metabolite, estimate.2)
      z<-inner_join(temp1, temp2, by='metabolite') 
      x_name <- paste0(activity_var[[i]]) 
      y_name <- paste0(activity_var[[d]])
      sample_name <- paste0(x)
        corr_plots(z, x_name=x_name,y_name=y_name, sample=sample_name, savefolder="by_activity")
    }
  }
}

# # remove duplicate files
# setwd("C:/Users/wattsel/OneDrive - National Institutes of Health/iDATA/Results/correlations/plots/by_activity")
# to_be_deleted <- list.files(pattern ="ap_methrs_ln_sd_ap_methrs_ln_sd|ag_soj3x_methrs_ln_sd_ag_soj3x_methrs_ln_sd|
# a24_total_methrs_ln_sd_a24_total_methrs_ln_sd|a24_sed_hrs_ln_sd_a24_sed_hrs_ln_sd|
# ag_soj3x_sed_hrs_ln_sd_ag_soj3x_sed_hrs_ln_sd|ap_sed_hrs_ln_sd_ap_sed_hrs_ln_sd|
# a24_light_hrs_ln_sd_a24_light_hrs_ln_sd|ag_soj3x_light_hrs_ln_sd_ag_soj3x_light_hrs_ln_sd|
# ap_light_hrs_ln_sd_ap_light_hrs_ln_sd|a24_mod_hrs_ln_sd_a24_mod_hrs_ln_sd|ag_soj3x_mod_hrs_ln_sd_ag_soj3x_mod_hrs_ln_sd|
# a24_vig_hrs_ln_sd_a24_vig_hrs_ln_sd|
# ag_soj3x_vig_hrs_ln_sd_ag_soj3x_vig_hrs_ln_sd|a24_mvpa_hrs_ln_sd_a24_mvpa_hrs_ln_sd|
# ag_soj3x_mvpa_hrs_ln_sd_ag_soj3x_mvpa_hrs_ln_sd|ap_mvpa_hrs_ln_sd_ap_mvpa_hrs_ln_sd")
# file.remove(to_be_deleted)


# BETWEEN SAMPLE TYPES
dta_list_plus <- c("urine",'fmv', "serum")
for (x in activity_var){
  for (i in 1:(length(dta_list) - 1)){
    for (d in (i + 1):length(dta_list)) {
      temp1<-full_dta %>% filter(str_detect(act_var, paste0(x))) %>% filter(str_detect(sample, dta_list[[i]])) %>% dplyr::select(metabolite, estimate.2)
      temp2<-full_dta %>% filter(str_detect(act_var, paste0(x))) %>% filter(str_detect(sample, dta_list[[d]])) %>% dplyr::select(metabolite, estimate.2)
      z<-inner_join(temp1, temp2, by='metabolite') 
      x_name <- paste0(dta_list[[i]]) 
      y_name <- paste0(dta_list[[d]])
      sample_name <- paste0(x)
      corr_plots(z, x_name=x_name,y_name=y_name, sample=sample_name, savefolder="by_sample")
    }
  }
}

# # remove duplicate files
# setwd("C:/Users/wattsel/OneDrive - National Institutes of Health/iDATA/Results/correlations/plots/by_sample")
# to_be_deleted <- list.files(pattern = "urine_urine|fmv_fmv|serum_serum")
# file.remove(to_be_deleted)
