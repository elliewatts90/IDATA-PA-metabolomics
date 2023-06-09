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

# # Function to calculate Q-het - this allows up to 4 comparison groups but can add more if needed
# r's are the correlation coefficients, and n's are the participant numbers for each correlation - r1 should correspond to n1 etc.
q_het <- function(r1, r2, r3=NA, r4=NA, n1, n2, n3=NA, n4=NA ) {
  # group1
  fish_z_1 <- 0.5*(log(1+r1)-log(1-r1))
  varience_1 <- (1/sqrt(n1-3))^2
  n_pool_1 <- fish_z_1/varience_1
  denom_pool_1 <- 1/varience_1
  # grp 2
  fish_z_2 <- 0.5*(log(1+r2)-log(1-r2))
  varience_2 <- (1/sqrt(n2-3))^2
  n_pool_2 <- fish_z_2/varience_2
  denom_pool_2 <- 1/varience_2
  # grp 3
  fish_z_3 <- 0.5*(log(1+r3)-log(1-r3))
  varience_3 <- (1/sqrt(n3-3))^2
  n_pool_3 <- fish_z_3/varience_3
  denom_pool_3 <- 1/varience_3
  # grp 4
  fish_z_4 <- 0.5*(log(1+r4)-log(1-r4))
  varience_4 <- (1/sqrt(n4-3))^2
  n_pool_4 <- fish_z_4/varience_4
  denom_pool_4 <- 1/varience_4
  
  # summed 
  numerator <-  sum(n_pool_1, n_pool_2, n_pool_3, n_pool_4, na.rm=TRUE)
  denominator <- sum(denom_pool_1, denom_pool_2, denom_pool_3, denom_pool_4, na.rm=TRUE)
  sum_n_d <- numerator / denominator
  
  # fish z - pooled
  fish_pooled_1 <- ((fish_z_1 - sum_n_d)^2)/varience_1
  fish_pooled_2 <- ((fish_z_2 - sum_n_d)^2)/varience_2
  fish_pooled_3 <- ((fish_z_3 - sum_n_d)^2)/varience_3
  fish_pooled_4 <- ((fish_z_4 - sum_n_d)^2)/varience_4
  
  # count number of comparison groups
  list<-c(fish_pooled_1, fish_pooled_2, fish_pooled_3, fish_pooled_4)
  ncompare <- sum(!is.na(list))
  
  # Cochrans q and p-value
  cochrans_q <- sum(fish_pooled_1, fish_pooled_2, fish_pooled_3, fish_pooled_4, na.rm=TRUE)
  p.value<-pchisq(cochrans_q, df=(ncompare - 1), lower.tail = FALSE)
  return(p.value)
}


######## INSTALL PACKAGES
package_list<-c("dplyr", "data.table",  "tibble", "psych", "stats")
install_load(package_list)


##########################   
# i) Heterogeneity by sex 
###########################  
dta_list <- c("serum","urine",'fmv')
PA_list <- c("act", "pal" )
for (d in PA_list) {
  for (x in dta_list) {

sex_data<-fread(paste("C:/Users/wattsel/OneDrive - National Institutes of Health/iDATA/Results/correlations/",d,"_corr_sex",x,".csv", sep= "")) 
sex_data$phet <- NA
    for (i in 1:nrow(sex_data)) {
      temp<- sex_data[i,]
      p<-q_het(r1 = temp$corr_female, r2=temp$corr_male, n1=temp$nfem, n2=temp$nmale)
      sex_data$phet[i] <-p
    }
  fwrite(data, paste("C:/Users/wattsel/OneDrive - National Institutes of Health/iDATA/Results/correlations/",d,"_corr_sex",x,"_qhet.csv", sep= "")) 
  }
}


######################################  
# ii) Heterogeneity by sampling type
######################################

PA_list <- c("act", "pal")
for (d in PA_list) {

fmv<-fread(paste("C:/Users/wattsel/OneDrive - National Institutes of Health/iDATA/Results/correlations/corr_",d,"_fmv.csv", sep = "")) %>% dplyr::select(c(metabolite, estimate.3 )) %>%
  rename(corr_fmv = estimate.3) 
fmv$npart<-560

urine_24<-fread(paste("C:/Users/wattsel/OneDrive - National Institutes of Health/iDATA/Results/correlations/corr_",d,"_urine.csv", sep = "")) %>% dplyr::select(c(metabolite, estimate.3 ))%>%
  rename(corr_24urine = estimate.3)
urine_24$npart<-560

serum<-fread(paste("C:/Users/wattsel/OneDrive - National Institutes of Health/iDATA/Results/correlations/corr_",d,"_serum.csv", sep = "")) %>% dplyr::select(c(metabolite, estimate.3 ))%>%
  rename(corr_serum = estimate.3)
serum$npart<-560

# Merge by metabolite
samples_combined <- merge(fmv, urine_24, by="metabolite", all.x=TRUE, all.y=TRUE) %>%
  merge(serum, by="metabolite", all.x=TRUE, all.y=TRUE) 
samples_combined$missings <- rowSums(is.na(samples_combined)) 
samples_combined <-  samples_combined %>% mutate(sample_count = case_when(missings == 0 ~ 3,
                                                                          missings == 2 ~ 2, 
                                                                          missings == 4 ~ 1)) %>% mutate_at(vars(2:7), as.numeric) 

# Heterogeneity by sample type
samples_combined$phet<- NA
for (i in 1:nrow(samples_combined)) {
  temp<-samples_combined[i,] 
  p <- q_het(r1 = temp$corr_fmv , r2= temp$corr_24urine , r3= temp$corr_serum, n1= temp$npart.x, n2=temp$npart.y , n3=temp$npart)
  samples_combined$phet[i]<-p
}
samples_combined$phet[samples_combined$sample_count == 1 ] <- NA
fwrite(samples_combined,paste("C:/Users/wattsel/OneDrive - National Institutes of Health/iDATA/Results/correlations/corr_",d,"_all_biosamples.csv", sep="" ))
}

########################################################  
# iii) Heterogeneity by serum and 24-hr urine only
########################################################

#### JUST FOR SERUM AND 24 HR URINE
PA_list <- c("act", "pal")
for (d in PA_list) {
  
  urine_24<-fread(paste("C:/Users/wattsel/OneDrive - National Institutes of Health/iDATA/Results/correlations/corr_",d,"_urine.csv", sep = "")) %>% dplyr::select(c(metabolite, estimate.3 ))%>%
    rename(corr_24urine = estimate.3)
  urine_24$npart<-560
  
  serum<-fread(paste("C:/Users/wattsel/OneDrive - National Institutes of Health/iDATA/Results/correlations/corr_",d,"_serum.csv", sep = "")) %>% dplyr::select(c(metabolite, estimate.3 ))%>%
    rename(corr_serum = estimate.3)
  serum$npart<-560
  
bl_urine <- merge(urine_24,serum , by="metabolite") %>% mutate_at(vars(2:5), as.numeric) 
bl_urine$phet<- NA
for (i in 1:nrow(bl_urine)) {
  temp<-bl_urine[i,] 
  p <- q_het(r1=temp$corr_24urine , r2=temp$corr_serum , n1=temp$npart.x, n2=temp$npart.y)
  bl_urine$phet[i]<-p
}

bl_urine$phet[bl_urine$sample_count == 1 ] <- NA
fwrite(samples_combined,paste("C:/Users/wattsel/OneDrive - National Institutes of Health/iDATA/Results/correlations/corr_",d,"_24hr_serum.csv", sep="" ))
}

