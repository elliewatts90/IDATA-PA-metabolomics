#### CORRELATIONS OF INDIVIDUAL METABOLITES WITH RMR, TEE AND PHYS - last edited June 9th 2023
# Iterates over serum, 24hr ur and FMV sample datasets and each type of avtivity. 
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

# Modified write xlsx function
write.xlsx.custom <- function(x, file, sheetName="Sheet1",
                              col.names=TRUE, row.names=TRUE, append=FALSE, showNA=FALSE)
{
  if (!is.data.frame(x))
    x <- data.frame(x)    # just because the error message is too ugly
  
  iOffset <- jOffset <- 0
  if (col.names)
    iOffset <- 1
  if (row.names)
    jOffset <- 1
  
  if (append && file.exists(file)){
    wb <- loadWorkbook(file)
  } else {
    ext <- gsub(".*\\.(.*)$", "\\1", basename(file))
    wb  <- createWorkbook(type=ext)
  }  
  sheet <- createSheet(wb, sheetName)
  
  noRows <- nrow(x) + iOffset
  noCols <- ncol(x) + jOffset
  if (col.names){
    rows  <- createRow(sheet, 1)                  # create top row
    cells <- createCell(rows, colIndex=1:noCols)  # create cells
    mapply(setCellValue, cells[1,(1+jOffset):noCols], colnames(x))
  }
  if (row.names)             # add rownames to data x                   
    x <- cbind(rownames=rownames(x), x)
  
  if(nrow(x) > 0) {
    colIndex <- seq_len(ncol(x))
    rowIndex <- seq_len(nrow(x)) + iOffset
    
    .write_block(wb, sheet, x, rowIndex, colIndex, showNA)
  }
  saveWorkbook(wb, file)
  
  invisible()
}

.write_block <- function(wb, sheet, y, rowIndex=seq_len(nrow(y)),
                         colIndex=seq_len(ncol(y)), showNA=TRUE)
{
  rows  <- createRow(sheet, rowIndex)      # create rows
  cells <- createCell(rows, colIndex)      # create cells
  
  for (ic in seq_len(ncol(y)))
    mapply(setCellValue, cells[seq_len(nrow(cells)), colIndex[ic]], y[,ic], FALSE, showNA)
  
  # Date and POSIXct classes need to be formatted
  indDT <- which(sapply(y, function(x) inherits(x, "Date")))
  if (length(indDT) > 0) {
    dateFormat <- CellStyle(wb) + DataFormat(getOption("xlsx.date.format"))
    for (ic in indDT){
      lapply(cells[seq_len(nrow(cells)),colIndex[ic]], setCellStyle, dateFormat)
    }
  }
  
  indDT <- which(sapply(y, function(x) inherits(x, "POSIXct")))
  if (length(indDT) > 0) {
    datetimeFormat <- CellStyle(wb) + DataFormat(getOption("xlsx.datetime.format"))
    for (ic in indDT){
      lapply(cells[seq_len(nrow(cells)),colIndex[ic]], setCellStyle, datetimeFormat)
    }
  }
  
}


######## INSTALL PACKAGES
options(java.parameters = "-Xmx8000m")

package_list<-c("plyr","dplyr", "data.table", "ggplot2", "ppcor", "fastDummies", "janitor", "tibble", "psych", "stringr", "xlsx", "tidyverse")
install_load(package_list)

# fix plyr/dplyr issue
if("dplyr" %in% (.packages())){
  detach("package:dplyr", unload=TRUE) 
  detach("package:plyr", unload=TRUE) 
} 
library(plyr)
library(dplyr)


##### CORRELATIONS FUNCTION
correlations_fun <- function(df, column_name){
  column_name = enquo(column_name)
  data<-df %>% dplyr::select(!!column_name) 
  colnames(data)<- "var" 
  # Minimally adjusted data
  corr_res<-sapply(1:ncol(correlations), function(x) { 
    pcor.test(data$var,correlations[,x],list(age$AGE,  sex$sex, smok$smok_curr, raceG$raceG), method="pearson") }) %>%
    as.data.frame() %>%
    slice(1:2) %>%
    t() %>% as.data.frame()
  corr_res$metabolite<-corr_names
  corr_res <- as.data.frame(lapply(corr_res, unlist)) %>% dplyr::select(c(metabolite, estimate, p.value)) %>%
    dplyr::rename(estimate.1 = estimate, p.value.1 = p.value)
  corr_res$p.fdr.1 <-  p.adjust(corr_res$p.value.1, method = "BH") 
  corr_res$model.1<-"minimally adj"
  corr_res<- corr_res %>% relocate(metabolite, model.1)
  # Adjusted for percentage bodyfat 
  corr_res_body<-sapply(1:ncol(correlations), function(x) { 
    pcor.test(data$var,correlations[,x],list(age$AGE,  sex$sex, smok$smok_curr, raceG$raceG, pbf$pct_fat_avg_imp), method="pearson") }) %>%
    as.data.frame() %>%
    slice(1:2) %>%
    t() %>% as.data.frame()
  corr_res_body <- as.data.frame(lapply(corr_res_body, unlist)) %>% dplyr::select(c(estimate, p.value))%>%
    dplyr::rename(estimate.2 = estimate, p.value.2 = p.value)
  corr_res_body$model.2<-"adj body fat"
  corr_res_body$p.fdr.2 <-  p.adjust(corr_res_body$p.value.2, method = "BH") 
  corr_res_body<- corr_res_body %>% relocate(model.2)

  complete_data<-cbind(corr_res, corr_res_body)
  return(complete_data)
}

# cbind full function
cbind.fill <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

#  Loop through serum, 24-hr urine and 1st morning void urine to estimate partial correlations
# N.B. Need to remove files where previously created due to append when I write out the results 
file.remove("C:/Users/wattsel/OneDrive - National Institutes of Health/iDATA/Results/correlations/corr_all.csv")
file.remove("C:/Users/wattsel/OneDrive - National Institutes of Health/iDATA/Results/correlations/corr_all.xlsx")
file.remove("C:/Users/wattsel/OneDrive - National Institutes of Health/iDATA/Results/correlations/corr_all_fdr.xlsx")
file.remove("C:/Users/wattsel/OneDrive - National Institutes of Health/iDATA/Results/correlations/met_freq_counts_by_measurement.xlsx")
file.remove("C:/Users/wattsel/OneDrive - National Institutes of Health/iDATA/Results/correlations/met_freq_counts_by_PA_type.xlsx")
file.remove("C:/Users/wattsel/OneDrive - National Institutes of Health/iDATA/Results/correlations/met_name_by_sample_type.xlsx")
file.remove("C:/Users/wattsel/OneDrive - National Institutes of Health/iDATA/Results/correlations/corr_all_ln_sd.csv")


dta_list <- c("serum","urine",'fmv')
  for (x in dta_list) {
    dat<-fread(paste("PA_avg_met_avg_",x,".csv", sep="")) 
    

# LOOP THROUGH FOR EACH ACTIVITY 

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

temp <- data.frame()

 for (i in activity_var) {

   # drop participants with missing data
   dat <- dat %>% filter(!is.na(!!rlang::sym(i)))
      
   # drop participants with missing body comp data
   dat <- dat %>% filter(!is.na(pct_fat_avg_imp)) 
  
   # dat<-fread("PA_avg_met_avg_urine.csv")
   dat$raceG<-as.numeric(dat$raceG)
   dat$sex<-as.numeric(dat$sex)
   dat$smok_curr<-as.numeric(dat$smok_curr)
   
   #### 1 ) CORRELATIONS ADJUSTED FOR AGE SEX RACE ONLY
   correlations<- dat %>%
     dplyr::select((matches("[A-Z]", ignore.case = FALSE)| starts_with("_") )) %>% 
     dplyr::select(-c(AGE, ageG, raceG, COTININE ))
   
   corr_names<-colnames(correlations)
   correlations<-as.matrix(correlations)
   
   raceG<- dat %>% dplyr::select(raceG)
   sex<- dat %>% dplyr::select(sex)
   age<- dat %>% dplyr::select(AGE)
   smok<-dat %>% dplyr::select(smok_curr)
   # leanmass <- dat %>%  dplyr::select(ffmc_avg)
   pbf <- dat %>%  dplyr::select(pct_fat_avg_imp)
  
   # CORRELATIONS 
  res <-correlations_fun(df= dat, column_name= paste(i))
  res$act_var<-paste(i)
  res$sample<-paste(x)
  temp<-rbind(temp, res)
 
}

# save out full dataset

fwrite(temp, "C:/Users/wattsel/OneDrive - National Institutes of Health/iDATA/Results/correlations/corr_all_ln_sd.csv", append = TRUE)

  }

### GENERATE COUNTS FOR TABLES ####

# Count for each metabolite significant hit
# SPLIT BY SERUM, FMV, 24-hr 
dta_list <- c("serum","urine",'fmv')
counts_res <- data.frame()
  for (x in dta_list) {
    sample_type <- fread("C:/Users/wattsel/OneDrive - National Institutes of Health/iDATA/Results/correlations/corr_all.csv") %>%
    filter(sample==paste(x)) %>% 
    filter(p.fdr.2 <0.05) %>% select(c(metabolite, estimate.2, p.fdr.2)) 
    count <- table(sample_type$metabolite) %>% as.data.frame() %>% arrange(desc(Freq))
    count$sample_type <- paste0(x)
    counts_res <- cbind.fill(counts_res, count )
  }
counts_res <- counts_res %>% as.data.frame()
fwrite(counts_res, "C:/Users/wattsel/OneDrive - National Institutes of Health/iDATA/Results/correlations/met_freq_counts.csv")

# BY MEASUREMENT TYPE 
dta_list <- c("serum","urine",'fmv')
measure_list <- c("dlw", "a24", "ag_soj3x" , "ap_")
  for (x in dta_list) {
    counts_res <- data.frame()
    measurement_res <- fread("C:/Users/wattsel/OneDrive - National Institutes of Health/iDATA/Results/correlations/corr_all.csv") %>%
      filter(sample==paste(x)) %>%
      filter(p.fdr.2 <0.05) %>% select(c(metabolite, estimate.2, p.fdr.2, act_var)) 
        for (d in measure_list) {
        measurement_res_exposure <- measurement_res %>% filter(str_detect(act_var, paste(d)) )
        count <- table(measurement_res_exposure$metabolite) %>% as.data.frame() %>% arrange(desc(Freq))
        count$sample_type <- paste0(x,"_",d)
        counts_res <- cbind.fill(counts_res, count)
        }
      counts_res <- counts_res %>% as.data.frame()
    write.xlsx.custom(counts_res, "C:/Users/wattsel/OneDrive - National Institutes of Health/iDATA/Results/correlations/met_freq_counts_by_measurement.xlsx", sheetName = paste0(x), append = TRUE)
      
  }


## BY TYPE OF ACTIVITY
type_list <- c("tee", "paee", "pal" , "sed",  "light", "mod", "vig", "mvpa", "methrs" )
dta_list <- c("serum","urine",'fmv')
for (x in dta_list) {
      counts_res <- data.frame()
      measurement_res <- fread("C:/Users/wattsel/OneDrive - National Institutes of Health/iDATA/Results/correlations/corr_all.csv") %>%
      filter(sample==paste(x)) %>%
      filter(p.fdr.2 <0.05) %>% select(c(metabolite, estimate.2, p.fdr.2, act_var))
    for (d in type_list) {
      measurement_res_exposure <- measurement_res %>% filter(str_detect(act_var, paste(d)) )
      if (nrow(measurement_res_exposure) > 0) {
        count <- table(measurement_res_exposure$metabolite) %>% as.data.frame() %>% arrange(desc(Freq))
        count$sample_type <- paste0(d)
        counts_res <- cbind.fill(counts_res, count )
      }
    }
  counts_res <- counts_res %>% as.data.frame()
  write.xlsx.custom(counts_res, "C:/Users/wattsel/OneDrive - National Institutes of Health/iDATA/Results/correlations/met_freq_counts_by_PA_type.xlsx", sheetName = paste0(x), append = TRUE)
}


# list most commonly occurring metabolites by significance 
dta_list <- c("serum","urine",'fmv')
  for (x in dta_list) {
      measurement_res <- fread("C:/Users/wattsel/OneDrive - National Institutes of Health/iDATA/Results/correlations/corr_all.csv") %>%
        filter(sample==paste(x)) %>%filter(p.fdr.2 <0.05) %>% select(c(metabolite, estimate.2, p.fdr.2, act_var)) %>%
        group_by(metabolite) %>% mutate(Activity = sort(act_var)) %>% 
       mutate(rid = row_number()) %>%
      pivot_wider(id_cols = metabolite, names_from = rid, values_from = c(act_var, estimate.2), names_sep = '') 
    nm1 <- names(sample_type)[colSums(is.na(sample_type)) >0]
    sample_type <- sample_type %>% arrange_at(vars(nm1), funs(is.na(.))) 
    write.xlsx.custom(as.data.frame(sample_type), "C:/Users/wattsel/OneDrive - National Institutes of Health/iDATA/Results/correlations/met_name_by_sample_type.xlsx", 
                      sheetName = paste0(x), col.names = TRUE, row.names=FALSE, append = TRUE)
    }
  

