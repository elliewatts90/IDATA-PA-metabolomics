#########################
### R MASTER SCRIPT #####
#########################

## Created ELW June 9th 2023  FOR THE IDATA STUDY ANALYSIS ##
## THIS SCRIPT PULLS IN THE IDATA PREP AND OTHER SUBSEQUENT ANALYSES

setwd("C:/Users/wattsel/OneDrive - National Institutes of Health/iDATA/Rscripts")

# 1) DATA PREPARING - N.B. CHUCK RAN THE DATA PREP FOR THE PA ASSESSMENTS USING SAS 
# SCRIPT AVAILABLE FROM "C:\Users\wattsel\OneDrive - National Institutes of Health\iDATA\Data\Chuck PA data\chuck edited script\chuck script EW 04192023.R"
# CODE MERGES METABOLITE AND PHENOTYPE DATA AND DOES SOME BASIC VARIABLE CODING
# THIS ASSUMES CHUCK'S CODE HAS ALREADY BEEN RUN
source("chuck_dta_merge.R")

# 2) PEARSON PARTIAL CORRELATIONS
source("partial_correlations_log_transformed_06092023.R")

# 3) TESTING FOR HETEROGENEITY BY PEARSON PARTIAL CORRELATIONS
# I can't find the results files that this script pulls in are old so this won't run
# You will probably need you to rejiggle a bit to suit your files, but the function may be useful to you 
# The formula was provided by Steve and then I turned it into an R function to read in correlation coefficients and test for heterogeneity
source("qhet_corr.R")

# 4) PLOTS FOR CORRELATIONS
# Reads in output from the partial correlations csv files and outputs correlation plots for different types of activities 
source("correlation_plots.R")


# 5) LASSO 
# I created this function to run LASSO models both with and without adjustment covariates, however Grace and her student Sang Kyu have developed an alternative script, 
# Might be worth getting in touch with Grace if we're running things in tandem 
source("lasso_06092023.R")