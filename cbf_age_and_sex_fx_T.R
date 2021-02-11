######################################
### Initial construction 1/14/2021 ###
######################################

####################
##### Summary ######
####################

#input: asc files, and pnc demographics, cnb, clinical files
#output: 10242 length vector csvs with T and/or p values for vertex-wide regression
#uses: goes vertex by vertex and does regression (coupling by age, sex, cognition, etc), pulls out T and p from these values and sticks it in a vector. The vector can then be used for visualization in matlab
#dependencies: R (3.6.3 is my current default in pmacs)

####################
###  Libraries   ###
####################

library(mgcv)
library(dplyr)
library(ggplot2)

#####################################################################################
####              Makes the 831x10242 matrices, both left and right              ####
#####################################################################################
#set home directory, switch this depending on whether running from PMACS or from home directory
#homedir <- "/Users/eballer/BBL/imco/pmacs/PMACS_remote/"
homedir <- "/project/imco/"

# read in demos
cbf_subjDemos <- read.csv(paste0(homedir, "/baller/subjectLists/n831_alff_cbf_finalSample_imageOrder.csv"))

#some verification preprocessing
cbf_subjDemos$sex <- as.factor(cbf_subjDemos$sex)
cbf_subjDemos$race <- as.factor(cbf_subjDemos$race)
cbf_subjDemos$race2 <- as.factor(cbf_subjDemos$race2)

#add osex category for use in gam later
cbf_subjDemos$osex <- ordered(cbf_subjDemos$sex)

#add psych bifactor scores
psych <- read.csv(paste0(homedir, "/pnc/clinical/n1601_goassess_itemwise_bifactor_scores_20161219.csv"), header = TRUE)

#remove 4factorv2 from title
names(psych) <-gsub("_4factorv2", "", names(psych))

#merge
cbf_subjDemos <- merge(cbf_subjDemos, psych, by = "bblid")

#drop the scanid.y
cbf_subjDemos <- subset(cbf_subjDemos, select = -scanid.y)

#rename scanid.x to scanid
names(cbf_subjDemos) <- gsub("scanid.x", "scanid", names(cbf_subjDemos)) 

#make list of bblid/scanid
bblid_scanid <- paste0(cbf_subjDemos$bblid, "_", cbf_subjDemos$datexscanid)

#####################
##### Left Side #####
#####################

#initiate matrix for storage
numrows <- dim(cbf_subjDemos)[1]
lh_matrix <- matrix(nrow = numrows, ncol = 10242)

#go through each subject, grab 5th column in asc, transpose and stick in matrix
# output is 831 x 10242 matrix
for (subj in 1:numrows) {
  
  bblid <- cbf_subjDemos$bblid[subj]
  datexscanid <- cbf_subjDemos$datexscanid[subj]
  file_path <- paste0(homedir, "/surfaceMaps/cbf_from_chead/", bblid, "_", datexscanid, "_lh_fs5_surf.asc")
  cbf_data <- read.table(file_path, stringsAsFactors = FALSE)
  lh_matrix[subj,] <- t(cbf_data$V5)
  
}

#append with demographics
cbf_subjDemos_with_lh <- cbind(cbf_subjDemos, lh_matrix)

#write output
write.table(cbf_subjDemos_with_lh, file = paste0(homedir, "/baller/processed_data/cbf_matrices/cbf_subjDemos_with_lh_", numrows, "x10242.csv"), sep = ",")

#####################
#### Right Side #####
#####################

#initiate matrix for storage
numrows <- dim(cbf_subjDemos)[1]
rh_matrix <- matrix(nrow = numrows, ncol = 10242)

#go through each subject, grab 5th column in asc, transpose and stick in matrix
# output is 831 x 10242 matrix
for (subj in 1:numrows) {
  
  bblid <- cbf_subjDemos$bblid[subj]
  datexscanid <- cbf_subjDemos$datexscanid[subj]
  file_path <- paste0(homedir, "/surfaceMaps/cbf_from_chead/", bblid, "_", datexscanid, "_rh_fs5_surf.asc")
  cbf_data <- read.table(file_path, stringsAsFactors = FALSE)
  rh_matrix[subj,] <- t(cbf_data$V5)
  
}

#append with demographics
cbf_subjDemos_with_rh <- cbind(cbf_subjDemos, rh_matrix)

#write output
write.table(cbf_subjDemos_with_rh, file = paste0(homedir, "/baller/processed_data/cbf_matrices/cbf_subjDemos_with_rh_", numrows, "x10242.csv"), sep = ",")


####-----------------------------------------------------------------------------####
####---------------------------End of Part 1- Making matrices---_----------------####
####-----------------------------------------------------------------------------####

#####################################################################################
####                   Run regression, both left and right                       ####
#####################################################################################

#make easier to reference names
lh_cbf_asl <- cbf_subjDemos_with_lh #can also read directly from files if you'd like
rh_cbf_asl <- cbf_subjDemos_with_rh


#####################################################
#                       lm/gams                     #
#####################################################

#initialize vectors for models

hemis <- c("lh", "rh") #hemispheres
models <- c("age", "sex", "age_sex") #models of interest
coeffs <- c("p", "t") #p or t value
corrs <- c("uncor", "fdr") #correction

for (hemi in hemis){
  for (model in models) {
    for (coeff in coeffs) {
      for (corr in corrs) {
        vector_init_cmd <- paste0(hemi, "_gam_", model, "_", coeff, "_", corrs, " <- vector(length = 10242)")
        print(vector_init_cmd)
        eval(parse(text=as.name(vector_init_cmd)))
      }
    }
  }
}


#make linear models as well
lm_models <- c("age")
for (hemi in hemis) {
  for (model in lm_models) {
    for (coeff in coeffs) {
      for (corr in corrs) {
        vector_init_cmd <- paste0(hemi, "_lm_", model, "_", coeff, "_", corrs, "<- vector(length= 10242)")
        eval(parse(text=as.name(vector_init_cmd)))
      }
    }
  }
}


#######################
######## Left #########
#######################

#get # of items in df for calculation of column)
numcolumns <- dim(lh_cbf_asl)[2]
#run gams models and store info in respective vectors
for (i in 1:10242) {
  curcol = (numcolumns - 10242 + i) # will start you counting at the right part of the df
  
  #gams
  age_sex_model <- gam(lh_cbf_asl[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion +
                         osex + s(ageAtScan1, k = 4, fx = T), data=lh_cbf_asl)
  
  age_sex_intx_model <- gam(lh_cbf_asl[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion +
                              osex + s(ageAtScan1, k = 4, fx = T) + s(ageAtScan1, by = osex, k = 4, fx = T), data=lh_cbf_asl)
  
  #lm_for_directionality
  age_lm_model <- lm(lh_cbf_asl[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion + ageAtScan1, data=lh_cbf_asl)
  
  #put pvalue in it's appropriate lm
  lh_gam_age_p_uncor[i] <- summary(age_sex_model)$s.table[1,4] #smooth term for ageAtScan1
  lh_gam_sex_p_uncor[i] <- summary(age_sex_model)$p.table[4,4] #linear term
  lh_gam_age_sex_p_uncor[i] <- summary(age_sex_intx_model)$s.table[2,4] #smooth term for interaction, this was changed
  
  #lm to assess directionality
  lh_lm_age_p_uncor[i] <- summary(age_lm_model)$coeff[4,4]
  
  #pull tvalue into its appropriate lm
  lh_gam_age_t_uncor[i] <- summary(age_sex_model)$s.table[1,3] #smooth term for ageAtScan1
  lh_gam_sex_t_uncor[i] <- summary(age_sex_model)$p.table[4,3] #linear term
  lh_gam_age_sex_t_uncor[i] <- summary(age_sex_intx_model)$s.table[2,3] #smooth term for interaction
  
  #lm to assess directionality
  lh_lm_age_t_uncor[i] <- summary(age_lm_model)$coeff[4,3]
}

#####################
###### RIGHT ########
#####################
#get # of items in df for calculation of column)
numcolumns <- dim(rh_cbf_asl)[2]
#run gams models and store info in respective vectors
for (i in 1:10242) {
  curcol = (numcolumns - 10242 + i) # will start you counting at the right part of the df
  
  #gams
  age_sex_model <- gam(rh_cbf_asl[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion +
                         osex + s(ageAtScan1, k = 4, fx = T), data=rh_cbf_asl)
  
  age_sex_intx_model <- gam(rh_cbf_asl[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion +
                              osex + s(ageAtScan1, k = 4, fx = T) + s(ageAtScan1, by = osex, k = 4, fx = T), data=rh_cbf_asl)
  
  #lm_for_directionality
  age_lm_model <- lm(rh_cbf_asl[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion + ageAtScan1, data=rh_cbf_asl)
  
  #put pvalue in it's appropriate lm
  rh_gam_age_p_uncor[i] <- summary(age_sex_model)$s.table[1,4] #smooth term for ageAtScan1
  rh_gam_sex_p_uncor[i] <- summary(age_sex_model)$p.table[4,4] #linear term
  rh_gam_age_sex_p_uncor[i] <- summary(age_sex_intx_model)$s.table[2,4] #smooth term for interaction, this was changed
  
  #lm to assess directionality
  rh_lm_age_p_uncor[i] <- summary(age_lm_model)$coeff[4,4]
  
  #pull tvalue into its appropriate lm
  rh_gam_age_t_uncor[i] <- summary(age_sex_model)$s.table[1,3] #smooth term for ageAtScan1
  rh_gam_sex_t_uncor[i] <- summary(age_sex_model)$p.table[4,3] #linear term
  rh_gam_age_sex_t_uncor[i] <- summary(age_sex_intx_model)$s.table[2,3] #smooth term for interaction
  
  #lm to assess directionality
  rh_lm_age_t_uncor[i] <- summary(age_lm_model)$coeff[4,3]
}

#################################################################################
#################################################################################

#####################################################
#                     results                       #
#####################################################

#### FDR correction ####
for (hemi in hemis) {
  for (model in models) {
    hemi_model_p_unc <- paste0(hemi, "_gam_", model, "_p_uncor") 
    hemi_model_p_fdr <- paste0(hemi, "_gam_", model, "_p_fdr")    
    
    print(hemi_model_p_unc)
    
    #correct p values
    pfdr <- eval(substitute(p.adjust(i, method="fdr"), list(i = as.name(hemi_model_p_unc))))
    
    #figure out which values are < 0.05 and add to pfdr matrix
    pfdr <- as.data.frame(pfdr)
    pfdr$sig <- ifelse(pfdr<0.05, 1, 0)
    pfdr$sig_noNA <- ifelse(is.na(pfdr$sig), 0, pfdr$sig)
    names(pfdr) <- c("pfdr", "sig05", "sig05_noNA")
    hemi_model_p_fdr <- as.data.frame(pfdr[,1]) #sig05
    
    
    
    #multiply T values by fdr vector to get the list of Ts that are fdr corrected
    hemi_model_t_unc <- paste0(hemi, "_gam_", model, "_t_uncor")
    hemi_model_t_fdr <- paste0(hemi, "_gam_", model, "_t_fdr"
    )
    t_df <- eval(substitute(as.data.frame(i), list(i = as.name(hemi_model_t_unc))))
    names(t_df) <- c("tval")
    t_df$tfdr <- pfdr[,3] * t_df$tval
    hemi_model_t_fdr <- as.data.frame(t_df[,2])
    
    
    #######################
    #### write tables #####
    #######################
    
    ## uncorrected ##
    
    ### p
    filename <- paste0(homedir, "/baller/results/cbf_age_and_sex/", hemi_model_p_unc, ".csv")
    write_table_command <- paste0("write.table(x = ", hemi_model_p_unc, ", file = \"", filename,"\", row.names = FALSE, col.names = FALSE)")
    eval(parse(text=write_table_command))
    
    ### t
    
    filename <- paste0(homedir, "/baller/results/cbf_age_and_sex/", hemi_model_t_unc, ".csv")
    write_table_command <- paste0("write.table(x = ", hemi_model_t_unc, ", file = \"", filename,"\", row.names = FALSE, col.names = FALSE)")
    eval(parse(text=write_table_command))
    
    ## corrected ##
    
    ### p
    filename <- paste0(homedir, "/baller/results/cbf_age_and_sex/", hemi, "_gam_", model, "_p_fdr05.csv")
    write_table_command <- paste0("write.table(x = ", hemi_model_p_fdr, ", file = \"", filename,"\", row.names = FALSE, col.names = FALSE)")
    eval(parse(text=write_table_command))
    
    ### t
    filename <- paste0(homedir, "/baller/results/cbf_age_and_sex/", hemi, "_gam_", model, "_t_fdr05.csv")
    write_table_command <- paste0("write.table(x = ", hemi_model_t_fdr, ", file = \"", filename,"\", row.names = FALSE, col.names = FALSE)")
    eval(parse(text=write_table_command))
    
  }
}

######  linear model alone #####

#### FDR correction ####
for (hemi in hemis) {
  for (model in lm_models) {
    hemi_model_p_unc <- paste0(hemi, "_lm_", model, "_p_uncor") 
    hemi_model_p_fdr <- paste0(hemi, "_lm_", model, "_p_fdr")    
    
    print(hemi_model_p_unc)
    
    #correct p values
    pfdr <- eval(substitute(p.adjust(i, method="fdr"), list(i = as.name(hemi_model_p_unc))))
    
    #figure out which values are < 0.05 and add to pfdr matrix
    pfdr <- as.data.frame(pfdr)
    pfdr$sig <- ifelse(pfdr<0.05, 1, 0)
    pfdr$sig_noNA <- ifelse(is.na(pfdr$sig), 0, pfdr$sig)
    names(pfdr) <- c("pfdr", "sig05", "sig05_noNA")
    hemi_model_p_fdr <- as.data.frame(pfdr[,1]) #sig05
    
    
    
    #multiply T values by fdr vector to get the list of Ts that are fdr corrected
    hemi_model_t_unc <- paste0(hemi, "_lm_", model, "_t_uncor")
    hemi_model_t_fdr <- paste0(hemi, "_lm_", model, "_t_fdr")
    t_df <- eval(substitute(as.data.frame(i), list(i = as.name(hemi_model_t_unc))))
    names(t_df) <- c("tval")
    t_df$tfdr <- pfdr[,3] * t_df$tval
    hemi_model_t_fdr <- as.data.frame(t_df[,2])
    
    
    #######################
    #### write tables #####
    #######################
    
    ## uncorrected ##
    
    ### p
    filename <- paste0(homedir, "/baller/results/cbf_age_and_sex/", hemi_model_p_unc, ".csv")
    write_table_command <- paste0("write.table(x = ", hemi_model_p_unc, ", file = \"", filename,"\", row.names = FALSE, col.names = FALSE)")
    eval(parse(text=write_table_command))
    
    ### t
    
    filename <- paste0(homedir, "/baller/results/cbf_age_and_sex/", hemi_model_t_unc, ".csv")
    write_table_command <- paste0("write.table(x = ", hemi_model_t_unc, ", file = \"", filename,"\", row.names = FALSE, col.names = FALSE)")
    eval(parse(text=write_table_command))
    
    ## corrected ##
    
    ### p
    filename <- paste0(homedir, "/baller/results/cbf_age_and_sex/", hemi, "_lm_", model, "_p_fdr05.csv")
    write_table_command <- paste0("write.table(x = ", hemi_model_p_fdr, ", file = \"", filename,"\", row.names = FALSE, col.names = FALSE)")
    eval(parse(text=write_table_command))
    
    ### t
    filename <- paste0(homedir, "/baller/results/cbf_age_and_sex/", hemi, "_lm_", model, "_t_fdr05.csv")
    write_table_command <- paste0("write.table(x = ", hemi_model_t_fdr, ", file = \"", filename,"\", row.names = FALSE, col.names = FALSE)")
    eval(parse(text=write_table_command))
  }
}