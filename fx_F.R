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

# read in demos
alffCbf_subjDemos <- read.csv("/project/imco/baller/subjectLists/n831_alff_cbf_finalSample_imageOrder.csv")

#some verification preprocessing
alffCbf_subjDemos$sex <- as.factor(alffCbf_subjDemos$sex)
alffCbf_subjDemos$race <- as.factor(alffCbf_subjDemos$race)
alffCbf_subjDemos$race2 <- as.factor(alffCbf_subjDemos$race2)

#add osex category for use in gam later
alffCbf_subjDemos$osex <- ordered(alffCbf_subjDemos$sex)

#add psych bifactor scores
psych <- read.csv("/project/imco/pnc/clinical/n1601_goassess_itemwise_bifactor_scores_20161219.csv", header = TRUE)

#remove 4factorv2 from title
names(psych) <-gsub("_4factorv2", "", names(psych))

#merge
alffCbf_subjDemos <- merge(alffCbf_subjDemos, psych, by = "bblid")

#cognitive data
cog <- read.csv("/project/imco/pnc/cnb/n1601_cnb_factor_scores_tymoore_20151006.csv")
cog <- subset(cog, select = c("bblid", "Overall_Efficiency", "Overall_Accuracy", "Overall_Speed"))

#merge
alffCbf_subjDemos <- merge(alffCbf_subjDemos, cog, by = "bblid")

#drop the scanid.y
alffCbf_subjDemos <- subset(alffCbf_subjDemos, select = -scanid.y)

#rename scanid.x to scanid
names(alffCbf_subjDemos) <- gsub("scanid.x", "scanid", names(alffCbf_subjDemos)) 

#make list of bblid/scanid
bblid_scanid <- paste0(alffCbf_subjDemos$bblid, "_", alffCbf_subjDemos$datexscanid)

#####################
##### Left Side #####
#####################

#initiate matrix for storage
lh_831_x_10242_matrix <- matrix(nrow = 831, ncol = 10242)

#go through each subject, grab 5th column in asc, transpose and stick in matrix
# output is 831 x 10242 matrix
for (subj in 1:831) {
  
  bblid <- alffCbf_subjDemos$bblid[subj]
  datexscanid <- alffCbf_subjDemos$datexscanid[subj]
  file_path <- paste0("/project/imco/couplingSurfaceMaps/alffCbf/lh/stat/", bblid, "_", datexscanid, "_lh.coupling_coef_alff_cbf.fwhm15.fsaverage5.asc")
  alffCbf_data <- read.table(file_path, stringsAsFactors = FALSE)
  lh_831_x_10242_matrix[subj,] <- t(alffCbf_data$V5)
  
}

#append with demographics
alffCbf_subjDemos_with_lh_831x10242 <- cbind(alffCbf_subjDemos, lh_831_x_10242_matrix)

#write output
write.table(alffCbf_subjDemos_with_lh_831x10242, file = "/project/imco/baller/results/alffCbf_subjDemos_with_lh_831x10242.csv", sep = ",")

#####################
#### Right Side #####
#####################
#initiate matrix for storage
rh_831_x_10242_matrix <- matrix(nrow = 831, ncol = 10242)

#go through each subject, grab 5th column in asc, transpose and stick in matrix
# output is 831 x 10242 matrix
for (subj in 1:831) {
  
  bblid <- alffCbf_subjDemos$bblid[subj]
  datexscanid <- alffCbf_subjDemos$datexscanid[subj]
  file_path <- paste0("/project/imco/couplingSurfaceMaps/alffCbf/rh/stat/", bblid, "_", datexscanid, "_rh.coupling_coef_alff_cbf.fwhm15.fsaverage5.asc")
  alffCbf_data <- read.table(file_path, stringsAsFactors = FALSE)
  rh_831_x_10242_matrix[subj,] <- t(alffCbf_data$V5)
  
}

#append with demographics
alffCbf_subjDemos_with_rh_831x10242 <- cbind(alffCbf_subjDemos, rh_831_x_10242_matrix)

#write output
write.table(alffCbf_subjDemos_with_rh_831x10242, file = "/project/imco/baller/results/alffCbf_subjDemos_with_rh_831x10242.csv", sep = ",")

####-----------------------------------------------------------------------------####
####---------------------------End of Part 1- Making matrices---_----------------####
####-----------------------------------------------------------------------------####

#####################################################################################
####                   Run regression, both left and right                       ####
#####################################################################################

#make easier to reference names
lh_cbf_asl <- alffCbf_subjDemos_with_lh_831x10242 #can also read directly from files if you'd like
rh_cbf_asl <- alffCbf_subjDemos_with_rh_831x10242


#####################################################
#                       lm/gams                     #
#####################################################

#initialize vectors for models

hemis <- c("lh", "rh") #hemispheres
models <- c("age", "sex", "age_sex", "accuracy", "speed", "efficiency", "mood", "psychopathology") #models of interest
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
lm_models <- c("age", "accuracy", "speed", "efficiency")
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
  age_sex_model <- gam(lh_cbf_asl[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion +
                         osex + s(ageAtScan1, k = 4, fx = F), data=lh_cbf_asl)
  
  age_sex_intx_model <- gam(lh_cbf_asl[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion +
                              osex + s(ageAtScan1, k = 4, fx = F) + s(ageAtScan1, by = osex, k = 4, fx = F), data=lh_cbf_asl)
  
  mood_model <- gam(lh_cbf_asl[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion +
                      osex + s(ageAtScan1, k = 4, fx = F) + mood, data=lh_cbf_asl)
  
  psychopathology_model <- gam(lh_cbf_asl[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion +
                                 osex + s(ageAtScan1, k = 4, fx = F) + overall_psychopathology, data=lh_cbf_asl)
  
  accuracy_model <- gam(lh_cbf_asl[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion +
                          osex + s(ageAtScan1, k = 4, fx = F) + Overall_Accuracy, data=lh_cbf_asl)
  
  speed_model <- gam(lh_cbf_asl[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion +
                          osex + s(ageAtScan1, k = 4, fx = F) + Overall_Speed, data=lh_cbf_asl)
  
  efficiency_model <- gam(lh_cbf_asl[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion +
                          osex + s(ageAtScan1, k = 4, fx = F) + Overall_Efficiency, data=lh_cbf_asl)
  
  
  age_lm_model <- lm(lh_cbf_asl[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion + ageAtScan1, 
                     data=lh_cbf_asl)
  
  accuracy_lm_model <- lm(lh_cbf_asl[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion + ageAtScan1 + Overall_Accuracy, 
                     data=lh_cbf_asl)
  
  speed_lm_model <- lm(lh_cbf_asl[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion + ageAtScan1 + Overall_Speed, 
                     data=lh_cbf_asl)
  
  efficiency_lm_model <- lm(lh_cbf_asl[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion + ageAtScan1 + Overall_Efficiency, 
                     data=lh_cbf_asl)
  
  #put pvalue in it's appropriate lm
  lh_gam_age_p_uncor[i] <- summary(age_sex_model)$s.table[1,4] #smooth term for ageAtScan1
  lh_gam_sex_p_uncor[i] <- summary(age_sex_model)$p.table[4,4] #linear term
  lh_gam_age_sex_p_uncor[i] <- summary(age_sex_intx_model)$s.table[2,4] #smooth term for interaction, this was changed
  
  lh_gam_mood_p_uncor[i] <- summary(mood_model)$p.table[5,4]
  lh_gam_psychopathology_p_uncor[i] <- summary(psychopathology_model)$p.table[5,4]
  lh_gam_accuracy_p_uncor[i] <- summary(accuracy_model)$p.table[5,4] #accuracy term
  lh_gam_speed_p_uncor[i] <- summary(speed_model)$p.table[5,4] #speed term
  lh_gam_efficiency_p_uncor[i] <- summary(efficiency_model)$p.table[5,4] #efficiency term
  
  #lm to assess directionality
  lh_lm_age_p_uncor[i] <- summary(age_lm_model)$coeff[4,4]
  lh_lm_accuracy_p_uncor[i] <- summary(accuracy_lm_model)$coeff[5,4]
  lh_lm_speed_p_uncor[i] <- summary(speed_lm_model)$coeff[5,4]
  lh_lm_efficiency_p_uncor[i] <- summary(efficiency_lm_model)$coeff[5,4]
  
  #pull tvalue into its appropriate lm
  lh_gam_age_t_uncor[i] <- summary(age_sex_model)$s.table[1,3] #smooth term for ageAtScan1
  lh_gam_sex_t_uncor[i] <- summary(age_sex_model)$p.table[4,3] #linear term
  lh_gam_age_sex_t_uncor[i] <- summary(age_sex_intx_model)$s.table[2,3] #smooth term for interaction
  
  lh_gam_mood_t_uncor[i] <- summary(mood_model)$p.table[5,3]
  lh_gam_psychopathology_t_uncor[i] <- summary(psychopathology_model)$p.table[5,3]
  lh_gam_accuracy_t_uncor[i] <- summary(accuracy_model)$p.table[5,3] #accuracy term
  lh_gam_speed_t_uncor[i] <- summary(speed_model)$p.table[5,3] #accuracy term
  lh_gam_efficiency_t_uncor[i] <- summary(efficiency_model)$p.table[5,3] #accuracy term
  
  #lm to assess directionality
  lh_lm_age_t_uncor[i] <- summary(age_lm_model)$coeff[4,3]
  lh_lm_accuracy_t_uncor[i] <- summary(accuracy_lm_model)$coeff[5,3]
  lh_lm_speed_t_uncor[i] <- summary(speed_lm_model)$coeff[5,3]
  lh_lm_efficiency_t_uncor[i] <- summary(efficiency_lm_model)$coeff[5,3]
}

#####################
###### RIGHT ########
#####################


#get # of items in df for calculation of column)
numcolumns <- dim(rh_cbf_asl)[2]
#run gams models and store info in respective vectors
for (i in 1:10242) {
  curcol = (numcolumns - 10242 + i) # will start you counting at the right part of the df
  age_sex_model <- gam(rh_cbf_asl[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion +
                         osex + s(ageAtScan1, k = 4, fx = F), data=rh_cbf_asl)
  
  age_sex_intx_model <- gam(rh_cbf_asl[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion +
                              osex + s(ageAtScan1, k = 4, fx = F) + s(ageAtScan1, by = osex, k = 4, fx = F), data=rh_cbf_asl)
  
  mood_model <- gam(rh_cbf_asl[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion +
                      osex + s(ageAtScan1, k = 4, fx = F) + mood, data=rh_cbf_asl)
  
  psychopathology_model <- gam(rh_cbf_asl[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion +
                                 osex + s(ageAtScan1, k = 4, fx = F) + overall_psychopathology, data=rh_cbf_asl)
  
  accuracy_model <- gam(rh_cbf_asl[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion +
                          osex + s(ageAtScan1, k = 4, fx = F) + Overall_Accuracy, data=rh_cbf_asl)
  
  speed_model <- gam(rh_cbf_asl[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion +
                       osex + s(ageAtScan1, k = 4, fx = F) + Overall_Speed, data=rh_cbf_asl)
  
  efficiency_model <- gam(rh_cbf_asl[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion +
                            osex + s(ageAtScan1, k = 4, fx = F) + Overall_Efficiency, data=rh_cbf_asl)
  
  age_lm_model <- lm(rh_cbf_asl[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion + ageAtScan1, data=rh_cbf_asl)
 

  accuracy_lm_model <- lm(rh_cbf_asl[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion + ageAtScan1 + Overall_Accuracy, 
                          data=rh_cbf_asl)
  
  speed_lm_model <- lm(rh_cbf_asl[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion + ageAtScan1 + Overall_Speed, 
                       data=rh_cbf_asl)
  
  efficiency_lm_model <- lm(rh_cbf_asl[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion + ageAtScan1 + Overall_Efficiency, 
                            data=rh_cbf_asl)
  
  #put pvalue in it's appropriate lm
  rh_gam_age_p_uncor[i] <- summary(age_sex_model)$s.table[1,4] #smooth term for ageAtScan1
  rh_gam_sex_p_uncor[i] <- summary(age_sex_model)$p.table[4,4] #linear term
  rh_gam_age_sex_p_uncor[i] <- summary(age_sex_intx_model)$s.table[2,4] #smooth term for interaction, this was changed
  
  rh_gam_mood_p_uncor[i] <- summary(mood_model)$p.table[5,4]
  rh_gam_psychopathology_p_uncor[i] <- summary(psychopathology_model)$p.table[5,4]
  rh_gam_accuracy_p_uncor[i] <- summary(accuracy_model)$p.table[5,4] #accuracy term
  rh_gam_speed_p_uncor[i] <- summary(speed_model)$p.table[5,4] #speed term
  rh_gam_efficiency_p_uncor[i] <- summary(efficiency_model)$p.table[5,4] #efficiency term
  
  #lm to assess directionality
  rh_lm_age_p_uncor[i] <- summary(age_lm_model)$coeff[4,4]
  rh_lm_accuracy_p_uncor[i] <- summary(accuracy_lm_model)$coeff[5,4]
  rh_lm_speed_p_uncor[i] <- summary(speed_lm_model)$coeff[5,4]
  rh_lm_efficiency_p_uncor[i] <- summary(efficiency_lm_model)$coeff[5,4]
  
  #pull tvalue into its appropriate lm
  rh_gam_age_t_uncor[i] <- summary(age_sex_model)$s.table[1,3] #smooth term for ageAtScan1
  rh_gam_sex_t_uncor[i] <- summary(age_sex_model)$p.table[4,3] #linear term
  rh_gam_age_sex_t_uncor[i] <- summary(age_sex_intx_model)$s.table[2,3] #smooth term for interaction
  
  rh_gam_mood_t_uncor[i] <- summary(mood_model)$p.table[5,3]
  rh_gam_psychopathology_t_uncor[i] <- summary(psychopathology_model)$p.table[5,3]
  rh_gam_accuracy_t_uncor[i] <- summary(accuracy_model)$p.table[5,3] #accuracy term
  rh_gam_speed_t_uncor[i] <- summary(speed_model)$p.table[5,3] #accuracy term
  rh_gam_efficiency_t_uncor[i] <- summary(efficiency_model)$p.table[5,3] #accuracy term
  
  #lm to assess directionality
  rh_lm_age_t_uncor[i] <- summary(age_lm_model)$coeff[4,3]
  rh_lm_accuracy_t_uncor[i] <- summary(accuracy_lm_model)$coeff[5,3]
  rh_lm_speed_t_uncor[i] <- summary(speed_lm_model)$coeff[5,3]
  rh_lm_efficiency_t_uncor[i] <- summary(efficiency_lm_model)$coeff[5,3]
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
    filename <- paste0("/project/imco/baller/processed_data/regression_matrices_n831/20210122_fx_F_csvs/", hemi_model_p_unc, ".csv")
    write_table_command <- paste0("write.table(x = ", hemi_model_p_unc, ", file = \"", filename,"\", row.names = FALSE, col.names = FALSE)")
    eval(parse(text=write_table_command))
    
    ### t
    
    filename <- paste0("/project/imco/baller/processed_data/regression_matrices_n831/20210122_fx_F_csvs/", hemi_model_t_unc, ".csv")
    write_table_command <- paste0("write.table(x = ", hemi_model_t_unc, ", file = \"", filename,"\", row.names = FALSE, col.names = FALSE)")
    eval(parse(text=write_table_command))
    
    ## corrected ##
    
    ### p
    filename <- paste0("/project/imco/baller/processed_data/regression_matrices_n831/20210122_fx_F_csvs/", hemi, "_gam_", model, "_p_fdr05.csv")
    write_table_command <- paste0("write.table(x = ", hemi_model_p_fdr, ", file = \"", filename,"\", row.names = FALSE, col.names = FALSE)")
    eval(parse(text=write_table_command))
    
    ### t
    filename <- paste0("/project/imco/baller/processed_data/regression_matrices_n831/20210122_fx_F_csvs/", hemi, "_gam_", model, "_t_fdr05.csv")
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
    filename <- paste0("/project/imco/baller/processed_data/regression_matrices_n831/20210122_fx_F_csvs/", hemi_model_p_unc, ".csv")
    write_table_command <- paste0("write.table(x = ", hemi_model_p_unc, ", file = \"", filename,"\", row.names = FALSE, col.names = FALSE)")
    eval(parse(text=write_table_command))
    
    ### t
    
    filename <- paste0("/project/imco/baller/processed_data/regression_matrices_n831/20210122_fx_F_csvs/", hemi_model_t_unc, ".csv")
    write_table_command <- paste0("write.table(x = ", hemi_model_t_unc, ", file = \"", filename,"\", row.names = FALSE, col.names = FALSE)")
    eval(parse(text=write_table_command))
    
    ## corrected ##
    
    ### p
    filename <- paste0("/project/imco/baller/processed_data/regression_matrices_n831/20210122_fx_F_csvs/", hemi, "_lm_", model, "_p_fdr05.csv")
    write_table_command <- paste0("write.table(x = ", hemi_model_p_fdr, ", file = \"", filename,"\", row.names = FALSE, col.names = FALSE)")
    eval(parse(text=write_table_command))
    
    ### t
    filename <- paste0("/project/imco/baller/processed_data/regression_matrices_n831/20210122_fx_F_csvs/", hemi, "_lm_", model, "_t_fdr05.csv")
    write_table_command <- paste0("write.table(x = ", hemi_model_t_fdr, ", file = \"", filename,"\", row.names = FALSE, col.names = FALSE)")
    eval(parse(text=write_table_command))
  }
}