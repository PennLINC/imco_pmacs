############################
### Summarize Mask Data ####
############################
###  Author: Erica Baller ##
############################
###  2/23/2020           ###
############################

#Pre: Takes in outputs from convert_Ts_to_masks_for_display.R (i.e. /Users/eballer/BBL/imco/pmacs/PMACS_remote/baller/results/coupling_accuracy/rh_pos_lm_sex_t_fdr05_Hill2010_evo_fsaverage5.csv_10242.csv)
#Post: Table that includes summaries of # of vertices within each mask that are in the direction as expected (i.e. mask has + value where direction is +)
#Uses: I wanted a useful way to summarize which mask had best correspondence to the coupling data. My assumption is that the mask that shares the most directionality 
      #with my coupling maps will most closely match on to function
#Dependencies: Any R will do. I used 3.2.5

#homedir
homedir <- '/Users/eballer/BBL/imco/pmacs/PMACS_remote/'
#homedir <- '/project/imco/'

workingdir <- '/baller/results/coupling_accuracy/'



######### functions ########
num_positives <- function(x){
  return(length(which(x>0)))
}

num_negatives <- function(x){
  return(length(which(x<0)))
}

num_nonzero <- function(x) {
  return(length(which(abs(x)>0)))
}
##### Read in the CSVs #####

#positive direction lm age, sex, exec_acc
gi_pos_age <- rbind(read.csv(paste0(homedir, workingdir, "lh_pos_lm_age_t_fdr05_GI_fsaverage5_10242.csv_10242.csv"), header = F),
                    read.csv(paste0(homedir, workingdir, "rh_pos_lm_age_t_fdr05_GI_fsaverage5_10242.csv_10242.csv"), header = F))
hill2010_pos_age <- rbind(read.csv(paste0(homedir, workingdir, "lh_pos_lm_age_t_fdr05_Hill2010_evo_fsaverage5.csv_10242.csv"), header = F),
                          read.csv(paste0(homedir, workingdir, "rh_pos_lm_age_t_fdr05_Hill2010_evo_fsaverage5.csv_10242.csv"), header = F))
mean_cbf_pos_age <- rbind(read.csv(paste0(homedir, workingdir, "lh_pos_lm_age_t_fdr05_MeanCBF.fsaverage5.csv_10242.csv"), header = F),
                          read.csv(paste0(homedir, workingdir, "rh_pos_lm_age_t_fdr05_MeanCBF.fsaverage5.csv_10242.csv"), header = F))
cmrglu_pos_age <- rbind(read.csv(paste0(homedir, workingdir, "lh_pos_lm_age_t_fdr05_CMRGlu_fsaverage5.csv_10242.csv"), header = F),
                        read.csv(paste0(homedir, workingdir, "rh_pos_lm_age_t_fdr05_CMRGlu_fsaverage5.csv_10242.csv"), header = F))
allometric_scaling_pos_age <- rbind(read.csv(paste0(homedir, workingdir, "lh_pos_lm_age_t_fdr05_AllometricScaling_fsaverage5.csv_10242.csv"), header = F),
                                    read.csv(paste0(homedir, workingdir, "rh_pos_lm_age_t_fdr05_AllometricScaling_fsaverage5.csv_10242.csv"), header = F))
  
gi_pos_sex <- rbind(read.csv(paste0(homedir, workingdir, "lh_pos_lm_sex_t_fdr05_GI_fsaverage5_10242.csv_10242.csv"), header = F),
                    read.csv(paste0(homedir, workingdir, "rh_pos_lm_sex_t_fdr05_GI_fsaverage5_10242.csv_10242.csv"), header = F))
hill2010_pos_sex <- rbind(read.csv(paste0(homedir, workingdir, "lh_pos_lm_sex_t_fdr05_Hill2010_evo_fsaverage5.csv_10242.csv"), header = F),
                          read.csv(paste0(homedir, workingdir, "rh_pos_lm_sex_t_fdr05_Hill2010_evo_fsaverage5.csv_10242.csv"), header = F))
mean_cbf_pos_sex <- rbind(read.csv(paste0(homedir, workingdir, "lh_pos_lm_sex_t_fdr05_MeanCBF.fsaverage5.csv_10242.csv"), header = F),
                          read.csv(paste0(homedir, workingdir, "rh_pos_lm_sex_t_fdr05_MeanCBF.fsaverage5.csv_10242.csv"), header = F))
cmrglu_pos_sex <- rbind(read.csv(paste0(homedir, workingdir, "lh_pos_lm_sex_t_fdr05_CMRGlu_fsaverage5.csv_10242.csv"), header = F),
                        read.csv(paste0(homedir, workingdir, "rh_pos_lm_sex_t_fdr05_CMRGlu_fsaverage5.csv_10242.csv"), header = F))
allometric_scaling_pos_sex <- rbind(read.csv(paste0(homedir, workingdir, "lh_pos_lm_sex_t_fdr05_AllometricScaling_fsaverage5.csv_10242.csv"), header = F),
                                    read.csv(paste0(homedir, workingdir, "rh_pos_lm_sex_t_fdr05_AllometricScaling_fsaverage5.csv_10242.csv"), header = F))

gi_pos_exec_accuracy <- rbind(read.csv(paste0(homedir, workingdir, "lh_pos_lm_exec_accuracy_t_fdr05_GI_fsaverage5_10242.csv_10242.csv"), header = F),
                    read.csv(paste0(homedir, workingdir, "rh_pos_lm_exec_accuracy_t_fdr05_GI_fsaverage5_10242.csv_10242.csv"), header = F))
hill2010_pos_exec_accuracy <- rbind(read.csv(paste0(homedir, workingdir, "lh_pos_lm_exec_accuracy_t_fdr05_Hill2010_evo_fsaverage5.csv_10242.csv"), header = F),
                          read.csv(paste0(homedir, workingdir, "rh_pos_lm_exec_accuracy_t_fdr05_Hill2010_evo_fsaverage5.csv_10242.csv"), header = F))
mean_cbf_pos_exec_accuracy <- rbind(read.csv(paste0(homedir, workingdir, "lh_pos_lm_exec_accuracy_t_fdr05_MeanCBF.fsaverage5.csv_10242.csv"), header = F),
                          read.csv(paste0(homedir, workingdir, "rh_pos_lm_exec_accuracy_t_fdr05_MeanCBF.fsaverage5.csv_10242.csv"), header = F))
cmrglu_pos_exec_accuracy <- rbind(read.csv(paste0(homedir, workingdir, "lh_pos_lm_exec_accuracy_t_fdr05_CMRGlu_fsaverage5.csv_10242.csv"), header = F),
                        read.csv(paste0(homedir, workingdir, "rh_pos_lm_exec_accuracy_t_fdr05_CMRGlu_fsaverage5.csv_10242.csv"), header = F))
allometric_scaling_pos_exec_accuracy <- rbind(read.csv(paste0(homedir, workingdir, "lh_pos_lm_exec_accuracy_t_fdr05_AllometricScaling_fsaverage5.csv_10242.csv"), header = F),
                                    read.csv(paste0(homedir, workingdir, "rh_pos_lm_exec_accuracy_t_fdr05_AllometricScaling_fsaverage5.csv_10242.csv"), header = F))

#make data_frames
positive_mask_age <- data.frame(gi_pos_age,hill2010_pos_age,mean_cbf_pos_age,cmrglu_pos_age, allometric_scaling_pos_age)
names(positive_mask_age) <- c("gi", "hill2010", "mean_cbf", "cmrglu", "allometric_scaling")
for (i in names(positive_mask_age)) {
  eval(parse(text=paste0("positive_mask_age$",i,"[10243] <- num_positives(positive_mask_age$",i,")")))
  eval(parse(text=paste0("positive_mask_age$",i,"[10244] <- num_nonzero(positive_mask_age$",i,")")))
}

positive_mask_sex <- data.frame(gi_pos_sex,hill2010_pos_sex,mean_cbf_pos_sex,cmrglu_pos_sex, allometric_scaling_pos_sex)
names(positive_mask_sex) <- c("gi", "hill2010", "mean_cbf", "cmrglu", "allometric_scaling")
for (i in names(positive_mask_sex)) {
  eval(parse(text=paste0("positive_mask_sex$",i,"[10243] <- num_positives(positive_mask_sex$",i,")")))
  eval(parse(text=paste0("positive_mask_sex$",i,"[10244] <- num_nonzero(positive_mask_sex$",i,")")))
}


positive_mask_exec_accuracy <- data.frame(gi_pos_exec_accuracy,hill2010_pos_exec_accuracy,mean_cbf_pos_exec_accuracy,cmrglu_pos_exec_accuracy, allometric_scaling_pos_exec_accuracy)
names(positive_mask_exec_accuracy) <- c("gi", "hill2010", "mean_cbf", "cmrglu", "allometric_scaling")
for (i in names(positive_mask_exec_accuracy)) {
  eval(parse(text=paste0("positive_mask_exec_accuracy$",i,"[10243] <- num_positives(positive_mask_exec_accuracy$",i,")")))
  eval(parse(text=paste0("positive_mask_exec_accuracy$",i,"[10244] <- num_nonzero(positive_mask_exec_accuracy$",i,")")))
  
}


#negative direction lm age, sex, exec_acc
gi_neg_age <- rbind(read.csv(paste0(homedir, workingdir, "lh_neg_lm_age_t_fdr05_GI_fsaverage5_10242.csv_10242.csv"), header = F),
                    read.csv(paste0(homedir, workingdir, "rh_neg_lm_age_t_fdr05_GI_fsaverage5_10242.csv_10242.csv"), header = F))
hill2010_neg_age <- rbind(read.csv(paste0(homedir, workingdir, "lh_neg_lm_age_t_fdr05_Hill2010_evo_fsaverage5.csv_10242.csv"), header = F),
                          read.csv(paste0(homedir, workingdir, "rh_neg_lm_age_t_fdr05_Hill2010_evo_fsaverage5.csv_10242.csv"), header = F))
mean_cbf_neg_age <- rbind(read.csv(paste0(homedir, workingdir, "lh_neg_lm_age_t_fdr05_MeanCBF.fsaverage5.csv_10242.csv"), header = F),
                          read.csv(paste0(homedir, workingdir, "rh_neg_lm_age_t_fdr05_MeanCBF.fsaverage5.csv_10242.csv"), header = F))
cmrglu_neg_age <- rbind(read.csv(paste0(homedir, workingdir, "lh_neg_lm_age_t_fdr05_CMRGlu_fsaverage5.csv_10242.csv"), header = F),
                        read.csv(paste0(homedir, workingdir, "rh_neg_lm_age_t_fdr05_CMRGlu_fsaverage5.csv_10242.csv"), header = F))
allometric_scaling_neg_age <- rbind(read.csv(paste0(homedir, workingdir, "lh_neg_lm_age_t_fdr05_AllometricScaling_fsaverage5.csv_10242.csv"), header = F),
                                    read.csv(paste0(homedir, workingdir, "rh_neg_lm_age_t_fdr05_AllometricScaling_fsaverage5.csv_10242.csv"), header = F))

gi_neg_sex <- rbind(read.csv(paste0(homedir, workingdir, "lh_neg_lm_sex_t_fdr05_GI_fsaverage5_10242.csv_10242.csv"), header = F),
                    read.csv(paste0(homedir, workingdir, "rh_neg_lm_sex_t_fdr05_GI_fsaverage5_10242.csv_10242.csv"), header = F))
hill2010_neg_sex <- rbind(read.csv(paste0(homedir, workingdir, "lh_neg_lm_sex_t_fdr05_Hill2010_evo_fsaverage5.csv_10242.csv"), header = F),
                          read.csv(paste0(homedir, workingdir, "rh_neg_lm_sex_t_fdr05_Hill2010_evo_fsaverage5.csv_10242.csv"), header = F))
mean_cbf_neg_sex <- rbind(read.csv(paste0(homedir, workingdir, "lh_neg_lm_sex_t_fdr05_MeanCBF.fsaverage5.csv_10242.csv"), header = F),
                          read.csv(paste0(homedir, workingdir, "rh_neg_lm_sex_t_fdr05_MeanCBF.fsaverage5.csv_10242.csv"), header = F))
cmrglu_neg_sex <- rbind(read.csv(paste0(homedir, workingdir, "lh_neg_lm_sex_t_fdr05_CMRGlu_fsaverage5.csv_10242.csv"), header = F),
                        read.csv(paste0(homedir, workingdir, "rh_neg_lm_sex_t_fdr05_CMRGlu_fsaverage5.csv_10242.csv"), header = F))
allometric_scaling_neg_sex <- rbind(read.csv(paste0(homedir, workingdir, "lh_neg_lm_sex_t_fdr05_AllometricScaling_fsaverage5.csv_10242.csv"), header = F),
                                    read.csv(paste0(homedir, workingdir, "rh_neg_lm_sex_t_fdr05_AllometricScaling_fsaverage5.csv_10242.csv"), header = F))

gi_neg_exec_accuracy <- rbind(read.csv(paste0(homedir, workingdir, "lh_neg_lm_exec_accuracy_t_fdr05_GI_fsaverage5_10242.csv_10242.csv"), header = F),
                              read.csv(paste0(homedir, workingdir, "rh_neg_lm_exec_accuracy_t_fdr05_GI_fsaverage5_10242.csv_10242.csv"), header = F))
hill2010_neg_exec_accuracy <- rbind(read.csv(paste0(homedir, workingdir, "lh_neg_lm_exec_accuracy_t_fdr05_Hill2010_evo_fsaverage5.csv_10242.csv"), header = F),
                                    read.csv(paste0(homedir, workingdir, "rh_neg_lm_exec_accuracy_t_fdr05_Hill2010_evo_fsaverage5.csv_10242.csv"), header = F))
mean_cbf_neg_exec_accuracy <- rbind(read.csv(paste0(homedir, workingdir, "lh_neg_lm_exec_accuracy_t_fdr05_MeanCBF.fsaverage5.csv_10242.csv"), header = F),
                                    read.csv(paste0(homedir, workingdir, "rh_neg_lm_exec_accuracy_t_fdr05_MeanCBF.fsaverage5.csv_10242.csv"), header = F))
cmrglu_neg_exec_accuracy <- rbind(read.csv(paste0(homedir, workingdir, "lh_neg_lm_exec_accuracy_t_fdr05_CMRGlu_fsaverage5.csv_10242.csv"), header = F),
                                  read.csv(paste0(homedir, workingdir, "rh_neg_lm_exec_accuracy_t_fdr05_CMRGlu_fsaverage5.csv_10242.csv"), header = F))
allometric_scaling_neg_exec_accuracy <- rbind(read.csv(paste0(homedir, workingdir, "lh_neg_lm_exec_accuracy_t_fdr05_AllometricScaling_fsaverage5.csv_10242.csv"), header = F),
                                              read.csv(paste0(homedir, workingdir, "rh_neg_lm_exec_accuracy_t_fdr05_AllometricScaling_fsaverage5.csv_10242.csv"), header = F))

#make data_frames
negative_mask_age <- data.frame(gi_neg_age,hill2010_neg_age,mean_cbf_neg_age,cmrglu_neg_age, allometric_scaling_neg_age)
names(negative_mask_age) <- c("gi", "hill2010", "mean_cbf", "cmrglu", "allometric_scaling")
for (i in names(negative_mask_age)) {
  eval(parse(text=paste0("negative_mask_age$",i,"[10243] <- num_negatives(negative_mask_age$",i,")")))
  eval(parse(text=paste0("negative_mask_age$",i,"[10244] <- num_nonzero(negative_mask_age$",i,")")))
}

negative_mask_sex <- data.frame(gi_neg_sex,hill2010_neg_sex,mean_cbf_neg_sex,cmrglu_neg_sex, allometric_scaling_neg_sex)
names(negative_mask_sex) <- c("gi", "hill2010", "mean_cbf", "cmrglu", "allometric_scaling")
for (i in names(negative_mask_sex)) {
  eval(parse(text=paste0("negative_mask_sex$",i,"[10243] <- num_negatives(negative_mask_sex$",i,")")))
  eval(parse(text=paste0("negative_mask_sex$",i,"[10244] <- num_nonzero(negative_mask_sex$",i,")")))
  
}


negative_mask_exec_accuracy <- data.frame(gi_neg_exec_accuracy,hill2010_neg_exec_accuracy,mean_cbf_neg_exec_accuracy,cmrglu_neg_exec_accuracy, allometric_scaling_neg_exec_accuracy)
names(negative_mask_exec_accuracy) <- c("gi", "hill2010", "mean_cbf", "cmrglu", "allometric_scaling")
for (i in names(negative_mask_exec_accuracy)) {
  eval(parse(text=paste0("negative_mask_exec_accuracy$",i,"[10243] <- num_negatives(negative_mask_exec_accuracy$",i,")")))
  eval(parse(text=paste0("negative_mask_exec_accuracy$",i,"[10244] <- num_nonzero(negative_mask_exec_accuracy$",i,")")))
}

all_together_df <- data.frame(rbind(positive_mask_age[10243,],
                           positive_mask_sex[10243,], 
                           positive_mask_exec_accuracy[10243,], 
                           negative_mask_age[10243,], 
                           negative_mask_sex[10243,],
                           negative_mask_exec_accuracy[10243,]))

all_together_nonzero_df <- data.frame(rbind((positive_mask_age[10244,]),
                                            (positive_mask_sex[10244,]), 
                                            (positive_mask_exec_accuracy[10244,]), 
                                            (negative_mask_age[10244,]),
                                            (negative_mask_sex[10244,]), 
                                            (negative_mask_exec_accuracy[10244,])))

#all_together_df<- rbind(all_together_df, all_together_nonzero_df[,1])
row.names(all_together_df) <- c("positive_mask_age", "positive_mask_sex", "positive_mask_exec_accuracy",
                                "negative_mask_age", "negative_mask_sex", "negative_mask_exec_accuracy")#, "total_num_vertices_in_mask")

row.names(all_together_nonzero_df) <- c("positive_mask_age", "positive_mask_sex", "positive_mask_exec_accuracy",
                                "negative_mask_age", "negative_mask_sex", "negative_mask_exec_accuracy")

#calculating the percentage -i.e. number of voxels going in the right direction

all_together_df_percentages <- (all_together_df/all_together_nonzero_df)*100

write.csv(all_together_df, paste0(homedir, "baller/results/coupling_accuracy/mask_summaries_correct_directionality.csv"))
write.csv(all_together_df_percentages, paste0(homedir, "baller/results/coupling_accuracy/mask_summaries_correct_directionality_percentages.csv"))
