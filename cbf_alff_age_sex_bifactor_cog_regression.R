library(mgcv)
library(dplyr)
library(ggplot2)

source("~/BBL/from_chead/ballerDepHeterogen/ballerDepHeterogenScripts/Hydra_functions.R")

#####################################################
# read in cognition, demographic and bifactor stuff #
#####################################################
set.seed(1)

demographics <- read.csv("/Users/eballer/BBL/from_chead/ballerDepHeterogen/data/n9498_demographics_go1_20161212.csv", header = TRUE, sep = ",") 

cnb_summary <- read.csv("/Users/eballer/BBL/from_chead/ballerDepHeterogen/data/n9498_cnb_zscores_fr_20170202.csv", header = TRUE, sep = ",") 

bifactor_summary <- read.csv("/Users/eballer/BBL/from_chead/ballerDepHeterogen/data/n9498_goassess_itemwise_bifactor_scores_20161219.csv", header = TRUE, sep = ",") 

T1_summary <- read.csv("/Users/eballer/BBL/from_chead/ballerDepHeterogen/data/neuroimaging/t1struct/n1601_t1QaData_20170306.csv")

rest_summary <- read.csv("/Users/eballer/BBL/from_chead/ballerDepHeterogen/data/neuroimaging/rest/n1601_RestQAData_20170714.csv")

asl_summary <- read.csv("/Users/eballer/BBL/from_chead/ballerDepHeterogen/data/neuroimaging/neuroimaging/n1601_PcaslQaData_20170403.csv")

#####################################################
#                 merge summary scores              #
#####################################################
#demo and cnb
df_summaries <- merge(demographics, cnb_summary, by = "bblid")

# add bifactor
df_summaries <- merge(df_summaries, bifactor_summary, by = "bblid")

#subset T1 to reduce size of df and merge
t1_subset <- subset(x = T1_summary, select = c("bblid", "scanid", "t1Exclude"))
df_summaries <- merge(df_summaries, t1_subset, by = "bblid")

#subset rest and merge
rest_subset <- subset(x = rest_summary, select = c("bblid", "scanid", "restRelMeanRMSMotion", "restExclude"))
df_summaries <- merge(df_summaries, rest_subset, by = c("bblid", "scanid"))

#subset asl and merge
asl_subset <- subset(x = asl_summary, select = c("bblid", "scanid", "pcaslRelMeanRMSMotion", "pcaslExclude"))
df_summaries <- merge(df_summaries, asl_subset, by = c("bblid", "scanid"))

df_summaries <- df_summaries[which(df_summaries$pcaslExclude == 0),]
df_summaries <- df_summaries[which(df_summaries$restExclude == 0),]

#####################################################
#             read in big matrices                  #
#####################################################
#files used for demographics:
#/data/joy/BBL/studies/pnc/n1601_dataFreeze/cnb/n1601_cnb_factor_scores_tymoore_20151006.csv
#/data/joy/BBL/studies/pnc/n1601_dataFreeze/clinical/n1601_goassess_itemwise_bifactor_scores_20161219.csv

#these files are based on LB's stuff which I believe parsed incorrectly, I have commented out
#lh_cbf_asl <- read.csv("/Users/eballer/BBL/imco/data/from_chead/n831_lh_cbf_asl_with_demo.csv", sep = ",")
#rh_cbf_asl <- read.csv("/Users/eballer/BBL/imco/data/from_chead/n831_rh_cbf_asl_with_demo.csv", sep = ",")

#### As of 20210112, I have made these myself using a c-ish looking bit of code in R ####
lh_cbf_asl <- read.csv("/Users/eballer/BBL/imco/data/from_chead/alffCbf_subjDemos_with_lh_831x10242.csv", sep = ",")
rh_cbf_asl <- read.csv("/Users/eballer/BBL/imco/data/from_chead/alffCbf_subjDemos_with_rh_831x10242.csv", sep = ",")

#some preprocessing
lh_cbf_asl <- subset(lh_cbf_asl, select = -scanid.y)
names(lh_cbf_asl)[3] <- "scanid"
lh_cbf_asl$sex <- as.factor(lh_cbf_asl$sex)
lh_cbf_asl$osex <- ordered(lh_cbf_asl$sex) #ordered factor for gam, males = 1, females = 2

rh_cbf_asl <- subset(rh_cbf_asl, select = -scanid.y)
names(rh_cbf_asl)[3] <- "scanid"
rh_cbf_asl$sex <- as.factor(rh_cbf_asl$sex)
rh_cbf_asl$osex <- ordered(rh_cbf_asl$sex) #ordered factor for gam, males = 1, females = 2

#####################################################
#                       merge                       #
#####################################################

#####################################################
#                       lm/gams                     #
#####################################################

# LEFT

#set vector to save info in - all are 10242 - # vertices in model
lh_gam_age <- vector(length = 10242)
lh_gam_sex <- vector(length = 10242)
lh_gam_age_sex <- vector(length = 10242)
lh_gam_mood <- vector(length = 10242)
lh_gam_psychopathology <- vector(length = 10242)
lh_gam_accuracy <- vector(length = 10242)


#get # of items in df for calculation of column)
numcolumns <- dim(lh_cbf_asl)[2]
#run gams models and store info in respective vectors
for (i in 1:10242) {
  curcol = (numcolumns - 10243 + i) # will start you counting at the right part of the df
  age_sex_model <- gam(lh_cbf_asl[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion +
                         osex + s(ageAtScan1, k = 4, fx = T) + 
                         s(ageAtScan1, by = osex, k = 4, fx = T) + 
                         mood + overall_psychopathology + Overall_Accuracy, 
                       data=lh_cbf_asl)
  
  #put pvalue in it's appropriate lm
  lh_gam_age[i] <- summary(age_sex_model)$s.table[1,4] #smooth term for ageAtScan1
  lh_gam_sex[i] <- summary(age_sex_model)$p.table[4,4] #linear term
  lh_gam_age_sex[i] <- summary(age_sex_model)$s.table[2,4] #smooth term for interaction
  
  lh_gam_mood[i] <- summary(age_sex_model)$p.table[5,4]
  lh_gam_psychopathology[i] <- summary(age_sex_model)$p.table[6,4]
  lh_gam_accuracy[i] <- summary(age_sex_model)$p.table[7,4] #accuracy term
  
}

#####################################################
#                     results                       #
#####################################################

lh_models <- c("lh_gam_mood", "lh_gam_age", "lh_gam_sex", "lh_gam_age_sex", "lh_gam_accuracy", "lh_gam_accuracy", "lh_gam_psychopathology")
for (model in lh_models){
  pfdr_anova <- eval(substitute(p.adjust(i, method="fdr"), list(i = as.name(model))))
  pfdr_anova <- as.data.frame(pfdr_anova)
  pfdr_round_anova <- round(pfdr_anova[pfdr_anova<0.05],3)
  print(paste0(model, " - num corrected: ", length(pfdr_round_anova)))
  
  #general uncorrected
  filename <- paste0("/Users/eballer/BBL/imco/results/coupling_", model, "_20210112.csv")
  write_table_command <- paste0("write.table(x = ", model, ", file = \"", filename,"\", row.names = FALSE, col.names = FALSE)")
  eval(parse(text=write_table_command))
  
  #corrected
  filename <- paste0("/Users/eballer/BBL/imco/results/coupling_", model, "fdr05_20210112.csv")
  write_table_command <- paste0("write.table(x = ", pfdr_anova, ", file = \"", filename,"\", row.names = FALSE, col.names = FALSE)")
  eval(parse(text=write_table_command))
  
}

########################################
### Save vectors for adam ##############
########################################
write.table(x = lh_gam_age, file = "/Users/eballer/BBL/imco/results/coupling_age_for_adam_lh_20210112.csv", row.names = FALSE, col.names=FALSE)
#####################
###### RIGHT ########
#####################


#set vector to save info in - all are 10242 - # vertices in model
rh_gam_age <- vector(length = 10242)
rh_gam_sex <- vector(length = 10242)
rh_gam_age_sex <- vector(length = 10242)
rh_gam_mood <- vector(length = 10242)
rh_gam_psychopathology <- vector(length = 10242)
rh_gam_accuracy <- vector(length = 10242)


#get # of items in df for calculation of column)
numcolumns <- dim(rh_cbf_asl)[2]
#run gams models and store info in respective vectors
for (i in 1:10242) {
  curcol = (numcolumns - 10243 + i) # will start you counting at the right part of the df
  age_sex_model <- gam(rh_cbf_asl[,curcol] ~ pcaslRelMeanRMSMotion + restRelMeanRMSMotion +
                         osex + s(ageAtScan1, k = 4, fx = T) + 
                         s(ageAtScan1, by = osex, k = 4, fx = T) + 
                         mood + overall_psychopathology + Overall_Accuracy, 
                       data=rh_cbf_asl)
  
  
  #put pvalue in it's appropriate lm
  
  rh_gam_sex[i] <- summary(age_sex_model)$p.table[4,4] #linear term
  rh_gam_age[i] <- summary(age_sex_model)$s.table[1,4] #smooth term for ageAtScan1
  rh_gam_age_sex[i] <- summary(age_sex_model)$s.table[2,4] #smooth term for interaction
  
  rh_gam_mood[i] <- summary(age_sex_model)$p.table[5,4]
  rh_gam_psychopathology[i] <- summary(age_sex_model)$p.table[6,4]
  rh_gam_accuracy[i] <- summary(age_sex_model)$p.table[7,4] #accuracy term
}

#####################################################
#                     results                       #
#####################################################

rh_models <- c("rh_gam_mood", "rh_gam_age", "rh_gam_sex", "rh_gam_age_sex", "rh_gam_accuracy", "rh_gam_psychopathology")
for (model in rh_models){
  pfdr_anova <- eval(substitute(p.adjust(i, method="fdr"), list(i = as.name(model))))
  pfdr_anova <- as.data.frame(pfdr_anova)
  pfdr_round_anova <- round(pfdr_anova[pfdr_anova<0.05],3)
  print(paste0(model, " - num corrected: ", length(pfdr_round_anova)))
  #save vectors
  
  #uncorrected
  filename <- paste0("/Users/eballer/BBL/imco/results/coupling_", model, "_20210112.csv")
  write_table_command <- paste0("write.table(x = ", model, ", file = \"", filename,"\", row.names = FALSE, col.names = FALSE)")
  eval(parse(text=write_table_command))
  
  #corrected
  filename <- paste0("/Users/eballer/BBL/imco/results/coupling_", model, "fdr05_20210112.csv")
  write_table_command <- paste0("write.table(x = ", pfdr_anova, ", file = \"", filename,"\", row.names = FALSE, col.names = FALSE)")
  eval(parse(text=write_table_command))
  
}


########################################
### Save vectors for adam ##############
########################################
write.table(x = rh_gam_age, file = "/Users/eballer/BBL/imco/results/coupling_age_for_adam_rh_20210112.csv", row.names = FALSE, col.names=FALSE)