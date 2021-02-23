#### Makes the 831x10242 matrices, both left and right, for use in later analyses

# read in demos
alffCbf_subjDemos <- read.csv("/data/jux/BBL/projects/coupling/subjectsLists/n831_rest_cbf_finalSample_imageOrder.csv")
alffCbf_subjDemos$sex <- as.factor(alffCbf_subjDemos$sex)

##### EB added 10/08/2020
#add psych bifactor scores
psych <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/clinical/n1601_goassess_itemwise_bifactor_scores_20161219.csv", header = TRUE)
#remove 4factorv2 from title
names(psych) <-gsub("_4factorv2", "", names(psych))

alffCbf_subjDemos <- merge(alffCbf_subjDemos, psych, by = "bblid")

#cognitive data
cog <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/cnb/n1601_cnb_factor_scores_tymoore_20151006.csv")
cog <- subset(cog, select = c("bblid", "Overall_Efficiency", "Overall_Accuracy", "Overall_Speed"))

alffCbf_subjDemos <- merge(alffCbf_subjDemos, cog, by = "bblid")

#make list of bblid/scanid
bblid_scanid <- paste0(alffCbf_subjDemos$bblid, "_", alffCbf_subjDemos$datexscanid)


##### Left Side #####
#initiate matrix for storage
lh_831_x_10242_matrix <- matrix(nrow = 831, ncol = 10242)

#go through each subject, grab 5th column in asc, transpose and stick in matrix
# output is 831 x 10242 matrix
for (subj in 1:831) {
  #print(bblid_scanid[subj])
  file_path <- paste0("/data/jux/BBL/projects/coupling/couplingSurfaceMaps/alffCbf/lh/stat/", bblid_scanid[subj], "_lh.coupling_coef_alff_cbf.fwhm15.fsaverage5.asc")
  #print(file_path)
  alffCbf_data <- read.table(file_path, stringsAsFactors = FALSE)
  lh_831_x_10242_matrix[subj,] <- t(alffCbf_data$V5)
}
#append with demographics
alffCbf_subjDemos_with_lh_831x10242 <- cbind(alffCbf_subjDemos, lh_831_x_10242_matrix)

#write output
write.table(alffCbf_subjDemos_with_lh_831x10242, file = "/data/jux/BBL/projects/coupling/coupling_test_eb_20200918/results/alffCbf_subjDemos_with_lh_831x10242.csv", sep = ",")


#### Right Side #######
#initiate matrix for storage
rh_831_x_10242_matrix <- matrix(nrow = 831, ncol = 10242)

#go through each subject, grab 5th column in asc, transpose and stick in matrix
# output is 831 x 10242 matrix
for (subj in 1:831) {
  #print(bblid_scanid[subj])
  file_path <- paste0("/data/jux/BBL/projects/coupling/couplingSurfaceMaps/alffCbf/rh/stat/", bblid_scanid[subj], "_rh.coupling_coef_alff_cbf.fwhm15.fsaverage5.asc")
  #print(file_path)
  alffCbf_data <- read.table(file_path, stringsAsFactors = FALSE)
  rh_831_x_10242_matrix[subj,] <- t(alffCbf_data$V5)
}
#append with demographics
alffCbf_subjDemos_with_rh_831x10242 <- cbind(alffCbf_subjDemos, rh_831_x_10242_matrix)

#write output
write.table(alffCbf_subjDemos_with_rh_831x10242, file = "/data/jux/BBL/projects/coupling/coupling_test_eb_20200918/results/alffCbf_subjDemos_with_rh_831x10242.csv", sep = ",")
