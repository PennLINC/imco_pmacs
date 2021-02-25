################################
## 2/19/2021 Convert Parcels ###
################################

################################
#### Author: Erica Baller ######
################################

### This script emerged out of the desire to take vertex output and display it within different networks/parcels

#pre: input: 
#1) a vector(n=10242) of T values, with 0s indicating vertices you don't want to include
#- the default is to make maps for ALL interesting analyses (coupling, alff, cbf; lm and gam; uncor and fdr05 corrected). 
#If you would like to change it, please do so in the following section by commenting out:
#analyses <- c("coupling_accuracy", "cbf_accuracy", "alff_accuracy") 
#models <- c("gam_age", "lm_age", "lm_sex", "lm_accuracy", "lm_exec_accuracy")
#corrs <- c("uncor", "fdr05") #correction
#2) mask vectors
#a) lh and rh 10242 matrices with the vertex #s
#b) the ID number matrix (length = 7 for Yeo7)
#c) the names of the networks (length = 7 for Yeo7)
#3) this is optional second mask, I am using glycolytic index so I don't have to rewrite the script

#post: 
#1) a vector corresponding to the parcel value for the regions you want to display
#2) Optional - a vector corresponding to second map

#uses: 
#1) Will take the input vector, convert all non-zeros to 1s. 
#will also do this for positive and negative only, for better visualization
#2) will multiple this vector with the parcel assignment (or second mask)
#3) Will save output in chosen directory

#dependencies
#any R should do
#I am using PMACS, and R 3.2.5

### get arguments if needed
#args = commandArgs(trailingOnly=TRUE)

#### set # parcels in case I want to do 7 or 17 or something else in the future
parcel_type = "Yeo" 
parcel_num = 7 
input_parcel_array_length = 10242

##### Alternative second mask
make_second_mask_flag = TRUE #i.e. I want to make an additional mask(s)
second_mask_path = "/baller/processed_data/zaixu_maps/fsaverage5/"
masks = c("GI_fsaverage5_10242.csv", "AllometricScaling_fsaverage5.csv", "CMRGlu_fsaverage5.csv", "Hill2010_evo_fsaverage5.csv","MeanCBF.fsaverage5.csv")
mask_length = 10242

## set abs and relative paths. Must toggle before running locally/on cluster
#homedir <- "/Users/eballer/BBL/imco/pmacs/PMACS_remote"
homedir <- "/project/imco"


#will loop through each of these
analyses <- c("coupling_accuracy", "cbf_accuracy", "alff_accuracy") 
models <- c("gam_age", "lm_age", "lm_sex", "lm_accuracy", "lm_exec_accuracy")
corrs <- c("uncor", "fdr05") #correction



for (analysis in analyses) {
  for (model in models) {
    for (corr in corrs) {  
      
      ### set results path
      stat_path <- paste0("/", analysis, "/")
      print(stat_path)
      result_path <- paste0(model, "_t_", corr)
      print(result_path)
      
      ### set paths
      ## input
      lh_stat_map <- read.csv(paste0(homedir, "/baller/results/", stat_path, "lh_", result_path, ".csv"), header = F)
      rh_stat_map <- read.csv(paste0(homedir, "/baller/results/", stat_path, "rh_", result_path, ".csv"), header = F)
      
      parcelID <- read.csv(paste0(homedir, "/baller/processed_data/yeo_network_data/NetworkIDnumbers", parcel_type, parcel_num, ".csv"), header = F)
      parcelName <- t(read.csv(paste0(homedir, "/baller/processed_data/yeo_network_data/NetworkNames", parcel_type, parcel_num, ".csv"), header = F))
      
      lh_parcel_nums <- read.csv(paste0(homedir, "/baller/processed_data/yeo_network_data/lh_", input_parcel_array_length, "_vertex_nums_", parcel_type, parcel_num, ".csv"), header = F)
      rh_parcel_nums <- read.csv(paste0(homedir, "/baller/processed_data/yeo_network_data/rh_", input_parcel_array_length, "_vertex_nums_", parcel_type, parcel_num, ".csv"), header = F)
      
      ## output
      lh_outdir <- paste0(homedir, "/baller/results/", stat_path, "lh_", result_path, "_", parcel_type, parcel_num, ".csv")
      rh_outdir <- paste0(homedir, "/baller/results/", stat_path, "rh_", result_path, "_", parcel_type, parcel_num, ".csv")
      
      lh_outdir_pos <- paste0(homedir, "/baller/results/", stat_path, "lh_pos_", result_path, "_", parcel_type, parcel_num, ".csv")
      rh_outdir_pos <- paste0(homedir, "/baller/results/", stat_path, "rh_pos_", result_path, "_", parcel_type, parcel_num, ".csv")
      
      lh_outdir_neg <- paste0(homedir, "/baller/results/", stat_path, "lh_neg_", result_path, "_", parcel_type, parcel_num, ".csv")
      rh_outdir_neg <- paste0(homedir, "/baller/results/", stat_path, "rh_neg_", result_path, "_", parcel_type, parcel_num, ".csv")
      
      ### print output so you know what is going on
      
      print(paste0("Converting ", homedir, "/baller/results/", stat_path, "lh_", result_path, ".csv to ", parcel_type, parcel_num ))
      print("Of note, these are the parcel names and associated numbers")
      refnum_networknum_name <- cbind(parcelName, parcelID, 1:dim(parcelID)[1])
      names(refnum_networknum_name) <- c("Network_name", "Net_number", "Mapping_for_matlab_PBP")
      print(refnum_networknum_name)
      
      
      #convert stat map to boolean
      lh_stat_boolean <- ifelse(lh_stat_map == 0, 0, 1)
      rh_stat_boolean <- ifelse(rh_stat_map == 0, 0, 1)
      
      #positive and negative vectors
      lh_stat_boolean_pos <- ifelse(lh_stat_map > 0, 1, 0)
      rh_stat_boolean_pos <- ifelse(rh_stat_map > 0, 1, 0)
      
      lh_stat_boolean_neg <- ifelse(lh_stat_map < 0, 1, 0)
      rh_stat_boolean_neg <- ifelse(rh_stat_map < 0, 1, 0)
      
      #make a column of numbers for mapping
      parcelID$network_num <- c(1:dim(parcelID)[1])
      
      #add extra row to parcelID, not clear why this didn't come from Yeo labels, maybe cerebellum?... 8 will equal 65793
      # comment this out if not using yeo 
      parcelID<- rbind(parcelID, c(65793, 8))
      
      #make vector for lh and rh with mapping
      lh_numerical_map <- lh_parcel_nums
      rh_numerical_map <- rh_parcel_nums
      
      #foreach vertex, which contains a bunch of numbers, match it to the appropriate column, and take the network num (i.e. yeo 2, which would correspond to Motor), associated with it
      lh_numerical_map[] <- lapply(lh_parcel_nums, function(x) parcelID$network_num[match(x, parcelID$V1)])
      rh_numerical_map[] <- lapply(rh_parcel_nums, function(x) parcelID$network_num[match(x, parcelID$V1)])
      
      #multiply
      
      lh_stat_booleanxnetwork <- lh_stat_boolean * lh_numerical_map
      rh_stat_booleanxnetwork <- rh_stat_boolean * rh_numerical_map 
      
      lh_stat_booleanxnetwork_pos <- lh_stat_boolean_pos * lh_numerical_map
      rh_stat_booleanxnetwork_pos <- rh_stat_boolean_pos * rh_numerical_map 
      
      lh_stat_booleanxnetwork_neg <- lh_stat_boolean_neg * lh_numerical_map
      rh_stat_booleanxnetwork_neg <- rh_stat_boolean_neg * rh_numerical_map 
      
      
      #write output
      write.table(x = lh_stat_booleanxnetwork, file = lh_outdir, quote = F, row.names = F, col.names = F)
      write.table(x = rh_stat_booleanxnetwork, file = rh_outdir, quote = F, row.names = F, col.names = F)
      
      #write output
      write.table(x = lh_stat_booleanxnetwork_pos, file = lh_outdir_pos, quote = F, row.names = F, col.names = F)
      write.table(x = rh_stat_booleanxnetwork_pos, file = rh_outdir_pos, quote = F, row.names = F, col.names = F)
      
      #write output
      write.table(x = lh_stat_booleanxnetwork_neg, file = lh_outdir_neg, quote = F, row.names = F, col.names = F)
      write.table(x = rh_stat_booleanxnetwork_neg, file = rh_outdir_neg, quote = F, row.names = F, col.names = F)
      
      
      ##############################################
      ######### Optional Second Mask Code ##########
      ##############################################
      
      ##### Alternative second mask
      
      if (make_second_mask_flag == TRUE){
        for (mask in masks) {
          lh_mask_nums <- read.csv(paste0(homedir, second_mask_path, "lh.", mask), header = F)
          rh_mask_nums <- read.csv(paste0(homedir, second_mask_path, "rh.", mask), header = F)
          print(paste0(homedir, second_mask_path, "lh.", mask), header = F)
          print(paste0(homedir, second_mask_path, "rh.", mask), header = F)
        
          
          #output
          lh_outdir <- paste0(homedir, "/baller/results/", stat_path, "lh_", result_path, "_", mask, "_", mask_length, ".csv")
          rh_outdir <- paste0(homedir, "/baller/results/", stat_path, "rh_", result_path, "_", mask, "_", mask_length, ".csv")
          
          lh_outdir_pos <- paste0(homedir, "/baller/results/", stat_path, "lh_pos_", result_path, "_", mask, "_", mask_length, ".csv")
          rh_outdir_pos <- paste0(homedir, "/baller/results/", stat_path, "rh_pos_", result_path, "_", mask, "_", mask_length, ".csv")
          
          lh_outdir_neg <- paste0(homedir, "/baller/results/", stat_path, "lh_neg_", result_path, "_", mask, "_", mask_length, ".csv")
          rh_outdir_neg <- paste0(homedir, "/baller/results/", stat_path, "rh_neg_", result_path, "_", mask, "_", mask_length, ".csv")
          
          ## multiply
          lh_stat_booleanxnetwork <- lh_stat_boolean * lh_mask_nums
          rh_stat_booleanxnetwork <- rh_stat_boolean * rh_mask_nums
          
          lh_stat_booleanxnetwork_pos <- lh_stat_boolean_pos * lh_mask_nums
          rh_stat_booleanxnetwork_pos <- rh_stat_boolean_pos * rh_mask_nums 
          
          lh_stat_booleanxnetwork_neg <- lh_stat_boolean_neg * lh_mask_nums
          rh_stat_booleanxnetwork_neg <- rh_stat_boolean_neg * rh_mask_nums
          
          #write output
          write.table(x = lh_stat_booleanxnetwork, file = lh_outdir, quote = F, row.names = F, col.names = F)
          write.table(x = rh_stat_booleanxnetwork, file = rh_outdir, quote = F, row.names = F, col.names = F)
          
          #write output
          write.table(x = lh_stat_booleanxnetwork_pos, file = lh_outdir_pos, quote = F, row.names = F, col.names = F)
          write.table(x = rh_stat_booleanxnetwork_pos, file = rh_outdir_pos, quote = F, row.names = F, col.names = F)
          
          #write output
          write.table(x = lh_stat_booleanxnetwork_neg, file = lh_outdir_neg, quote = F, row.names = F, col.names = F)
          write.table(x = rh_stat_booleanxnetwork_neg, file = rh_outdir_neg, quote = F, row.names = F, col.names = F)
        }
        
      }
    }
  }
}

