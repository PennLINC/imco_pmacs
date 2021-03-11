#########################################
## Correlations of coupling with masks ##
#########################################

### Author: Erica Baller
## Date: 2/26/2021

### This script emerged out of the desire to take vertex output and display it within different networks/parcels

#pre: input: 
#1) mean coupling map, r and l side
#2) masks, GI, Hill, CBF, etc

#post - matrix of correlations for each of the maps

#uses - We want to see which of the masks best correlates with mean coupling. To do this, will correlate the uncorrected mean map with the GI map, Hill2010, and all the others

#dependencies: Any R will do, I used 3.2.5

# set home directory
homedir = '/Users/eballer/BBL/imco/pmacs/PMACS_remote/'
#homedir = '/project/imco'

#mean coupling maps
mean_coupling_map <- rbind(read.csv(paste0(homedir, "/baller/results/mean_maps/n831_lh.coupling_coef_alff_cbf.fwhm15.fsaverage5.csv"), header = F),
                           read.csv(paste0(homedir, "/baller/results/mean_maps/n831_rh.coupling_coef_alff_cbf.fwhm15.fsaverage5.csv"), header = F))

#gi maps
gi_map <- rbind(read.csv(paste0(homedir, "/baller/processed_data/zaixu_maps/fsaverage5/lh.GI_fsaverage5_10242.csv"), header = F),
                read.csv(paste0(homedir, "/baller/processed_data/zaixu_maps/fsaverage5/rh.GI_fsaverage5_10242.csv"), header = F))

#hill2010
hill_map <- rbind(read.csv(paste0(homedir, "/baller/processed_data/zaixu_maps/fsaverage5/lh.Hill2010_evo_fsaverage5.csv"), header = F),
                  read.csv(paste0(homedir, "/baller/processed_data/zaixu_maps/fsaverage5/rh.Hill2010_evo_fsaverage5.csv"), header = F))

#allometric scaling
as_map <- rbind(read.csv(paste0(homedir, "/baller/processed_data/zaixu_maps/fsaverage5/lh.AllometricScaling_fsaverage5.csv"), header = F),
                read.csv(paste0(homedir, "/baller/processed_data/zaixu_maps/fsaverage5/rh.AllometricScaling_fsaverage5.csv"), header = F))

#cmrglu
cmrglu_map <- rbind(read.csv(paste0(homedir, "/baller/processed_data/zaixu_maps/fsaverage5/lh.CMRGlu_fsaverage5.csv"), header = F),
                    read.csv(paste0(homedir, "/baller/processed_data/zaixu_maps/fsaverage5/rh.CMRGlu_fsaverage5.csv"), header = F))

#meancbf
meancbf_map <- rbind(read.csv(paste0(homedir, "/baller/processed_data/zaixu_maps/fsaverage5/lh.MeanCBF.fsaverage5.csv"), header = F),
                     read.csv(paste0(homedir, "/baller/processed_data/zaixu_maps/fsaverage5/rh.MeanCBF.fsaverage5.csv"), header = F))

maps <- c('mean_coupling_map','gi_map', 'hill_map', 'as_map', 'cmrglu_map', 'meancbf_map')

correlations <- data.frame(matrix(0, nrow = 1, ncol = 6))
names(correlations) <- maps
row.names(correlations) <- "mean_coupling_map"

i = 1
for (map in maps) {
    map_to_corr <- map
    corr_cmd <- paste0("cor(x = mean_coupling_map, ", "y =", map_to_corr, ")")
    corr_results <- eval(parse(text = as.character(corr_cmd)))
    print(paste0("mean_coupling_map:", map_to_corr, "--> R = ", corr_results))
    correlations[1,i] = corr_results
    i <- i + 1

}

write.csv(correlations, file = paste0(homedir, "/baller/results/mean_coupling_x_mask_correlations/mean_coupling_map_x_mask_correlations.csv"), quote = F)
