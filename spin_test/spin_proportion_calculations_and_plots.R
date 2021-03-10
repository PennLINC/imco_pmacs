##############################################
########### Spin Test Distribution for #######
##########Intermodal Coupling Paper ##########

####### Author: Erica Baller
#### Date: 3/9/2021

#######
##pre: right and left 10242 x 1000 matrices from matlab SpinPermuFS, yeo R & L assignments
##post: 2 7 x 1000 matrices (r & l) that contain the proportion of vertices within a network divided by the total number of vertices, and plots
## uses: Takes output of spin test, and calcualted the number of vertices within each of yeo's 7 networks out of the number of total possible vertices within the network
    #### 1) Read in the yeo network assignments and calculate total number of vertices per network
    #### 2) Multiply the yeo networks x the matrices (so every value is 1 -7 if they were within the mask, -1--7 if they were medial wall, and 0 otherwise)
    #### 3) Foreach permutation (r and l separately), and for each network, calculate the (# of vertices with a 1) divided(/) by the (number of total vertices within network minus number of negative vertices
    #### 4) Store
    #### 5) Plot

### dependencies: ggplot2, bigmemory, vroom


#library(bigmemory.sri)
library(ggplot2)
library(tidyr)



#################
### set home directory
#homedir <- "/Users/eballer/BBL/imco/pmacs/PMACS_remote"
homedir <- "/project/imco"

source(paste0(homedir, "/baller/scripts/imco_functions.R"))

#initialize
hemis <- c("lh", "rh")
permNum <- 1000
yeo_num <- 7
model = "gam_sex"
################
### Read in matrices 

#{lh and rh}_gam_sex_t_fdr05 -> actual results
lh_gam_sex_t_fdr05_results <- read.table(paste0(homedir, "/baller/results/coupling_accuracy/lh_", model, "_t_fdr05_Yeo7_1_0_-1.csv"))
rh_gam_sex_t_fdr05_results <- read.table(paste0(homedir, "/baller/results/coupling_accuracy/rh_", model, "_t_fdr05_Yeo7_1_0_-1.csv"))

#spins
lh_spin <- t(read.table(paste0(homedir, "/baller/results/coupling_accuracy/spin_test_results/lh_spintest_", model,"_output.csv"), sep = ","))
rh_spin <- t(read.table(paste0(homedir, "/baller/results/coupling_accuracy/spin_test_results/rh_spintest_", model,"_output.csv"), sep = ","))
                                        
#bring together, with original values as first column
lh_act_results_and_spin <- cbind(lh_gam_sex_t_fdr05_results, lh_spin)
rh_act_results_and_spin <- cbind(rh_gam_sex_t_fdr05_results, rh_spin)

#grab list of yeo 7 networks in fsaverage5 space
yeo_networks <- get_parcel_mapping_yeo(yeo_num)

#separate into right and left
lh_yeo_network <- yeo_networks[[1]]
rh_yeo_network <- yeo_networks[[2]]

#count up number of vertices per network
lh_yeo_network_count_table <- table(lh_yeo_network)
rh_yeo_network_count_table <- table(rh_yeo_network)

#multiply yeo network x spin test
lh_spinxyeo <- lh_act_results_and_spin*lh_yeo_network
rh_spinxyeo <- rh_act_results_and_spin*rh_yeo_network

#proportions
#go through each hemisphere, go through each perm, and go through each network

lh_hemi_spin_proportions <- data.frame(matrix(nrow = yeo_num, ncol = (permNum + 1)))
rh_hemi_spin_proportions <- data.frame(matrix(nrow = yeo_num, ncol = (permNum + 1)))
for (hemi in hemis){
 # print(paste0("hemi", hemi))
  for (perm in 1:(permNum + 1)){
    #print(paste0("perm: ", perm))
    for (network in 1:yeo_num){
      #print(paste0("network: ", network))
      #to evaluate
      
      #number of vertices within network that are fdr corrected
      num_pos_to_parse<- paste0("length(which(", hemi, "_spinxyeo[", perm, "] == ", network, "))")
      #print(num_pos_to_parse)
      num_vertices_in_spin <- eval(parse(text = as.character(num_pos_to_parse)))
      #print(paste0("num vertices in spin", num_vertices_in_spin))
      
      #number of vertices within network that are negative (i.e., medial wall)
      num_neg_to_parse <- paste0("length(which(", hemi, "_spinxyeo[", perm, "] == -", network, "))")
     # print(num_neg_to_parse)
      num_neg <- eval(parse(text = as.character(num_neg_to_parse)))
      #print(paste0("num neg: ", num_neg))
      
      #total number of vertices in normal network
      total_possible_to_parse <- paste0(hemi, "_yeo_network_count_table[", network, "]")
      #print(total_possible_to_parse)
      total_possible <- eval(parse(text = as.character(total_possible_to_parse)))
      
      #print(paste0("total possible:", total_possible))
      #proportion of vertices within network , with denominator being total possible by # in medial wall
      proportion_potential_vertices <- num_vertices_in_spin/(total_possible - num_neg)
      
      #print(paste0("prop:", proportion_potential_vertices))
      
      #store in matrix
      storing_to_parse <- paste0(hemi, "_hemi_spin_proportions[", network, ",", perm, "] = ", proportion_potential_vertices)
      #print(storing_to_parse)
      eval(parse(text = as.character(storing_to_parse)))
    }
  }
}

write.table(lh_hemi_spin_proportions, file = paste0(homedir, "/baller/results/coupling_accuracy//spin_test_results/lh_spin_test_", model, "_proportions.csv"), sep = ",", col.names = F, row.names = F)
write.table(rh_hemi_spin_proportions, file = paste0(homedir, "/baller/results/coupling_accuracy//spin_test_results/rh_spin_test_", model, "_proportions.csv"), sep = ",", col.names = F, row.names = F)
#then plot
