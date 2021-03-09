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


library(bigmemory.sri)
library(ggplot2)
library(tidyr)



#################
### set home directory
homedir <- "/Users/eballer/BBL/imco/pmacs/PMACS_remote"
#homedir <- "/project/imco"

source(paste0(homedir, "/baller/scripts/imco_functions.R"))

#initialize
hemi <- c("lh", "rh")
permNum <- 1000
yeo_num <- 7

################
### Read in matrices 

lh_spin <- read.csv("lh.csv")
rh_spin <- read.csv("rh.csv")

#grab list of yeo 7 networks in fsaverage5 space
yeo_networks <- get_parcel_mapping_yeo(yeo_num)

#separate into right and left
lh_yeo_network <- yeo_networks[[1]]
rh_yeo_network <- yeo_networks[[2]]

#count up number of vertices per network
lh_yeo_network_count_table <- table(lh_yeo_network)
rh_yeo_network_count_table <- table(rh_yeo_network)

#multiply yeo network x spin test
lh_spinxyeo <- lh_spin*lh_yeo_network
rh_spinxyeo <- rh_spin*rh_yeo_network

#pseudocode
#go through each hemisphere, go through each perm, and go through each network
lh_hemi_spin_proportions <- matrix(nrow = yeo_num, ncol = (permNum + 1))
rh_hemi_spin_proportions <- matrix(nrow = yeo_num, ncol = (permNum + 1))
for (hemi in hemis){
  for (perm in 1:permNum) {
    for (network in 1:yeo_num){
      #to evaluate
      eval_to_parse<- paste0("length(which(", hemi, "_spinxyeo == ", network, "\"")
      num_vertices_in_spin <- eval(parse(text = as.character(eval_to_parse)))
      num_neg_to_parse <- paste0("length(which(", hemi, "_spinxyeo == -", network, "\"")
      num_neg <- eval(parse(text = as.character(num_neg_to_parse)))
      total_possible <- eval(parse(text = as.character(paste0(hemi, "yeo_network_count_table[", yeo_num, "]"))))
      proportion_potential_vertices <- num_vertices_in_spin/(total_possible - num_neg)
      eval(parse(text = as.character(paste0(hemi, "_hemi_spin_proportions[", network, "][", perm, "] = ", proportion_potential_vertices))))
    }
  }
}

#then plot