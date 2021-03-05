require('visreg')
require('mgcv')
require('tableone')
require('dplyr')

make_demographics_table<- function(data_frame) {
  
  #subset demographics
  listVars <- c("Age", "Sex", "Race", "Maternal Ed") #Race 1 = caucasian, Maternal Ed = years, age = years
  demo <- data.frame(data_frame$ageAtScan1, data_frame$sex, data_frame$race2, data_frame$medu1)
  names(demo) <- c(listVars)
  
  #Change categorical values to have names
  demo$Race <- ifelse(demo$Race == 1, "Caucasian", "Non-caucasian")
  demo$Sex <- ifelse(demo$Sex == 1, "Male", "Female")
 
  #Define Categorical Variables
  cat_variables <- c("Sex", "Race")
  title <- paste0("IMCO Demographics, n = ", dim(demo)[1])
  
  #create demographics table
  demo_table <- CreateTableOne(vars = listVars, data = demo, factorVars = cat_variables)
  print(demo_table, showAllLevels = TRUE)
}

get_parcel_mapping_yeo <- function(parcel_num){
  
  #pre: input parcel #, either 7 or 17
  #post: list lh and rh yeo networks that map onto code #s
  #uses: easy way to translate the weird numerical maps in fsaverage 5 space into something we are more familiar with
  #dependencies: Any R will do, I used 3.2.5
  
  ## Set Yeo info
  #### set # parcels in case I want to do 7 or 17 or something else in the future
  parcel_type = "Yeo" 
  parcel_num = parcel_num 
  input_parcel_array_length = 10242
  
  # read in yeo fsaverage5 vectors
  parcelID <- read.csv(paste0(homedir, "/baller/processed_data/yeo_network_data/NetworkIDnumbers", parcel_type, parcel_num, ".csv"), header = F)
  parcelName <- t(read.csv(paste0(homedir, "/baller/processed_data/yeo_network_data/NetworkNames", parcel_type, parcel_num, ".csv"), header = F))
  
  lh_parcel_nums <- read.csv(paste0(homedir, "/baller/processed_data/yeo_network_data/lh_", input_parcel_array_length, "_vertex_nums_", parcel_type, parcel_num, ".csv"), header = F)
  rh_parcel_nums <- read.csv(paste0(homedir, "/baller/processed_data/yeo_network_data/rh_", input_parcel_array_length, "_vertex_nums_", parcel_type, parcel_num, ".csv"), header = F)
  
  # map Yeo numbers to parcels
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

  lh_and_rh_numerical_map_list <- list(lh_numerical_map$V1, rh_numerical_map$V1)
  return(lh_and_rh_numerical_map)
  
}