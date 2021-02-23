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

