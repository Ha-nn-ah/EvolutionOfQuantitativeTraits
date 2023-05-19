## script for merging the simulation output of ISM.cpp

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

# generation of simulation
Dgen <- 1

# import files
MEAN <- read.table(file = "mean.tsv",  sep = '\t', header = FALSE)
VAR <- read.table(file = "variance.tsv",  sep = '\t', header = FALSE)
MUT <- read.table(file = "mutants.tsv", sep = '\t', header = FALSE)
FIX <- read.table(file = "fixed.tsv", sep = '\t', header = FALSE)
LOSS <- read.table(file = "extinct.tsv", sep = '\t', header = FALSE)

# create dataframe
gen <- c()
mean_mean <- c()
mean_var <- c()
variance_mean <- c()
variance_var <- c()
mutloci <- c()
fixloci <- c()
lossloci <- c()
segloci <- c()
genNew <- 0
for (i in 1:(length(MEAN)-1)){
  genNew <- genNew + Dgen
  gen <- c(gen, genNew)
  mean_mean <- c(mean_mean, mean(MEAN[,i]))
  mean_var <- c(mean_var, var(MEAN[,i]))
  variance_mean <- c(variance_mean, mean(VAR[,i]))
  variance_var <- c(variance_var, var(VAR[,i]))
  mutloci <- c(mutloci, mean(MUT[,i]))
  fixloci <- c(fixloci, mean(FIX[,i]))
  lossloci <- c(lossloci, mean(LOSS[,i]))
  segloci <- c(segloci, mean(MUT[,i]) - (mean(FIX[,i]) + mean(LOSS[,i])))
}

data <- data.frame(gen, mean_mean, variance_mean, mutloci, fixloci, lossloci, segloci, mean_var, variance_var)
  
# export file
write.table(data, file = "simulation.tsv", quote=FALSE, sep='\t', row.names = FALSE)

