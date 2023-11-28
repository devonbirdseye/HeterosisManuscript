set.seed(586)

#Path to input data
PathToData <- "Data/"
#Path to output plots
PathToPlots <- "Plots/"
#Path to output tables
PathToTables <- "Tables/"
#Path to output datasets
PathToDatasets <- "Datasets/"

source("Scripts/Functions.R")

source("Scripts/01_ImportData.R")

source("Scripts/02_ProcessProteinData.R")
source("Scripts/03_ProcessRNAdata.R")

source("Scripts/PlotFigures.R")
