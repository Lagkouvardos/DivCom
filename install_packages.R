#' This script installs all required libraries automatically
#' Please install libraries manually if it was not possible to install an library automatically
#' Missing libraries are listed in missing_packages.txt
#' To install an library please use following two command:
#' install.packages("name of the missing library")
#' library("name of the missing library)

##################################################################################
######             Set parameters in this section manually                  ######
##################################################################################

#' Please set the directory of the present script as the working folder
#' Note: the path is denoted by forward slash "/"
setwd("D:/path/to/DivCom")


######                  NO CHANGES ARE NEEDED BELOW THIS LINE               ######
##################################################################################
######                             Main Script                              ###### 
##################################################################################
###################       Load all required libraries     ########################


# Check if required packages are already installed, and install if missing
packages <-c("ade4","ape","caTools","cluster","cowplot","data.table","dplyr","factoextra",
             "fpc","ggplot2","ggpubr","ggtree","graphics","grid","gtable","gridExtra","lattice","mclust",
             "GUniFrac","permute","phangorn","RColorBrewer","stats","tidyr","tools","vegan") 

# Function to check whether the package is installed
InsPack <- function(pack)
{
  if ((pack %in% installed.packages()) == FALSE) {
    install.packages(pack,repos = "http://cloud.r-project.org/",dependencies = TRUE)
  } 
}

if (("ade4" %in% installed.packages()) == FALSE) {
  install.packages("ade4",repos ="http://cloud.r-project.org/",dependencies = TRUE)
} 

# Applying the installation on the list of packages
lapply(packages, InsPack)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (("ggtree" %in% installed.packages()) == FALSE) {
BiocManager::install("ggtree",update=FALSE,force = TRUE)
}

# Make the libraries
lib <- lapply(packages, require, character.only = TRUE)
not_installed <- which(lib==FALSE)
missing_packages <- lapply(not_installed, function(x) print(packages[x]))

if (length(not_installed)==0) {
  # Adding log file in analysis
  sink(file = "missing_packages.txt")
  cat ("**********************************************************","\n")
  cat ("All the required packages have been succesfully installed.","\n")
  cat ("**********************************************************","\n")
  sink()
} else {
# Adding log file in analysis
sink(file = "missing_packages.txt")
cat ("******************************************","\n")
cat ("Please install following packages manually","\n")
cat (as.character(missing_packages))
sink()
}


#################################################################################
######                           End of Script                             ######
#################################################################################                           
