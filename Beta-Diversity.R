
#' Script: Beta-Diversity
#' #'This script was last modified on 05/01/2022

#' Author: Ilias Lagkouvardos
#' 
#'##############################################################################################
#'######################### Beta-Diversity Script Overview #####################################
#'############################################################################################## 
#'
#' 1.Calculate beta-diversity for microbial communities
#' based on permutational multivariate analysis of variances (PERMANOVA) using multiple distance matrices
#' computed from phylogenetic distances between observed organisms
#' 2.Produce the plots of three different validation indices
#'
#' Input:
#' 1. Set the path to the directory where the file is stored 
#' 2. Write the name of the OTU table without taxonomy information
#' 3. Write if the OTUs table is normalized or not
#' 4. Write the name of the mapping file that includes the samples groups
#' 5. Write if you will provide distance matrix or phylogenetic tree
#' 6. Write the name of the OTU tree or the phylogenetic tree
#' 7. Write the name of the mapping file of the variable (sample group) used for comparison
#' 8. Write the name of groups that will be used for as reference groups
#' 9. Write the name of groups that will be used for as reference groups
#' 
#' 
#'
#' Output: 
#' The script generates three graphical outputs (pdf), one text file and a newick tree
#' 1. A phylogram with colour-coded group clustering
#' 2. MDS and NMDS plots showing information about beta-diversity across all sample groups
#' 3. MDS and NMDS plots of all pairwise comparisons
#' 4. The distance matrix
#' 5. Plots showing the optimal number of clusters  
#' 6. PDFs that suggest the optimal number of clusters for every group
#' 7. Dendogram for all samples in a newick tree file
#'
#' Concept:
#' A distance matrix is calculated based on the generalized UniFrac approach
#' (Chen J, et al. Associating microbiome composition with environmental covariates using generalized UniFrac distances. 2012)
#' Samples are clustered based on the distance matrix using the Ward's hierarchical clustering method
#' To determine similarities between samples, a multivariate analysis is applied
#' and sample distribution is illustrated by means of MDS and NMDS (non-metric) plots
#' The Calinski-Harabasz (CH), silhouette Index, Within Sum of Squares (WSS) plot, and Prediction Strength plot 
#' are used to assess the optimal number of clusters the dataset was most robustly partitioned into  


##################################################################################
######             Set parameters in this section manually                  ######
##################################################################################

#' Please set the directory of the script as the working folder (e.g D:/studyname/NGS-Data/Rhea/beta-diversity/)
#' Note: the path is denoted by forward slash "/"
setwd("D:/path/to/DivCom/")     #<--- CHANGE ACCORDINGLY !!!

#' Please give the name of the OTUs or ASVs table
input_otu = "OTUs-Table.tab"              #<--- CHANGE ACCORDINGLY !!!


#' Please insert "YES" if the OTUs/ASVs table is normalized otherwise, insert "NO"
#' In case the table is not normalized, the OTUs/ASVs table will be normalized based on the minimum sum of reads of the samples.
normalized = "NO"

#' Choose if you will insert a distances matrix or the phylogenetic tree of the OTUs/ASVs sequences
#' There are two options: "distances matrix" or "tree"
#' 1) Please insert "distances matrix" if you will provide a distances matrix
#' 2) Please insert "tree" if you will provide a phylogenetic tree (In this case, the Generalized unifrac distances will be calculated)
tree_or_matrix = "tree"


#' Please insert the name of the distances matrix or the anme of the phylogenetic tree 
#' -> -> !!! In case you will choose the "mean" or "median" option you HAVE TO provide a phylogenetic tree !!! <- <-
input_tree_or_matrix = "OTUs-NJTree.tre" 


#' Please give the name of the mapping file which contains the labels of the samples
#' !! CAUTION: The rows of the mapping file should have the same sample names as the OTUs table !!
input_meta = "mapping_file.tab"              #<--- CHANGE ACCORDINGLY !!!


#' Please provide the name or the number of the column (of the mapping file) based on which the samples will be partitioned into groups
mapping_column = "Condition"                    #<--- CHANGE ACCORDINGLY !!!

#' Please place in the vector one or more names which will be used to identify the samples that composing 
#'         the REFERENCE group (e.g reference_name <- c("group_a","group_b"))
#' -> -> !!! CAUTION: You should provide at least one name!!! <- <-
Reference_name = c("NI")

#' Please provide the names of the test groups 
#' There are two options: a User-defined vector or "None"
#' 1) User-defined vector --> Form a vector with one or more elements referring to the name of the groups (e.g test_name <- c("group_a","group_b"))
#' 2)     c()             --> In case you insert an empty vector,there won't be any test group. 
#'                            Only the indeces of the reference samples will be calculated (Not recommended)
Test_name =c("IBD")




#-------------------- Additional parameters ----------------------#


#' Turn on sample labeling
#' 0 = Samples are not labeled in the MDS/NMDS plots
#' 1 = All Samples are labeled in the MDS/NMDS plots
label_samples = 0

#' Determine which sample labeled should appear
#' Write the name of samples (in quotation marks), which should appear in the MDS/NMDS plots, in the vector (c) below
#' If more than one sample should be plotted, please separate their IDs by comma (e.g. c("sample1","sample2"))
label_id =c("")

#' De-Novo Clustering will be performed for the number of samples or maximal for the set limit
#' Default Limit is 10
kmers_limit=10



#################################################################################################################################################################
################################################# DO NOT CHANGE ANYTHING BELOW THIS LINE #################################################################################
#################################################################################################################################################################


############################################################################################################################################################
############################################################ Main Script ###################################################################################
############################################################################################################################################################

# Save the variables of the active working environment
ls <- c("meta_file","otu_file","unifract_dist_all",ls())


###################       Load all required libraries     ########################



# Check if required packages are already installed, and install if missing
packages <-c("ade4","GUniFrac","phangorn","factoextra","cluster","fpc","dplyr","graphics","vegan","stats","data.table","caTools","gridExtra","grid","gtable","tools") 

# Function to check whether the package is installed
InsPack <- function(pack)
{
  if ((pack %in% installed.packages()) == FALSE) {
    install.packages(pack,repos ="http://cloud.r-project.org/",dependencies = TRUE)
  } 
}

if (("ade4" %in% installed.packages()) == FALSE) {
  install.packages("ade4",repos ="http://cloud.r-project.org/",dependencies = TRUE)
} 

# Applying the installation on the list of packages
lapply(packages, InsPack)

# Make the libraries
lib <- lapply(packages, require, character.only = TRUE)

# Check if it was possible to install all required libraries
flag <- all(as.logical(lib))


###################################################################################
#################            Necessary functions                  #################
###################################################################################

#------------------------ Function that generates MDS plots-----------------------#

# distance -> the given distance matrix 
#   groups -> a vector that contains the indexes of the groups
#     kp   -> Number of clusters 

sclass <- function(dist,all_groups, kp) {
  
  
  mds<- cmdscale(dist,eig=T, x.ret=T)
  mds.variation<- round(mds$eig/sum(mds$eig)*100,1) # axes variance calculation
  plot_color<-rainbow(length(levels(all_groups)))[all_groups]
  all_groups_comp <- all_groups[!is.na(all_groups)]
  all_groups_comp<-factor(all_groups_comp,levels(all_groups_comp)[unique(all_groups_comp)])
  s.class(
    mds$points, col = unique(plot_color), cpoint =
      2, fac = all_groups_comp )
  graphics:: title (main=paste0("MDS graph k= ",kp ), xlab = paste("MDS1:",mds.variation[1],"%"), ylab=paste("MDS2:",mds.variation[2],"%"))
  
}

##################  Create all the necessary paths and directories ###############

# Save the current path in one variable (will be used later)
OriginalPath <- getwd()

groups_name <- c(Reference_name,Test_name)


# Create the directory where all output files will be saved 
results_folder <- paste0( "Beta-Diversity_",mapping_column,"(", paste(groups_name,collapse="-"),")")
results_path <- paste0(OriginalPath,"/",results_folder)
dir.create(results_folder)

# Create the directory where all output files (optimal number of clusters) will be saved
optimal_number <- "Optimal number of Clusters"
optimal_number_path <- paste0(results_path,"/",optimal_number)
dir.create(paste0(results_path,"/",optimal_number))

# Create the directory where all output files (Beta-diversity) will be saved
beta_diversity <- "Beta Diversity"
beta_diversity_path <- paste0(results_path,"/",beta_diversity)
dir.create(paste0(results_path,"/",beta_diversity))



##################################################################################
##############          Read all required input files          ###################
##################################################################################

# Check the tree_or_matrix and input_tree_or_matrix are consistent to each other
if ( file_ext(input_tree_or_matrix)=="nwk" | file_ext(input_tree_or_matrix)=="tre"){tree_or_matrix = "tree"}

#------------------- Meta_file--------------------------#

# Load the mapping file containing individual sample information (sample names in the first column)
meta_file <- data.frame(read.table (file = input_meta, check.names = FALSE, header = TRUE, dec = ".", sep = "\t", row.names = 1, comment.char = ""))

if (ncol(meta_file)==0 | nrow(meta_file)==0){
  # Load the mapping file containing individual sample information (sample names in the first column)
  meta_file <- data.frame(read.table (file = input_meta, check.names = FALSE, header = TRUE, dec = ".", sep = ",", row.names = 1, comment.char = ""))
}


# Clean table from empty lines
meta_file <- data.frame(meta_file[!apply(is.na(meta_file) | meta_file=="",1,all),])

# Order meta file
if (ncol(meta_file)==1){
meta_file <- data.frame(meta_file[order(meta_file[,mapping_column]),, drop = FALSE])
colnames(meta_file) <- mapping_column
meta_file <- data.frame(meta_file[meta_file[,mapping_column] %in% groups_name,,drop=FALSE])
} else {meta_file <- data.frame(meta_file[order(meta_file[,mapping_column]),])
meta_file <- data.frame(meta_file[meta_file[,mapping_column] %in% groups_name,])
}


# Convert the selected column (mapping_column) of meta_file to factor
meta_file[,mapping_column] <- as.factor(meta_file[,mapping_column])


#-------------------------------------------------------#


#------------------------ OTUs table--------------------------#

# Load the tab-delimited file containing the values to be checked (row names in the first column)
otu_table <-  read.table (input_otu,check.names = FALSE,header = TRUE,dec = ".",sep = "\t", row.names = 1,comment.char = "")

if (ncol(otu_table)==0 | nrow(otu_table)==0){
  # Load the tab-delimited file containing the values to be checked (row names in the first column)
  otu_table <-  read.table (input_otu,check.names = FALSE,header = TRUE,dec = ".",sep = ",", row.names = 1,comment.char = "")
}


# Load the tab-delimited file containing the values to be checked (row names in the first column)
otu_table <-  read.table (input_otu,check.names = FALSE,header = TRUE,dec = ".",sep = "\t", row.names = 1,comment.char = "")

# Clean table from empty lines
otu_table <- otu_table[!apply(is.na(otu_table) | otu_table=="",1,all),, drop = FALSE]

# Check if exists a column with taxonomy information
taxonomy <- otu_table %>% select_if(is.factor)
if (ncol(taxonomy)==0) {
  taxonomy <-  otu_table %>% select_if(is.character)
}
# Delete the taxonomy column (if exists)
if (ncol(taxonomy)!=0) {
  otu_table[,colnames(taxonomy)] <- NULL
}

# Check if the mapping file and the OTUs table have the same sample names
{ if(all(rownames(meta_file) %in% rownames(otu_table))==FALSE & all(rownames(meta_file) %in% colnames(otu_table))==FALSE)
  stop ("Some rows of the mapping file are not present present in the OTUS file. \n Fix that problem and try again!")
  
  #--Normalize the OTUs table (if the user ask for)--#
  if (normalized=="YES") {
    
    # Check if Samples are in rows or columns
    if (any(colnames(otu_table) %in% rownames(meta_file))==TRUE) {
      
      # Keep only those rows that appear in the mapping file
      otu_table <- otu_table[,rownames(meta_file)]
      
      # Transpose OTU-table and convert format to a data frame
      otu_table<- data.frame(t(otu_table))
      
      # Create the otu_file
      otu_file <- otu_table[rownames(meta_file),] 
      
    } else {
      
      # Keep only those rows that appear in the mapping file
      otu_table <- otu_table[rownames(meta_file),]
      
      # Convert format to a data frame
      otu_table<- data.frame(otu_table)
      
      # Create the otu_file
      otu_file <- otu_table[rownames(meta_file),] 
    }
    
  } else {
    
    if (any(colnames(otu_table) %in% rownames(meta_file))==TRUE) {
      
      # Calculate the minimum sum of all columns/samples
      min_sum <- min(colSums(otu_table))
      
      # Divide each value by the sum of the sample and multiply by the minimal sample sum
      otu_file <- t(min_sum * t(otu_table) / colSums(otu_table))
      
      # Clean table from empty lines
      otu_file <- otu_file[!apply(is.na(otu_file) | otu_file =="",1,all),]
      
      # Transpose OTU-table and convert format to a data frame
      otu_file <- data.frame(t(otu_file))
      
      # Create the otu_file
      otu_file <- otu_file[rownames(meta_file),]
      
    } else {
      
      # Transpose the OTU-table
      otu_table <- data.frame(t(otu_table))
      
      # Calculate the minimum sum of all columns/samples
      min_sum <- min(colSums(otu_table))
      
      # Divide each value by the sum of the sample and multiply by the minimal sample sum
      otu_file <- t(min_sum * t(otu_table) / colSums(otu_table))
      
      # Clean table from empty lines
      otu_file <- otu_file[!apply(is.na(otu_file) | otu_file =="",1,all),]
      
      # Transpose the OTU-table and convert format to a data frame
      otu_file <- data.frame(t(otu_file))
      
      # Create the otu_file
      otu_file <- otu_file[rownames(meta_file),]
      
    }
  } 
#-------------------------------------------------------#



##############################################################################
############               Calculate beta-diversity               ############
##############################################################################
  
  if (tree_or_matrix=="tree"){ 
    
    #---------------- Tree -----------------------#
    # Load the phylogenetic tree calculated from the OTU sequences 
    tree_file <- read.tree(input_tree_or_matrix)
    
    # Remove single quotes from the tips of the tree
    tree_file$tip.label <- gsub("'", "", tree_file$tip.label)
    
    # Root the OTU tree at midpoint 
    rooted_tree <- midpoint(tree_file)
    #--------------------------------------------#
    
    # Calculate the UniFrac distance matrix for comparing microbial communities
    unifracs <- GUniFrac(otu_file, rooted_tree, alpha = c(0.0,0.5,1.0))$unifracs
    # Weight on abundant lineages so the distance is not dominated by highly abundant lineages with 0.5 having the best power
    unifract_dist_all <- unifracs[, , "d_0.5"]
  } else {
    # Read the distances matrix
    unifract_dist_all <- read.table (input_tree_or_matrix,check.names = FALSE,header = TRUE,dec = ".",sep = "\t", row.names = 1,comment.char = "")

  }

  
# OTUs table and mapping file transformation in order to include only the chosen by the user groups
otu_file <- otu_file[meta_file[,mapping_column] %in% groups_name,]

# Select the rows of the mapping file with the selected samples
meta_file <- data.frame(meta_file[meta_file[,mapping_column] %in% groups_name,,drop=FALSE])

# Calculate the UniFrac distance matrix only for the selected groups
unifract_dist <- unifract_dist_all[rownames(meta_file) , rownames(meta_file)]

# Add an extra column in the mapping file
meta_file[,(ncol(meta_file)+1)] <- rep(1,nrow(meta_file))
for (i in 1:nrow(meta_file)){
  if (meta_file[i,mapping_column] %in% Reference_name) {meta_file[i,ncol(meta_file)] <- paste(Reference_name,collapse="+")
  }else {meta_file[i,ncol(meta_file)] <-as.character(meta_file[i,mapping_column])}
} 

# Convert that column into factor
meta_file[,ncol(meta_file)] <- as.factor(meta_file[,ncol(meta_file)])

# Select metadata group based on the pre-set group name
all_groups <- as.factor(meta_file[,ncol(meta_file)])




################ Generate tree #######################

# Change working directory
setwd(beta_diversity_path)

# Save the UniFrac output as distance object
all_dist_matrix <- as.dist(unifract_dist)

# Apply a hierarchical cluster analysis on the distance matrix based on the Ward's method
all_fit <- hclust(all_dist_matrix, method = "ward.D2")

# Generates a tree from the hierarchically generated object
tree <- as.phylo(all_fit)
my_tree_file_name <- paste("phylogram","(",paste(groups_name,collapse="-"),").pdf",sep="")
plot_color<-rainbow(length(levels(all_groups)))[all_groups]

# Save the generated phylogram in a pdf file
pdf(my_tree_file_name)

# The tree is visualized as a Phylogram color-coded by the selected group name
plot(tree, type = "phylogram",use.edge.length = TRUE, tip.color = (plot_color), label.offset = 0.01)
print.phylo(tree)
axisPhylo()
tiplabels(pch = 16, col = plot_color)
dev.off()

#################            Build NMDS plot           ########################

# Generated figures are saved in a pdf file 
file_name <-paste("Beta-diversity","(",paste(groups_name,collapse="-"),").pdf",sep="")
pdf(file_name)

# Calculate the significance of variance to compare multivariate sample means (including two or more dependent variables)
# Omit cases where there isn't data for the sample (NA)

all_groups_comp <- all_groups[!is.na(all_groups)]
unifract_dist_comp <- unifract_dist[!is.na(all_groups), !is.na(all_groups)]
if (nlevels(all_groups_comp)>1){
adonis<-adonis(as.dist(unifract_dist_comp) ~ all_groups_comp)}
all_groups_comp<-factor(all_groups_comp,levels(all_groups_comp)[unique(all_groups_comp)])

# Calculate and display the MDS plot (Multidimensional Scaling plot)
if(nlevels(all_groups_comp)>1){
s.class(
  cmdscale(unifract_dist_comp, k = 2), col = unique(plot_color), cpoint =
    2, fac = all_groups_comp, sub = paste("MDS plot of Microbial Profiles\n(p-value ",adonis[[1]][6][[1]][1],")",sep="")
)} else {s.class(
  cmdscale(unifract_dist_comp, k = 2), col = unique(plot_color), cpoint =
    2, fac = all_groups_comp)}
if (label_samples==1) {
  lab_samples <- row.names(cmdscale(unifract_dist_comp, k = 2))
  ifelse (label_id != "",lab_samples <- replace(lab_samples, !(lab_samples %in% label_id), ""), lab_samples)
  text(cmdscale(unifract_dist_comp, k = 2),labels=lab_samples,cex=0.7,adj=c(-.1,-.8))
}

# Calculate and display the NMDS plot (Non-metric Multidimensional Scaling plot)
meta <- metaMDS(unifract_dist_comp,k = 2)
if (nlevels(all_groups_comp)>1){
  s.class(
  meta$points, col = unique(plot_color), cpoint = 2, fac = all_groups_comp,
  sub = paste("metaNMDS plot of Microbial Profiles\n(p-value ",adonis[[1]][6][[1]][1],")",sep="")
  )} else {s.class(
    meta$points, col = unique(plot_color), cpoint = 2, fac = all_groups_comp  )}
if (label_samples==1){
  lab_samples <- row.names(meta$points)
  ifelse (label_id != "",lab_samples <- replace(lab_samples, !(lab_samples %in% label_id), ""), lab_samples)
  text(meta$points,labels=lab_samples,cex=0.7,adj=c(-.1,-.8))
}

#close the pdf file
dev.off()

###############          NMDS for pairwise analysis        ###################

# This plot is only generated if there are more than two groups included in the comparison
# Calculate the pairwise significance of variance for group pairs
# Get all groups contained in the mapping file
unique_groups <- levels(all_groups_comp)
if (dim(table(unique_groups)) > 2) {
  
  # Initialise vector and lists
  pVal = NULL
  pairedMatrixList <- list(NULL)
  pair_1_list <- NULL
  pair_2_list <- NULL
  
  for (i in 1:length(combn(unique_groups,2)[1,])) {
    
    # Combine all possible pairs of groups
    pair_1 <- combn(unique_groups,2)[1,i]
    pair_2 <- combn(unique_groups,2)[2,i]
    
    # Save pairs information in a vector
    pair_1_list[i] <- pair_1
    pair_2_list[i] <- pair_2
    
    # Generate a subset of all samples within the mapping file related to one of the two groups
    inc_groups <-
      rownames(subset(meta_file, meta_file[,ncol(meta_file)] == pair_1
                      |
                        meta_file[,ncol(meta_file)] == pair_2))
    
    # Convert UniFrac distance matrix to data frame
    paired_dist <- as.data.frame(unifract_dist_comp)
    
    # Save all row names of the mapping file
    row_names <- rownames(paired_dist)
    
    # Add row names to the distance matrix
    paired_dist <- cbind(row_names,paired_dist)
    
    # Generate distance matrix with samples of the compared groups (column-wise)
    paired_dist <- paired_dist[sapply(paired_dist[,1], function(x) all(x %in% inc_groups)),]
    
    # Remove first column with unnecessary group information
    paired_dist[,1] <- NULL
    paired_dist <- rbind(row_names,paired_dist)
    
    # Generate distance matrix with samples of the compared group (row-wise)
    paired_dist <- paired_dist[,sapply(paired_dist[1,], function(x) all(x %in% inc_groups))]
    
    # Remove first row with unnecessary group information 
    paired_dist <- paired_dist[-1,]
    
    # Convert generated distance matrix to data type matrix (needed by multivariate analysis)
    paired_matrix <- as.matrix(paired_dist)
    class(paired_matrix) <- "numeric"
    
    # Save paired matrix in list
    pairedMatrixList[[i]] <- paired_matrix
    
    # Applies multivariate analysis to a pair out of the selected groups
    adonis <- adonis2(paired_matrix ~ all_groups_comp[all_groups_comp == pair_1 |
                                                       all_groups_comp == pair_2])
    
    # List p-values
    pVal[i] <- adonis[[1]][6][[1]][1]
    
  }
  
  # Adjust p-values for multiple testing according to Benjamini-Hochberg method
  pVal_BH <- p.adjust(pVal,method="BH", n=length(pVal))
  
  # Generated NMDS plots are stored in one pdf file called "pairwise-beta-diversity-nMDS.pdf"
  file_name <-paste("pairwise-beta-diversity-NMDS","(",paste(groups_name,collapse="-"),").pdf",sep="")
  pdf(file_name)
  
  for(i in 1:length(combn(unique_groups,2)[1,])){
    meta <- metaMDS(pairedMatrixList[[i]], k = 2)
    s.class(
      meta$points,
      col = rainbow(length(levels(all_groups_comp))), cpoint = 2,
      fac = as.factor(all_groups_comp[all_groups_comp == pair_1_list[i] |
                                        all_groups_comp == pair_2_list[i]]),
      sub = paste("NMDS plot of Microbial Profiles\n ",pair_1_list[i]," - ",pair_2_list[i], "\n(p-value ",pVal[i],","," corr. p-value ", pVal_BH[i],")",sep="")
    )
  }
  dev.off()
  
  # Generated MDS plots are stored in one pdf file called "pairwise-beta-diversity-MDS.pdf"
  file_name <- paste("pairwise-beta-diversity-MDS","(",paste(groups_name,collapse="-"),").pdf",sep="")
  pdf(file_name)
  
  for(i in 1:length(combn(unique_groups,2)[1,])){
    # Calculate and display the MDS plot (Multidimensional Scaling plot)
    s.class(
      cmdscale(pairedMatrixList[[i]], k = 2), col = rainbow(length(levels(all_groups_comp))), cpoint =
        2, fac = as.factor(all_groups_comp[all_groups_comp == pair_1_list[i] |
                                             all_groups_comp == pair_2_list[i]]), sub = paste("MDS plot of Microbial Profiles\n ",pair_1_list[i]," - ",pair_2_list[i], "\n(p-value ",pVal[i],","," corr. p-value ", pVal_BH[i],")",sep="")
    )
  }
  dev.off()                                     
  
}


##############################################################################################
#########################      Determine the number of clusters   ############################
##############################################################################################


for (i in 1:(length(Test_name)+1)){
  
if (i==1){
  name <- paste(Reference_name,collapse="+")
  group <- Reference_name
} else {
  name <-Test_name[i-1]
  group <- Test_name[i-1]
  
  
}
  
  # Create the directory where the outputs for every group will be saved 
  group_results_folder <- name
  group_results_path <- paste0(optimal_number_path,"/",group_results_folder)
  dir.create(paste0(optimal_number_path,"/",group_results_folder))
  
  # Change working directory
  setwd(group_results_path)
  
  
  
  # Variable where the Calinski-Harabasz values will be stored
  ch_nclusters=NULL
  
  # Variable where the silhouette values will be stored
  sil_nclusters=NULL
  
  # Create a PDF file where the MDS plot for the different k will be stored
  pdf(paste0("MDS plots-",name,".pdf"))
  
  if (dim(otu_file[meta_file[,mapping_column]==group,])[1]-1 <= kmers_limit) {
    kmers_limit=dim(otu_file[meta_file[,mapping_column]==group,])[1]-1
  }
  for (k in 1:kmers_limit) { 
    if (k==1) {
      ch_nclusters[k]=NA 
      sil_nclusters[k]=NA
      data_cluster=as.vector(pam(as.dist(unifract_dist[meta_file[,mapping_column]%in%group,meta_file[,mapping_column]%in% group]), k, diss=TRUE)$clustering)
      
      #MDS plot for k=1
      sclass(unifract_dist[meta_file[,mapping_column]%in%group,meta_file[,mapping_column]%in% group],as.factor(data_cluster),k)
    } else {
      # Partitioning the data into k clusters (max k is number of samples within the dataset)
      data_cluster=as.vector(pam(as.dist(unifract_dist[meta_file[,mapping_column]%in%group,meta_file[,mapping_column]%in% group]), k, diss=TRUE)$clustering)
      
      # Calculate Calinski-Harabasz and silhouette Index 
      index=cluster.stats(as.dist(unifract_dist[meta_file[,mapping_column]%in% group,meta_file[,mapping_column]%in% group]),data_cluster)
      ch_nclusters[k] <- index[["ch"]]
      sil_nclusters[k] <-index[["avg.silwidth"]]
      print(k)
      
      #  MDS plot for k 2-kmers_limit
      sclass(unifract_dist[meta_file[,mapping_column]%in% group,meta_file[,mapping_column]%in% group],as.factor(data_cluster),k)
      
      
    }
  }
  dev.off()
  
  # Calculate Within sum of squares
  wss <- fviz_nbclust(otu_file[meta_file[,mapping_column]%in% group,] ,diss=unifract_dist[meta_file[,mapping_column]%in% group,meta_file[,mapping_column]%in% group], k.max = kmers_limit  ,cluster::pam, method =  "wss")
  
 
  
#################################################################################
######                        Write Output Files                           ######
#################################################################################

  
  # Generated plots showing the optimal number of clusters
  pdf(paste0("Indices-",name,".pdf"))
  
  # Plot Calinski-Harabasz plot
  plot(ch_nclusters, type="h", xlab="k clusters", ylab="CH index",main=paste0("'",name,"'"," Optimal number of clusters (CH)"))
  # Plot silhouette Index plot
  plot(sil_nclusters, type="h", xlab="k clusters", ylab="Average silhouette width",main=paste0("'",name,"'"," Optimal number of clusters(Silhouette)"))
  # Plot WSS plot
  plot(1:kmers_limit,wss[["data"]][["y"]], type="b", xlab="k clusters", ylab="Within sum of squares",main=paste0("'",name,"'"," Optimal number of clusters(WSS)"),pch=19)
 
  
  dev.off()
  
}


# Change working directory
setwd(results_path)
# Write the distance matrix table in a file 
write.table( unifract_dist, paste0("distance-matrix-gunif","(",paste(groups_name,collapse="-"),").tab"), sep = "\t", col.names = NA, quote = FALSE)

# Change working directory
setwd(beta_diversity_path)
# Write samples' tree in Newwick format
write.tree(tree,"samples-Tree.nwk",tree.names = FALSE)

# Return to the original path
setwd(OriginalPath)

# Graphical output files are generated in the main part of the script
if(!flag) { stop("
    It was not possible to install all required R libraries properly.
                 Please check the installation of all required libraries manually.\n
                 Required libaries:ade4, GUniFrac, phangorn")
}
}

# Remove unnecessary variables
rm(list=ls()[!(ls()%in%ls)])
# Clear the memory
invisible(gc())

#################################################################################
######                           End of Script                             ######
#################################################################################
