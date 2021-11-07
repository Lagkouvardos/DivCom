#'
#' The script performs de novo clustering of the groups and then calculate the distances from 
#' the most representative points of the each cluster. 
#' Then, conducts statistical analysis and produces plots and tables. 
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
#' 9. Write the desired number of clusters for the reference group
#' 9. Write the name of groups that will be used for as reference groups
#'10. Write the desired number of clusters for the test groups
#'11  Write the names of the columns used for chi-square comparisons
#'12. Write if the central point of each cluster will be the medoid, the median or mean point
#'13. Write the preffered type of plot 
#' 
#'
#' Output: 
#' The script generates three graphical outputs (pdf), one text file and a newick tree
#' 1. A Distances Based Analysis report
#' 2. A De Novo Analysis report
#' 3. A folder with the p-values tables for the different tests 
#' 4. A folder with the plots presenting in the reports
#' 5. A folder with various tables containing information about the distances 
#' 6. The mapping file with an extra column 
#' 
#'
#' Concept:
#' A common problem in data analysis is to efficiently comparing a group of control samples with a set of test samples.
#' The relation existing between different groups can be highly affected by their substructure.
#' Instead of comparing the groups as entireties, the program performs de novo clustering(with PAM algorithm) 
#' to both the control and test groups.
#' Then, finds the most representative points of the reference  dataset and calculates the distances(Generalized Unifrac) 
#' of the test groups from these points.
#' To determine  the level of similarity between the groups, the script performs statistical analysis and produces various plots.
#' 


###########################################################################################################################################################
###################################################### Set parameters in this section manually ############################################################
###########################################################################################################################################################


#' Please set the directory of the script as the working folder (e.g D:/studyname/NGS-Data/Rhea/beta-diversity/)
#' Note: the path is denoted by forward slash "/"
setwd("C:/...../..../Distances")  #<--- CHANGE ACCORDINGLY !!!


#' Please give the name of the OTUs or ASVs table
input_otu = "SOTUs-Table.tab" 


#' Please insert "YES" if the OTUs/ASVs table is normalized otherwise, insert "NO"
#' In case the table is not normalized, the OTUs/ASVs table will be normalized based on the minimum sum of reads of the samples.
normalized = "NO"


#' Choose if you will insert a distances matrix or the phylogenetic tree of the OTUs/ASVs sequences
#' There are two options: "distances matrix" or "tree"
#' 1) Please insert "distances matrix" if you will provide a distances matrix
#' 2) Please insert "tree" if you will provide a phylogenetic tree (In this case, the Generalized unifrac distances will be calculated)
tree_or_matrix =  "tree" 


#' Please insert the name of the distances matrix or the anme of the phylogenetic tree 
#' -> -> !!! In case you will choose the "mean" or "median" option you HAVE TO provide a phylogenetic tree !!! <- <-
input_tree_or_matrix= "SOTUs-NJTree.tre"  


#' Please give the name of the mapping file which contains the labels of the samples
#' !! CAUTION: The rows of the mapping file should have the same sample names as the OTUs table !!
input_meta = "mapping_file.tab"


#' Please provide the name or the number of the column (of the mapping file) based on which the samples will be partitioned into groups
mapping_column = "Type4"


#' Please place in the vector one or more names which will be used to identify the samples that composing 
#'         the REFERENCE group (e.g reference_name <- c("group_a","group_b"))
#' -> -> !!! CAUTION: You should provide at least one name!!! <- <-
reference_name = c("NI")

#' Please insert the desired number of clusters for the reference group  
#' -> -> CAUTION: You have to provide this information <- <-
reference_clusters = 2



#' Please provide the names of the test groups 
#' There are two options: a User-defined vector or "None"
#' 1) User-defined vector --> Form a vector with one or more elements referring to the name of the groups (e.g test_name <- c("group_a","group_b"))         
#' 2)     "None"          --> There won't be any test group. Only the distances of the reference samples will be calculated (Not recommended)
Test_name = c("IBD.B","IBD.3M")


#' OPTIONAL-  Please insert the desired number of clusters for every test group (e.g test_clusters= c(3,2))
#' -> -> !! CAUTION: The number of clusters should be more than one !! <- <-
#' If you do not wish to perform de novo clustering to the test groups just insert an empty vector  (e.g test_clusters= c())
test_clusters = c(3,3)


#' OPTIONAL- Please insert the names of the columns of the mapping file you wish to analyze against the de novo clusters 
#' -> -> !! CAUTION: The columns must contain categorical variables; if that is not the case, the analysis will NOT be conducted!! <- <- 
exploratory_columns = c("Type5","Category","Iron")


#' OPTIONAL- Please choose if the medoid, mean or median points of every cluster will be used as the most central(representative) points
#' There are three options: "medoid", "mean", "median"
#' 1) Insert "medoid" if you you wish to use the medoids as the most representative points. !! DEFAULT OPTION !!
 
#' 2) Insert  "mean"  if you you wish to use the means as the most representative points.
#'  -> -> !! WARNING: In this case you should provide the phylogenetic tree of the OTUs/ASVs sequences !! <- <-
  
#' 3) Insert "median" if you you wish the medians to be the most representative points.
#'  -> -> !! WARNING: In this case you should provide the phylogenetic tree constructed of the OTUs/ASVs sequences !! <- <-
central_point = "medoid" 


#' OPTIONAL- Please select the type of the output plots. 
#' There are three options: "Boxplots", "Point plots", "Violin plots"
#' 1) "Boxplots"     -> This is the default type of plot.
#' 2) "Point plots"  -> All the samples will be presented as individual points. This option is preferable when the samples are low in number.
#' 3) "Violin plots" -> This option displays in a more detailed view the distribution of the values.
plot_type = "Boxplots"


#################################################################################################################################################################
################################################# DO NOT CHANGE ANYTHING BELOW THIS LINE #################################################################################
#################################################################################################################################################################


############################################################################################################################################################
############################################################ Main Script ###################################################################################
############################################################################################################################################################


##################  Create all the necessary paths and directories ###############

# Save the current path in one variable (will be used later)
OriginalPath   <- getwd()

# Create the directory where all output files will be saved 
outputs_folder <- paste0(mapping_column,"-",paste(reference_name,collapse="+"),"-",reference_clusters," (",paste(Test_name,collapse="+"),")")
outputs_path   <- paste0(OriginalPath,"/",outputs_folder)
dir.create(outputs_folder)

# Path and Directory where the results will be stored
results_folder <- paste0("results","-",paste(reference_name,collapse="+"),"-",reference_clusters)
results_path   <- paste0(outputs_path,"/",results_folder)
dir.create(paste0(outputs_path,"/",results_folder))

# Path and Directory where the p-values will be stored
pvalues_folder <- paste0("P-values","-",paste(reference_name,collapse="+"),"-",reference_clusters)
pvalues_path   <- paste0(results_path,"/",pvalues_folder)
dir.create(paste0(results_path,"/",pvalues_folder))

# Path and Directory where the plots will be stored
plots_folder   <- paste0("Plots","-",paste(reference_name,collapse="+"),"-",reference_clusters)
plots_path     <- paste0(results_path,"/",plots_folder)
dir.create(paste0(results_path,"/",plots_folder))

# Path and Directory where the rest tables will be stored
tables_folder  <- paste0("Tables","-",paste(reference_name,collapse="+"),"-",reference_clusters)
tables_path    <- paste0(results_path,"/",tables_folder)
dir.create(paste0(results_path,"/",tables_folder))


############################################################################################################################################################
#############################################  Load all required libraries #################################################################################
############################################################################################################################################################


# Check if required packages are already installed, and install if missing
packages <- c("ade4","GUniFrac","phangorn","cluster","fpc","ggplot2","gridExtra","grid","gtable","stats","cowplot",
             "graphics","vegan", "dplyr","data.table","tidyr","caTools","RColorBrewer","ggpubr","ggtree","tools") 

# Function to check whether the package is installed
InsPack <- function(pack)
{
  if ((pack %in% installed.packages()) == FALSE) {
    install.packages(pack,repos ="http://cloud.r-project.org/")
  } 
}

if (("ade4" %in% installed.packages()) == FALSE) {
  install.packages("ade4",repos ="http://cloud.r-project.org/")
} 

# Applying the installation on the list of packages
lapply(packages, InsPack)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree",update=FALSE,force = TRUE)

# Make the libraries
lib <- lapply(packages, require, character.only = TRUE)

# Check if it was possible to install all required libraries
flag <- all(as.logical(lib))


###########################################################################################################################################################
#################################################  Functions that will be used in main Script. ############################################################
###########################################################################################################################################################

#----------------------------------------------------------------------------#
#************** Fuction that saves plots in png format **********************#
#----------------------------------------------------------------------------#

# name -> name of the output file
# plot -> name of the plot

png_plot <- function(name,plot) {
  
  # Set the working directory
  setwd(plots_path)
  
  # Set the parameters of the output file
  png(name,width = 595,height=842)
  
  # Save the plot in png format
  print(plot)
  dev.off()
  
  # Return to the output folder
  setwd(outputs_path)
}


#----------------------------------------------------------------------------#
#************** Fuction that saves tables in tab format *********************#
#----------------------------------------------------------------------------#

# name   -> name of the output file
# table  -> name of the table
# folder -> Define where the tables will be stored (In the p-values of the table folder)

write_table <- function(name, table,folder){
  
  # Check where the file will be stored
  if(folder=="pvalues"){setwd(pvalues_path)} else if (folder=="table"){setwd(tables_path)} 
  
  # Write the table
  write.table(table, name,col.names =NA, row.names = TRUE,quote = FALSE)
  
  # Return to the output folder
  setwd(outputs_path)
  
}

##############################################################################
#*************** Necessary functions for the plots **************************#
##############################################################################



#------------ Theme of Gtables ---------------#

# Set the theme to change text for plotting (ggplot - Gtable)
mytheme <- gridExtra::ttheme_default(
  
  # Adjust settings for the text inside table
  core = list(fg_params = list(cex = 0.8)),
  
  # Adjust the test for column and row header
  colhead = list(fg_params = list(cex = 0.85)),
  rowhead = list(fg_params = list(cex = 0.85))
)

#----------------------------------------------------------------------------#
#************************ Fuction for gTables *******************************#
#----------------------------------------------------------------------------#

# title -> title of the table
# table -> input table that will be converted into tableGrob
# setDT -> if NO the table will remain data frame otherwise will be converted into data table
# size  -> size of the title 

gtableGrob <- function(title,table,setDT,size){
  
  if (setDT=="NO"){
    # Create the Grob object
    pvaltable <- tableGrob(table,rows = NULL ,theme = mytheme)
    } else {
      # Create the Grob object
      pvaltable <- tableGrob(setDT(table, keep.rownames = "")[],rows = NULL,theme = mytheme)
    }
  
  # Title of tables 
  title <- textGrob(title,gp = gpar(fontsize = size))
  padding <- unit(2,"mm")
  # Adjust the grob placement accordingly
  pvaltable <- gtable_add_rows(pvaltable, heights = grobHeight(title) + padding,pos = 0)
  pvaltable <- gtable_add_grob(pvaltable, title, 1, 1, 1, ncol(pvaltable),clip = "off")
}

#------------------------------------------------------------------------------#
#*********** Function that calculates the Prevalence of the groups ************#
#------------------------------------------------------------------------------#


# dist   -> a data frame, in the first column will have the group indentities, and in another column with distances
# column -> the number of the column that contain the distances

prevalence <- function(dist, column) {
  # Find the factors of the first column
  Levels <- as.factor(unique(as.character(dist[,1])))
  
  # Add a new column where the prevalence of non-Na elements will be placed
  dist$label <- c(rep(0,nrow(dist)))
  # Convert the label column into characters
  dist$label <- as.character(dist$label)
  
  #--- Calculate the prevalence of non-Na elements of each group and assign the prevalence value of a group to each sample---#
  
  #Place label vector in the input_table
  
  for (l in 1:length(Levels)) {
    # Calculate the prevalence 
    label <- paste0(Levels[l],"(",length(which(!is.na(dist[dist[,1]==Levels[l],column]))),"/",length(dist[dist[,1]==Levels[l],column]),")")
    dist[which(dist[,1]== Levels[l]),"label"] =label
  }
  # Convert the label column into factor
  dist$label <- as.factor(dist$label)
}

#------------------------------------------------------------------------------#
#********** Functions that calculates the p-values of the Wilcoxon test *******#
#------------------------------------------------------------------------------#

#********** 1st function *******#

# dist_matrix -> matrix containg the distances and the labels
#   report    -> "DeBaAn" or "DeNoAn"  
#    page     -> In which page of the final report this table will be placed

pvalues_function <- function(dist_matrix,report,page){
  
  # Find all the available pairwise combinations
  combinations <- combs(levels(as.factor(dist_matrix[,1])), 2)
  
  # pvalues -> data frame where the p-values will be placed
  pvalues <- data.frame()
  
  
  # index indicating the number of row 
  row <- 0
  
  for (g in 1:nrow(combinations)){
    # Add a row
    row <- row+1
    # Performing Wilcoxon Rank Sum Test
    wilcox <- tryCatch (round(wilcox.test(dist_matrix[dist_matrix[,1]==combinations[g,1],2],dist_matrix[dist_matrix[,1]==combinations[g,2],2])$p.value,4), error = function(i) {wilcox <- NA})
    # In the first column add the names of the two groups
    pvalues[row,1]  <- paste0(combinations[g,1],"-",combinations[g,2])
    # In the second column add the p-value
    pvalues[row,2] <- wilcox
  }
  # Performing fdr correction
  p.adjust <- round(p.adjust(pvalues[,2], method = "BH"),4)
  # Add the adjusted p-values in the pvalues matrix
  pvalues <- data.frame(pvalues,p.adjust)
  
  # Name the columns
  colnames(pvalues) <- c("Groups","p-value","Adj. p-value")
  
  # Save all the p-values in table format
  write_table(paste0(report," page ",page," Wilcoxon p-values.tab"),pvalues,"pvalues")
  
  # Convert p-values column into numeric
  pvalues[,2] <- as.numeric(pvalues[,2])
  # Order the p-values
  pvalues <- pvalues[order(pvalues[,2]),]
  # Choose the first 10 most significant p-values
  if (nrow(pvalues)>10){ pvalues2 <- pvalues[1:10,] } else { pvalues2 <- pvalues }
}


#********** 2nd function *******#


#      dist_matrix   -> matrix containg the distances and the labels
# reference_clusters -> Number of reference clusters
#      report        -> "DeBaAn" or "DeNoAn"
#       page         -> In which page of the final report this table will be placed

pvalues_function2 <- function(dist_matrix,reference_clusters,report,page){
  
  # pvalues -> data frame where the p-values will be placed
  pvalues <- data.frame()
  
  # index indicating the number of row 
  row <- 0
  
  for (cl in 1:reference_clusters){
    
    # matrix with the distances of every cluster
    wilcoxon_dist_matrix <- dist_matrix[dist_matrix[,3]==cl,]
    
    if(length(as.character(unique(wilcoxon_dist_matrix[,1])))>1){
      # Find all the available pairwise combinations
      combinations <- combs(as.character(unique(wilcoxon_dist_matrix[,1])), 2)
      
      for (g in 1:nrow(combinations)){
        # Add a row
        row <- row+1
        # Performing Wilcoxon Rank Sum Test
        wilcox <- tryCatch (round(wilcox.test(wilcoxon_dist_matrix[wilcoxon_dist_matrix[,1]==combinations[g,1],2],wilcoxon_dist_matrix[wilcoxon_dist_matrix[,1]==combinations[g,2],2])$p.value,4), error = function(i) {wilcox <- NA})
        # In the first column add the names of the two groups
        pvalues[row,1]  <- paste0(combinations[g,1],"-",combinations[g,2])
        # In the second column add the p-value
        pvalues[row,2] <- wilcox
      }
    }
  }
  
  # Performing fdr correction
  p.adjust <- round(p.adjust(pvalues[,2], method = "BH"),4)
  # Add the adjusted p-values in the pvalues matrix
  pvalues <- data.frame(pvalues,p.adjust)
  
  
  # Name the columns
  colnames(pvalues) <- c("Groups","p-value","Adj. p-value")
  
  # Save all the p-values in table format
  write_table(paste0(report," page ",page," Wilcoxon p-values.tab"),pvalues,"pvalues")
  
  # Convert p-values column into numeric
  pvalues[,2] <- as.numeric(pvalues[,2])
  # Order the p-values
  pvalues <- pvalues[order(pvalues[,2]),]
  # Choose the first 10 most significant p-values
  if (nrow(pvalues)>10){ pvalues2 <- pvalues[1:10,] } else { pvalues2 <- pvalues }
  
}

###########################################################################################################
########################## Functions for the visualization ################################################
###########################################################################################################


#------------------------------------------------------------------------------#
#****************** Function that generates cladorgams ************************#
#------------------------------------------------------------------------------#

# distance -> The given distances matrix 
# groups   -> Vector that contains the indexes of the groups 
# colours  -> Vector with the colours of each group

tree <- function(distance,groups,colours) {
  
  # Save the rownames of the distances matrix
  all_groups <- c(rownames(distance))
  # Save the UniFrac output as distance object
  all_dist_matrix <- as.dist(distance)
  
  # Apply a hierarchical cluster analysis on the distances matrix based on the Ward's method
  all_fit <- hclust(all_dist_matrix, method = "ward.D2")
  
  # Matrix containing the samples and the corresponding groups the belong
  Samples_matrix <- data.frame(all_groups,groups)
  # Name the rows of the Samples_matrix
  rownames(Samples_matrix) <- Samples_matrix[,1]
  
  # Generates a tree from the hierarchically generated object
  tree <- as.phylo(all_fit)
  
  # Generate the circular cladogram
  ggtree_plot <- ggtree(tree,  layout='circular')
  # Add the annotation ring
  gheatmap( ggtree_plot, Samples_matrix[,2,drop=F] , offset=0, width=0.1, colnames_position="bottom", colnames_angle=90, colnames_offset_y = 0.25,hjust=10)+
    theme(plot.title=element_text( hjust=0.5, face='bold',size=18))+
    scale_fill_manual(breaks=c(levels(as.factor(groups))), values=colours,name="Groups")
  
}


#------------------------------------------------------------------------------#
#*************** Function that generates MDS and NMDS plots *******************#
#------------------------------------------------------------------------------#

#   distance  -> the given distances matrix 
#    groups   -> a vector that contains the indexes of the groups 
# individuals -> a vector containing Samples of the distance matrix that will be represented on MDS plot as individual points  
#    colours  -> Vector with the colours of each group


sclass <- function(distance,groups,individuals=NULL,col) {
  
  # Performing classical multidimensional scaling
  mds <- cmdscale(distance,eig=T, x.ret=T)
  # Computing the variance of each axis
  mds.variation <- round(mds$eig/sum(mds$eig)*100,1) # axes variance calculation
  # Create a color palette
  plot_color <- rainbow(length(levels(groups)))[groups]
  all_groups_comp <- groups[!is.na(groups)]
  
  # Performing analysis of variance using the distance matrix(only if individuals=NULL)
  if(is.null(individuals)==TRUE & nlevels(groups) > 1) {
    # PERMANOVA test
    adonis <- adonis(distance ~ all_groups_comp)
    # Create Subtitle
    sub <-  paste("MDS plot of Microbial Profiles\n(p-value ",adonis[[1]][6][[1]][1],")",sep="")
  } else {
    sub <- c("")
  }
  # Find the factors of all_groups_comp vector
  all_groups_comp <- factor(all_groups_comp,levels(all_groups_comp)[unique(all_groups_comp)])
  if (is.null(individuals) == FALSE) {
    # Create the colour vector
    colors <- as.character(brewer.pal(n = (length(Test_name)+1), name = "Set2")[2:(nlevels(as.factor(individuals))+1)])
    # find the levels of the individual points
    plot_levels <- levels(individuals)
    # Create a color palette for the individual points
    colo <- c()
    for (i in 1:length(plot_levels)) {colo <- c(colo,rep(colors[i],length(which(as.factor(mapping_cluster[which(input_table[,1] %in% levels)])==plot_levels[i]))))
    }
  }
  
  # Display the MDS plot (Multidimensional Scaling plot)
  if (is.null(individuals)==TRUE){
    s.class(
      mds$points, col = col, cpoint =
        2, fac = all_groups_comp, grid = T, sub = sub)
    graphics:: title (main=paste0("MDS graph " ), xlab = paste("MDS1:",mds.variation[1],"%"), ylab=paste("MDS2:",mds.variation[2],"%"))
  } else {
    s.class(
      mds$points[1:length(groups),], col = col, cpoint =
        1.5, grid = T,fac = all_groups_comp[1:length(groups)],xlim = c(1.4*min(mds[["points"]][,1]),1.2*max(mds[["points"]][,1])),ylim = c(1.2*min(mds[["points"]][,2]),1.5*max(mds[["points"]][,2])))
    points(mds$points[(length(groups)+1):(length(groups)+length(individuals)),1],mds$points[(length(groups)+1):(length(groups)+length(individuals)),2],pch = 19,col=colo,cex=0.70)
    graphics:: title (main=paste0(" " ), xlab = paste("MDS1:",mds.variation[1],"%"), ylab=paste("MDS2:",mds.variation[2],"%"))
    legend(1.3*min(mds[["points"]][,1]),1.1*max(mds[["points"]][,2]),c(levels(as.factor(groups)),levels(as.factor(individuals))),
           fill=c(rep(brewer.pal(n = (length(Test_name)+1), name = "Set2")[1],reference_clusters),colors),cex=0.8)
  }
}


#------------------------------------------------------------------------------#
#********** Functions that generate boxplots, violinplots, pointplots *********#
#------------------------------------------------------------------------------#
  
# name      -> The name will be used as the title for the plots
# plot_df   -> data frame containing the distances and the names
# label     -> name of the column of plot_df containing the names of the boxplots
# abundance -> name of the column of plot_df containing the distances
# colours   -> Vector with the colours of each group


# Create a ggplot object of the plotted layout (including axis labels and scaling)-Boxplots
# Generate plots 
my_boxplot <- function(name,plot_df,label,abundance ,colour=NULL){
  g <-ggplot(plot_df,aes(x = label ,y = abundance,fill=label))
  if(is.null(colour)==TRUE) {
    my_boxplot <- g + stat_boxplot(geom = "errorbar", width = 0.25) + geom_boxplot(varwidth = FALSE,width = 0.7)+scale_fill_manual(values=rainbow(length(unique(plot_df[,2]))))
  } else {
    my_boxplot <- g + stat_boxplot(geom = "errorbar", width = 0.25) + geom_boxplot(varwidth = FALSE,width = 0.7)+scale_fill_manual(values=colour) 
  }
  my_boxplot <-my_boxplot + ggtitle(name) + guides(fill = "none") + ylab("Distances") + xlab("") + theme_bw() + theme(axis.text.x = element_text(colour = "grey20",size = 12,angle = 45,hjust = 1,vjust = 1,face = "plain"),
                                                                                                                      axis.text.y = element_text(colour = "grey20",size = 7,angle = 0,hjust = 1,vjust = 0,face = "plain"),axis.title.y = element_text(colour = "grey20",size = 14,angle = 90,hjust = .5,vjust = .5,face = "plain"),
                                                                                                                      plot.title = element_text(colour = "grey22",size = 15,hjust = .5,vjust = .5,face = "bold")) + theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))
}

# Create a ggplot object of the plotted layout (including axis labels and scaling)-violinplots
# Generate plots
my_violinplot <- function(name,plot_df,label,abundance,colour=NULL) {
  # Create a ggplot object of the plotted layout (including axis labels and scaling)
  g <-ggplot(plot_df,aes(x = label ,y = abundance,color=label))
  if(is.null(colour)==TRUE) {
    my_violinplot <- g + geom_violin(width = 0.7) + geom_boxplot(width = 0.1)+scale_fill_manual(values=rainbow(length(unique(plot_df[,2]))))} else {
      my_violinplot <- g + geom_violin(width = 0.7) + geom_boxplot(width = 0.1)+scale_fill_manual(values=colour)
    }
  my_violinplot <-my_violinplot + ggtitle(name) + ylab("Distances") + guides(fill = "none") + xlab("") + theme_bw() +theme(axis.text.x = element_text(colour = "grey20",size = 12,angle = 45,hjust = 1,vjust = 1,face = "plain"),
                                                                                                                           axis.text.y = element_text(colour = "grey20",size = 7,angle = 0,hjust = 1,vjust = 0,face = "plain"),axis.title.y = element_text(colour = "grey20",size = 14,angle = 90,hjust = .5,vjust = .5,face = "plain"),
                                                                                                                           plot.title = element_text(colour = "grey22",size = 15,hjust = .5,vjust = .5,face = "bold")) + theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+ theme(legend.position="none")
}

# Create a ggplot object of the plotted layout (including axis labels and scaling)-pointplots
# Generate plots
my_point_boxplot <- function(name,plot_df,label,abundance,colour=NULL){
  # Create a ggplot object of the plotted layout (including axis labels and scaling)
  g <-ggplot(plot_df,aes(x = label ,y = abundance,fill=label))
  if(is.null(colour)==TRUE) {
    my_point_boxplot <- g + geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.7)+scale_fill_manual(values=rainbow(length(unique(plot_df[,2]))))} else {
      my_point_boxplot <- g + geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.7)+scale_fill_manual(values=colour)
    }
  my_point_boxplot <- my_point_boxplot + ggtitle(name) + theme_bw() + ylab("Distances") +stat_summary(fun = median, fun.min = median, fun.max = median,geom = "crossbar", width = 0.5) +
    guides(fill = "none") + xlab("") + theme(legend.position = "none") + guides(colour = "none") + theme( axis.text.x = element_text(colour = "grey20",size = 12,angle = 45,hjust = 1,vjust =1,face = "plain"),
                                                                                                          axis.text.y = element_text(colour = "grey20",size = 7,angle = 0,hjust = 1,vjust = 0,face = "plain"),axis.title.y = element_text(colour = "grey20",size = 14,angle = 90,hjust = .5,vjust = .5,face = "plain"),
                                                                                                          plot.title = element_text(colour = "grey22",size = 15,hjust = .5,vjust = .5,face = "bold")) + theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))
}


###########################################################################################################################################################
########################################################          Data Loading and      ###################################################################
########################################################           pre-processing       ###################################################################
###########################################################################################################################################################



############################################################
#################### Data Loading ##########################
############################################################


#------------------- Meta_file--------------------------#

# Load the mapping file containing individual samples information (sample names should be in the first column of the file)
meta_file <- data.frame(read.table (file = input_meta, check.names = FALSE, header = TRUE, dec = ".", sep = "\t", row.names = 1, comment.char = ""))

# Clean table from empty lines
#meta_file <- data.frame(meta_file[!apply(is.na(meta_file) | meta_file=="",1,all),])

# Order the rows of the meta file
if (ncol(meta_file)==1){
  meta_file <- data.frame(meta_file[order(meta_file[,mapping_column ]),, drop = FALSE])
  colnames(meta_file) <- mapping_column 
  meta_file <- data.frame(meta_file[meta_file[,mapping_column] %in% c(reference_name ,Test_name),,drop=FALSE])
} else {
  meta_file <- data.frame(meta_file[order(meta_file[,mapping_column]),])
  meta_file <- data.frame(meta_file[meta_file[,mapping_column] %in% c(reference_name ,Test_name),])
}


# Convert the selected column (mapping_column) of meta_file into factor
meta_file[,mapping_column] <- as.factor(meta_file[,mapping_column])

#-------------------------------------------------------#




#------------------------ OTUs table--------------------------#

# Load the tab-delimited file containing the values to be checked (row names should be in the first column of the file)
otu_table <-  read.table (input_otu,check.names = FALSE,header = TRUE,dec = ".",sep = "\t", row.names = 1,comment.char = "")

# Clean table from empty lines
otu_table <- otu_table[!apply(is.na(otu_table) | otu_table=="",1,all),]

# Check if a column with taxonomy information exists
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
  
  #---- Normalization of the OTUs table (if the user ask for)----#
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
      
      # Reorder the OTUs table
      otu_table <- otu_table[,rownames(meta_file)]
      
      # Calculate the minimum sum of all columns/samples
      min_sum <- min(colSums(otu_table))
      
      # Divide each value by the sum of the sample and multiply by the minimal sample sum
      otu_file <- t(min_sum * t(otu_table) / colSums(otu_table))
      
      # Clean table from empty lines
      otu_file <- otu_file[!apply(is.na(otu_file) | otu_file =="",1,all),]
      
      # Transpose OTU-table and convert format into a data frame
      otu_file <- data.frame(t(otu_file))
      
      # Create the otu_file
      otu_file <- otu_file[rownames(meta_file),]
      
    } else {
      
      # Transpose the OTU-table
      otu_table <- data.frame(t(otu_table))
      
      # Reorder the OTUs table
      otu_table <- otu_table[,rownames(meta_file)]
      
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
  
  
  
  #---------------- Tree -----------------------#
  
  if (tree_or_matrix=="tree") {
    # Load the phylogenetic tree calculated from the OTU sequences 
    tree_file <- read.tree(input_tree_or_matrix)
    #tree_file <- keep.tip(tree_file,colnames(otu_file))
    
    
    # Root the OTU tree at midpoint 
    rooted_tree <- midpoint(tree_file)
  }
  #--------------------------------------------#
  
  
  
  
  ########################################################
  ################ Data pre-processing ###################
  ########################################################
  
  
  # Forming the levels vector 
  levels <- levels(meta_file[,mapping_column])
  levels <- levels[!(levels %in% reference_name)]
  
  
  # index= 0 -> De novo clustering to the test groups will be performed
  # index= 1 -> De novo clustering to the test groups will not be performed (The remaing procedures will be performed normally)
  index = 0
  
  
   #------------- Perform all the necessary checks ---------------# 
  
  # Check if has been provided the information about the desired number of clusters for the test groups
  if (length(test_clusters) == 0 | Test_name == "None" | any(test_clusters == 1)) {test_clusters <- c()
  index <- 1}
  
  
  #  Check if the names of the test groups have been provided, if not Test_name = "None"
  if (length(Test_name) == 0 | length(levels) == 0) {Test_name <- "None"
  test_clusters <- c()
  index <- 1}
  
  
  # Check if the number of the given test clusters is equal to the number of the test groups
  if (length(test_clusters) != length(levels) & Test_name == "All") {test_clusters <- c()
  index <- 1}
  
  # Check if the number of the given test clusters is equal to the number of the test groups
  if (length(test_clusters )!= length(Test_name)) {test_clusters <- c()
  index <- 1}
  
  # Check if a phylogenetic tree has been provided when the user chose the "mean" or "median" option
  if (central_point != "medoid" & tree_or_matrix != "tree") {central_point <- "medoid"}
  
  # Check if the tree_or_matrix and input_tree_or_matrix are consistent to each other
  if ( file_ext(input_tree_or_matrix) == "nwk" | file_ext(input_tree_or_matrix) == "tre"){tree_or_matrix = "tree"}
  
  # Check if the the names of the exploratory_column are present in the mapping file
  if (all(exploratory_columns %in% colnames(meta_file)) == TRUE & length(exploratory_columns) >= 1){
    chi_square <- 1
  } else { chi_square <- 0 
  }
 #--------------------------------------------------------------#
  
  
  #------- OTUs and mapping file preparation ---------#
  
  # OTUs table and mapping file transformation in order to include only the chosen by the user groups
  if (Test_name == "None") {
    levels <- c()
    otu_file <- rbind(otu_file[meta_file[,mapping_column] %in% reference_name,])
    meta_file <- meta_file[meta_file[,mapping_column] %in% reference_name,,drop=FALSE]
  } else {
    levels <- Test_name
    otu_file <- rbind(otu_file[meta_file[,mapping_column] %in% reference_name,],otu_file[meta_file[,mapping_column] %in% levels,])
    meta_file <- rbind(meta_file[meta_file[,mapping_column] %in% reference_name,,drop=FALSE],meta_file[meta_file[,mapping_column] %in% levels,,drop=FALSE])
  }
  
  #-------------------------------------------------#
  
  
  #------------- Additional Checking ---------------#
  # Check if it is meaningful to continue the procedure
  { if (sum(meta_file[,mapping_column] %in% reference_name) < reference_clusters) 
    stop("The number of the reference samples are less than the number of clusters!")
    if (all(is.na(otu_file)))
      stop("The mapping file and the OTUs table does not have same rows names")
    if (any(reference_name%in%Test_name))
      stop("You have entered the same names as reference and test groups.\n Each group can be either reference or test not  both!")
    if (all(c(Test_name,reference_name)%in%levels(factor(meta_file[,mapping_column]))) == FALSE)
      stop("The names of the Reference and Test groups does not exist in the mapping file!")
    
    # Check if empty lines exist in otu file. If yes then these lines will be removed from otu file and mapping file
    if(anyNA(otu_file) == TRUE) {
      cat("There is inconsistency between rownames of mapping file and rownames of OTUs file.\nWARNING all the empty lines will be removed from OTUs table and mapping file!\n")
      # Clear the empty lines from otu_file
      otu_file <- otu_file[!apply(is.na(otu_file) | otu_file =="",1,all),]
      # Clear the empty lines from meta_file
      meta_file <- meta_file[rownames(otu_file),]
      # Convert the mapping_column into factor
      meta_file[,mapping_column] <- factor(meta_file[,mapping_column])
      # Form the new levels vector
      levels <- levels(meta_file[,mapping_column])
      levels <- levels[!(levels %in% reference_name)]
    }
    #---------------------------------------------#
    
    
    #---------------- Calculate Distances --------------#
    
    if (tree_or_matrix == "tree"){ 
      # Calculate the UniFrac distance matrix for comparing microbial communities
      unifracs <- GUniFrac(otu_file, rooted_tree, alpha = c(0.0,0.5,1.0))$unifracs
      # Weight on abundant lineages so the distance is not dominated by highly abundant lineages with 0.5 having the best power
      unifract_dist <- unifracs[, , "d_0.5"]
    } else {
      # Read the distances matrix
      unifracs <- read.table (input_tree_or_matrix,check.names = FALSE,header = TRUE,dec = ".",sep = "\t", row.names = 1,comment.char = "")
      # Select only the samples of the mapping file
      unifract_dist <- unifracs[rownames(meta_file) , rownames(meta_file)]
    }
    
    # Number of test groups whose the distances will be computed
    test_groups_number <- length(levels)
    
    
    
    #########################################################################################################################################################
    ################################################################     PAM clustering        ##############################################################
    #########################################################################################################################################################
    
    
    ####################################################
    #------Pam Clustering For the Reference Group------#
    ####################################################
    
    # Check if the number of the reference samples is equal to the number of clusters
    if (sum(meta_file[,mapping_column] %in% reference_name) == reference_clusters) { 
      # Create a matrix that includes only the centroids
      cluster_centroid <- otu_file[meta_file[,mapping_column] %in% reference_name,]
      ref_clusters <- as.factor(rep(1:reference_clusters))
      names(ref_clusters) <- rownames(otu_file[meta_file[,mapping_column] %in% reference_name,])
    } else {
      pam <- pam(unifract_dist[meta_file[,mapping_column] %in% reference_name,meta_file[,mapping_column] %in% reference_name], reference_clusters, diss=TRUE)
      
      # Vector with the indexes of the clusters
      ref_clusters <- as.factor(as.vector(pam$clustering))
      # Naming the elements of the cluster vector
      names(ref_clusters) <- rownames(otu_file[meta_file[,mapping_column] %in% reference_name,])
      
      
      # Calculate the new centroids If the user selected the "mean" or "median" option
      if (central_point != "medoid"){
        
        # Create a matrix that includes only the centroids
        cluster_centroid <- c()
        for (i in 1:reference_clusters) { 
          # Create a matrix that includes only the sample of each cluster will be stored
          cl<-otu_file[pam$clustering==i,]
          if (central_point=="median"){
            # Calculate the new centroids 
            cluster_centroid <-  rbind(cluster_centroid,apply(cl,2,median))
          } else if (central_point=="mean") {
            # Calculate the new centroids 
            cluster_centroid <-  rbind(cluster_centroid,apply(cl,2,mean))
          }
        }
      } else {
        # Create a matrix that includes only the centroids
        cluster_centroid <- otu_file[pam[["medoids"]],]
      }
    }
    #-------------#
    
    
    ########################################################################
    #-------Pam Clustering For Test Groups (1 cluster for each  group)-----#
    ########################################################################
    
    # Each test group will we clustered into 1 cluster
    clusters <- ref_clusters
    # Selecting a medoid for every Test group
    if (Test_name != "None") {
      for (i in 1:length(levels)) {
        # If the test group has only element, this one will be used as the medoid. Otherwise pam clustering will be performed (1 cluster)
        if (length(which(meta_file[,mapping_column]==levels[i])) == 1) {
          # Add the new medoid to the cluster_centroid matrix
          cluster_centroid <- rbind(cluster_centroid,otu_file[which(meta_file[,mapping_column] == levels[i]),])
        } else {
          if (central_point == "medoid"){
            # Perform Pam algorithm
            m <- otu_file[pam(unifract_dist[meta_file[,mapping_column] == levels[i],meta_file[,mapping_column] == levels[i]], 1, diss=TRUE)[["medoids"]],]
            # Add the new medoid to the cluster_centroid matrix
            cluster_centroid <- rbind(cluster_centroid,m)
          } else {
            if (central_point == "median"){
              
              # Add the new median point to the cluster_centroid matrix
              cluster_centroid <-  rbind.data.frame(cluster_centroid,apply(otu_file[which(meta_file[,mapping_column] == levels[i]),],2,median))
            } else if (central_point=="mean"){
              # Add the new mean point to the cluster_centroid matrix
              cluster_centroid <-  rbind.data.frame(cluster_centroid,apply(otu_file[which(meta_file[,mapping_column] == levels[i]),],2,mean))
            }
          }
        }
        # Create a vector with the clusters indices
        test_groups_clusters <- rep(reference_clusters+i,sum(meta_file[,mapping_column] == levels[i]))
        # Naming the elements of the test_groups_clusters vector
        names(test_groups_clusters) <- rownames(otu_file[meta_file[,mapping_column] == levels[i],])
        # Convert cluster vector to factor
        clusters <- as.factor(c(clusters,test_groups_clusters))
      }
    }
    
    if (central_point != "medoid"){
      # Name the rows of the cluster_centroid(IN case of "mean" or "median")
      rownames(cluster_centroid) <- c(paste0("mean",1:nrow(cluster_centroid)))
    }
    
    #------------#
    
    
    # Total number of clusters that have formed from the clustering of the reference and test groups
    # cluster_centroid -> matrix contaning the medoids of the reference and test groups
    clusters_number <- nrow(cluster_centroid)
    
    # Naming the elements of the cluster vector
    # clusters -> vector containing the clusters of the reference and test groups
    clusters <- clusters[rownames(meta_file)]
    
    #--------------------------------------#
    
    
    
    ####################################################################################
    #--------------- Performing De novo clustering to test groups ---------------------#
    ####################################################################################
    
    
    # Vector where the test clusters will be placed 
    # denovo_clusters -> vector with the clusters of the reference and test groups (de novo clustering has been performed to both the reference and test groups)
    denovo_clusters <- c()
    
    if (Test_name != "None" & index == 0) { 
      for (i in 1:length(levels)) {
        if (length(which(meta_file[,mapping_column] == levels[i])) > test_clusters[i]){
          # Perform de novo clustering for every test group
          pam2 <- pam(unifract_dist[meta_file[,mapping_column] == levels[i],meta_file[,mapping_column] == levels[i]], test_clusters[i], diss=TRUE)
          # Store the clusters
          test_groups_clusters2 <- as.factor(as.vector(pam2$clustering))
          # Name denovo_clusters
          names(test_groups_clusters2) <- rownames(otu_file[meta_file[,mapping_column] == levels[i],])
          denovo_clusters <- c(denovo_clusters,test_groups_clusters2)
          
        } else {index <- 1}
      }
    }
    if (index == 0){
      # De novo clusters of the reference and test groups
      denovo_clusters <- c(ref_clusters,denovo_clusters)
      # Name denovo_clusters
      denovo_clusters <- denovo_clusters[rownames(meta_file)]
    }
    
    # Add in the cluster_centroid the mean or median points for every test cluster (that derived from the denovo clustering)
    if (Test_name != "None" & index == 0) { 
      for (i in 1:length(levels)) {
        if (central_point != "medoid"){
          for (j in 1:test_clusters[i]) { 
            # Create a matrix that includes only the sample of each cluster will be stored
            cl <- otu_file[names(which(denovo_clusters[meta_file[,mapping_column] == levels[i]] == j)),]
            if (central_point == "median"){
              # Add the new median point to the cluster_centroid matrix
              cluster_centroid <-  rbind(cluster_centroid,apply(cl,2,median))
              # Name the new row of the cluster_centroid
              rownames(cluster_centroid)[nrow(cluster_centroid)] <- paste(levels[i],j)
            } else if (central_point == "mean") {
              # Add the new mean point to the cluster_centroid matrix
              cluster_centroid <-  rbind(cluster_centroid,apply(cl,2,mean))
              # Name the new row of the cluster_centroid
              rownames(cluster_centroid)[nrow(cluster_centroid)] <- paste(levels[i],j)
              
            }
          }
        }
      }
    }
    
    ##########################################################################################################################################################
    #######################################################           DISTANCES       ########################################################################
    ##########################################################################################################################################################
    
    
    if(central_point != "medoid"){
      # Combine the OTUs table with the samples of the cluster_cenrtoid 
      otu_file2 <- rbind.data.frame(cluster_centroid,otu_file)
      # Calculate the UniFrac distance matrix for comparing microbial communities
      unifracs2 <- GUniFrac(otu_file2, rooted_tree, alpha = c(0.0,0.5,1.0))$unifracs
      # Weight on abundant lineages so the distance is not dominated by highly abundant lineages with 0.5 having the best power
      unifract_dist2 <- unifracs2[, , "d_0.5"]
      
      
      # Initialization of the reference distance matrix
      distances_matrix <- data.frame(matrix(ncol = clusters_number, nrow = nrow(otu_file)))
      # Name the rows of the distances_matrix
      rownames(distances_matrix) <- rownames(otu_file)
      # Placing in the distances_matrix all the distances from each sample to the closest central point 
      for (i in (1:nrow(otu_file))) {
        for (j in 1:clusters_number) {
          dist <- unifract_dist2[rownames(otu_file)[i],rownames(cluster_centroid)[j]]
          distances_matrix[i,j] <- dist
        }
      }
    } else {
      # Initialization of the reference distance matrix
      distances_matrix <- data.frame(matrix(ncol = clusters_number, nrow = nrow(otu_file)))
      # Name the rows of the distances_matrix
      rownames(distances_matrix) <- rownames(otu_file)
      
      
      
      # Placing in the distances_matrix all the distances from each sample to the closest central point 
      for (i in (1:nrow(otu_file))) {
        for (j in 1:clusters_number) {
          dist <- unifract_dist[rownames(otu_file)[i],rownames(cluster_centroid)[j]]
          distances_matrix[i,j] <- dist
        }
      }
    }
    
    # Add cluster vector as a new column in the distances_matrix
    distances_matrix[,(clusters_number+1)] <- c(clusters) 
    
    
    if (clusters_number == 1){
      
      # Find the minimum distance of each row
      minimum_distance <- as.numeric(c(distances_matrix[,1]))
      # Add in the distances_matrix the closest cluster
      distances_matrix[,(clusters_number+2)] <- c(rep(1,length(minimum_distance)))
      
      # Insert that information in the distances_matrix as a new column
      distances_matrix[,(clusters_number+3)] <- c(rep(1,length(minimum_distance)))
      
      # Find the distance to the closest reference central point for every sample
      minimum_reference_distance <- as.numeric(c(distances_matrix[,1]))
      
      # Name the columns of the distances_matrix
      colnames(distances_matrix) <- c(paste("Distance", central_point,1:clusters_number,sep="_"),"Pam cluster","Nearest cluster","Nearest reference cluster")
      
    } else { 
      
      # Find the minimum distance of each row 
      minimum_distance <- as.numeric(apply(distances_matrix[,1:(ncol(distances_matrix)-1)],1,function(x) min(x[1:clusters_number])))
      # Add in the distances_matrix the closest cluster
      for (i in 1:length(minimum_distance)) {
        distances_matrix[i,(clusters_number+2)] <- which(distances_matrix[i,1:clusters_number] == min(distances_matrix[i,1:clusters_number]))[1]
      }
      
      # Find closest reference central point for each Sample
      minimum_reference_cluster <- apply(distances_matrix[,1:(ncol(distances_matrix)-2)],1,function(x) which(x[1:(clusters_number-test_groups_number)] == min(x[1:(clusters_number-test_groups_number)])))
      
      # Insert that information in the distances_matrix as a new column
      distances_matrix[,(clusters_number+3)] <- minimum_reference_cluster
      
      # Find the distance to the closest reference central point for every sample
      minimum_reference_distance <- apply(distances_matrix[,1:(ncol(distances_matrix)-3)],1,function(x) min(x[1:(clusters_number-test_groups_number)]))
      
      # Name the columns of the distances_matrix
      colnames(distances_matrix) <- c(paste("Distance",central_point,1:clusters_number,sep="_"),"Pam cluster","Nearest cluster","Nearest reference cluster")
      
    }
    
    #########################################################################################################################################################
    ############################################      Forming aucillary variables that will     ##############################################################
    ############################################          be used in the following steps       ##############################################################
    #########################################################################################################################################################
    
    
    ###################################################################
    ##################### Make all the necessary Changes ##############
    ###################################################################
    
    
    ###############  Change names of clusters indexes in order to be more readable and informative  #################
    
    # Convert clusters vector to numeric to avoid errors
    clusters <- as.numeric(clusters)
    # Prepare a new vector for the mapping file table(will be the first column)
    mapping_cluster <- as.numeric(clusters)
    
    #----- Change the names in the cluster vector ------#
    
    # Change the names of the reference samples
    for (i in 1:reference_clusters) {
      clusters[clusters == i] <- paste(paste(reference_name,collapse="+"),i)
    }
    
    # Change the names of the test samples
    for (i in reference_clusters+1:(reference_clusters+1+length(levels))) {
      clusters[clusters == i] <- levels[i-reference_clusters]
    }
    
    #----- Change the names in mapping_cluster vector -----#
    
    # Change the names of the reference samples
    if (length(reference_name) == 1){
      for (i in 1:reference_clusters) {
        mapping_cluster[mapping_cluster == i] <- paste(reference_name,i)
      }
    } else {
      for (i in 1:reference_clusters) {
        mapping_cluster[mapping_cluster == i] <- paste(paste(reference_name,collapse="+"),i)
      }
    }
    
    # Change the names of the test samples
    for (i in reference_clusters+1:(reference_clusters+1+length(levels))) {
      mapping_cluster[mapping_cluster == i] <- levels[i-reference_clusters]
    }
    
    #------ Change the names of the samples that derived from the de novo clustering (denovo_clusters) --------#
    
    if (index == 0) {
      for (i in 1:length(levels)) {
        # Change the names of the test samples
        denovo_clusters[meta_file[,mapping_column] == levels[i]] <- paste(levels[i],denovo_clusters[meta_file[,mapping_column] == levels[i]])
      }
      # Change the names of the reference samples
      denovo_clusters[meta_file[,mapping_column] %in% reference_name] <- mapping_cluster[meta_file[,mapping_column] %in% reference_name]
    }
    
    #---------------------------------#
    
    
    
    
    ###################   Mapping file #######################
    
    
    # Form the name of the first column of the mapping file
    column_name <- paste0(mapping_column,"-",paste(reference_name,collapse="+"),"-",reference_clusters)
    # Create the new mapping file
    mapping_file <- cbind(column_name=mapping_cluster,meta_file)
    # Name the columns of the mapping file 
    colnames(mapping_file)[1] <- column_name
    # Create an extra column with the clusters derived from the de novo clustering of the reference ad test groups
    if (index == 0 & Test_name != "None"){ 
      mapping_file <- cbind(column_name=denovo_clusters, mapping_file)
      colnames(mapping_file)[1] <- "De novo clusters"
    }
    
    #-------------#
    
    
    
    ####################  central points matrix #################
    
    
    # Change the row names of Medoid in order to be more readable
    rownames(cluster_centroid) <- paste0(rownames(cluster_centroid),"(",levels(factor(clusters,levels=unique(clusters))),")")
    # Add the taxonomy column
    cluster_centroid <- cbind(t(cluster_centroid),taxonomy)
    
    
    
    ##########################################################################################
    ##########################   Summarized Distance Matrix   ################################
    ##########################################################################################
    
    
    
    # Form a vector containing distances from each sample to the self Medoid (clusters originate from the initial pam clustering)
    self_dist_initial_pam <- as.numeric(apply(distances_matrix , 1 , function(x) x[as.numeric(x[(clusters_number+1)])]))
    
    
    #----------- Perform PAM clustering to reference group (1 cluster) ----------------#
    
    # Selecting a single Medoid for the reference group by performing PAM clustering
    self_pam <- pam(unifract_dist[meta_file[,mapping_column]%in%reference_name,meta_file[,mapping_column]%in% reference_name], 1, diss=TRUE)[["medoids"]]
    
    # Form a vector containing distances to self Medoid (in this case reference group has only one Medoid)
    if (Test_name != "None"){
      self_distance <- c(as.numeric(unifract_dist[self_pam,1:sum(meta_file[,mapping_column] %in% reference_name)]),
                         as.numeric(apply(distances_matrix[(sum(meta_file[,mapping_column] %in% reference_name)+1):length(clusters),],1,function(x) x[as.numeric(x[(clusters_number+1)])])))
      } else{
                           self_distance <- c(as.numeric(unifract_dist[self_pam,1:sum(meta_file[,mapping_column] %in% reference_name)]))
                         }
    #----------------------------------#
    
    
    # Formation of the summarized_distance_matrix. 
    summarized_distance_matrix <- data.frame(clusters,self_dist_initial_pam,self_distance,minimum_distance,as.numeric(minimum_reference_distance))
    # Convert first column into factor
    summarized_distance_matrix[,1] <- as.factor(summarized_distance_matrix[,1])
    
    # Naming the rows of summarized_distance_matrix
    rownames(summarized_distance_matrix) <- rownames(otu_file)
    
    #Naming the columns of summarized_distance_matrix
    colnames(summarized_distance_matrix) <- c("Clusters",paste0("Distances From Self ",central_point,"(",reference_clusters,"clusters",")"),paste0("Distances from Self ", central_point,"(1 cluster)"),"Minimum_distance","Closest Reference Cluster")
    
    
    # Replace zeros with NAs
    summarized_distance_matrix[summarized_distance_matrix == 0] <- NA
    
    
    
    
    ####################################################################################################### 
    ##############################  plot_matrix (will be used later for the plots) ########################
    #######################################################################################################  
    
    
    # plot_matrix -> matrix that will be used in the following plots
    if (index == 0 & Test_name != "None") {
      denovo_clusters <- denovo_clusters[rownames(distances_matrix)]
      plot_matrix <- data.frame(as.factor(c(denovo_clusters)),distances_matrix)
      } else {plot_matrix <- data.frame(as.factor(rep(0,nrow(distances_matrix))),distances_matrix)
      }
    # Add a column with the names of the clusters
    plot_matrix[meta_file[,mapping_column] %in% reference_name,(ncol(plot_matrix)+1)] <- paste(paste(reference_name,collapse="+"),plot_matrix[meta_file[,mapping_column]%in%reference_name,"Nearest.reference.cluster" ])
    
    # Change the names of the distances-based (DB) clusters
    if (Test_name != "None") {
      for (i in 1:length(levels)) {
        plot_matrix[meta_file[,mapping_column] == levels[i],ncol(plot_matrix)] <- paste0(levels[i]," ",plot_matrix[meta_file[,mapping_column] == levels[i],"Nearest.reference.cluster"],"(DB)")}
    }
    
    
    
    ##########################################################################################
    ##########################    Input table(will be used for the     ####################### 
    ##########################   visualization and the Chi square test)  #####################
    ##########################################################################################
    
    
    #------------- Input table(will be used for the visualization and the Chi square test) ------------------#
    
    # Form the data frame that will be used for the visualization 
    input_table <- cbind.data.frame(meta_file[,mapping_column],clusters,self_dist_initial_pam,self_distance,as.numeric(minimum_reference_distance))
    # Name the columns of the input_table
    colnames(input_table) <- c("Category","Clusters",paste("Distances From Self Medoid","(",reference_clusters,"clusters",")"),"Distances From Self Medoid (1 cluster)","Closest Reference Cluster")
    # Name the rows of the input_table
    rownames(input_table) <- rownames(meta_file)
    
    # Replace zeros with NAs
    input_table[input_table == 0] <- NA
    
    # If the reference_name variable has more than 1 elements then these inputs will be combined and form the new name
    input_table[,1] <- as.character(input_table[,1])
    if (length(reference_name) >= 2) {
      input_table[1:length(which(input_table[,1] %in% reference_name)),1] <- c(rep(paste(reference_name,collapse="+"),length(which(input_table[,1] %in% reference_name))))
    }
    
    
    #---------------------#
    
    
    
    ###############################################################  
    ###################    Colours Matrix  ########################
    ###############################################################
    
    #-------------------------------------------------------------------#
    # Create the data frame with the colours of the different groups
    # These colours will be used for all the following plots
    #-------------------------------------------------------------------#
    
    
    # Colours for the reference clusters
    colour_groups <- levels(factor(plot_matrix[meta_file[,mapping_column]%in%reference_name,ncol(plot_matrix)], levels=unique(plot_matrix[meta_file[,mapping_column]%in%reference_name,ncol(plot_matrix)])))
    # Choose the colours from the 'Set2' palette
    colour_vector <- c(rep(brewer.pal(n = (length(Test_name)+1), name = "Set2")[1],length(colour_groups)))
    
    
    
    
    # Add the colours for the test groups
    if (Test_name != "None"){
      for (i in 1:length(levels)){
        # Distances-based groups
        colour_groups <- c(colour_groups,levels(factor(plot_matrix[meta_file[,mapping_column]%in%levels[i],ncol(plot_matrix)], levels=unique(plot_matrix[meta_file[,mapping_column]%in%levels[i],ncol(plot_matrix)])))) 
        if (index == 0){
          # De novo clustering groups
          colour_groups <- c(colour_groups,levels(factor(plot_matrix[meta_file[,mapping_column]%in%levels[i],1], levels=unique(plot_matrix[meta_file[,mapping_column]%in%levels[i],1]))))}
        # colours for the distances-based groups
        colour_vector <- c(colour_vector,c(rep(brewer.pal(n = (length(Test_name)+1), name = "Set2")[i+1],length(levels(factor(plot_matrix[meta_file[,mapping_column]%in%levels[i],ncol(plot_matrix)], levels=unique(plot_matrix[meta_file[,mapping_column]%in%levels[i],ncol(plot_matrix)])))))))
        if (index == 0){
          # colours for the de novo clustering groups
          colour_vector <- c(colour_vector,c(rep(brewer.pal(n = (length(Test_name)+1), name = "Set2")[i+1],length(levels(factor(plot_matrix[meta_file[,mapping_column]%in%levels[i],1], levels=unique(plot_matrix[meta_file[,mapping_column]%in%levels[i],1])))))))}
      }
      
      for(i in 1:length(levels)) {
        # Add the test groups
        colour_groups <- c(colour_groups,levels[i])
        # Colour for the test groups
        colour_vector <- c(colour_vector,brewer.pal(n = (length(Test_name)+1), name = "Set2")[i+1])
      }
    }
    
    
    # Add the Reference group (1 group)
    colour_groups <- c(colour_groups,paste(reference_name,collapse="+"))
    # Colour for the Reference group (1 group)
    colour_vector <- c(colour_vector,brewer.pal(n = (length(Test_name)+1), name = "Set2")[1])
    
    
    # Add the Reference groups when exist more than one 
    if (length(reference_name) > 1){
      for(i in 1:length(reference_name)){
        colour_groups <- c(colour_groups,reference_name[i])
        colour_vector <- c(colour_vector,brewer.pal(n = (length(Test_name)+1), name = "Set2")[1])
      }
    }
    
    
    # Form the final data frame with the colours
    colour_matrix <- data.frame(colour_groups,colour_vector)
    rownames(colour_matrix) <- colour_matrix[,1]
    
    
    
    
    
    #########################################################################################################################################################
    ##############################################                   Descriptive Statistics  and           ##################################################
    ##############################################             Statistical analysis of the clusters        ##################################################
    #########################################################################################################################################################
    
    
    
    
    ###########################################################
    ################### Descriptive Statistics ################
    ###########################################################
    
    
    
    #------------------ Reference Group ---------------------#
    
    
    # Index whether the Descriptive statistics will be calculated 
    stats <- "NO"
    
    # Statistics will be calculated in case that the reference sample are more than the clusters
    if (sum(meta_file[,mapping_column] %in% reference_name) > 3*reference_clusters) {
      
      # Index whether the Descriptive statistics will be calculated 
      stats <- "YES"
      
      # Ordering all the distances
      ordered_distance <- data.frame(minimum_reference_distance[order(minimum_reference_distance)])
      
      # Replace zeros with NA in minimum_reference_distance
      minimum_reference_distance[as.numeric(minimum_reference_distance) == 0] <- NA
      
      # Calculate the mean of the Reference Samples
      mean <- mean(as.numeric(minimum_reference_distance[meta_file[,mapping_column] %in% reference_name]),na.rm = TRUE)
      # Calculate the median of the Reference Samples
      median <- median(as.numeric(minimum_reference_distance[meta_file[,mapping_column] %in% reference_name]),na.rm = TRUE)
      # Calculate the range of the Reference Samples
      range <- as.numeric(max(minimum_reference_distance[meta_file[,mapping_column] %in% reference_name],na.rm = TRUE)) - as.numeric(min (minimum_reference_distance[meta_file[,mapping_column] %in% reference_name],na.rm = TRUE))
      # Calculate the maximum distance Samples for every reference cluster
      max <- as.numeric(max(minimum_reference_distance[meta_file[,mapping_column] %in% reference_name],na.rm = TRUE))
      
      # Test whether the distances of Reference Sample from reference medoid are normally distributed 
      shapino <- shapiro.test(as.numeric(minimum_reference_distance[meta_file[,mapping_column] %in% reference_name]))[["p.value"]]
      # if the distances are normally distributed then calculate standard deviation 
      if (shapino > 0.05) {
        sd <- sd(as.numeric(minimum_reference_distance[meta_file[,mapping_column] %in% reference_name]),na.rm = TRUE)
        # Forming the table contains the statistics about Reference samples
        reference_statistic <- cbind(max,mean,median,range,sd)
        colnames(reference_statistic) <- c("Max distance","Mean","Median","Range","SD")
        rownames(reference_statistic) <- paste(reference_name,collapse="+")
      } else {
        # Forming the table contains the statistics about Reference samples without the standard deviation
        reference_statistic <- cbind(max,mean,median,range)
        colnames(reference_statistic) <- c("Max distance","Mean","Median","Range")
        rownames(reference_statistic) <- paste(reference_name,collapse="+")
      }
      
      
      #------------------------- Statistic for each reference cluster separately ---------------------#
      
      # Preparing the vectors
      cluster_means <- c()
      cluster_median <- c()
      cluster_range <- c()
      cluster_max <- c()
      
      for (i in 1:reference_clusters) {
        
        # Calculate the mean of the Reference Samples for every reference cluster
        cluster_means[i] <- mean(as.numeric(minimum_reference_distance[summarized_distance_matrix[,1] == paste(paste(reference_name,collapse="+"),i)]),na.rm = TRUE)
        # Calculate the median of the Reference Samples for every reference cluster
        cluster_median[i] <- median(as.numeric(minimum_reference_distance[summarized_distance_matrix[,1] == paste(paste(reference_name,collapse="+"),i)]),na.rm = TRUE)
        # Calculate the Range of the Reference Samples for every reference cluster
        cluster_range[i] <- as.numeric(max(minimum_reference_distance[summarized_distance_matrix[,1] == paste(paste(reference_name,collapse="+"),i)],na.rm = TRUE)) - as.numeric(min(minimum_reference_distance[summarized_distance_matrix[,1] == paste(paste(reference_name,collapse="+"),i)],na.rm = TRUE))
        # Calculate the maximum distance Samples for every reference cluster
        cluster_max[i] <- as.numeric(max(minimum_reference_distance[summarized_distance_matrix[,1] == paste(paste(reference_name,collapse="+"),i)],na.rm = TRUE))
      }
      
      
      # Forming the table that contains the statistics of every reference cluster(maximum distance, Mean, Median, range of each cluster)
      cluster_statistics <- cbind(cluster_max,cluster_means,cluster_median,cluster_range)
      colnames(cluster_statistics) <- c("Max distance","Mean","Median","Range")
      rownames(cluster_statistics) <- c(paste(paste(reference_name,collapse="+"),1:reference_clusters))
      if (shapino > 0.05) {
        cluster_sd <- c(rep(NA,reference_clusters))
        cluster_statistics <- cbind(cluster_statistics,cluster_sd )
        colnames(cluster_statistics) <- c("Max distance","Mean","Median","Range","SD")
      }
      
      # Print the reference_statistic and cluster_statistics tables
      write_table("DeNoAn- page 4 Table 1.tab",data.frame(rbind(reference_statistic,cluster_statistics)),"table")
      
      
      # Table containing the Descriptive statistics of the reference clusters
      reference_statistis_table <- gtableGrob("Descriptive statistics of the reference clusters (Table 1)",data.frame(rbind(reference_statistic,cluster_statistics)),setDT="YES",size=12)
      
      #-----------------------------#
      
      
      
      #------------------------- Statistic for each test cluster separately ---------------------#
      
      if (index == 0 & Test_name!="None") {
        
        # Preparing the vectors
        all_cluster_means <- c()
        all_cluster_median <- c()
        all_cluster_range <- c()
        all_cluster_max <- c()
        all_names <- c()
        
        names(self_distance) <- rownames(meta_file)
        for (i in 1:length(levels)){
          for(j in 1:test_clusters[i]){
            
            # Calculate the mean of the Reference Samples for every reference cluster
            all_cluster_means <- c(all_cluster_means,mean(self_distance[denovo_clusters == paste(levels[i],j)],na.rm = TRUE))
            # Calculate the median of the Reference Samples for every reference cluster
            all_cluster_median <- c(all_cluster_median, median(self_distance[denovo_clusters == paste(levels[i],j)],na.rm = TRUE))
            # Calculate the Range of the Reference Samples for every reference cluster
            all_cluster_range <-  c(all_cluster_range, max(self_distance[denovo_clusters == paste(levels[i],j)],na.rm = TRUE)-min(self_distance[denovo_clusters==paste(levels[i],j)],na.rm = TRUE))
            # Calculate the maximum distance Samples for every reference cluster
            all_cluster_max <- c(all_cluster_max, max(self_distance[denovo_clusters == paste(levels[i],j)],na.rm = TRUE))
            # Save the names of the clusters
            all_names <- c(all_names,as.character(paste(levels[i],j)))
          }
        }
        
        # Forming the table that contains the statistics of every reference cluster(maximum distance, Mean, Median, range of each cluster)
        all_cluster_statistics <- cbind(all_cluster_max,all_cluster_means,all_cluster_median,all_cluster_range)
        colnames(all_cluster_statistics) <- c("Max distance","Mean","Median","Range")
        rownames(all_cluster_statistics) <- all_names
        
        # Print the all_cluster_statistics table
        write_table("DeNoAn- page 4 Table 2.tab",data.frame(all_cluster_statistics),"table")
        
        # Table containing the Descriptive statistics of the test clusters
        statistic_table_test_clusters <- gtableGrob("Descriptive statistics of the test clusters (Table 2)",data.frame(all_cluster_statistics),setDT="YES",size=12)
        
      }
      
      #-----------------------------#
      
      
      
      #----------------------------- Order cluster distances ---------------------------#
      
      # Ordering the distances of Reference samples (only) based on the reference clusters they belong
      ordered_ref_clusters <- c()
      for (i in 1:reference_clusters) {
        # Subset only the Reference Samples 
        subset_matrix <- subset(distances_matrix,distances_matrix$`Pam cluster`==i)
        # Order the distance of the above matrix 
        subset_matrix  <- subset_matrix[order(subset_matrix[,i]),][,c(i,(clusters_number+1))]
        colnames(subset_matrix) <- c("Distances","Cluster")
        # Forming the final table with the ordered distances of the Reference Samples
        ordered_ref_clusters <- rbind(ordered_ref_clusters,subset_matrix)
      }
      
      # Ordering the distances of all samples based on the reference clusters they belong
      if (Test_name != "None") {
        ordered_clusters <- c()
        for (i in 1:reference_clusters) {
          # Subset all the Samples 
          subset_matrix2 <- subset(distances_matrix,distances_matrix$`Nearest reference cluster`==i)
          # Order the distance of the above matrix 
          subset_matrix2  <- subset_matrix2[order(subset_matrix2[,i]),][,c(i,(clusters_number+3))]
          colnames(subset_matrix2) <- c("Distances","Cluster")
          # Forming the final table with the orderd distances of the Reference Samples
          ordered_clusters <- rbind(ordered_clusters,subset_matrix2)
        }
      }
      
      #---------------------------------#
      
    }  
    
    
    
    # Create a new vector that contains "YES" and "NO" 
    # YES indicates that the sample compared to all the other reference samples(except outliers) of the same cluster, is the most distant object 
    # NO will indicate the opposite 
    if (Test_name != "None") {
      maximum <- c()
      
      # find the outliers and exclude them from the vector and then calculate the maximum distance of each reference cluster
      for (i in 1:reference_clusters) {
        all <- subset(distances_matrix[meta_file[,mapping_column] %in% reference_name,i],distances_matrix[meta_file[,mapping_column] %in% reference_name,(clusters_number+1)] == i)
        outliers <-  boxplot.stats(subset(distances_matrix[meta_file[,mapping_column] %in% reference_name,i],distances_matrix[meta_file[,mapping_column] %in% reference_name,(clusters_number+1)] == i))$out
        maximum[i] <- max(all[!(all %in% outliers)])
      }
      # Compare the distances with the number found above
      ref_dist <- apply(distances_matrix[meta_file[,mapping_column] != reference_name,],1,function(x) {a <- as.numeric(x[(clusters_number+3)])
      if (min(x[1:reference_clusters])<maximum[a]) {"NO"} else {"YES"}})
      
      # For all the reference sample "NO" will be placed
      reference <- c(rep("NO", sum(meta_file[,mapping_column] %in% reference_name)))
      names(reference) <- rownames(distances_matrix[meta_file[,mapping_column] %in% reference_name,])
      temporary_vector <- c(reference,ref_dist)
      temporary_vector <- temporary_vector[rownames(distances_matrix)]
      # Add "YES" and "NO" as a new column to the distances_matrix
      distances_matrix[,(clusters_number+4)] <- temporary_vector
      # Name the columns of the distances_matrix
      colnames(distances_matrix) <- c(paste("Distance Medoid",1:clusters_number,sep="_"),"Pam cluster","Nearest cluster","Nearest reference cluster","Most distant object" )
    }
    
    
    
    
    median_distances <- data.frame()
    for(j in 1:reference_clusters){
      median_distances[1,j] <- median(distances_matrix[rownames(meta_file[meta_file[,mapping_column] %in% reference_name,]),j])
      rownames(median_distances) <- paste(reference_name,collapse="+")
    }
    
    if (length(Test_name) > 0 & Test_name != "None"){
      for(i in 1:length(Test_name)){
        for(j in 1:reference_clusters){
          median_distances[1+i,j] <- median(distances_matrix[rownames(meta_file[meta_file[,mapping_column] %in% Test_name[i],]),j])
        }
      }
      rownames(median_distances) <- c(paste(reference_name,collapse="+"),Test_name)
    }
    
    min_dist <- apply(median_distances,1,min)
    for (i in 1:nrow(median_distances)){
      median_distances[i,(reference_clusters+1)] <- which(median_distances[1,] == min_dist[1])
    }
    
    colnames(median_distances) <- c(paste(paste(reference_name,collapse="+"),1:reference_clusters),"Closest Ref Cluster")
    
    # Write the Chi_square_matrix matrix
    write_table("DiBaAn-page 7 Median Distances.tab",median_distances,"table")
    
    # Table containing the observations per Group
    median_table <- gtableGrob("Median distances from reference clusters",median_distances,setDT="YES",size=15)
    
    
    
    ###################################################################################
    ##################### Statistical analysis of the clusters  #######################
    #####################              and groups               #######################
    ###################################################################################
    
    
    
    
    ###################### Chi square test ###########################
    
    
    if (length(reference_name) >= 2) {
      reference_name <- c(paste(reference_name,collapse="+"))
    }
    
    if (reference_clusters > 1) {
      # Convert independent variable into factor to avoid errors
      input_table[,1] <- as.factor(input_table[,1]) 
      input_table[,1] <- factor(input_table[,1],levels=unique(input_table[,1]))
      
      # Groups that will be used for the Chi square test
      chi_levels <- levels(input_table[,1])
      
      # Initialization of the Chi_square_matrix
      Chi_square_matrix <- data.frame(matrix(ncol = (length(chi_levels)) , nrow = reference_clusters))
      # Naming the row and the columns of Chi_square_matrix
      colnames(Chi_square_matrix) <- c(chi_levels)
      rownames(Chi_square_matrix) <- c(paste(paste(reference_name,collapse="+"),1:reference_clusters,sep=" "))
      
      # Initialization of the chi_names vector
      # Form a matrix where the number of observations for each cluster will be added
      for (l in 1:length(chi_levels)){
        for (c in 1:reference_clusters){
          Chi_square_matrix[c,(l)] <-  length(which(distances_matrix[input_table[,1]==chi_levels[l],"Nearest reference cluster"] == c))
        }
      }
      
      #-----####-------####--------###--- Chi square test ---#####----####----####------####----------#
      
      
      #--------------------- Chi square independence testing for reference Cluster(expected vs observed) --------------------------#
      
      # Conduct Chi square test
      exp_p_value <- data.frame(round(chisq.test(Chi_square_matrix[,reference_name],p=c(rep(1/reference_clusters,reference_clusters)))$p.value,4))
      
      # Forming the data frame containing the p values
      colnames(exp_p_value) <- c("p value")
      # Naming the rows of the exp_p_value data frame
      rownames(exp_p_value) <- c("expected vs observed")
      
      
      #--------------------- Chi square goodness of fit testing for Groups --------------------------#
      if (Test_name != "None") {
        combs <- combs(levels(input_table[,1]), 2)
        # Save the names of that combinations
        group_names <- apply(combs,1,function(x)paste0(x[1],"-",x[2]))
        
        
        # P_value -> a vector where the p values will be stored
        chi_p_value <- c(); test <- c()
        for(chi in 1:nrow(combs)) { 
          # Performed the Chi square test
          chi_p_value[chi] <- round(chisq.test(Chi_square_matrix[,combs[chi,2]],p=Chi_square_matrix[,combs[chi,1]],rescale.p = TRUE)$p.value,4)
          }
        
        
        # Forming the data frame containing the p values
        p_value <- data.frame(group_names)
        # Naming the rows of the p_value data frame
        p_value[,2] <- chi_p_value
        # Naming the columns of p_value data frame
        colnames(p_value) <- c("Groups","p value")
        
        
        # Convert p-values column into numeric
        p_value[,2] <- as.numeric(p_value[,2])
        # Order the p-values
        p_value <- p_value[order(p_value[,2]),]
        # Choose the first 10 most significant p-values
        if (nrow(p_value)>10) {p_value2 <- p_value[1:10,]} else {p_value2 <- p_value}
      }
      
      
      #--------------------------------- Prepare the tables for the pdf------------------------#
      
      if (chi_square == 1){
        # Write the Chi_square_matrix matrix
        write_table("DiBaAn-page 9 Table 1.tab",Chi_square_matrix,"table")
        
      }else{
        # Write the Chi_square_matrix matrix
        write_table("DiBaAn-page 8 Table 1.tab",Chi_square_matrix,"table")
      }
      
      # Table containing the observations per Group
      DB_matrix_table <- gtableGrob("Distances based (DB) clustering (Table 1) ",Chi_square_matrix,setDT="YES",size=12)
      
      
      
      ######################################################################################################################################################
      ##########################################################  Visualization ############################################################################  
      ######################################################################################################################################################
      
      
      #---------------------------------- Calculate percentages of each cluster -------------------------------#
      
      
      # Convert Chi_square_matrix into data frame
      Chi_square_matrix <- data.frame(Chi_square_matrix)
      
      
      # Forming the data frame where the percentages will be saved
      percentage <- data.frame(c(paste("Cluster",1:reference_clusters,sep="")))
      if (ncol(Chi_square_matrix) == 1){
        colsum <- sum (Chi_square_matrix[1:nrow(Chi_square_matrix),1])
      } else {
        colsum <- colSums(Chi_square_matrix[,1:ncol(Chi_square_matrix)])
      }
      
      for(p in 1:ncol(Chi_square_matrix)) {
        # Calculating the percentages
        percentage <- data.frame(percentage,round((Chi_square_matrix[,(p)]/colsum[p])*100,1))}
      # Name the columns
      colnames(percentage) <- c("Clusters",chi_levels)
      
      
      # Transformation of the percentages data frame in order to be plotted
      barplot_percentage <- data.frame(pivot_longer(percentage, cols = 2:ncol(percentage), names_to = "GROUP"))
      # Naming the columns
      colnames(barplot_percentage) <- c("Clusters","GROUP","Percentages")
      
      # Create a vector with the colours for the barplot
      barplot_percentage[,"GROUP"] <- as.factor(barplot_percentage[,"GROUP"])
      color <- colour_matrix[levels(barplot_percentage[,"GROUP"]),2]
      
      # Barplot of the percentages
      barplot <- ggplot(data=barplot_percentage, aes(x=Clusters, y=Percentages, fill=GROUP)) +
        geom_bar(stat="identity", position=position_dodge())+
        geom_text(aes(label=Percentages), vjust=1.6,position = position_dodge(0.9), size=3.5)+
        ylab("Percentage (%)")+
        scale_fill_manual(breaks=c(levels(as.factor(barplot_percentage[,"GROUP"]))), values=color)+
        theme_classic()
      
      
      # Table containing the percentages of each cluster
      percentage_table <- gtableGrob("",percentage,setDT="NO",size=18)
      
      
      
      #--------------------------------- Prepare the tables for the pdf------------------------#
      
      if (Test_name!="None") {
        
        if (chi_square == 1){
          
          # Write the Chi_square_matrix
          write_table("DiBaAn-page 8 Table 2.tab",p_value,"table")
          
        }else{
          
          # Write the Chi_square_matrix
          write_table("DiBaAn-page 7 Table 2.tab",p_value,"table")
        }
        
        
        # Table containing p values of groups independence testing
        chisquare_pvalues <- gtableGrob("P-values of Chi square test\n Goodness of fit test- (Table 2)",p_value2,setDT="NO",size=12)
        
      }
    }
    
    
    
    
    if (Test_name != "None") {
      
      ################################################################################################################################
      #-------------------------- Boxplots presenting the differnces between the refernence clusters  -------------------------------#
      #--------------------------     and the samples that are closer to each of these clusters       -------------------------------#
      #--------------------------                            (plot 1)                                 -------------------------------#
      ################################################################################################################################
      
      
      # Subset the Reference samples from the plot_matrix
      cluster_dist_ref <- plot_matrix[rownames(input_table[which(input_table[,1] %in% reference_name),]),]
      # Subset the Test samples from the plot_matrix
      cluster_dist_test <- plot_matrix[rownames(input_table[which(input_table[,1] %in% levels),]),]
      
      
      # dist_plot_1 -> matrix with the clusters that have been formed based on the the distances from the closest reference medoid (DB) 
      dist_plot_1 <- c()
      
      for (j in 1:reference_clusters) {
        # dist1 -> aucillary matrix with the distances of the reference samples from the closest reference medoid
        dist1 <- cluster_dist_ref[cluster_dist_ref$Nearest.reference.cluster == j,c(ncol(plot_matrix),(j+1))]
        # Create a column with the index of the closest reference cluster
        dist1[,3] <- c(rep(j,nrow(dist1)))
        # Name the columns
        colnames(dist1) <- c("names","distances","cluster")
        
        # dist3 -> aucillary matrix with the distances of the reference samples and test samples
        dist3 <- c()
        for (k in 1:length(levels)){
          # Form a temporary matrix with the distances of each test group
          temporary_matrix <- plot_matrix[input_table[,1] == levels[k],]
          # dist2 -> aucillary matrix with the distances of the test samples from the closest reference medoid
          dist2 <- temporary_matrix[temporary_matrix$Nearest.reference.cluster == j,c(ncol(plot_matrix),(j+1))]
          # Create a column with the index of the closest reference cluster
          dist2[,3] <- c(rep(j,nrow(dist2)))
          # Name the columns
          colnames(dist2) <- c("names","distances","cluster")
          dist3 <- rbind(dist3,dist2)
        }
        
        # Final dist_plot_1 matrix
        dist_plot_1 <- rbind(dist_plot_1,dist1,dist3)
      }
      
      
      # Convert first column into character
      dist_plot_1[,1] <- as.character(dist_plot_1[,1])
      
      # Replace zeros with NAs
      dist_plot_1[dist_plot_1 == 0] <- NA
      
      # Calculate and add the prevalence of the groups
      dist_plot_1 <- data.frame(dist_plot_1,label=prevalence(dist_plot_1,column=2))
      
      # Convert label column into factor
      dist_plot_1$label <- factor(dist_plot_1$label,levels = unique(dist_plot_1$label))
      
      # Create a new data frame that will be used for the graphs
      plot_df_a <- cbind.data.frame(Group=as.factor(dist_plot_1[,1]),abundance = as.numeric(dist_plot_1[,2]),label = factor(dist_plot_1[,4]))
      
      # Store the colours that will be used in the plot
      color <- colour_matrix[levels(factor(plot_df_a[,1],levels = unique(plot_df_a[,1]))),2]
      
      
      #------------------- Wilcoxon Rank Sum Test p.values  ----------------------#
      
      
      # Calculate the p-values of the Wilcoxon Rank sum test 
      pvalues2 <- pvalues_function2(dist_plot_1,3,"DiBaAn-",10)
      
      
      # gtable of the p-values
      pvaltable <- gtableGrob("Wilcoxon Rank Sum Test - pairwise",pvalues2,setDT="NO",size=9)
      
      #------------------#
      
      
      #------------------------ Create the Plots -------------------------#
      
      # The name that will be used as the title for the plots 
      my_name <- paste0("Closest Reference Cluster")
      
      
      # Save the boxplot object in the list
      list_box <- my_boxplot(my_name,plot_df_a,label,abundance,colour=color)
      
      # Save the boxplot with points object in the list
      list_point <- my_point_boxplot(my_name,plot_df_a,label,abundance,color)
      
      # Save the violin object in the list
      list_violin <- my_violinplot(my_name,plot_df_a,label,abundance,color)
      
      
      
      #################################################################################################
      #------------------------- Comparative boxplots of Groups for self-medoids  --------------------# 
      #-------------------------          and minimum reference distance          --------------------#
      #-------------------------                  (plot 2)                        --------------------#
      #################################################################################################
      
      
      # In the end, plots based on Category and Cluster columns will be generated 
      
      
      #---------------------- Making all the necessary objects -------------------- #
      
      # Making an object for lists of boxplot
      list_box2 <- list()
      
      # Making an object for lists of point boxplot
      list_point2 <- list()
      
      # Making an object for lists of violin boxplot
      list_violin2 <- list()
      
      # Making a list for pairwise p-value table of Wilcoxon test
      pvaltable2 <- list()
      
      
      
      #-----------------------------------------------------------------------------#  
      
      # Convert independent variables (1st column) into factor in order to avoid errors
      input_table[,1] <- as.factor(input_table[,1]) 
      input_table[,1] <- factor(input_table[,1],levels=unique(input_table[,1]))
      
      
      # Find all the individual groups of the independent variable
      Levels <- levels(input_table[,1]) 
      
      
      for (j in 4:ncol(input_table)) {
        
        # Delete label column
        input_table$label <- NULL
        
        # Calculate and add the prevalence of the groups
        input_table <- data.frame(input_table,label=prevalence(input_table,column=j))
        
        # Create a new data frame that will be used for the graphs
        plot_df <- cbind.data.frame(Group=as.factor(input_table[,1]),abundance = as.numeric(input_table[,j]),label = factor(input_table[,6], levels = unique(input_table[,6])))
        
        # Store the colours that will be used in the plot
        color <- colour_matrix[levels(factor(plot_df[,1],levels = unique(plot_df[,1]))),2]
        
        
        #------------------- Wilcoxon Rank Sum Test p.values ----------------------#
        
        if (j == 4){pages <- 4} else {pages <- 6}
        
        # Calculate the p-values of the Wilcoxon Rank sum test 
        pvalues2 <- pvalues_function(plot_df,"DiBaAn-",pages)
        
        
        # Create a gtable containing text grobs 
        pvaltable2[[j-3]] <- list()
        pvaltable2[[j-3]] <- gtableGrob("Wilcoxon Rank Sum Test - pairwise",pvalues2,setDT="NO",size=9)
        
        #--------------------#
        
        
        
        #------------------------ Plots -------------------------#
        
        # The name that will be used as the title for the plots 
        my_name <- c("Category","Clusters",paste("Distances From Self Medoid","(",reference_clusters,"clusters",")"),paste0("Distances From Self-",central_point, " (1 cluster)"),"Closest Reference Cluster")[j]
        
        
        
        # Save the boxplot object in the list
        plot_df$label <- factor(plot_df$label,levels = unique(plot_df$label))
        list_box2[[j-3]] <- list()
        list_box2[[j-3]] <- my_boxplot(my_name,plot_df,label,abundance,color)
        
        
        # Save the boxplot with points object in the list
        list_point2[[j-3]] <- list()
        list_point2[[j-3]] <- my_point_boxplot(my_name,plot_df,label,abundance,color)
        
        # Save the violin object in the list
        list_violin2[[j-3]] <- list()
        list_violin2[[j-3]] <- my_violinplot(my_name,plot_df,label,abundance,color)
        
        
      }
      
      #-----------------#
      
      
      
      
      if (index == 0) {
        
        
        #########################################################################################################################
        #--------------------------     Compare the compactness of the cluster that divived from -------------------------------#
        #--------------------------     the de novo clustering of the reference and test groups  -------------------------------#
        #--------------------------                          (plot 3)                            -------------------------------#
        #########################################################################################################################
        
        
        
        # Form the vector with the distances from the closest self medoid
        # de_novo_distances -> Vector where the distances will be placed
        de_novo_distances <- c()
        
        
        # Test subclusters
        l <- levels(as.factor(denovo_clusters[meta_file[,mapping_column] %in% Test_name]))
        for(i in 1:length(l)){
          if (central_point == "medoid"){
            # Perform PAM clustering for every test cluster
            de_novo_pam <-pam(unifract_dist[names(denovo_clusters[denovo_clusters == l[i]]),names(denovo_clusters[denovo_clusters == l[i]])], 1, diss=TRUE)[["medoids"]]
            # Distances from the self medoids
            pam_distances <-as.numeric(unifract_dist[de_novo_pam,c(names(denovo_clusters[denovo_clusters == l[i]]))])
          } else {
            pam_distances <-as.numeric(unifract_dist2[l[i],c(names(denovo_clusters[denovo_clusters == l[i]]))])
            
          }
          # Naming the distances vector
          names(pam_distances) <- names(denovo_clusters[denovo_clusters == l[i]])
          de_novo_distances <- c(de_novo_distances,pam_distances)
        }
        
        # Add the reference subclusters
        denovo_distances <- c(minimum_reference_distance[input_table[,1]%in%reference_name],de_novo_distances)
        denovo_distances <- denovo_distances[names(denovo_clusters)]
        
        
        
        # Form data frame with the de novo clusters and the distances from the closest self-medoid
        dist_plot_3 <- data.frame(denovo_clusters,as.numeric(denovo_distances))
        
        
        # Replace zeros with NAs
        dist_plot_3[dist_plot_3 == 0] <- NA
        
        # Calculate and add the prevalence of the groups
        dist_plot_3 <- data.frame(dist_plot_3,label=prevalence(dist_plot_3,column=2))
        
        # Convert label column into factor
        dist_plot_3$label <- factor(dist_plot_3$label,levels = unique(dist_plot_3$label))
        
        # Create a new data frame that will be used for the graphs
        plot_df <- cbind.data.frame(Group=as.factor(dist_plot_3[,1]),abundance = as.numeric(dist_plot_3[,2]),label = factor(dist_plot_3[,3], levels=unique(dist_plot_3[,3])))
        
        # Store the colours that will be used in the plot
        color <- colour_matrix[levels(factor(plot_df[,1],levels = unique(plot_df[,1]))),2]
        
        
        
        
        #------------------- Wilcoxon Rank Sum Test p.values  ----------------------#
        
        # Calculate the p-values of the Wilcoxon Rank sum test 
        pvalues2 <- pvalues_function(dist_plot_3,"DeNoAn-",5)
        
        # gtable of the p-values
        pvaltable3 <- gtableGrob("Wilcoxon Rank Sum Test - pairwise",pvalues2,setDT="NO",size=9)
        
        #---------------------------#
        
        
        #------------------------  Create the Plots -------------------------#
        
        # The name that will be used as the title for the plots 
        my_name <- paste0("Distances From Self-",central_point)
        
        # Save the boxplot object in the list
        list_box3 <- my_boxplot(my_name,plot_df,label,abundance,color)
        
        # Save the boxplot with points object in the list
        list_point3<- my_point_boxplot(my_name,plot_df,label,abundance,color)
        
        # Save the violin object in the list
        list_violin3 <- my_violinplot(my_name,plot_df,label,abundance,color)
        
        
        
        
        
        
        
        #########################################################################################################################
        #--------------------------    Compare the distances of the de novo based clusters from     ----------------------------#
        #--------------------------    the closest(based on the median distance) reference group   -----------------------------#
        #--------------------------                           (plot 4)                             -----------------------------#
        #########################################################################################################################
        
        
        
        
        ###################### Median and observations matrices ####################################
        
        # median_matrix -> data frame with the median distance of each test subcluster from the reference clusters
        median_matrix <- data.frame()
        
        # matrix_median ->
        matrix_median <- c()
        
        # max_observations -> 
        max_observations <- c()
        
        # observations_matrtix -> matrix with the observations
        observations_matrtix <- data.frame()
        
        for (i in 1:length(levels)){
          # Subset from the plot_matrix the samples of each test group
          temporary_matrix <- plot_matrix[meta_file[,mapping_column] == levels[i],]
          # Vector with the factors of the temporary_matrix
          levels_matrix <- as.factor(unique(temporary_matrix[,1]))
          
          for(j in 1:length(levels_matrix)){
            for (l in 1:reference_clusters){
              # Find the median distance of ecsh test subcluster from every reference cluster
              matrix_median[l] <- round(median(temporary_matrix[temporary_matrix[,1] == levels_matrix[j],(l+1)]),3)
              # Find the number of the samples of every test cluster that are closer to each reference cluster
              max_observations[l] <- sum(temporary_matrix[temporary_matrix[,1] == levels_matrix[j],"Nearest.reference.cluster"] == l, na.rm = TRUE)
            }
            # Find the minimun median distance
            matrix_median[(reference_clusters+1)] <- which(matrix_median==min(matrix_median))
            median_matrix <- rbind(median_matrix,matrix_median)
            # Name the rows of the median_matrix
            rownames(median_matrix)[nrow(median_matrix)] <-as.character(levels_matrix[j])
            # Name the columns of the median_matrix
            colnames(median_matrix) <- c(paste(paste(reference_name,collapse="+"),1:reference_clusters,sep=" "),"Closest Ref cluster")
            
            # Create the matrix with the observations
            observations_matrtix <- rbind(observations_matrtix,max_observations)
            # Name the rows of the observations_matrtix
            rownames(observations_matrtix)[nrow(observations_matrtix)] <-as.character(levels_matrix[j])
            # Name the columns of the observations_matrtix
            colnames(observations_matrtix) <- c(paste(paste(reference_name,collapse="+"),1:reference_clusters,sep=" "))
          }
        }
        
        # Final matrix containing the median distances of the test clusters from the reference clusters
        median_matrix <- data.frame(rownames(median_matrix),median_matrix)
        # Naming the columns of the median_matrix
        colnames(median_matrix) <- c("",c(paste(paste(reference_name,collapse="+"),1:reference_clusters,sep=" "),"Closest Ref cluster"))
        
        
        
        
        # dist_plot_4 ->  matrix showing the closest test cluster for every reference cluster based on the median distances
        dist_plot_4 <- c()
        
        for (j in 1:reference_clusters) {
          # Forming for every cluster the appropriate distance matrix
          which <- which(median_matrix[,ncol(median_matrix)] == j)
          
          # dist1 -> aucillary matrix for the test groups showing the closest reference cluster
          dist1 <- cluster_dist_test[cluster_dist_test[,1]%in% rownames(median_matrix)[which],c(1,(j+1))]
          dist1[,3] <- c(rep(j,nrow(dist1)))
          colnames(dist1) <- c("names","distances","cluster")
          
          # dist1 -> aucillary matrix for the test groups
          dist2 <- cluster_dist_ref[cluster_dist_ref$Nearest.reference.cluster == j,c(1,(j+1))]
          dist2[,1] <- c(rep(paste(paste(reference_name,collapse="+"),j),nrow(dist2)))
          dist2[,3] <- c(rep(j,nrow(dist2)))
          
          colnames(dist2) <- c("names","distances","cluster")
          
          #Form the final dist_plot_4 matrix
          dist_plot_4 <- rbind(dist_plot_4,dist2,dist1)
        }
        
        # Convert first column into character
        dist_plot_4[,1] <- as.character(dist_plot_4[,1])
        
        # Replace zeros with NAs
        dist_plot_4[dist_plot_4 == 0] <- NA
        
        # Calculate and add the prevalence of the groups
        dist_plot_4 <- data.frame(dist_plot_4,label=prevalence(dist_plot_4,column=2))
        
        # Convert label column into factor
        dist_plot_4$label <- factor(dist_plot_4$label,levels = unique(dist_plot_4$label))
        
        
        
        # Create a new data frame that will be used for the graphs
        plot_df_b <- cbind.data.frame(Group=as.factor(dist_plot_4[,1]),abundance = as.numeric(dist_plot_4[,2]),label = factor(dist_plot_4[,4]))
        
        # Store the colours that will be used in the plot
        color <- colour_matrix[levels(factor(plot_df_b[,1],levels = unique(plot_df_b[,1]))),2]
        
        
        #------------------- Wilcoxon Rank Sum Test p.values  ----------------------#
        
        # Calculate the p-values of the Wilcoxon Rank sum test 
        pvalues2 <- pvalues_function2(dist_plot_4,3,"DeNoAn-",9)
        
        
        # gtable of the p-values
        pvaltable4 <- gtableGrob("Wilcoxon Rank Sum Test - pairwise",pvalues2,setDT="NO",size=9)
        
        #--------------------#
        
        
        
        #------------------------ Create the Plots -------------------------#
        
        # The name that will be used as the title for the plots 
        my_name <- paste0("Closest reference cluster")
        
        # Save the boxplot object in the list
        list_box4 <- my_boxplot(my_name,plot_df_b,label,abundance,color)
        
        # Save the boxplot with points object in the list
        list_point4 <- my_point_boxplot(my_name,plot_df_b,label,abundance,color)
        
        # Save the violin object in the list
        list_violin4 <- my_violinplot(my_name,plot_df_b,label,abundance,color)
        
        
        
        #--------Create a list of text--------------#
        
        if (chi_square == 1){
          write_table(paste("DeNoAn-page", (length(exploratory_columns)*3+7), "Table.tab"),observations_matrtix,"table")
        }else{
          write_table("DeNoAn-page 6 Table.tab",observations_matrtix,"table")}
        
        # Create a gtable containing text grobs of the observations_matrtix (will be used later)
        observation <- gtableGrob("Closest Reference cluster",observations_matrtix,setDT="YES",size=15)
        
        if(chi_square == 1){
          write_table(paste("DeNoAn-page",(length(exploratory_columns)*3+7) , "Table.tab"),median_matrix,"table")
        }else{
          write_table("DeNoAn-page 8 Table.tab",median_matrix,"table")}
        
        # Create a gtable containing text grobs of the median_matrix (will be used later)
        median_matrix_table <- gtableGrob("Median distances from reference clusters",median_matrix,setDT="NO",size=15)
        
        
        
        
        ######################################################################################################################
        #-------------------------- Compare the distances of the de novo based clusters from  -------------------------------#
        #--------------------------              the closest reference medoid                 -------------------------------#
        #--------------------------                         (plot 5)                          -------------------------------#
        ######################################################################################################################      
        
        
        # Form the matrix with the distances of the samples from the closest reference medois
        reference_distances <- summarized_distance_matrix[,c(1,5)]
        # Reorder the rows of the reference_distances matrix
        reference_distances <- reference_distances[names(denovo_clusters),]
        # Add in the first column the de novo based clusters
        reference_distances[,1] <- as.factor(c(denovo_clusters))
        # Naming the rows 
        rownames(reference_distances) <- rownames(summarized_distance_matrix)
        
        # Subset the Reference samples from the distances matrix
        cluster_dist_ref_b <- summarized_distance_matrix[rownames(input_table[which(input_table[,1] %in% reference_name),]),c(1,5)]
        
        
        
        # Forming for every cluster the appropriate distance matrix
        dist_plot_5 <- rbind(cluster_dist_ref_b ,reference_distances[rownames(input_table[which(input_table[,1 ]%in% levels),]),])
        
        # Convert first column into character
        dist_plot_5[,1] <- as.character(dist_plot_5[,1])
        
        # Calculate and add the prevalence of the groups
        dist_plot_5 <- data.frame(dist_plot_5,label=prevalence(dist_plot_5,column=2))
        
        #-----#
        
        # Convert label column into factor
        dist_plot_5$label <- factor(dist_plot_5$label,levels = unique(dist_plot_5$label))
        
        
        # Create a new data frame that will be used for the graphs
        plot_df <- cbind.data.frame(Group=as.factor(dist_plot_5[,1]),abundance = as.numeric(dist_plot_5[,2]),label = factor(dist_plot_5[,3]))
        
        # Store the colours that will be used in the plot
        color <- colour_matrix[levels(factor(plot_df[,1],levels = unique(plot_df[,1]))),2]
        
        
        #------------------- Wilcoxon Rank Sum Test p.values  ----------------------#
        
        
        # Calculate the p-values of the Wilcoxon Rank sum test 
        pvalues2 <- pvalues_function(dist_plot_5,"DeNoAn-",7)
        
        # Create a gtable containing text grobs 
        pvaltable5 <- gtableGrob("Wilcoxon Rank Sum Test - pairwise",pvalues2,setDT="NO",size=9)
        
        #------------------#
        
        
        #------------------------ Plots -------------------------#
        
        # The name that will be used as the title for the plots 
        my_name <- paste("Closest Reference",central_point)
        
        # Save the boxplot object in the list
        list_box5 <- my_boxplot(my_name,plot_df,label,abundance,color)
        
        # Save the boxplot with points object in the list
        list_point5 <- my_point_boxplot(my_name,plot_df,label,abundance,color)
        
        # Save the violin object in the list
        list_violin5 <- my_violinplot(my_name,plot_df,label,abundance,color)
        
      }
    }
    
    
    #########################################################################################################################################################
    ###############################################           Write Output Files in Separate folder        ##################################################
    #########################################################################################################################################################
    
    
    
    
    ########################################################
    #################  Print the plots   ###################
    ########################################################
    
    
    
    setwd(results_path)
    
    
    
    if (reference_clusters>1) {
      
      # Print the MDS plot only for the reference clusters
      png_plot(paste0("MDS reference clusters k=",reference_clusters,".png"),
               sclass(unifract_dist[input_table[,1] %in% reference_name,input_table[,1]%in%reference_name],factor(mapping_cluster[which(input_table[,1] %in% reference_name)], levels=unique(mapping_cluster[which(input_table[,1] %in% reference_name)])),
                      NULL, colour_matrix[levels(factor(mapping_cluster[which(input_table[,1] %in% reference_name)],levels = unique(mapping_cluster[which(input_table[,1] %in% reference_name)]))),2]))
      
    }
    
    if (Test_name != "None"){
      
      setwd(plots_path)
      
      # Print the MDS plot where all the test samples are represented as individual points
      pdf("DiBaAn-page 5 MDS Plot.pdf")
      
      layout(matrix(c(1, 1, 1,
                      2, 2, 2,
                      3, 3, 3), nrow=3, byrow=TRUE),heights=c(0.3,3,0.3) )
      
      par(mar = c(0.1, 0.1, 0.1, 0.1))
      
      plot.new()
      
      
      sclass(unifract_dist,as.factor(mapping_cluster[which(input_table[,1] %in% reference_name)]),as.factor(mapping_cluster[which(input_table[,1] %in% levels)]),
             colour_matrix[levels(factor(mapping_cluster[which(input_table[,1] %in% reference_name)],levels = unique(mapping_cluster[which(input_table[,1] %in% reference_name)]))),2] )
      par(mar = c(0.1, 0.1, 0.1, 0.1))
      
      plot.new()
      
      
      dev.off()
      
      
      
      # Return to results path
      setwd(outputs_path)
    }
    
    
    
    if (Test_name != "None") {
      
      
      ########################################################################################
      #################################      REPORTS         #################################
      ########################################################################################
      
      
      
      
      #################################################
      ########## Distsances-based report ##############
      #################################################
      
      
      # Initialize the list that  will later contain the plots
      boxplot_list <- list()
      pointplot_list <- list()
      violinplot_list <- list()
      
      # Open a PDF file to store the outputs
      pdf("Distances-based report.pdf")
      
      
      
      ######### Cover page #############
      
      # Index indicating the page of the report
      page <- 1
      
      # Set the layout of the cover page
      layout( matrix(c(1,1), ncol=1))
      # Default margins
      par(mar= c(5.1, 4.1, 8, 2.1))
      
      # Create an empty plot 
      plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      # Title of the plot
      title(main=paste0("Distances Based Analysis\n","(DiBaAn-",central_point,")"),cex.main = 1.75,font.main = 2)
      # Description of this PDF
      text(x = 0.5, y = 0.5, paste0("This report includes the results derived from the analysis of the distances\n",
                                    "of the test samples from their closest reference ", central_point,".\n",
                                    "The Reference samples were clustered into ",reference_clusters, " distinct groups based on\n",
                                    "the results of the PAM algorithm. The remaining samples of the test \n groups ",
                                    "were compared with the ",central_point,"s of the reference clusters\n",
                                    "\n\n","* All the elements of this report (plots and tables) will be found in the \n results folder with the prefix DiBaAn-"),cex = 1.1, col = "black",font = 3)
      
      # Default margins
      par(mar= c(5.1, 4.1, 4.1, 2.1))
      
      
      
      ########## Page 2 ############
      
      # Index indicating the page of the report
      page <- page+1
      
      # Set the layout of the cover page (3 rows, 1 column)
      layout(matrix(c(1, 1, 1,
                      2, 2, 2,
                      3, 3, 3), nrow=3, byrow=TRUE),heights=c(1,3,0.3))
      
      # Create an empty plot
      plot.new()
      # Title of the plot
      title(main="MDS Plot",cex.main = 1.7,font.main = 2)
      # Text of the plot
      text(x=0.5, y = 0.45, paste0("MDS plot presenting the reference Samples (", paste(reference_name,collapse="+"), ") \n",
                                   "and the test samples (", paste(Test_name,collapse="-"), ") " ),cex = 1.5, col = "black")
      
      # MDS plot where the refernce groups has been clustered and all the test samples are represented as individual points
      sclass(unifract_dist,factor(input_table[,1],levels=unique(input_table[,1])),NULL,colour_matrix[levels(factor(input_table[,1],levels=unique(input_table[,1]))),2] )
      
      # Set narrow margings for the last plot
      par(mar = c(0.1, 0.1, 0.1, 0.1))
      # Create an empty plot
      plot.new()
      
      # Default margins
      par(mar= c(5.1, 4.1, 4.1, 2.1))
      
      # Print the MDS plot in png format
      png_plot(paste("DiBaAn-page",page,"MDS plot.png"),
               sclass(unifract_dist,factor(input_table[,1],levels=unique(input_table[,1])),NULL,colour_matrix[levels(factor(input_table[,1],levels=unique(input_table[,1]))),2] ))
      
      
      
      
      ########## Page 3 ###########
      
      # Index indicating the page of the report
      page <- page+1
      
      
      # Title of the 3rd page
      text <- "Phylogram"
      # Create a text grob for the title
      tgrob1 <- text_grob(text, face = "bold.italic", color = "Black",size=16)
      
      # Description of the Phylogram
      text2 <- paste0("Phylogram of the reference (",paste(reference_name,collapse="+"),") and test (", paste(Test_name,collapse="-"), ")\n" , "samples")
      # Create a text grob for the Description
      tgrob2 <- text_grob(text2, face = "italic", color = "Black",size=12)
      
      # Tree of the reference end test groups
      tree_plot <- tree(unifract_dist,factor(meta_file[,mapping_column],levels = unique(meta_file[,mapping_column])),
                        colour_matrix[levels(factor(meta_file[,mapping_column],levels = unique(meta_file[,mapping_column]))),2] )
      
      
      # Arrange in the grid the above elements
      grid.arrange(tgrob1,tgrob2,tree_plot,nrow = 3,heights=c(0.5,0.5,2.5))
      
      # Print the Phylogram in png format
      png_plot(paste("DiBaAn-page",page,"Phylogram.png"),tree_plot)
      
      
      ############ Page 4 ##############
      
      # Index indicating the page of the report
      page <- page+1
      
      
      # Description of the list_box2[[1]]
      text <- paste0(plot_type," displaying the distribution of the distances around their own self-",central_point,".\n",
                     "This plot can provide additional information about the compactness of the groups.")
      # Create a text grob for the Description of the list_box2[[1]]
      tgrob <- text_grob(text, face = "italic", color = "Black",size=12)
      
      
      # Choose the selected type of plot
      if (plot_type == "Boxplots") {plot <-list_box2[[1]] } else if (plot_type == "Point plots") {plot <-list_point2[[1]] } else {plot <-list_violin2[[1]]}
      
      # Arrange in the grid the above elements
      grid.arrange(arrangeGrob( plot,pvaltable2[[1]],nrow = 1,ncol=2,widths = c(1.1,1)),tgrob,nrow = 2,as.table = T,heights=c(2.5,1))
      
      # Print the plot in png format
      png_plot(paste("DiBaAn-page",page,plot_type,".png"),plot)
      
      
      
      ########## Page 5 ############
      
      # Index indicating the page of the report
      page <- page+1
      
      
      # Set the layout of the cover page (3 rows, 1 column)
      layout(matrix(c(1, 1, 1,
                      2, 2, 2,
                      3, 3, 3), nrow=3, byrow=TRUE),heights=c(1.2,2.9,0.3))
      
      # Create an empty plot
      plot.new()
      # Title of the plot
      title(main="MDS Plot",cex.main = 1.7,font.main = 2)
      # Text of the plot
      text(x=0.5, y = 0.45, paste0("MDS plot presenting the reference Samples (", paste(reference_name,collapse="+"), ") \n"," clustered into 3 Groups\n",
                                   "and the test samples (", paste(Test_name,collapse="-"), ") as individual points" ),cex = 1.5, col = "black")
      
      # MDS plot where the refernce groups has been clustered and all the test samples are represented as individual points
      sclass(unifract_dist,as.factor(mapping_cluster[which(input_table[,1] %in% reference_name)]),as.factor(mapping_cluster[which(input_table[,1] %in% levels)]),
             colour_matrix[levels(factor(mapping_cluster[which(input_table[,1] %in% reference_name)],levels = unique(mapping_cluster[which(input_table[,1] %in% reference_name)]))),2] )
      
      # Set narrow margings for the last plot
      par(mar = c(0.1, 0.1, 0.1, 0.1))
      # Create an empty plot
      plot.new()
      
      # Default margins
      par(mar= c(5.1, 4.1, 4.1, 2.1))
      
      
      ########## Page 6 ############
      
      # Index indicating the page of the report
      page <- page+1
      
      
      # Description of the list_box2[[2]]
      text <- paste0(plot_type," displaying the distribution of the distances around the closest referenece ",central_point,".\n",
                     "For each sample of the test groups (",paste(Test_name,collapse="-"), ") the nearest reference ",central_point," \nwas found ", "and then these distances were calculated.")
      
      # Create a text grob for the Description of the list_box2[[2]]
      tgrob <- text_grob(text, face = "italic", color = "Black",size=12)
      
      # Choose the selected type of plot
      if (plot_type == "Boxplots") {plot <-list_box2[[2]] } else if (plot_type == "Point plots") {plot <-list_point2[[2]] } else {plot <-list_violin2[[2]]}
      
      
      # Arrange in the grid the above elements
      grid.arrange(arrangeGrob( plot,pvaltable2[[2]],nrow = 1,ncol=2,widths = c(1.1,1)),tgrob,nrow = 2,as.table = T,heights=c(2.5,1))
      
      # Print the plot in png format
      png_plot(paste("DiBaAn-page",page,plot_type,".png"),plot)
      
      
      
      ############ Page 7 ##########
      
      # Index indicating the page of the report
      page <- page+1
      
      
      # Title of the 7th page 
      text1 <- "Distances from the closest reference cluster"
      # Create a text grob for the title
      tgrob <- text_grob(text1,face = "bold.italic", color = "Black",size=15)
      
      # Description of the matrix of the median distances
      text2 <- paste("For each group the median distance of the sample from \n every reference",central_point, "was calculated. \n Then, the closest reference cluster was found.\n")
      # Create a text grob for the median matrix
      tgrob2 <- text_grob(text2,face = "italic", color = "Black",size=13)
      
      
      # Arrange in the grid the above elements
      grid.arrange(tgrob,median_table,tgrob2,nrow=3,heights=c(0.4,1,1.5))
      
      
      if (reference_clusters > 1) {
        
        ############ Page 8-undefined ##########
        
        if (chi_square == 1){
          
          # Index indicating the page of the report
          page <- page+1
          
          # Tile of the page
          text1 <- "Chi-square test \nGoodness of Fit analysis"
          
          # Create a text grob for the title
          tgrob <- text_grob(text1,face = "bold.italic", color = "Black",size=20)
          
          # Text of the footnote
          text2 <- paste("*For the next",(length(exploratory_columns)*length(Test_name)+length(exploratory_columns)+1),"pages")
          
          # Create a text grob for the footnote
          tgrob2 <- text_grob(text2,face = "italic", color = "Black",size=12)
          
          # Arrange in the grid the above elements
          grid.arrange(tgrob,tgrob2, nrow=2,heights=c(2,0.5))
          
        }
        
        # Index indicating the page of the report
        page <- page+1
        
        # Tile of the 7th page
        text1 <- paste("Distances from the closest reference",central_point)
        
        # Create a text grob for the title
        tgrob <- text_grob(text1,face = "bold.italic", color = "Black",size=15)
        
        # Description of the two tables
        text2 <- paste(" De novo clustering has been applied to the reference dataset.\n Following, the test samples have been clustered based on the distances \n from the closest reference", central_point, "(matrix on the left).
                   \n The matrix on the right presents the p-values of the Chi square test for every \n pair of groups." )
        
        # Create a text grob for the description
        tgrob2 <- text_grob(text2, face = "italic", color = "Black",size=12)
        
        # Arrange in the grid the above elements
        grid.arrange(tgrob,arrangeGrob(DB_matrix_table,chisquare_pvalues,nrow=1),tgrob2,nrow=3,heights=c(0.4,1,1))
        
        
        if (chi_square == 1){
          for (k in 1:length(exploratory_columns)){
            
            # Index indicating the page of the report
            page <- page+1
            
            # Check if the column is character of factor
            if (is.character(meta_file[,exploratory_columns[k]])==TRUE |is.factor(meta_file[,exploratory_columns[k]])==TRUE ){
              
              # Tile of the page
              text1 <- paste0("'",exploratory_columns[k],"' Column")
              
              # Create a text grob for the title
              tgrob <- text_grob(text1,face = "bold.italic", color = "Black",size=20)
              
              # Main text of the page
              text <- paste0("The relationship between the reference clusters and the  '",exploratory_columns[k],"' column  \nwill be examined in the following pages.\n\n",
                             "The distribution of the samples across the reference clusters is presented in Table(s) \nof the first line.\n\n",
                             "Table 2 presents the p-values of the Chi-square test of the observed distribution \nof the counts against the expected distribution if the counts were \nuniformly distributed.\n\n",
                             "Table 3 contains the p-values of the pairwise comparisons.\n\n",
                             "The tables contain only the 10 most significant p-values. The rest of the p-values\n
                           have been printed in the results folder.")
              
              # Create a text grob for the text
              tgrob2 <- text_grob(text,face = "italic", color = "Black",size=13)
              
              # Create a text grob for the footnote
              text3 <- "*If the number of counts is limited or zero counts occurs, then the Chi-square approximation \n may be incorrect." 
              
              # Create a text grob for the footnote
              tgrob3 <- text_grob(text3,face = "italic", color = "Black",size=11)
              
              # Arrange in the grid the above elements
              grid.arrange(tgrob,tgrob2,tgrob3 ,nrow=3)
              
              
              for(i in 1:length(Test_name)){
                
                # Index indicating the page of the report
                page <- page+1
                
                # Subset the names of the reference group
                names_ref <-rownames(input_table[which(input_table[,1] %in% reference_name),])
                # Subset the names of each test group
                names_test <- rownames(input_table[input_table[,1] == Test_name[i],])
                
                # Extract the factors of the reference group
                factor_ref <- factor(meta_file[names_ref,exploratory_columns[k]],levels=unique(meta_file[names_ref,exploratory_columns[k]]))
                # Extract the factors of each test group
                factor_test <- factor(meta_file[names_test,exploratory_columns[k]],levels=unique(meta_file[names_test,exploratory_columns[k]]))
                
                if (nlevels(factor_ref) > 0 | nlevels(factor_test) > 0){
                  # Check if one or two matrices will be formed 
                  if (length(which(levels(factor_ref) %in% levels(factor_test) == TRUE)) > 1 & nlevels(factor_ref) > 1 & nlevels(factor_test) > 1){
                    
                    # Form the matrix where the counts of the observations will be placed
                    chi_ref <- data.frame(matrix(ncol = (nlevels(factor_ref)) , nrow = reference_clusters))
                    # Naming the rows of the chi_ref
                    rownames(chi_ref) <- paste(reference_name,1:reference_clusters)
                    # Naming the columns of the chi_ref
                    colnames(chi_ref) <- levels(factor_ref)
                    
                    # Form the matrix where the counts of the observations will be placed
                    chi_test <- data.frame(matrix(ncol = (nlevels(factor_test)) , nrow = reference_clusters))
                    # Naming the rows of the chi_test
                    rownames(chi_test) <- paste(reference_name,1:reference_clusters)
                    # Naming the columns of the chi_ref
                    colnames(chi_test) <- levels(factor_test)
                    
                    for (ii in 1:reference_clusters){
                      # Subset the names of the reference group for each cluster
                      names2_ref <- rownames(mapping_file[mapping_file[,1]==paste(reference_name,ii),])
                      for (iii in 1:nlevels(factor_ref)){
                        # Place the counts in the chi_ref
                        chi_ref[ii,iii] <- length(which(meta_file[names2_ref,exploratory_columns[k]]==levels(factor_ref)[iii]))
                      }
                    }
                    
                    # Replace 0 with 0.01
                    chi_ref[chi_ref == 0] <- 0.01
                    
                    for (ii in 1:reference_clusters){
                      
                      # Subset the names of the test group for each cluster
                      names2_test <- distances_matrix[distances_matrix[,"Nearest reference cluster"] == ii,]
                      names2_test <- rownames(names2_test[rownames(names2_test)%in% names_test,])
                      
                      for (iii in 1:nlevels(factor_test)){
                        # Place the counts in the chi_test matrix
                        chi_test[ii,iii] <- length(which(meta_file[names2_test,exploratory_columns[k]]==levels(factor_test)[iii]))
                      }
                    }
                    
                    # Replace 0 with 0.01
                    chi_test[chi_test == 0] <- 0.01
                    
                    
                    # Performing Chi-square test(Goodness of fit) for each column of chi_ref
                    exp_col_ref <- apply(chi_ref,2,function(x) chisq.test(x)$p.value)
                    
                    # Performing Chi-square test(Goodness of fit) for each row of chi_ref
                    exp <-apply(chi_ref,2,sum)/reference_clusters
                    exp_r_ref <- apply(chi_ref,1,function(x) chisq.test(x,p=exp,rescale.p = TRUE)$p.value)
                    
                    # Performing Chi-square test(Goodness of fit) for each column of chi_test
                    exp_col_test <- apply(chi_test,2,function(x) chisq.test(x)$p.value)
                    
                    # Performing Chi-square test(Goodness of fit) for each row of chi_test
                    exp <-apply(chi_test,2,sum)/reference_clusters
                    exp_r_test <- apply(chi_test,1,function(x) chisq.test(x,p=exp,rescale.p = TRUE)$p.value)
                    
                    # Performing fdr correction
                    adjustb <- round(p.adjust(c(exp_col_ref,exp_r_ref,exp_col_test,exp_r_test), method = "BH"),4)
                    # Add the adjusted p-values in the expected matrix
                    expected <- data.frame(cbind(c(exp_col_ref,exp_r_ref,exp_col_test,exp_r_test),adjustb))
                    # Naming the columns of p_value data frame
                    colnames(expected) <- c( "P-values","Adj p-values")
                    # Naming the rows of p_value data frame
                    rownames(expected) <- c(paste(names(exp_col_ref),"(Ref)"),paste(names(exp_r_ref),"(Ref)"),paste0(names(exp_col_test),"(",Test_name[i],")"),paste(names(exp_r_test),"(",Test_name[i],")"))
                    
                    
                    # Order the p-values
                    expected <- expected[order(expected[,2]),]
                    # Choose the first 10 most significant p-values
                    if (nrow(expected)>10) {expected2 <- expected[1:10,]} else {expected2 <- expected}
                    
                    
                    # Find all the available pairwise combinations for the reference group
                    combs <- combs(levels(factor_ref), 2)
                    # Save the names of that combinations
                    group_names <- apply(combs,1,function(x)paste0(x[1],"(Ref)-",x[2],"(Ref)"))
                    
                    
                    # P_value -> a vector where the p values will be stored
                    chi_p_value_ref <- c(); test <- c()
                    for(chii in 1:nrow(combs)) { 
                      # Performed the Chi square test
                      chi_p_value_ref[chii] <- round(chisq.test(chi_ref[,combs[chii,2]],p=chi_ref[,combs[chii,1]],rescale.p = TRUE)$p.value,4) 
                      }
                    
                    #------------ Rows pairwise comparisons ----------#
                    # Find all the available pairwise combinations for the reference group
                    combs <- combs(rownames(chi_ref), 2)
                    # Save the names of that combinations
                    group_names_rows <- apply(combs,1,function(x)paste0(x[1],"(Ref)-",x[2],"(Ref)"))
                    
                    
                    # P_value -> a vector where the p values will be stored
                    chi_p_value_ref_rows <- c(); test <- c()
                    for(chii in 1:nrow(combs)) { 
                      # Performed the Chi square test
                      chi_p_value_ref_rows[chii] <- round(chisq.test(as.numeric(chi_matrix[combs[chii,2],]),p=as.numeric(chi_matrix[combs[chii,1],]),rescale.p = TRUE)$p.value,4) }
                    
                    # Replace 0.01 with 0
                    chi_ref[chi_ref == 0.01] <- 0
                    
                    
                    # Forming the data frame containing the p values
                    p_value_ref <- data.frame(group_names)
                    # Naming the rows of the p_value data frame
                    p_value_ref[,2] <- chi_p_value_ref
                    # Naming the columns of p_value data frame
                    colnames(p_value_ref) <- c("Groups","p value")
                    
                    
                    # Forming the data frame containing the p values
                    p_value_ref_rows <- data.frame(group_names_rows)
                    # Naming the rows of the p_value data frame
                    p_value_ref_rows[,2] <- chi_p_value_ref_rows
                    # Naming the columns of p_value data frame
                    colnames(p_value_ref_rows) <- c("Groups","p value")
                    
                    # Find all the available pairwise combinations for the test groups
                    combs <- combs(levels(factor_test), 2)
                    # Save the names of that combinations
                    group_names <- apply(combs,1,function(x)paste0(x[1],"(",Test_name[i],")-",x[2],"(",Test_name[i],")"))
                    
                    
                    # P_value -> a vector where the p values will be stored
                    chi_p_value_test <- c(); test <- c()
                    for(chii in 1:nrow(combs)) { 
                      # Performed the Chi square test
                      chi_p_value_test[chii] <- round(chisq.test(chi_test[,combs[chii,2]],p=chi_test[,combs[chii,1]],rescale.p = TRUE)$p.value,4)
                      }
                    
                    
                    #------------ Rows pairwise comparisons ----------#
                    
                    # Find all the available pairwise combinations for the test groups
                    combs <- combs(rownames(chi_test), 2)
                    # Save the names of that combinations
                    group_names_rows <- apply(combs,1,function(x)paste0(x[1],"(",Test_name[i],")-",x[2],"(",Test_name[i],")"))
                    
                    
                    # P_value -> a vector where the p values will be stored
                    chi_p_value_test_rows <- c(); test <- c()
                    for(chii in 1:nrow(combs)) { 
                      # Performed the Chi square test
                      chi_p_value_test_rows[chii] <- round(chisq.test(as.numeric(chi_matrix[combs[chii,2],]),p=as.numeric(chi_matrix[combs[chii,1],]),rescale.p = TRUE)$p.value,4) }
                    
                    # Replace 0.01 with 0
                    chi_test[chi_test == 0.01] <- 0
                    
                    
                    
                    # Forming the data frame containing the p values
                    p_value_test <- data.frame(group_names)
                    # Naming the rows of the p_value data frame
                    p_value_test[,2] <- chi_p_value_test
                    # Naming the columns of p_value data frame
                    colnames(p_value_test) <- c("Groups","p value")
                    
                    
                    # Forming the data frame containing the p values
                    p_value_test_rows <- data.frame(group_names_rows)
                    # Naming the rows of the p_value data frame
                    p_value_test_rows[,2] <- chi_p_value_test_rows
                    # Naming the columns of p_value data frame
                    colnames(p_value_test_rows) <- c("Groups","p value")
                    
                    
                    
                    # Form a data frame where the p-values for the coompn factors will be placed 
                    common <- data.frame()
                    # row index initiation
                    row <- 0
                    for(l in 1:nlevels(factor_ref)){
                      
                      if (length(which(levels(factor_ref)[l] == levels(factor_test)) == 1)){
                        # Add a row
                        row <- row+1
                        # Add the name of the row
                        common[row,1] <- paste0(levels(factor_ref)[l],"(Ref) -",levels(factor_ref)[l],"(",Test_name[i],")")
                        # Perform Chi-square test
                        common[row,2] <-  round(chisq.test(chi_test[,levels(factor_ref)[l]],p=chi_ref[,levels(factor_ref)[l]],rescale.p = TRUE)$p.value,4)
                        # Naming the columns
                        colnames(common) <- c("Groups","p value")
                      }
                    }
                    
                    # rbind all the p-values
                    p_value <- rbind(p_value_ref,p_value_test,p_value_ref_rows,p_value_test_rows,common)
                    
                    # Convert p-values column into numeric
                    p_value[,2] <- as.numeric(p_value[,2])
                    # Perform fdr correction
                    p.adjust <- round(p.adjust(p_value[,2], method = "BH"),4)
                    # Add the adjusted p-values in the pvalues matrix
                    p_value[,3] <- p.adjust
                    # Naming the columns of p_value data frame
                    colnames(p_value) <- c("Groups","p value","Adj p-values")
                    # Order the p-values
                    p_value <- p_value[order(p_value[,3]),]
                    # Choose the first 10 most significant p-values
                    if (nrow(p_value)>10) {p_value2 <- p_value[1:10,]} else {p_value2 <- p_value}
                    
                    # Write Reference counts table at result folder
                    write_table(paste0("DeBaAn-page ",page,reference_name," Group (Table 1).tab"),chi_ref,"table")
                    # Write test counts table at result folder
                    write_table(paste0("DeBaAn-page ",page,Test_name[i]," Group (Table 2).tab"),chi_test,"table")
                    # Write Expected-observed p-values at result folder
                    write_table(paste0("DeBaAn-page ",page, " Expected-observed (Table 2).tab"),expected,"pvalues")
                    # Write pairwise p-values at result folder
                    write_table(paste0("DeBaAn-page ",page, " Pairwise comparisons (Table 3).tab"),p_value,"pvalues")
                    
                    if ((sum(nchar(colnames(chi_ref)))+sum(nchar(colnames(chi_test)))) < 80){
                      # Tile of the page
                      text1 <- paste0(exploratory_columns[k],"-",Test_name[i])
                      
                      # Create a text grob for the title
                      tgrob <- text_grob(text1,face = "bold.italic", color = "Black",size=15)
                      
                      # Tile of the p-value tables
                      text1 <- "P-values of the Chi-square test (Goodness of Fit)"
                      
                      # Create a text grob for the title
                      tgrob2 <- text_grob(text1,face = "bold.italic", color = "Black",size=12)
                      
                      # Arrange in the grid the above elements
                      grid.arrange(tgrob,arrangeGrob(gtableGrob(paste(reference_name,"Group (Table 1)"),chi_ref,setDT="YES",size=12),gtableGrob(paste(Test_name[i],"Group (Table 2)"),chi_test,setDT="YES",size=12),nrow=1)
                                   , arrangeGrob(tgrob2,arrangeGrob(gtableGrob("Expected-observed (Table 2)",expected2,setDT="YES",size=12),gtableGrob("Pairwise comparisons (Table 3)",p_value2,setDT="NO",size=12),nrow=1),nrow=2,heights=c(0.3,1))
                                   , ncol = 1, nrow = 3,heights=c(0.4,1,2))
                    } else {
                      text1 <- paste(exploratory_columns[k],"-",Test_name[i])
                      
                      # Create a text grob for the title
                      tgrob <- text_grob(text1,face = "bold.italic", color = "Black",size=15)
                      
                      # Tile of the p-value tables
                      text2 <-"The tables are too big to be printed!!\nYou can find them in .tab format in the results folder."
                      
                      # Create a text grob for the title
                      tgrob2 <- text_grob(text2,face = "bold.italic", color = "Black",size=12)
                      
                      grid.arrange(tgrob,tgrob2,nrow=2)
                    }
                    
                    
                  } else {
                    
                    # Form the matrix where the counts of the observations will be placed
                    chi_matrix <- data.frame(matrix(ncol = (nlevels(factor_ref)+nlevels(factor_test)) , nrow = reference_clusters))
                    # Naming the rows of the chi matrix
                    rownames(chi_matrix) <- paste(reference_name,1:reference_clusters)
                    if (nlevels(factor_ref) == 0){
                      # Naming the columns of the chi matrix
                      colnames(chi_matrix) <- c(paste(paste0(levels(factor(factor_test,levels=unique(factor_test))),"(",Test_name[i],")")))
                    }else if (nlevels(factor_test) == 0){
                      # Naming the columns of the chi matrix
                      colnames(chi_matrix) <- c(paste(levels(factor(factor_ref,levels=unique(factor_ref))),"(Ref)"))
                    } else {
                      # Naming the columns of the chi matrix
                      colnames(chi_matrix) <- c(paste(levels(factor(factor_ref,levels=unique(factor_ref))),"(Ref)"),paste0(levels(factor(factor_test,levels=unique(factor_test))),"(",Test_name[i],")"))
                    }
                    
                    if (nlevels(factor_ref) > 0){
                      for (ii in 1:reference_clusters){
                        # Subset the names of the reference group
                        names2 <- rownames(mapping_file[mapping_file[,1] == paste(reference_name,ii),])
                        for (iii in 1:nlevels(factor_ref)){
                          # Place the counts of the reference group in the chi matrix
                          chi_matrix[ii,iii] <- length(which(meta_file[names2,exploratory_columns[k]]==levels(factor_ref)[iii]))
                        }
                      }
                    }
                    
                    if (nlevels(factor_test) > 0){
                      for (ii in 1:reference_clusters){
                        # Subset the names of the test group for each cluster
                        names2_test <- distances_matrix[distances_matrix[,"Nearest reference cluster"] == ii,]
                        names2_test <- rownames(names2_test[rownames(names2_test) %in% names_test,])
                        
                        
                        for (iii in 1:nlevels(factor_test)){
                          # Place the counts of the test group in the chi matrix
                          chi_matrix[ii,(iii+nlevels(factor_ref))] <- length(which(meta_file[names2_test,exploratory_columns[k]] == levels(factor_test)[iii]))
                        }
                      }
                    }
                    
                    # Replace 0 with 0.01
                    chi_matrix[chi_matrix == 0] <- 0.01
                    
                    # Performing Chi-square test(Goodness of fit) for each column
                    exp_col <- apply(chi_matrix,2,function(x) chisq.test(x)$p.value)
                    
                    if (ncol(chi_matrix) > 1){
                      # Performing Chi-square test(Goodness of fit) for each row
                      exp <-apply(chi_matrix,2,sum)/reference_clusters 
                      exp_r <- apply(chi_matrix,1,function(x) chisq.test(x,p=exp,rescale.p = TRUE)$p.value)
                      
                      
                      # Performing fdr correction
                      adjustb <- round(p.adjust(c(exp_col,exp_r), method = "BH"),4)
                      # Add the adjusted p-values in the pvalues matrix
                      expected <- data.frame(cbind(c(exp_col,exp_r),adjustb))
                      colnames(expected) <- c( "P-values","Adj p-values")
                      
                    }else{
                      
                      # Performing fdr correction
                      adjustb <- round(p.adjust(c(exp_col), method = "BH"),4)
                      # Add the adjusted p-values in the pvalues matrix
                      expected <- data.frame(cbind(exp_col,adjustb))
                      # Naming the rows of p_value data frame
                      rownames(expected) <- c(names(exp_col))
                    }
                    # Naming the columns of expected data frame
                    colnames(expected) <- c( "P-values","Adj p-values")
                    
                    
                    # Order the p-values
                    expected <- expected[order(expected[,2]),]
                    # Choose the first 10 most significant p-values
                    if (nrow(expected) > 10) {expected2 <- expected[1:10,]} else {expected2 <- expected}
                    
                    if(ncol(chi_matrix) > 1){
                      # Find all the available pairwise combinations
                      combs <- combs(colnames(chi_matrix),2)
                      # Save the names of that combinations
                      group_names <- apply(combs,1,function(x)paste0(x[1],"-",x[2]))
                      
                      
                      # P_value -> a vector where the p values will be stored
                      chi_p_value <- c(); test <- c()
                      for(chii in 1:nrow(combs)) { 
                        # Performed the Chi square test
                        chi_p_value[chii] <- round(chisq.test(chi_matrix[,combs[chii,2]],p=chi_matrix[,combs[chii,1]],rescale.p = TRUE)$p.value,4) }
                      
                      #---------- Rows pairwise comparisons ----------#
                      
                      # Find all the available pairwise combinations
                      combs <- combs(rownames(chi_matrix),2)
                      # Save the names of that combinations
                      group_names_rows <- apply(combs,1,function(x)paste0(x[1],"-",x[2]))
                      
                      
                      # P_value -> a vector where the p values will be stored
                      chi_p_value_rows <- c(); test <- c()
                      for(chii in 1:nrow(combs)) { 
                        # Performed the Chi square test
                        chi_p_value_rows[chii] <- round(chisq.test(as.numeric(chi_matrix[combs[chii,2],]),p=as.numeric(chi_matrix[combs[chii,1],]),rescale.p = TRUE)$p.value,4) }
                      
                      # Replace 0.01 with 0
                      chi_matrix[chi_matrix == 0.01] <- 0
                      
                      
                      # Forming the data frame containing the p values
                      p_value <- data.frame(group_names)
                      # Naming the rows of the p_value data frame
                      p_value[,2] <- chi_p_value
                      # Naming the columns of p_value data frame
                      colnames(p_value) <- c("Groups","p value")
                      
                      
                      # Forming the data frame containing the p values
                      p_value_rows <- data.frame(group_names_rows)
                      # Naming the rows of the p_value data frame
                      p_value_rows[,2] <- chi_p_value_rows
                      # Naming the columns of p_value data frame
                      colnames(p_value_rows) <- c("Groups","p value")
                      
                      # rbind all the p-values
                      p_value <- rbind(p_value,p_value_rows)
                      
                      
                      # Convert p-values column into numeric
                      p_value[,2] <- as.numeric(p_value[,2])
                      # Performing fdr correction
                      p.adjust <- round(p.adjust(p_value[,2], method = "BH"),4)  
                      # Add the adjusted p-values in the pvalues matrix
                      p_value[,3] <- p.adjust
                      # Naming the columns of p_value data frame
                      colnames(p_value) <- c("Groups","p value","Adj p-values")
                      # Order the p-values
                      p_value <- p_value[order(p_value[,3]),]
                      # Choose the first 10 most significant p-values
                      if (nrow(p_value) > 10) {p_value2 <- p_value[1:10,]} else {p_value2 <- p_value}
                      
                    } else {
                      
                      #-------------- Rows pairwise comparisons -------------#
                      # Find all the available pairwise combinations
                      combs <- combs(rownames(chi_matrix),2)
                      # Save the names of that combinations
                      group_names_rows <- apply(combs,1,function(x)paste0(x[1],"-",x[2]))
                      
                      
                      # P_value -> a vector where the p values will be stored
                      chi_p_value_rows <- c(); test <- c()
                      for(chii in 1:nrow(combs)) { 
                        # Performed the Chi square test
                        chi_p_value_rows[chii] <-round(chisq.test(as.numeric(chi_matrix[combs[chii,2],]),p=as.numeric(chi_matrix[combs[chii,1],]),rescale.p = TRUE)$p.value,4) }
                      
                      # Replace 0.01 with 0
                      chi_matrix[chi_matrix == 0.01] <- 0
                      
                      # Forming the data frame containing the p values
                      p_value <- data.frame(group_names_rows)
                      # Naming the rows of the p_value data frame
                      p_value[,2] <- chi_p_value_rows
                      # Naming the columns of p_value data frame
                      colnames(p_value) <- c("Groups","p value")
                      
                      # Convert p-values column into numeric
                      p_value[,2] <- as.numeric(p_value[,2])
                      # Performing fdr correction
                      p.adjust <- round(p.adjust(p_value[,2], method = "BH"),4)  
                      # Add the adjusted p-values in the pvalues matrix
                      p_value[,3] <- p.adjust
                      # Naming the columns of p_value data frame
                      colnames(p_value) <- c("Groups","p value","Adj p-values")
                      # Order the p-values
                      p_value <- p_value[order(p_value[,3]),]
                      # Choose the first 10 most significant p-values
                      if (nrow(p_value) > 10){p_value2 <- p_value[1:10,]} else {p_value2 <- p_value}
                    }
                    
                    # Write counts table at results folder
                    write_table(paste0("DeBaAn-page ",page, " Counts per Cluster (Table 1).tab"),chi_matrix,"table")
                    # Write Expected-observed p-values at results folder
                    write_table(paste0("DeBaAn-page ",page, " Expected-observed (Table 2).tab"),expected,"pvalues")
                    # Write pairwise p-values at results folder
                    write_table(paste0("DeBaAn-page ",page, " Pairwise comparisons (Table 3).tab"),p_value,"pvalues")
                    
                    if (sum(nchar(colnames(chi_matrix))) < 80) {
                      
                      # Tile of the page
                      text1 <- paste(exploratory_columns[k],"-",Test_name[i])
                      
                      # Create a text grob for the title
                      tgrob <- text_grob(text1,face = "bold.italic", color = "Black",size=15)
                      
                      # Tile of the p-value tables
                      text2 <- "P-values of the Chi-square test (Goodness of Fit)"
                      
                      # Create a text grob for the title
                      tgrob2 <- text_grob(text2,face = "bold.italic", color = "Black",size=12)
                      
                      # Arrange in the grid the above elements
                      grid.arrange(tgrob,gtableGrob("Counts per Cluster (Table 1)",chi_matrix,setDT="YES",size=12),
                                   arrangeGrob(tgrob2,arrangeGrob(gtableGrob("Expected-observed (Table 2)",expected2,setDT="YES",size=12),gtableGrob("Pairwise comparisons (Table 3)",p_value2,setDT="NO",size=12),nrow=1),nrow=2,heights=c(0.3,1), ncol = 1), 
                                   nrow = 3,heights=c(0.4,1,2))
                    } else{
                      text1 <- paste(exploratory_columns[k],"-",Test_name[i])
                      
                      # Create a text grob for the title
                      tgrob <- text_grob(text1,face = "bold.italic", color = "Black",size=15)
                      
                      # Tile of the p-value tables
                      text2 <-  "The tables are too big to be printed!!\nYou can find them in .tab format in the results folder."
                      
                      # Create a text grob for the title
                      tgrob2 <- text_grob(text2,face = "bold.italic", color = "Black",size=12)
                      
                      grid.arrange(tgrob,tgrob2,nrow=2)
                      
                    }
                  }
                  
                }else {
                  # Tile of the page
                  text1 <- paste0("'",exploratory_columns[k],"' Column")
                  
                  # Create a text grob for the title
                  tgrob <- text_grob(text1,face = "bold.italic", color = "Black",size=20)
                  
                  # Main text of the page
                  text2 <- paste0("It was not possible to perform the chi-square test for the '",Test_name[i],"'. There are no Factors")
                  
                  # Create a text grob for the text
                  tgrob2 <- text_grob(text2,face = "italic", color = "Black",size=12)
                  
                  
                  # Arrange in the grid the above elements
                  grid.arrange(tgrob,tgrob2 ,nrow=2)
                }
              }
              
            } else {
              
              
              # Tile of the page
              text1 <- paste0("'",exploratory_columns[k],"' Column")
              
              # Create a text grob for the title
              tgrob <- text_grob(text1,face = "bold.italic", color = "Black",size=20)
              
              # Main text of the page
              text2 <- paste0("It was not possible to perform the chi-square test for the '",exploratory_columns[k],"'. The variable must be categorical ")
              
              # Create a text grob for the text
              tgrob2 <- text_grob(text2,face = "italic", color = "Black",size=12)
              
              
              # Arrange in the grid the above elements
              grid.arrange(tgrob,tgrob2 ,nrow=2)
            }
          }
        }
        
        ############ Page (Percentages of each cluster) ##########
        
        # Index indicating the page of the report
        page <- page+1
        
        # Tile of the page
        text1 <- "Percentages of each cluster"
        
        # Create a text grob for the title
        tgrob <- text_grob(text1,face = "bold.italic", color = "Black",size=15)
        
        # Arrange in the grid the above elements
        grid.arrange(tgrob,percentage_table, barplot, ncol = 1, nrow = 3,heights=c(0.4,1.2,1.5))
        
        # Print the Barplot in png format
        png_plot(paste("DiBaAn-page",page, "Barplot.png"),barplot)
        
        
      }
      
      
      ############ Page (MDS reference clusters and the distances-based clusters) ########
      
      # Index indicating the page of the report
      page <- page+1
      
      
      # Set the layout of the cover page (3 rows, 1 column)
      layout(matrix(c(1, 1, 1,
                      2, 2, 2,
                      3, 3, 3), nrow=3, byrow=TRUE),heights=c(0.3,3,1) )
      
      # Set narrow margings for the last plot
      par(mar = c(0.1, 0.1, 3.5, 0.1))
      # Create an empty plot
      plot.new()
      # Title for the empty plot
      title(main="MDS Plot",cex.main = 1.7)
      
      # MDS plot presenting the reference clusters and the distances-based clusters of the test groups
      sclass(unifract_dist,factor(plot_matrix[,ncol(plot_matrix)],levels = unique(plot_matrix[,ncol(plot_matrix)])),NULL,colour_matrix[colour_matrix[,1]%in%factor(plot_matrix[,ncol(plot_matrix)],levels = unique(plot_matrix[,ncol(plot_matrix)])),2])
      
      # Create an empty plot
      plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      # Description of the MDS plot
      text(x=0.5, y = 0.8, paste0("MDS plot presenting the reference Samples (", paste(reference_name,collapse="+"), ") clustered into 3 Groups\n", "and the test samples (" ,paste(Test_name,collapse="-"), ")\n",
                                  "clustered into groups based on their distances from the closest reference ",central_point,"." ), cex = 1.5, col = "black")
      
      # Default margins
      par(mar= c(5.1, 4.1, 4.1, 2.1))
      
      # Print the MDS plot in png format
      png_plot(paste("DiBaAn-page",page, "MDS Plot.png"),sclass(unifract_dist,factor(plot_matrix[,ncol(plot_matrix)],levels = unique(plot_matrix[,ncol(plot_matrix)])),NULL,colour_matrix[colour_matrix[,1]%in%factor(plot_matrix[,ncol(plot_matrix)],levels = unique(plot_matrix[,ncol(plot_matrix)])),2])
      )
      
      ######### Page (boxplots BD clustering) ########
      
      # Index indicating the page of the report
      page <- page+1
      
      
      # Description of list_box
      text <- paste0("The test samples have been clustered based on their distances from the\n", "closest reference ",central_point,".\n",
                     "The ",plot_type, " present the distribution of distances of each of these samples from \n",
                     "their nearest reference cluster\n",
                     "The p-values for every reference cluster and its nearest test clusters are displayed.\n",
                     "On the next pages each reference cluster will be examined seperatly.")
      
      # Create a text grob for the description of the list_box
      tgrob <- text_grob(text, face = "italic", color = "Black",size=12)
      
      # Choose the selected type of plot
      if (plot_type == "Boxplots") {plot <-list_box } else if (plot_type == "Point plots") {plot <-list_point } else {plot <-list_violin}
      
      # Arrange in the grid the above elements
      grid.arrange(arrangeGrob( plot,pvaltable,nrow = 1,ncol=2,widths = c(1,1)),tgrob,nrow = 2,as.table = T,heights=c(2.5,1))
      
      # Print the plot in png format
      png_plot(paste0("DiBaAn-page ",page," ",plot_type,".png"),plot)
      
      
      ############# Pages (individual cluster analysis) ############
      
      # Create the boxplots of every reference cluster
      for (cluster in 1:reference_clusters){
        # Subset the samples 
        report_plot <- plot_df_a[dist_plot_1[,3] == cluster,]
        # Name the rows
        rownames(report_plot) <- rownames(dist_plot_1[dist_plot_1[,3] == cluster,])
        # Convert label column into factor
        report_plot$label <- factor(report_plot$label,levels = unique(report_plot$label))
        
        # Store the colours that will be used in the plot
        color <- colour_matrix[levels(factor(report_plot[,1],levels = unique(report_plot[,1]))),2]
        
        # Save the boxplot object in the list
        my_name <- "Closest Reference Cluster"
        boxplot_list[[cluster]] <- list()
        boxplot_list[[cluster]] <- my_boxplot(my_name,report_plot,label,abundance,color)
        
        pointplot_list[[cluster]] <- list()
        pointplot_list[[cluster]] <- my_point_boxplot(my_name,report_plot,label,abundance,color)
        
        violinplot_list[[cluster]] <- list()
        violinplot_list[[cluster]] <- my_violinplot(my_name,report_plot,label,abundance,color)
        
      }
      
      for (cluster in 1:reference_clusters) {
        
        # Index indicating the page of the report
        page <- page+1
        
        
        if (nlevels(as.factor(plot_matrix[plot_matrix$Nearest.reference.cluster==cluster,ncol(plot_matrix)])) == 1) {pvaluesb <- data.frame("NA")
        colnames(pvaluesb) <- "pvalues"
        } else {
          # Available combination for the pairwise comparisons
          combinations <- combs(levels(as.factor(plot_matrix[plot_matrix$Nearest.reference.cluster==cluster,ncol(plot_matrix)])), 2)
          
          # Vector containing the p-values
          pvaluesb <- c()
          
          # Calculate the PERMANOVA p-values
          for (g in 1:nrow(combinations)){ 
            # Choose the groups for the pairwise comparisons
            groups <- rownames(plot_matrix[plot_matrix[,ncol(plot_matrix)] %in% combinations[g,],])
            # Calculate pairwise PERMANOVA p-values
            adonis <- adonis(unifract_dist[c(groups),c(groups)]~plot_matrix[c(groups),ncol(plot_matrix)])[[1]][6][[1]][1]
            # Form the matrix with the p-values
            pvaluesb <- rbind(pvaluesb,c(paste0(combinations[g,1],"-",combinations[g,2]),adonis))
          }
          
          # Performing fdr correction
          p.adjustb <- round(p.adjust(pvaluesb[,2], method = "BH"),4)
          
          # Add the adjusted p-values in the pvalues matrix
          pvaluesb <- data.frame(pvaluesb,p.adjustb)
          
          # Name the columns
          colnames(pvaluesb) <- c("Groups","p-value","Adj. p-value")
        }
        
        if (nlevels(as.factor(plot_matrix[plot_matrix$Nearest.reference.cluster == cluster,ncol(plot_matrix)])) != 1){
          # Convert p-values column into numeric
          pvaluesb[,2] <- as.numeric(pvaluesb[,2])
          # Order the p-values
          pvaluesb <- pvaluesb[order(pvaluesb[,2]),]
          # Choose the first 10 most significant p-values
          if (nrow(pvaluesb)>10){ pvaluesb2 <- pvaluesb[1:10,] } else { pvaluesb2 <- pvaluesb }
        }else {
          pvaluesb2 <- pvaluesb
          }
        
        # Create a gtable containing text grobs 
        pvaluetable <- gtableGrob("PERMANOVA test - pairwise",pvaluesb2,setDT="NO",size=9)
        
        
        # Default margins
        par(mar= c(5.1, 4.1, 4.1, 2.1))
        
        
        if (length(which(plot_matrix$Nearest.reference.cluster == cluster)) > 2){
          if (length(Test_name) <= 2){
            # Set the layout of the cover page (3 rows, 1 column)
            layout(matrix(c(1, 1, 1,
                            2, 2, 2,
                            3, 3, 3), nrow=3, byrow=TRUE),heights=c(1,3,0.3) )
          } else {
            # Set the layout of the cover page (3 rows, 1 column)
            layout(matrix(c(1, 1, 1,
                            2, 2, 2,
                            3, 3, 3), nrow=3, byrow=TRUE),heights=c(1.5,3,0.3) )
          }  
          
          
          # Create an empty plot
          plot.new()
          # Title of the page
          title(main=paste("Reference Cluster",cluster),cex.main = 1.7,font.main = 2)
          
          
          # Text that will we used as description for every reference cluster
          text <- c()
          
          for(group in 1:length((Test_name))) { 
            # Subset the samples of each test group
            temporary <- plot_matrix[rownames(input_table[which(input_table[,1] == levels[group]),]),]
            # Write how many samples of each test group is closer to every reference cluster
            paste2 <- paste(nrow(temporary[temporary$Nearest.reference.cluster == cluster,]),"of the", nrow(plot_matrix[rownames(input_table[which(input_table[,1] == levels[group]),]),]),
                            "samples of the", levels[group], " group are closer to the reference cluster",cluster,"\n")
            # Final text
            text <- paste(text,paste2)}
          
          # Information about the number of the test samples that are closer to each reference cluster
          text(x=0.5, y = 0.25, text,  cex = 1.3, col = "black",font=3)
          
          
          # MDS plot with the reference cluster and the test sample that are closer to it
          sclass(unifract_dist[plot_matrix$Nearest.reference.cluster == cluster,plot_matrix$Nearest.reference.cluster == cluster],factor(plot_matrix[plot_matrix$Nearest.reference.cluster == cluster,ncol(plot_matrix)],levels = unique(plot_matrix[plot_matrix$Nearest.reference.cluster == cluster,ncol(plot_matrix)])),
                 NULL, colour_matrix[levels(factor(plot_matrix[plot_matrix$Nearest.reference.cluster == cluster,ncol(plot_matrix)],levels = unique(plot_matrix[plot_matrix$Nearest.reference.cluster == cluster,ncol(plot_matrix)]))),2])
          
          # Change the margins in order to be created blank space at the end of the page
          par(mar = c(0.1, 0.1, 0.1, 0.1))
          
          # Create an empty plot
          plot.new()
          
          
          # Print the MDS plot in png format
          png_plot(paste0("DiBaAn-page ",page, " MDS plot.png"),sclass(unifract_dist[plot_matrix$Nearest.reference.cluster==cluster,plot_matrix$Nearest.reference.cluster==cluster],factor(plot_matrix[plot_matrix$Nearest.reference.cluster==cluster,ncol(plot_matrix)],levels = unique(plot_matrix[plot_matrix$Nearest.reference.cluster==cluster,ncol(plot_matrix)])),
                                                                       NULL, colour_matrix[levels(factor(plot_matrix[plot_matrix$Nearest.reference.cluster==cluster,ncol(plot_matrix)],levels = unique(plot_matrix[plot_matrix$Nearest.reference.cluster==cluster,ncol(plot_matrix)]))),2]))
          
          
          # Default margins
          par(mar= c(5.1, 4.1, 4.1, 2.1))
          
          
          # Index indicating the page of the report
          page <- page+1
          
          
          # tree with the samples of the reference cluster and the test samples that are closer to it
          tree_plot <- tree(unifract_dist[plot_matrix$Nearest.reference.cluster == cluster,plot_matrix$Nearest.reference.cluster==cluster],factor(plot_matrix[plot_matrix$Nearest.reference.cluster == cluster,ncol(plot_matrix)],levels = unique(plot_matrix[plot_matrix$Nearest.reference.cluster == cluster,ncol(plot_matrix)])),
                            colour_matrix[levels(factor(plot_matrix[plot_matrix$Nearest.reference.cluster == cluster,ncol(plot_matrix)],levels = unique(plot_matrix[plot_matrix$Nearest.reference.cluster == cluster,ncol(plot_matrix)]))),2])
          
          # Title of the plot
          text3 <- paste("Phylogram")
          # Create a text grob for the title
          tgrob3 <- text_grob(text3, face = "bold", color = "Black")
          
          # Choose the selected type of plot
          if (plot_type == "Boxplots") {plot <-boxplot_list[[cluster]] } else if (plot_type == "Point plots") {plot <- pointplot_list[[cluster]] } else {plot <-violinplot_list[[cluster]]}
          
          
          # Grid arrange all the above elements
          grid.arrange(arrangeGrob(plot,pvaluetable,ncol=2,nrow=1,widths = c(1,1)),tgrob3,tree_plot,nrow = 3,as.table = T,heights=c(1.55,0.1,1))
          
          # Write the PERMANOVA p-values
          write_table(paste0("DiBaAn-page ",page," PERMANOVA pvalues.tab"),pvaluesb,"pvalues")
          
          # Print the Phylogram  in png format
          png_plot(paste0("DiBaAn-page ",page," Phylogram.png"),tree_plot)
          
          # Print the plot in png format
          png_plot(paste0("DiBaAn-page ",page," ",plot_type,".png"),plot)
          
        } 
      }
      
      
      dev.off()
      
    }
    
    
    
    #################################################
    ########## De novo-based report #################
    #################################################
    
    if (index == 0 & Test_name != "None") {
      
      # Create the lists where the plots will be saved
      boxplot_list <- list()
      pointplot_list <- list()
      violinplot_list <- list()
      
      
      # Open the PDF
      pdf("De novo clustering report.pdf")
      
      
      
      ########### Cover Page ##########
      
      # Index indicating the page of the report
      page <- 1
      
      
      # Set the layout of the cover page
      layout( matrix(c(1,1), ncol=1) )
      # Default margins
      par(mar= c(5.1, 4.1, 8, 2.1))
      
      # Create an empty plot
      plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      # Title of the plot
      title(main=paste0("De novo clustering Analysis\n","(DeNoAn-",central_point,")"),cex.main = 1.8,font.main = 2)
      # Description of this PDF
      text(x = 0.5, y = 0.5, paste("This report includes results derived from the analysis that is based \n",
                                   "on the de novo clustering of both the reference and test samples\n",
                                   "Both the Reference and test samples were clustered into distinct\n groups ",
                                   "based on the results of the PAM algorithm.\n\n\n",
                                   "* All the elements of this report (plots and tables) will be found in the \n results folder with the prefix DeNoAn-"),cex = 1.1, col = "black",font = 3)
      # Default margins
      par(mar= c(5.1, 4.1, 4.1, 2.1))
      
      
      ############## Page 2 ##########
      
      # Index indicating the page of the report
      page <- page+1
      
      
      # Information about the refernce clusters
      clustering_text <- paste0("The reference dataset (",paste(reference_name,collapse="+"),") was clustered into ", reference_clusters ," groups:\n")
      # Vector where the propabilities will be placed
      propability <- c()
      # Vector where the number of samples will be placed
      number_samples <- c()
      
      for (i in 1:reference_clusters){
        clustering_text <- paste(clustering_text,paste(reference_name,collapse="+") ,i, "has", length(which(denovo_clusters==paste(paste(reference_name,collapse="+"),i))),"samples.\n")
        if (reference_clusters > 1) {
          number_samples <- c(number_samples,length(which(denovo_clusters==paste(paste(reference_name,collapse="+"),i))))
          propability <- c(propability,1/reference_clusters)
        }
      }
      
      if (reference_clusters > 1) {
        
        # Calculate the Chi-square p-values
        clustering_text <- paste(clustering_text,paste("The p-value of the Chi square (Goodness of fit) test is",round(chisq.test(number_samples,p=propability)$p.value,4)),"*\n")
      }
      
      # Information about the de novo clustering of test groups
      for(i in 1:length(levels)) {
        clustering_text <- paste(clustering_text,"\n", "The",levels[i] ,"group was clustered into",test_clusters[i],"groups:\n" )
        # Vector where the propabilities will be placed
        propability <- c()
        # Vector where the number of samples will be placed
        number_samples <- c()
        
        for(j in 1:test_clusters[i]) {
          clustering_text <- paste(clustering_text,levels[i],j,"has", length(which(denovo_clusters==paste(levels[i],j))),"samples.\n")
          if (test_clusters[i] > 1) {
            number_samples <- c(number_samples,length(which(denovo_clusters==paste(levels[i],j))))
            propability <- c(propability,1/test_clusters[i])
          }
        }
        
        if (test_clusters[i] > 1) {
          
          # Calculate the Chi-square p-values
          clustering_text <- paste(clustering_text,paste("The p-value of the Chi square (Goodness of fit) test is",round(chisq.test(number_samples,p=propability)$p.value,4)),"*\n")
        }
      }
      
      if (length(Test_name) <= 2){
        # Create a text grob for the clustering
        tgrob <- text_grob(clustering_text, face = "italic", color = "Black",size=13)
      } else {
        # Create a text grob for the clustering
        tgrob <- text_grob(clustering_text, face = "italic", color = "Black",size=9)
      }
      
      # Create a text grob for the ttle of the page
      tgrob2 <- text_grob("De novo clustering of the reference and \n test groups",face = "bold.italic", color = "Black",size=15)
      
      text <-  "* The Chi square (Goodness of fit) test checks if there is discrepancy between the observed number of observations \n of each cluster and the number of observation in case the samples were distributed uniformly."
      
      # Create a text grob for the clustering
      tgrob3 <- text_grob(text, face = "italic", color = "Black",size=9)
      
      # Arrange in the grid the above elements
      grid.arrange(tgrob2,tgrob,tgrob3,nrow=3,heights=c(0.5,3,0.5))
      
      
      ############ Page 3 #############
      
      # Index indicating the page of the report
      page <- page+1
      
      
      
      # Set the layout of the cover page (3 rows, 1 column)
      layout(matrix(c(1, 1, 1,
                      2, 2, 2,
                      3, 3, 3), nrow=3, byrow=TRUE),heights=c(1,3,0.3) )
      
      # Create an empty plot
      plot.new()
      # Title of the 3rd page
      title(main="MDS Plot",cex.main = 1.8,font.main = 2)
      
      # Description of the MDS plot
      text(x=0.5, y = 0.30, paste0("MDS plot presenting the groups which have formed from the de novo clustering of  \n the ", 
                                   "reference (", paste(reference_name,collapse="+"), ") and the test samples (", paste(Test_name,collapse="+"), ") \n"), cex = 1.5, col = "black")
      
      # MDS plot presenting the clusters of the de novo clustering of both the reference and test groups
      sclass(unifract_dist,factor(denovo_clusters,levels=unique(denovo_clusters)),
             NULL,colour_matrix[levels(factor(denovo_clusters,levels=unique(denovo_clusters))),2])
      
      
      # Change the margins in order to be created blank space at the end of the page
      par(mar = c(0.1, 0.1, 0.1, 0.1))
      # Create an empty plot
      plot.new()
      
      # Default margins
      par(mar= c(5.1, 4.1, 4.1, 2.1))
      
      # Print the MDS plot in png format
      png_plot(paste0("DeNoAn-page ",page, " MDS Plot.png"),sclass(unifract_dist,factor(denovo_clusters,levels=unique(denovo_clusters)),
                                                                   NULL,colour_matrix[levels(factor(denovo_clusters,levels=unique(denovo_clusters))),2]))
      
      
      
      ############### Page 4 ##########
      
      # Index indicating the page of the report
      page <- page+1
      
      
      # Title of the 4th page 
      text1 <- "Descriptive Statistics"
      
      # Create a text grob for the title
      tgrob <- text_grob(text1,face = "bold.italic", color = "Black",size=15)
      
      
      # Arrange in the grid the above elements
      grid.arrange(tgrob,reference_statistis_table,statistic_table_test_clusters,nrow=3,heights=c(0.4,1.3,1.6))
      
      
      
      ############### Page 5-(undefined) ##########
      
      if(chi_square == 1){
        
        # Index indicating the page of the report
        page <- page+1
        
        # Tile of the chi-square analysis
        text1 <- "Chi-square test \nGoodness of Fit analysis"
        
        # Create a text grob for the title
        tgrob <- text_grob(text1,face = "bold.italic", color = "Black",size=20)
        
        # Text of the footnote
        text2 <- paste("*For the next",(length(exploratory_columns)*length(Test_name)+length(exploratory_columns)),"pages")
        
        # Create a text grob for the footnote
        tgrob2 <- text_grob(text2,face = "italic", color = "Black",size=12)
        
        # Arrange in the grid the above elements
        grid.arrange(tgrob,tgrob2, nrow=2,heights=c(2,0.5))
        
        
        #------------------------------------------- Chi-Square Analysis ---------------------------------------#
        
        for (k in 1:length(exploratory_columns)){
          
          # Index indicating the page of the report
          page <- page+1
          
          # Check if the column is character of factor
          if (is.character(meta_file[,exploratory_columns[k]]) == TRUE | is.factor(meta_file[,exploratory_columns[k]]) == TRUE ){
            
            # Index indicating the page of the report
            page <- page+1
            
            # Tile of the page
            text1 <- paste0("'",exploratory_columns[k],"' Column")
            
            # Create a text grob for the title
            tgrob <- text_grob(text1,face = "bold.italic", color = "Black",size=20)
            
            # Main text of the page
            text2 <- paste0("The distribution of the samples across the clusters of the test group(s)\n is presented in Table 1.\n
                 Table 2 presents the p-values of the Chi-square test of the observed distribution \nof the counts across the clusters in comparison with the expected \ndistribution if the counts were uniformly distributed.\n
                 Table 3 contains the p-values of the pairwise comparisons.\n
                 The tables contain only the 10 most significant p-values. The rest of the p-values\n
                           have been printed in the results folder."
            )
            
            # Create a text grob for the text
            tgrob2 <- text_grob(text2,face = "italic", color = "Black",size=12)
            
            # Create a text grob for the footnote
            text3 <- "*If the number of counts is limited or zero counts occurs, then the Chi-square approximation \n may be incorrect."  
            
            # Create a text grob for the footnote
            tgrob3 <- text_grob(text3,face = "italic", color = "Black",size=11)
            
            # Arrange in the grid the above elements
            grid.arrange(tgrob,tgrob2,tgrob3 ,nrow=3)
            
            # Calculate the chi-square p-values for every test group
            for(i in 1:length(Test_name)){
              
              # Index indicating the page of the report
              page <- page+1
              
              # Subset only the names of the test Samples 
              test_names <- rownames(input_table[input_table[,1]==Test_name[i],])
              # Extract the factors of the column
              factors <- factor(meta_file[test_names,exploratory_columns[k]],levels=unique(meta_file[test_names,exploratory_columns[k]]))
              
              if (nlevels(factors) > 0){
                # Form the matrix where the counts of the observations will be placed
                chi_matrix <- data.frame(matrix(ncol = (nlevels(factors)) , nrow = test_clusters[i]))
                # Naming the rows of p_value data frame
                rownames(chi_matrix) <- paste(Test_name[i],1:test_clusters[i])
                # Naming the columns of p_value data frame
                colnames(chi_matrix) <- levels(factors)
                
                
                for (ii in 1:test_clusters[i]){
                  # Subset only the names of the test Samples for every cluster
                  names2 <- rownames(mapping_file[mapping_file[,1]==paste(Test_name[i],ii),])
                  for (iii in 1:nlevels(factors)){
                    # Place the counts in the Chi_matrix
                    chi_matrix[ii,iii] <- length(which(meta_file[names2,exploratory_columns[k]] == levels(factors)[iii]))
                  }
                }
                
                # Replace zero counts with 0.01
                chi_matrix[chi_matrix==0] <- 0.01
                
                
                # Performing Chi-square test(Goodness of fit) for each column
                exp_col <- apply(chi_matrix,2,function(x) chisq.test(x)$p.value)
                
                if (ncol(chi_matrix) > 1){
                  # Performing Chi-square test(Goodness of fit) for each row
                  exp <-apply(chi_matrix,2,sum)/test_clusters[i] 
                  exp_r <- apply(chi_matrix,1,function(x) chisq.test(x,p=exp,rescale.p = TRUE)$p.value)
                  
                  
                  # Performing fdr correction
                  adjustb <- round(p.adjust(c(exp_col,exp_r), method = "BH"),4)
                  # Add the adjusted p-values in the expected matrix
                  expected <- data.frame(cbind(c(exp_col,exp_r),adjustb))
                  # Naming the columns of p_value data frame
                  colnames(expected) <- c( "P-values","Adj p-values")
                  # Naming the rows of p_value data frame
                  rownames(expected) <- c(names(exp_col),names(exp_r))
                } else {
                  # Performing fdr correction
                  adjustb <- round(p.adjust(exp_col, method = "BH"),4)
                  # Add the adjusted p-values in the expected matrix
                  expected <- data.frame(cbind(exp_col,adjustb))
                  # Naming the columns of p_value data frame
                  colnames(expected) <- c( "P-values","Adj p-values")
                  # Naming the rows of p_value data frame
                  rownames(expected) <- c(names(exp_col))
                }
                
                
                # Order the p-values
                expected <- expected[order(expected[,2]),]
                # Choose the first 10 most significant p-values
                if (nrow(expected) > 10){expected2 <- expected[1:10,]} else { expected2 <- expected }
                
                if (ncol(chi_matrix) > 1) {
                  # Find all the available pairwise combinations
                  combs <- combs(levels(factors), 2)
                  # Save the names of that combinations
                  group_names <- apply(combs,1,function(x)paste0(x[1],"-",x[2]))
                  
                  
                  # chi_p_value -> a vector where the p values will be stored
                  chi_p_value <- c(); test <- c()
                  for(chi in 1:nrow(combs)) { 
                    # Performed the Chi square test
                    chi_p_value[chi] <- round(chisq.test(chi_matrix[,combs[chi,2]],p=chi_matrix[,combs[chi,1]],rescale.p = TRUE)$p.value,4)
                  }
                  
                  
                  #--- Rows pairwise comparisons ----#  
                  # Find all the available pairwise combinations
                  combs <- combs(rownames(chi_matrix), 2)
                  # Save the names of that combinations
                  group_names_rows <- apply(combs,1,function(x)paste0(x[1],"-",x[2]))
                  
                  
                  # chi_p_value -> a vector where the p values will be stored
                  chi_p_value_rows <- c(); test <- c()
                  for(chi in 1:nrow(combs)) { 
                    # Performed the Chi square test
                    chi_p_value_rows[chi] <- round(chisq.test(as.numeric(chi_matrix[combs[chi,2],]),p=as.numeric(chi_matrix[combs[chi,1],]),rescale.p = TRUE)$p.value,4)
                  }
                  
                  
                  # Replace 0.01 with 0
                  chi_matrix[chi_matrix == 0.01] <- 0
                  
                  # Forming the data frame containing the p values
                  p_value <- data.frame(group_names)
                  # Place the p-values
                  p_value[,2] <- chi_p_value
                  # Naming the columns of p_value data frame
                  colnames(p_value) <- c("Groups","p value")
                  
                  # Forming the data frame containing the p values
                  p_value_rows <- data.frame(group_names_rows)
                  # Naming the rows of the p_value data frame
                  p_value_rows[,2] <- chi_p_value_rows
                  # Naming the columns of p_value data frame
                  colnames(p_value_rows) <- c("Groups","p value")
                  
                  # rbind all the p-values
                  p_value <- rbind(p_value,p_value_rows)
                  
                  # Convert p-values column into numeric
                  p_value[,2] <- as.numeric(p_value[,2])
                  # Performing fdr correction
                  p.adjust <- round(p.adjust(p_value[,2], method = "BH"),4)
                  # Add the adjusted p-values in the pvalues matrix
                  p_value[,3] <- p.adjust
                  # Naming the columns of p_value data frame
                  colnames(p_value) <- c("Groups","p value","Adj p-values")
                  # Order the p-values
                  p_value <- p_value[order(p_value[,3]),]
                  # Choose the first 10 most significant p-values
                  if (nrow(p_value) > 10) {p_value2 <- p_value[1:10,]} else {p_value2 <- p_value}
                  
                } else {
                  
                  #--- Rows pairwise comparisons ----#  
                  # Find all the available pairwise combinations
                  combs <- combs(rownames(chi_matrix), 2)
                  # Save the names of that combinations
                  group_names <- apply(combs,1,function(x)paste0(x[1],"-",x[2]))
                  
                  
                  # chi_p_value -> a vector where the p values will be stored
                  chi_p_value_rows <- c(); test <- c()
                  for(chi in 1:nrow(combs)) { 
                    # Performed the Chi square test
                    chi_p_value_rows[chi] <- round(chisq.test(as.numeric(chi_matrix[combs[chi,2],]),p=as.numeric(chi_matrix[combs[chi,1],]),rescale.p = TRUE)$p.value,4)
                  }
                  
                  
                  # Replace 0.01 with 0
                  chi_matrix[chi_matrix == 0.01] <- 0
                  
                  # Forming the data frame containing the p values
                  p_value <- data.frame(group_names_rows)
                  # Naming the rows of the p_value data frame
                  p_value[,2] <- chi_p_value_rows
                  # Naming the columns of p_value data frame
                  colnames(p_value) <- c("Groups","p value")
                  # Convert p-values column into numeric
                  p_value[,2] <- as.numeric(p_value[,2])
                  # Performing fdr correction
                  p.adjust <- round(p.adjust(p_value[,2], method = "BH"),4)
                  # Add the adjusted p-values in the pvalues matrix
                  p_value[,3] <- p.adjust
                  # Naming the columns of p_value data frame
                  colnames(p_value) <- c("Groups","p value","Adj p-values")
                  # Order the p-values
                  p_value <- p_value[order(p_value[,3]),]
                  # Choose the first 10 most significant p-values
                  if (nrow(p_value)>10) {p_value2 <- p_value[1:10,]} else {p_value2 <- p_value}
                  
                }
                
                # Print the counts table in the results folder
                write_table(paste0("DeNoAn-page ",page, " Counts per cluster (Table 1).tab"),chi_matrix,"table")
                # Print the Expected-observed p-values table in the results folder
                write_table(paste0("DeNoAn-page ",page, " Expected-observed (Table 2).tab"),expected,"pvalues")
                # Print the pairwise p-values table in the results folder
                write_table(paste0("DeNoAn-page ",page, " Pairwise comparisons (Table 3).tab"),p_value,"pvalues")
                
                if (sum(nchar(colnames(chi_matrix))) < 62){
                  
                  # Tile of the page
                  text1 <- paste(exploratory_columns[k],"-",Test_name[i])
                  
                  # Create a text grob for the title
                  tgrob <- text_grob(text1,face = "bold.italic", color = "Black",size=15)
                  
                  # Tile of the p-value tables
                  text2 <- "P-values of the Chi-square test (Goodness of Fit)"
                  
                  # Create a text grob for the title
                  tgrob2 <- text_grob(text2,face = "bold.italic", color = "Black",size=12)
                  
                  # Arrange in the grid the above elements
                  grid.arrange(tgrob,gtableGrob("Counts per cluster (Table 1)",chi_matrix,setDT="YES",size=12), 
                               arrangeGrob(tgrob2,arrangeGrob(gtableGrob("Expected-observed (Table 2)",expected2,setDT="YES",size=12),gtableGrob("Pairwise comparisons (Table 3)",p_value2,setDT="NO",size=12),nrow=1),nrow=2,heights=c(0.3,1))
                               , nrow = 3,heights=c(0.4,1,2))
                } else {
                  text1 <- paste(exploratory_columns[k],"-",Test_name[i])
                  
                  # Create a text grob for the title
                  tgrob <- text_grob(text1,face = "bold.italic", color = "Black",size=15)
                  
                  # Tile of the p-value tables
                  text2 <- "The tables are too big to be printed!!\nYou can find them in .tab format in the results folder."
                  
                  # Create a text grob for the title
                  tgrob2 <- text_grob(text2,face = "bold.italic", color = "Black",size=12)
                  
                  grid.arrange(tgrob,tgrob2,nrow=2)
                }
                
                
                
              } else {
                # Tile of the page
                text1 <- paste0("'",exploratory_columns[k],"' Column")
                
                # Create a text grob for the title
                tgrob <- text_grob(text1,face = "bold.italic", color = "Black",size=20)
                
                # Main text of the page
                text2 <- paste0("It was not possible to perform the chi-square test for the '",Test_name[i],"'. There are no Factors")
                
                # Create a text grob for the text
                tgrob2 <- text_grob(text2,face = "italic", color = "Black",size=12)
                
                
                # Arrange in the grid the above elements
                grid.arrange(tgrob,tgrob2 ,nrow=2)
              }  
            }
            
          } else {
            
            
            # Tile of the page
            text1 <- paste0("'",exploratory_columns[k],"' Column")
            
            # Create a text grob for the title
            tgrob <- text_grob(text1,face = "bold.italic", color = "Black",size=20)
            
            # Main text of the page
            text2 <- paste0("It was not possible to perform the chi-square test for the '",exploratory_columns[k],"'. The variable must be categorical ")
            
            # Create a text grob for the text
            tgrob2 <- text_grob(text2,face = "italic", color = "Black",size=12)
            
            
            # Arrange in the grid the above elements
            grid.arrange(tgrob,tgrob2 ,nrow=2)
          }
        }
      }
      
      
      ############### Page (self-medoid) ############
      
      # Index indicating the page of the report
      page <- page+1
      
      
      # Description of the list_box3
      text <- paste0(plot_type," displaying the distribution of the distances around their own self-",central_point,".\n",
                     "The plots can provide information about the compactness of the groups.")
      
      
      # Create a text grob for the Description of the list_box3
      tgrob <- text_grob(text, face = "italic", color = "Black",size=12)
      
      # Choose the selected type of plot
      if (plot_type == "Boxplots") {plot <-list_box3 } else if (plot_type == "Point plots") {plot <- list_point3 } else {plot <- list_violin3}
      
      # Arrange in the grid the above elements
      grid.arrange(arrangeGrob(plot,pvaltable3,nrow = 1,ncol=2,widths = c(1.1,1)),tgrob,nrow = 2,as.table = T,heights=c(2.5,1))
      
      
      # Print the plot in png format
      png_plot(paste0("DeNoAn-page ",page," ",plot_type,".png"),plot)
      
      ############# Page (closest reference cluster) ##########
      
      # Index indicating the page of the report
      page <- page+1
      
      
      
      # Title of the 6th page 
      text1 <- paste("Distances From the Closest Reference",central_point)
      
      # Create a text grob for the title
      tgrob <- text_grob(text1,face = "bold.italic", color = "Black",size=15)
      
      # Description of the observations matrix
      text2 <- paste0("The closest reference ",central_point, " was found for every sample of the test groups.\n 
                   The above matrix presents the number of the sample of each test group \n from the closest reference cluster(",central_point,").\n
                   Based on this information in the next page the ",plot_type," of the \n distances from the nearest reference ",central_point," were formed.")
      
      
      # Create a text grob for the description of the observations matrix
      tgrob2 <- text_grob(text2,face = "italic", color = "Black",size=13)
      
      # Arrange in the grid the above elements
      grid.arrange(tgrob,observation,tgrob2,nrow=3,heights=c(0.4,1,1.5))
      
      
      
      ############### Page (Boxplot-closest reference medoid) ###########  
      
      # Index indicating the page of the report
      page <- page+1
      
      
      # Description of the list_box5
      text <- paste0(plot_type," displaying the distribution of the distances \n",
                     "from their nearest reference ",toupper(central_point),".\n",
                     "The columns represent the clusters that have been derived \n",
                     "from the de novo clustering of the reference and tests samples.\n",
                     "\n",
                     
                     "*Only the top-10 most significant p-values are displayed.\n"," The rest of them are at the results folder.")
      
      # Create a text grob for the list_box5
      tgrob <- text_grob(text, face = "italic", color = "Black",size=13)
      
      # Choose the selected type of plot
      if (plot_type == "Boxplots") {plot <-list_box5 } else if (plot_type == "Point plots") {plot <- list_point5 } else {plot <- list_violin5}
      
      
      # Arrange in the grid the above elements
      grid.arrange(arrangeGrob( plot,pvaltable5,nrow = 1,ncol=2,widths = c(1.1,1)),tgrob,nrow = 2,as.table = T,heights=c(2.5,1))
      
      # Print the plot in png format
      png_plot(paste0("DeNoAn-page ",page," ",plot_type,".png"),plot)
      
      
      ############### Page (median distances table) ############
      
      # Index indicating the page of the report
      page <- page+1
      
      
      # Title of the page 
      text1 <- "Distances from the closest reference cluster"
      # Create a text grob for the title
      tgrob <- text_grob(text1,face = "bold.italic", color = "Black",size=15)
      
      # Description of the matrix of the median distances
      text2 <- paste("For each group the median distance of the sample from \n every reference",central_point, "was calculated. \n Then, the closest reference cluster was found.\n
              Based on the information of the above matrix, the",plot_type,"of the distances from \n the closest reference GROUP were formed.\n
              * More cluster-related results will be found at the results folder.")
      # Create a text grob for the median matrix
      tgrob2 <- text_grob(text2,face = "italic", color = "Black",size=13)
      
      
      # Arrange in the grid the above elements
      grid.arrange(tgrob,median_matrix_table,tgrob2,nrow=3,heights=c(0.4,1,1.5))
      
      
      
      ############# Page (Boxplots closest reference cluster) ############
      
      # Index indicating the page of the report
      page <- page+1
      
      
      
      # Description of the list_box4
      text <- paste(plot_type,"displaying the distribution of the distances of the test samples \n",
                    "from their nearest reference cluster.\n",
                    "\n",
                    "On the next pages each reference cluster will be examined seperatly.")
      # Create a text grob for the list_box4
      tgrob <- text_grob(text, face = "italic", color = "Black",size=13)
      
      # Choose the selected type of plot
      if (plot_type == "Boxplots") {plot <-list_box4 } else if (plot_type == "Point plots") {plot <- list_point4 } else {plot <- list_violin4}
      
      
      # Arrange in the grid the above elements
      grid.arrange(arrangeGrob(plot,pvaltable4,nrow = 1,ncol=2,widths = c(1.1,1)),tgrob,nrow = 2,as.table = T,heights=c(2.5,1))
      
      # Print the plot in png format
      png_plot(paste0("DeNoAn-page ",page," ",plot_type,".png"),plot)
      
      
      ############ Pages (individual cluster analysis) ###############
      
      for (cluster in 1:reference_clusters){
        
        # Subset the samples
        report_plot <- data.frame(plot_df_b[dist_plot_4[,3] == cluster,])
        # Name the rows
        rownames(report_plot) <- rownames(dist_plot_4[dist_plot_4[,3]==cluster,])
        # Convert label column into factor
        report_plot$label <- factor(report_plot$label,levels = unique(report_plot$label))
        
        # Store the colours that will be used in the plot
        color <- colour_matrix[levels(factor(report_plot[,1],levels = unique(report_plot[,1]))),2]
        
        
        # Save the boxplot object in the list
        my_name <- "Closest Reference Cluster"
        boxplot_list[[cluster]] <- list()
        boxplot_list[[cluster]] <- my_boxplot(my_name,report_plot,label,abundance,color)
        
        pointplot_list[[cluster]] <- list()
        pointplot_list[[cluster]] <- my_point_boxplot(my_name,report_plot,label,abundance,color)
        
        violinplot_list[[cluster]] <- list()
        violinplot_list[[cluster]] <- my_violinplot(my_name,report_plot,label,abundance,color)
      }
      
      for (cluster in 1:reference_clusters) {
        if (length(rownames(dist_plot_4[dist_plot_4$cluster==cluster,]))>2){
          
          # Index indicating the page of the report
          page <- page+1
          
          if (nlevels(as.factor(dist_plot_4[dist_plot_4$cluster==cluster,1])) == 1) {
            pvaluesb <- data.frame("NA")
            colnames(pvaluesb) <- "pvalues"
          } else {
            # Available combination for the pairwise comparisons
            combinations <- combs(levels(as.factor(dist_plot_4[dist_plot_4$cluster == cluster,1])), 2)
            
            
            # Vector containing the p-values
            pvaluesb <- c()
            
            # Calculate the PERMANOVA p-values
            for (g in 1:nrow(combinations)){ 
              # Choose the groups for the pairwise comparison
              groups <- rownames(plot_matrix[plot_matrix[,1] %in% combinations[g,],])
              # Calculate pairwise PERMANOVA p-values
              adonis <- adonis(unifract_dist[c(groups),c(groups)]~plot_matrix[c(groups),1])[[1]][6][[1]][1]
              # Form the matrix with the p-values
              pvaluesb <- rbind(pvaluesb,c(paste0(combinations[g,1],"-",combinations[g,2]),adonis))
            }
            
            # Performing fdr correction
            p.adjustb <- round(p.adjust(pvaluesb[,2], method = "BH"),4)
            
            # Add the adjusted p-values in the pvalues matrix
            pvaluesb <- data.frame(pvaluesb,p.adjustb)
            
            # Name the columns
            colnames(pvaluesb) <- c("Groups","p-value","Adj. p-value")
          }
          
          if (nlevels(as.factor(dist_plot_4[dist_plot_4$cluster==cluster,1])) != 1){
            # Convert p-values column into numeric
            pvaluesb[,2] <- as.numeric(pvaluesb[,2])
            # Order the p-values
            pvaluesb <- pvaluesb[order(pvaluesb[,2]),]
            
            # Choose the first 10 most significant p-values
            if (nrow(pvaluesb) > 10){ pvaluesb2 <- pvaluesb[1:10,]} else {pvaluesb2 <- pvaluesb}
          } else { 
            pvaluesb2 <- pvaluesb
            }
          
          # Create a gtable containing text grobs 
          pvaluetable <- gtableGrob("PEMANOVA test - pairwise",pvaluesb2,setDT="NO",size=9)
          
          temporary <- rownames(median_matrix[which(median_matrix[,ncol(median_matrix)]==cluster),])
          paste2 <- paste0("According to the median distances from the reference ",central_point, ".\n", paste(paste(temporary,collapse=",")))
          
          
          
          # Set the layout of the cover page (3 rows, 1 column)
          layout(matrix(c(1, 1, 1,
                          2, 2, 2,
                          3, 3, 3), nrow=3, byrow=TRUE),heights=c(1.1,3,0.3) )
          
          
          # Create an empty plot
          plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
          # Title of the page
          title(main=paste("Reference Cluster",cluster),cex.main = 2)
          
          if (length(temporary) == 0){
            text(x = 0.5, y = 0.3, paste("None of the test groups are closer to Reference cluster ",cluster,"."),cex = 1.5, col = "black",font = 3)
          } else if (length(temporary)==1){
            text(x = 0.5, y = 0.3, paste0(paste2," is closer to Reference cluster ", cluster,"."),cex = 1.5, col = "black",font = 3)
          } else {
            text(x = 0.5, y = 0.3, paste0(paste2," are closer to Reference cluster ", cluster,"."),cex = 1.5, col = "black",font = 3)
          }
          
          
          # MDS presenting the every reference cluster with the closest test clusters
          sclass(unifract_dist[rownames(dist_plot_4[dist_plot_4$cluster==cluster,]),rownames(dist_plot_4[dist_plot_4$cluster==cluster,])],factor(denovo_clusters[rownames(dist_plot_4[dist_plot_4$cluster==cluster,])],levels=unique(denovo_clusters[rownames(dist_plot_4[dist_plot_4$cluster==cluster,])])),
                 NULL,colour_matrix[levels(factor(denovo_clusters[rownames(dist_plot_4[dist_plot_4$cluster==cluster,])],levels=unique(denovo_clusters[rownames(dist_plot_4[dist_plot_4$cluster==cluster,])]))),2])
          
          
          # Print the MDS plot in png format
          png_plot(paste0("DeNoAn-page ",page," MDS Plot.png"),
                   sclass(unifract_dist[rownames(dist_plot_4[dist_plot_4$cluster==cluster,]),rownames(dist_plot_4[dist_plot_4$cluster==cluster,])],factor(denovo_clusters[rownames(dist_plot_4[dist_plot_4$cluster==cluster,])],levels=unique(denovo_clusters[rownames(dist_plot_4[dist_plot_4$cluster==cluster,])])),
                          NULL,colour_matrix[levels(factor(denovo_clusters[rownames(dist_plot_4[dist_plot_4$cluster==cluster,])],levels=unique(denovo_clusters[rownames(dist_plot_4[dist_plot_4$cluster==cluster,])]))),2]))
          
          # Change the margins in order to be created blank space at the end of the page
          par(mar = c(0.1, 0.1, 0.1, 0.1))
          
          # Create an empty plot
          plot.new()
          
          # Defaults margins
          par(mar= c(5.1, 4.1, 4.1, 2.1))
          
          # Index indicating the page of the report
          page <- page+1
          
          # tree plot presenting the samples of each reference clusters and the samples of their closest test clusters
          tree_plot <- tree(unifract_dist[rownames(dist_plot_4[dist_plot_4$cluster==cluster,]),rownames(dist_plot_4[dist_plot_4$cluster==cluster,])],factor(denovo_clusters[rownames(dist_plot_4[dist_plot_4$cluster==cluster,])],levels=unique(denovo_clusters[rownames(dist_plot_4[dist_plot_4$cluster==cluster,])]))
                            ,colour_matrix[levels(factor(denovo_clusters[rownames(dist_plot_4[dist_plot_4$cluster==cluster,])],levels=unique(denovo_clusters[rownames(dist_plot_4[dist_plot_4$cluster==cluster,])]))),2])
          
          # Title of the cladogram
          text3 <- paste("Phylogram")
          
          # Create a text grob for the cladogram
          tgrob3 <- text_grob(text3, face = "bold", color = "Black")
          
          # Choose the selected type of plot
          if (plot_type == "Boxplots") {plot <-boxplot_list[[cluster]] } else if (plot_type == "Point plots") {plot <- pointplot_list[[cluster]] } else {plot <-violinplot_list[[cluster]]}
          
          
          # Arrange in the grid the above elements
          grid.arrange(arrangeGrob(plot,pvaluetable,ncol=2,nrow=1,widths = c(1,1)),tgrob3,tree_plot,nrow = 3,as.table = T,heights=c(1.5,0.2,1))
          
          # Write the PERMANOVA p-values
          write_table(paste0("DeNoAn-page ",page, " PERMANOVA pvalues.tab"),pvaluesb,"pvalues")
          
          # Print the Phylogram in png format
          png_plot(paste0("DeNoAn-page ",page,"Phylogram.png"),tree_plot)
          
          # Print theplot in png format
          png_plot(paste0("DeNoAn-page ",page," ",plot_type,".png"),plot)
          
        }
      }
      
      dev.off()
      
    } else if (index == 1 & Test_name == "None") {
      
      # Report in case there are no test groups
      
      #################################################
      ########## Cluster Analysis report ##############
      #################################################
      
      
      # Set working directory
      setwd(outputs_path)
      
      # Open the PDF
      pdf("Cluster Analysis.pdf")
      
      
      ########### Cover Page ##########
      
      # Index indicating the page of the report
      page <- 1
      
      # Set the layout of the cover page
      layout( matrix(c(1,1), ncol=1) )
      # Default margins
      par(mar= c(5.1, 4.1, 8, 2.1))
      
      # Create an empty plot
      plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      # Title of the plot
      title(main=paste0("Cluster Analysis \n ClAn-",central_point),cex.main = 1.8,font.main = 2)
      # Description of this PDF
      text(x = 0.5, y = 0.5, paste("This report contains results about the Reference group\n",
                                   "* All the elements of this report(plots and tables) will be found in the \n results folder with the prefix ClAn-"),cex = 1.3, col = "black",font = 3)
      
      # Default margins
      par(mar= c(5.1, 4.1, 4.1, 2.1))
      
      
      ############# Page 2 ############
      
      # Index indicating the page of the report
      page <- page+1
      
      # Title of the 2nd page
      text <- "Cladogram"
      # Create a text grob for the title
      tgrob1 <- text_grob(text, face = "bold.italic", color = "Black",size=16)
      
      # Description of the cladogram
      text2 <- paste0("Cladogram of the reference (",paste(reference_name,collapse="+"),") samples")
      # Create a text grob for the Description
      tgrob2 <- text_grob(text2, face = "italic", color = "Black",size=12)
      
      # Create a vector with the colours of the reference groups
      col <- brewer.pal(n = nlevels(factor(meta_file[,mapping_column],levels = unique(meta_file[,mapping_column]))) , name = "Set1")[1:nlevels(factor(meta_file[,mapping_column],levels = unique(meta_file[,mapping_column])))]
      
      # Tree of the reference end test groups
      tree_plot <- tree(unifract_dist,factor(meta_file[,mapping_column],levels = unique(meta_file[,mapping_column])), col )
      
      # Arrange in the grid the above elements
      grid.arrange(tgrob1,tgrob2,tree_plot,nrow = 3,heights=c(0.5,0.5,2.5))
      
      # Print the Phylogram in png format
      png_plot(paste0("ClAn-page ",page,"Phylogram.png"),tree_plot)
      
      ############# Page 3 ############
      
      if (reference_clusters > 1) {
        
        # Index indicating the page of the report
        page <- page+1
        
        
        # Information about the refernce clusters
        clustering_text <- paste0("The reference dataset (",paste(reference_name,collapse="+"),") was clustered into ", reference_clusters ," groups:\n")
        propability <- c()
        number_samples <- c()
        for (i in 1:reference_clusters){
          clustering_text <- paste(clustering_text,paste(reference_name,collapse="+") ,i, "has", length(which(clusters==paste(paste(reference_name,collapse="+"),i))),"samples.\n")
          if (reference_clusters > 1) {
            number_samples <- c(number_samples,length(which(clusters==paste(paste(reference_name,collapse="+"),i))))
            propability <- c(propability,1/reference_clusters)}
        }
        
        clustering_text <- paste(clustering_text,paste("The p-value of the Chi square (Goodness of fit) test is",round(chisq.test(number_samples,p=propability)$p.value,4)),"*\n")
        
        
        # Create a text grob for the clustering
        tgrob <- text_grob(clustering_text, face = "italic", color = "Black",size=13)
        
        # Create a text grob for the ttle of the page
        tgrob2 <- text_grob("De novo clustering of the reference group",face = "bold.italic", color = "Black",size=15)
        
        text <-  "* The Chi square (Goodness of fit) test checks if there is discrepancy betwwen the observed number of observations \n of each cluster and the number of observation if the samples were uniformly distributed "
        
        # Create a text grob for the clustering
        tgrob3 <- text_grob(text, face = "italic", color = "Black",size=9)
        
        # Arrange in the grid the above elements
        grid.arrange(tgrob2,tgrob,tgrob3,nrow=3,heights=c(0.5,3,0.5))
        
        
        ############# Page 4 ############
        
        # Index indicating the page of the report
        page <- page+1
        
        # Set the layout of the cover page (3 rows, 1 column)
        layout(matrix(c(1, 1, 1,
                        2, 2, 2,
                        3, 3, 3), nrow=3, byrow=TRUE),heights=c(1,3,0.3) )
        
        # Create an empty plot
        plot.new()
        # Title of the 3rd page
        title(main="MDS Plot",cex.main = 1.8,font.main = 2)
        
        # Description of the MDS plot
        text(x=0.5, y = 0.30, "MDS plot presenting the reference clusters.", cex = 1.5, col = "black")
        
        # MDS plot presenting the reference clusters
        sclass(unifract_dist[input_table[,1] %in% reference_name,input_table[,1]%in%reference_name],factor(mapping_cluster[which(input_table[,1] %in% reference_name)], levels=unique(mapping_cluster[which(input_table[,1] %in% reference_name)])),
               NULL, colour_matrix[levels(factor(mapping_cluster[which(input_table[,1] %in% reference_name)],levels = unique(mapping_cluster[which(input_table[,1] %in% reference_name)]))),2])
        
        
        # Change the margins in order to be created blank space at the end of the page
        par(mar = c(0.1, 0.1, 0.1, 0.1))
        # Create an empty plot
        plot.new()
        
        # Default margins
        par(mar= c(5.1, 4.1, 4.1, 2.1))
        
        # Print the MDS Plot in png format
        png_plot(paste0("ClAn-page ",page," MDS Plot.png"),
                 sclass(unifract_dist[input_table[,1] %in% reference_name,input_table[,1]%in%reference_name],factor(mapping_cluster[which(input_table[,1] %in% reference_name)], levels=unique(mapping_cluster[which(input_table[,1] %in% reference_name)])),
                        NULL, colour_matrix[levels(factor(mapping_cluster[which(input_table[,1] %in% reference_name)],levels = unique(mapping_cluster[which(input_table[,1] %in% reference_name)]))),2]))
        
        
        ############# Page 5 ############
        
        # Index indicating the page of the report
        page <- page+1
        
        # Tile of the 4th page
        text1 <- "Cluster Percentages"
        
        # Create a text grob for the title
        tgrob <- text_grob(text1,face = "bold.italic", color = "Black",size=15)
        
        # Arrange in the grid the above elements
        grid.arrange(tgrob,gtableGrob("Clusters",Chi_square_matrix,setDT="YES",size=12), barplot, ncol = 1, nrow = 3,heights=c(0.4,1.2,1.5))
        
      }
      
      dev.off()
      
    }
    
    #---------- Create the Description Analysis file ------------#
    
    description_analysis <- data.frame(c(input_otu,normalized,tree_or_matrix,input_tree_or_matrix,input_meta,mapping_column,reference_name,
                                         reference_clusters,paste(Test_name,collapse=","),paste(test_clusters,collapse=","),paste(exploratory_columns,collapse=",")
                                         ,central_point,plot_type))
    # Name the rows of the description_analysis
    rownames(description_analysis) <- c("OTUs/ASVs table:","Normalized:","Tree or Distances matrix:","Name of tree:","Name of the mapping file:",
                                        "Column of the mapping file:","Reference group(s):","Number of reference clusters:",
                                        "Names of the test group(s)","Number of test clusters:", "Exploratory columns:","Central points:","plot type:")
    # Name the rows of the description_analysis
    colnames(description_analysis) <- c("")
    
    
    #################################################################################
    ######################      Write Output tables   ###############################
    #################################################################################
    
    setwd(tables_path)
    
    # Write a table that contains all the Medoids' distances 
    write.table(distances_matrix, "Distances From Cluster Medoids.tab", sep = "\t",col.names =NA, row.names = TRUE,quote = FALSE)
    
    # Write a table that contains the ordered distance of all the samples
    if (Test_name != "None") {
      write.table(ordered_clusters, "Ordered cluster distances(All Samples).tab", sep = "\t",col.names =NA, row.names = TRUE,quote = FALSE)
    }
    
    # Return to results path
    setwd(outputs_path)
    
    # Write the new mapping file with the added column
    write.table(mapping_file, "mapping file.tab", sep = "\t",col.names =NA, row.names = TRUE,quote = FALSE)  
    
    # Write the Description Analysis file
    write.table(description_analysis, "Description Analysis.tab", sep = "\t",col.names =F, row.names = TRUE,quote = FALSE) 
    
    # Return to the original path
    setwd(OriginalPath)
    
    if(!flag) { stop("
        It was not possible to install all required R libraries properly.
                     Please check the installation of all required libraries manually.\n
                     Required libaries:ade4, GUniFrac, phangorn, cluster, fpc,ggplot2,gridExtra,grid,gtable,
                       stats,cowplot,graphics,vegan, dplyr,data.table,tidyr,caTools,RColorBrewer,tools")
    }
  }
}

############################################################################################################################################################
#############################################################       Script Ended         ###################################################################
############################################################################################################################################################
