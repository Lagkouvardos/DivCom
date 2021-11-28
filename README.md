![DivCom-Logo-final](https://user-images.githubusercontent.com/8244618/139091582-43c02470-9e6f-4711-bab7-1e04ad59c300.png)



## Introduction
When analyzing microbiome data, one of the main objectives is to effectively compare the microbial profiles of samples belonging to different groups. The visual representation of the dataset, usually in the form of a PCoA plot, in conjunction with the results obtained from a suitable statistical test, typically provides us with all the necessary information to draw conclusions about the overall similarity or dissimilarity of the observed microbial communities. However, this approach can lead us to misleading assumptions. To address those issues, we developed Divide and Compare (DivCom), an automated R-based tool for advanced beta diversity analysis. The main advantage of DivCom is that it does not rely only on the visual representation of the data or the results of the statistical tests. Instead, it takes into consideration the substructure of the data and uses beta diversity measures as distances metrics. Firstly, DivCom attempts to reveal the inner structure of the control group by dividing their samples into the appropriate number of subgroups. Then the program selects the most representative points of each cluster and calculates the distances of the profiles belonging to the test groups to those reference points. Finally, for each sample, the closest representative point is reported and stored. In this way, all the individual points are compared only with samples that are relative to them. As an extra layer to the DivCom analysis, the program also applies de novo clustering to the test groups.   

## Description

### Organization
DivCom consists of two scripts:

1.	Optimal Number of Clusters
2.	Distances

The first script that the user should run is named 'Optimal Number of Clusters'. Its purpose is to provide all the necessary information about the number of clusters of each group. This information will be used in the 'Distances' script, which is the main script of the program that applies the proposed methodology. Also, there is an extra folder that contains the template data of the DivCom. These data can be used for training and demonstration purposes. 

The pre-processing of the raw sequencing data through the [IMNGS](https://www.imngs.org/ "IMNGS download site") platform can provide the user with the necessary input files for the DivCom. Also, the  outputs derived from the [Rhea](https://github.com/Lagkouvardos/Rhea/ "Rhea download site") pipeline have the appropriate form to be used as inputs to the DivCom.

A detailed README file is provided for each of the scripts. In these files, the procedures followed by the program are described, and all the requirements are specified.


### Script structure
In order for the program to be as user-friendly as possible, the two scripts of DivCom follow the same structure. Each script contains the Commentary, the Initialization, and the Main section. 
The general rule that applies is that the lines concerning the Commentary and Initialization sections start with a hash sign followed by a back-tick (#`). These lines aim to inform the user about the details of the program and the actions that should be followed. The comments in the Main section start with the hash sign (#), and it is advised not to be changed by the user.


#### Commentary section
This section is always positioned at the beginning of each script and holds comment lines explaining the script's tasks and concept, the required input files, and the expected outputs. These lines are only for informational purposes, and they do not add any functionality to the script. Therefore, it is recommended that the user read this brief introduction and then proceed to the next stage of the program. An example of those comments can be seen below:


 
       #'This script was last modified on 21/11/2021
       #'Script Task: 
       #'
       #' The script performs de novo clustering of the groups and then calculate the distances from 
       #' the most representative points of each cluster. 
       #' Then, conducts statistical analysis and produces plots and tables. 
       #'
       #' Input: Please enter the following parameters
       #' 1. Set the path to the directory where the file is stored 
       #' 2. Write the name of the OTU table of interest in quotes
       #' 
       #'
       #' Output: 
       #' The script generates two reports in pdf format, 3 folder where the results are printed and a mapping file
       #' 1. A Distances Based Analysis report
       #' 2. A De Novo Analysis report
       #' 3. A folder with the p-values tables for the different tests 
       #' 
       #'
       #' Concept:
       #' A common problem in data analysis is to efficiently compare a group of control samples with a set of test samples.
       #' The relation existing between different groups can be highly affected by their substructure.
       #' Instead of comparing the groups as entireties, the program performs de novo clustering(with PAM algorithm) 
       #' to both the control and test groups.
       #' Then, finds the most representative points of the reference  dataset and calculates the distances(Generalized Unifrac) 
       #' of the test groups from these points.
       #' To determine  the level of similarity between the groups, the script performs statistical analysis 
       #' and produces various plots.


#### Initialization section
All the required and optional parameters necessary for executing the scripts are described in the Initiation phase. The user must fill in the required parameters; otherwise, it will not be possible for the program to be executed. The optional parameters provide an additional user-specified parameterization to the outputs. A detailed description is provided for each of the parameters. The first and most important setting that needs to be changed is the path where the input files are located. At the end of the Initialization phase, the program is ready to be executed; the only thing left for the user is to select all the code (Ctrl+A) and then run the script. Next, an example of the format of the Initiation section is presented:

	   ##################################################################################
	   ######             Set parameters in this section manually                  ######
	   ##################################################################################
	
       #' Please set the directory of the script as the working folder 
       #' Note: the path is denoted by forward slash "/"
       setwd("C:/...../..../Distances") 


       #' Please give the name of the OTUs or ASVs table
       input_otu = "OTUs-Table.tab" 


       #' Please insert "YES" if the OTUs/ASVs table is normalized otherwise, insert "NO"
       #' In case the table is not normalized, the OTUs/ASVs table will be normalized based on 
       #' the minimum sum of reads of the samples.
       normalized = "NO"
	
	   ###################    NO CHANGES ARE NEEDED BELOW THIS LINE    ####################
	

#### Main Section
The main section contains the actual code that applies the methodology of the scripts. The functions take as inputs the values that were declared during the Initialization phase. It is strongly recommended the lines of the main section not be changed as this can interfere with the rest of the code and affect the results and or cause fatal errors. An example of such a section can be seen below: 

	   ##################################################################################
	   ######                             Main Script                              ###### 
	   ##################################################################################
	
       # Load the mapping file containing individual sample information (sample names in the first column)
       meta_file <- data.frame(read.table (file = input_meta, check.names = FALSE, header = TRUE,
                               dec = ".", sep = "\t", row.names = 1, comment.char = ""))

	   if (ncol(meta_file)==0 | nrow(meta_file)==0){
	     # Load the mapping file containing individual sample information (sample names in the first column)
	     meta_file <- data.frame(read.table (file = input_meta, check.names = FALSE, header = TRUE,
	                                         dec = ".", sep = ",", row.names = 1, comment.char = ""))
	   }

	   # Order the rows of the meta file
	   if (ncol(meta_file)==1){
	     meta_file <- data.frame(meta_file[order(meta_file[,mapping_column ]),, drop = FALSE])
	     colnames(meta_file) <- mapping_column 
	     meta_file <- data.frame(meta_file[meta_file[,mapping_column] %in% c(reference_name ,Test_name),,drop=FALSE])
	   } else {
	     meta_file <- data.frame(meta_file[order(meta_file[,mapping_column]),])
	     meta_file <- data.frame(meta_file[meta_file[,mapping_column] %in% c(reference_name ,Test_name),])
	   }


### Requirements
In order to execute the scripts, it is required the [R language](https://www.r-project.org/ "R download site") to be locally installed. The use of the [R-Studio](https://www.rstudio.com/products/rstudio-desktop/ "R-studio download site") will simplify the procedure even for non-experienced users. During the first execution of the program, all the necessary packages will be installed, so at least for this first run, a stable internet connection is required. 


### Installation
DivCom does not require any installation process. 

In order to use DivCom, the following steps should be followed:
*	First, the user has to download the DivCom project from the Github repository.
*	Then decompress the files into the desired destination. 
*	For a smoother experience and reduced waiting times, we recommend running the install_packages.R script first. This script will automatically install all the packages required by DivCom scripts. Packages that, for any reason, cannot be installed automatically will be identified and noted in the corresponding report. 

Keeping one copy of the DivCom files in each study will help to avoid mixing up different analyses. We also recommend looking for new and improved versions of DivCom from time to time.

