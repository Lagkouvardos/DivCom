![DivCom-Logo-final](https://user-images.githubusercontent.com/8244618/139091582-43c02470-9e6f-4711-bab7-1e04ad59c300.png)



## Introduction
When analyzing microbiome data, one of the main objectives is to effectively compare the microbial profiles of samples belonging to different groups. Divide and Compare (DivCom) is an automated tool for advanced beta diversity analysis. DivCom reveals the inner structure of the groups by dividing their samples into the appropriate number of clusters and then compares the distances of every profile to the centers of these clusters. The advantage of this approach is that it takes into consideration the substructure of the dataset and uses beta diversity measures as distances metrics. This information can be used for determining the existing interrelation of the groups.

## Description


### Organization
DivCom consists of two scripts:

1.	Beta-Diversity
2.	DivCom

The first script is named 'Beta-Diversity' and its purpose is to provide all the necessary information about the number of clusters of each group. This information will be used in the 'DivCom' script, which is the main script of the program that applies the proposed methodology. 

A detailed README file is provided for each of the scripts. In these files, the procedures followed by the program are described, and all the requirements are specified.

Also, there is an extra folder that contains the template data of the DivCom. These data can be used for training and demonstration purposes. 

The Workflow of DivCom is presented on the following Figure:

<img src="https://user-images.githubusercontent.com/68451208/163047398-72f8bab1-62b2-4ebd-bfcb-569ef848cc32.jpg" width="250" height="450" align="center">


### Requirements and Installation
In order to execute the scripts, it is required the [R language](https://www.r-project.org/ "R download site") to be locally installed. The use of the [R-Studio](https://www.rstudio.com/products/rstudio-desktop/ "R-studio download site") will simplify the procedure even for non-experienced users. During the first execution of the program, all the necessary packages will be installed, so at least for this first run, a stable internet connection is required.

DivCom does not require any installation process.
In order to use DivCom, the following steps should be followed:

*	First, the user has to download the DivCom project from the Github repository.
*	Then decompress the files into the desired destination. 
*	For a smoother experience and reduced waiting times, we recommend running the install_packages.R script first. This script will automatically install all the packages required by DivCom scripts. Packages that, for any reason, cannot be installed automatically will be identified and noted in the corresponding report.

Note for Linux users: In case you experiencing problems installing the 'ggpubr' package first run the following line and then proceed to the scripts. 

	install.packages("https://cran.r-project.org/src/contrib/Archive/nloptr/nloptr_1.2.1.tar.gz", repos=NULL, type="source")

### Input Files
The user has to provide three mandatory files

* An OTUs or ASVs abundance table.
* A phylogenetic tree that corresponds to the OTUs or ASVs of the abundance table OR a dissimilarity matrix.
* A mapping file that contains the labels of the samples.

The pre-processing of the raw sequencing data through the [IMNGS platform](www.imngs.org/ "IMNGS site") can provide the user with the necessary input files for the DivCom. Also, the  outputs derived from the [Rhea pipeline](https://github.com/Lagkouvardos/Rhea/ "Rhea repository") have the appropriate form to be used as inputs to the DivCom.


### Script structure
In order for the program to be as user-friendly as possible, the two scripts of DivCom follow the same structure.

Each script contains: 
* The Commentary section
* The Initialization section
* The Main section 

The general rule that applies is that the lines concerning the Commentary and Initialization sections start with a hash sign followed by a back-tick (#`). These lines aim to inform the user about the details of the program and the actions that should be followed. The comments in the Main section start with the hash sign (#), and it is advised these lines not to be changed by the user.

Please follow the instractions and guidlines as they presented in the scripts.


Keeping one copy of the DivCom files in each study will help to avoid mixing up different analyses. We also recommend looking for new and improved versions of DivCom from time to time.


### Citation
If you use DivCom in your work, please cite the original publication:

 Intze E and Lagkouvardos I (2022) DivCom: A Tool for Systematic Partition of Groups of Microbial Profiles Into Intrinsic Subclusters and Distance-Based Subgroup Comparisons. Front. Bioinform. 2:864382. doi: 10.3389/fbinf.2022.864382
