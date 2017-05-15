# WWTP_Impact_on_Stream

# Introduction
This repository contains all of the scripts needed to replicate the results and output presented in:   

INSERT CITATION INFORMATION HERE

Please read this markdown file, including the NOTES section, in it's entirety before proceeding. 

All scripts have been written for and tested on Mac machines; some modification may be required for use with PC or Linux/Unix systems.

___  

# Replicate bioinformatics pipeline  
Follow the instructions in this section to carry out all of the bioinformatics processes up to and including the creation of the phyloseq object.

Several items are needed before starting: 
1) Raw reads (see NOTES)
2) SILVA files (see NOTES)
3) mapping file (in data directory)
All placed in thier correct locations (refer to scripts)

The file "BioinformaticsPathway.bash" is a bash script that will sequentially execute the required scripts in the "/R/" directory. This step is computationally heavy and will take approximately 40 hours to complete. It is suggested that users open up a second terminal window and excecute the "caffeinate" command to prevent the computer from going to sleep. 

___  

# Perform analysis peice-meal
Follow the instructions in this section to carry out analysis of the phyloseq object from the bioinformatics pipeline above. A copy of the phyloseq object is included in the data directory for users who do not want to build the phyloseq object from scratch.

The scripts for the results section generally follow this order:   
* qPCR.R  
* AlphaDiversity.R    
* RelAbund.R    
* Ordination_exploratory.R  
* DESeq2.R    
* DESeq2_Ind.R    
* DESeq2_Compare.R   
* DESeq2_2vs3.R    
* Ordination_db-RDA.R   
* BIOENV.R   

___  

#NOTES:  
**Change directory paths before starting**   
Prior to carrying out any replication of results, users should change the file paths found in the scripts, i.e. change "PATH/TO/DIR/" so that R is able to locate data and scripts as they are sourced.   

**Raw sequencing files**  
Raw sequencing files are currently hosted on github for reviewers/users convenience. Github may request the removal of this data. Raw sequence reads have been uploaded to the NCBI Sequence Read Archive:   
* SRA ID: SRP103534   
* Bioproject ID: PRJNA382371
Raw sequencing files should be placed within the "./data/rawreads/" directory

**SILVA**    
Users must obtain a SILVA liscense and copies of the SILVA training and species assignment datasets. These datasets should be placed in a directory called "assignTaxonomy" within the data/ directory.    
* The specific name/version of training set used for this paper: silva_nr_v123_train_set.fa.gz    
* The specific name/version of the species assignment set used for this paper: silva_species_assignment_v123.fa.gz   
* Refer to dada2's taxonomic training page: http://benjjneb.github.io/dada2/training.html