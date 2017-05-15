#----------------------------------------------------------
# Notes or References:
#
#----------------------------------------------------------
########
# Preliminary Items
########

# load utilities and functions
source("/PATH/TO/DIR/WWTP_Impact_on_Stream/R/utilities.R",
	verbose=TRUE,max.deparse.length=Inf)
source("/PATH/TO/DIR/WWTP_Impact_on_Stream/R/functions.R",
	verbose=TRUE,max.deparse.length=Inf)

# set temporary directories, if needed
temprds_path<-file.path(rds_path,"01")
if(!file_test("-d",temprds_path)) dir.create(temprds_path)

# prep environment and get data
.cran_packages<-c("ggplot2", "gridExtra","magrittr")
.bioc_packages<-c("dada2", "phyloseq", "DECIPHER", "phangorn")
sapply(c(.cran_packages,.bioc_packages),require,character.only=TRUE)

seqtab<-readRDS(file.path(temprds_path,"seqtab.RDS"))

# set seed
set.seed(20141004)

#----------------------------------------------------------
########
# assign taxonomy and create taxonomy table
########
# https://www.rdocumentation.org/packages/dada2/versions/1.0.3/topics/assignTaxonomy
# http://benjjneb.github.io/dada2/assign.html
# http://benjjneb.github.io/dada2/species.html
# http://benjjneb.github.io/dada2/training.html

########
# assign taxonomy - only goes to genus
taxtab.1to6<-assignTaxonomy(
	seqtab,
	refFasta=file.path(data_path,"assignTaxonomy/silva_nr_v123_train_set.fa.gz"),
	minBoot=80, # default=50
	verbose=TRUE
)
head(unname(taxtab.1to6))
saveRDS(taxtab.1to6,
	file.path(temprds_path,"taxtab.1to6.RDS"))

########
# merge taxonomy and genus-species classification information
taxtab<-addSpecies(
	taxtab.1to6,
	refFasta=file.path(data_path,"assignTaxonomy/silva_species_assignment_v123.fa.gz"),
	verbose=TRUE,
	allowMultiple=TRUE
)
head(unname(taxtab))
colnames(taxtab)
head(rownames(taxtab))
table(unname(taxtab)[,1],useNA="always")
table(unname(taxtab)[,2],useNA="always")
table(unname(taxtab)[,3],useNA="always")
table(unname(taxtab)[,4],useNA="always")
table(unname(taxtab)[,5],useNA="always")
table(unname(taxtab)[,6],useNA="always")
table(unname(taxtab)[,7],useNA="always")
saveRDS(taxtab,
	file.path(temprds_path,"taxtab.RDS"))

#----------------------------------------------------------
########
# clear temporary paths
########
rm(temprds_path)
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
########
# 
########