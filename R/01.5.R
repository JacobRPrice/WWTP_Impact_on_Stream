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

taxtab<-readRDS(file.path(temprds_path,"taxtab.RDS"))
samdf<-readRDS(file.path(temprds_path,"samdf.RDS"))
seqtab<-readRDS(file.path(temprds_path,"seqtab.RDS"))
fitGTR<-readRDS(file.path(temprds_path,"fitGTR.RDS"))

# set seed
set.seed(20141004)

#----------------------------------------------------------
########
# combine all produced objects into one PhyloSeq object
########

ps_orig <- phyloseq(
	tax_table(taxtab),
	sample_data(samdf),
	otu_table(seqtab, taxa_are_rows = FALSE),
	phy_tree(fitGTR$tree)
)
ps_orig
sample_names(ps_orig)
# save phyloseq object in temprds_path as a backup
saveRDS(ps_orig,file.path(temprds_path,"ps_orig_BAK.RDS"))
# save phyloseq object in data_path for use going forward
saveRDS(ps_orig,file.path(data_path,"ps_orig.RDS"))

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