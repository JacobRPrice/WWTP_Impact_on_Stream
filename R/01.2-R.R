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
tempfig_path<-file.path(figs_path,"01")
if(!file_test("-d",tempfig_path)) dir.create(tempfig_path)
temprds_path<-file.path(rds_path,"01")
if(!file_test("-d",temprds_path)) dir.create(temprds_path)

# prep environment and get data
.cran_packages<-c("ggplot2", "gridExtra","magrittr")
.bioc_packages<-c("dada2", "phyloseq", "DECIPHER", "phangorn")
sapply(c(.cran_packages,.bioc_packages),require,character.only=TRUE)

# set seed
set.seed(20141004)

#----------------------------------------------------------
########
# import derep
########
derepRs <- readRDS(file.path(temprds_path,"derepRs.RDS"))

#----------------------------------------------------------
########
# use dada2 to model substitution errors and distinguish sequencing errors from real biological variation. 
########
# https://www.rdocumentation.org/packages/dada2/versions/1.0.3/topics/dada
ddR <- dada(derepRs, err=NULL, selfConsist=TRUE, multithread=TRUE)
dada2:::checkConvergence(ddR[[1]])
saveRDS(ddR,file.path(temprds_path,"ddR.RDS"))
dadaRs <- dada(derepRs, err=ddR[[1]]$err_out, pool=TRUE, multithread=TRUE)
for (i in seq(1,length(dadaRs))) {
	print(names(dadaRs)[[i]])
	print(dadaRs[[i]])
}
saveRDS(dadaRs,file.path(temprds_path,"dadaRs.RDS"))

plot.err.r<-plotErrors(dadaRs,nominalQ=TRUE)
ggsave(file.path(tempfig_path,"dada.err.r.pdf"),plot.err.r,width=25,height=20,limitsize=FALSE)

#----------------------------------------------------------
########
# clear temporary paths
########
rm(tempfig_path)
rm(temprds_path)
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
########
# 
########