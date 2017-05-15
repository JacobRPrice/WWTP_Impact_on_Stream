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
derepFs <- readRDS(file.path(temprds_path,"derepFs.RDS"))

#----------------------------------------------------------
########
# use dada2 to model substitution errors and distinguish sequencing errors from real biological variation. 
########
# https://www.rdocumentation.org/packages/dada2/versions/1.0.3/topics/dada
ddF <- dada(derepFs, err=NULL, selfConsist=TRUE, multithread=TRUE)
dada2:::checkConvergence(ddF[[1]])
saveRDS(ddF,file.path(temprds_path,"ddF.RDS"))
dadaFs <- dada(derepFs, err=ddF[[1]]$err_out, pool=TRUE, multithread=TRUE)
for (i in seq(1,length(dadaFs))) {
	print(names(dadaFs)[[i]])
	print(dadaFs[[i]])
}
saveRDS(dadaFs,file.path(temprds_path,"dadaFs.RDS"))

plot.err.f<-plotErrors(dadaFs,nominalQ=TRUE)
ggsave(file.path(tempfig_path,"dada.err.f.pdf"),plot.err.f,width=25,height=20,limitsize=FALSE)

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