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
# import files
########
derepFs <- readRDS(file.path(temprds_path,"derepFs.RDS"))
dadaFs <- readRDS(file.path(temprds_path,"dadaFs.RDS"))
derepRs <- readRDS(file.path(temprds_path,"derepRs.RDS"))
dadaRs <- readRDS(file.path(temprds_path,"dadaRs.RDS"))

# verify that the error rates have been reasonably well-estimated
# inspect the fit between observed error rates and the fitted error rates
plot.err.f<-plotErrors(dadaFs,nominalQ=TRUE)
plot.err.r<-plotErrors(dadaRs,nominalQ=TRUE)
plot.err.model<-grid.arrange(plot.err.f,plot.err.r,nrow=1,ncol=2)
ggsave(file.path(tempfig_path,"dada.err.FULL.pdf"), plot.err.model,width=15,height=7,limitsize=FALSE)

#----------------------------------------------------------
########
# merge the inferred forward and reverse sequences, while removing paired sequnces that do not perfectly overlap as a final control against residual error
########
# https://www.rdocumentation.org/packages/dada2/versions/1.0.3/topics/mergePairs
mergers <- mergePairs(
	dadaFs, 
	derepFs, 
	dadaRs, 
	derepRs,
	minOverlap=20,
	verbose=TRUE
)
head(mergers[[1]])
saveRDS(mergers,file.path(temprds_path,"mergers.RDS"))

#----------------------------------------------------------
########
# Construct the sequence table and remove chimeras
########
# https://www.rdocumentation.org/packages/dada2/versions/1.0.3/topics/makeSequenceTable
names(mergers)
table(mergers[[1]]$nmatch)
seqtab.all <- makeSequenceTable(mergers) 
saveRDS(seqtab.all,file.path(temprds_path,"seqtab.all.RDS"))

########
# remove chimeras
########
# https://www.rdocumentation.org/packages/dada2/versions/1.0.3/topics/removeBimeraDenovo
seqtab <- removeBimeraDenovo(seqtab.all,verbose=TRUE)
saveRDS(seqtab,file.path(temprds_path,"seqtab.RDS"))

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