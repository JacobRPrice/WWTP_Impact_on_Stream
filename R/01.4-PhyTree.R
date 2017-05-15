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
# Construct a phylogenetic tree
# make tree using updated protocol (v2)
########

# https://www.rdocumentation.org/packages/dada2/versions/1.0.3/topics/getSequences
seqs <- getSequences(seqtab)
names(seqs) <- seqs 
alignment<-AlignSeqs(DNAStringSet(seqs),anchor=NA)
saveRDS(alignment,file.path(temprds_path,"alignment.RDS"))
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
saveRDS(phang.align,file.path(temprds_path,"phang.align.RDS"))
dm <- dist.ml(phang.align)
saveRDS(dm,file.path(temprds_path,"dm.RDS"))
treeNJ <- NJ(dm) 
saveRDS(treeNJ,file.path(temprds_path,"treeNJ.RDS"))
fit <- pml(treeNJ, data=phang.align) 
saveRDS(fit,file.path(temprds_path,"fit.RDS"))
fitGTR <- update(fit, k=4, inv=0.2)
saveRDS(fitGTR,file.path(temprds_path,"fitGTR.RDS"))
fitGTR <- optim.pml(
	fitGTR, 
	model="GTR", 
	optInv=TRUE, 
	optGamma=TRUE,
	rearrangement = "stochastic", 
	control = pml.control(trace = 0)
)
saveRDS(fitGTR,file.path(temprds_path,"fitGTR.RDS"))

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