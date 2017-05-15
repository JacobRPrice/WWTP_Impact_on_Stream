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

# prep environment and get data
pkgs.to.load<-c("phyloseq","ggplot2","DESeq2","vegan")
analysis_prep(pkgs.to.load)

# set seed
set.seed(20141004)

#----------------------------------------------------------
########
# prep data 
########
########
# apply VST to tax counts
########

psnorm<-ps

# prep dds object
dds<-phyloseq_to_deseq2(ps,~1)
dds
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds)
vst<-getVarianceStabilizedData(dds)
dim(vst)
dim(otu_table(ps))
vst[vst < 0.0] <- 0.0
otu_table(ps)<-otu_table(vst, taxa_are_rows=TRUE)

# pull out just the indicator species
targOrd<-c("Bacteroidales","Bifidobacteriales","Clostridiales")
ps
ps.ind<-subset_taxa(ps, Order %in% targOrd)
ps.ind

########
# prep/subset sample data
########
sd<-vegansd(ps)
dim(sd)
names(sd)
sd<-sd[,-c(1:61,68,79)]
dim(sd)

#----------------------------------------------------------
########
# full community
########
resf<-bioenv(
	comm=vegdist(t(vst),method="bray"),
	env=sd,
	parallel=2,
	upto=11,
	metric="euclidean"
	)
resf
# Best model has 2 parameters (max. 11 allowed):
# log_Cl log_SO4
# with correlation  0.6002505 
summary(resf)
bioenvdist(resf,which="best")

#----------------------------------------------------------
########
# indicator community
########
ind<-veganotu(ps.ind)

resi<-bioenv(
	comm=vegdist(ind,method="euclidean"),
	env=sd,
	parallel=2,
	upto=11,
	metric="euclidean"
	)
resi
# Best model has 4 parameters (max. 11 allowed):
# log_Cl log_NO2N log_Mg log_Si
# with correlation  0.6804091 
summary(resi)
bioenvdist(resi,which="best")

#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
########
# 
########