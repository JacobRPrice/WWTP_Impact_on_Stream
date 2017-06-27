#----------------------------------------------------------
# Notes or References:
#
#----------------------------------------------------------
########
# Preliminary Items
########

# load utilities and functions
source("/Users/jprice/Desktop/WWTP_Impact/R/utilities.R",
	verbose=TRUE,max.deparse.length=Inf)
source("/Users/jprice/Desktop/WWTP_Impact/R/functions.R",
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
# create sample data object
########
samdf<-read.csv(file.path(data_path,"map_stream.txt"),sep="\t",header=TRUE)

str(samdf, list.len=ncol(samdf))
# Make sure that data is being formatted correctly
samdf$site<-as.factor(samdf$site)
samdf$day<-as.factor(samdf$day)
samdf$RecordLabel<-as.character(samdf$RecordLabel)
samdf$Location<-factor(samdf$Location,labels=c("UG","Ambler","Below Confluence","Sandy"))
samdf$Position<-factor(samdf$Position,labels=c("Up","Down"))
samdf$stream<-factor(samdf$stream,labels=c("Wiss","Sandy"))
samdf$Date<-as.Date(samdf$Date,origin="1899-12-30")

# transformation of environmental data
samdf$log_F<-log10(samdf$F)
samdf$log_Cl<-log10(samdf$Cl)
samdf$log_NO2N<-log10(samdf$NO2N)
samdf$log_Br<-log10(samdf$Br)
samdf$log_NO3N<-log10(samdf$NO3N)
samdf$log_SO4<-log10(samdf$SO4)
samdf$log_PO4<-log10(samdf$PO4)
samdf$log_Ca<-log10(samdf$Ca)
samdf$log_Mg<-log10(samdf$Mg)
samdf$log_Na<-log10(samdf$Na)
samdf$log_K<-log10(samdf$K)
samdf$log_Fe<-log10(samdf$Fe)
samdf$log_TotDisP<-log10(samdf$TotDisP)
samdf$log_Si<-log10(samdf$Si)
samdf$log_Cu<-log10(samdf$Cu)
samdf$log_Mn<-log10(samdf$Mn)
samdf$log_Sr<-log10(samdf$Sr)

# Create SampleID var in samdf
samdf$SampleID<-rownames(seqtab)

# is the metadata in the correct order? 
cbind(
	rownames(seqtab),
	samdf$RecordLabel,
	samdf$SampleID,
	samdf$site,
	samdf$day,
	samdf$seqrun,
	c(
	"1.USUG.1","1.USUG.2",
	"2.DSUG.1","2.DSUG.2",
	"3.USAmb.1","3.USAmb.2",
	"4.DSAmb.1","4.DSAmb.2",
	"5.BC.1","5.BC.2",
	"6.SR.1","6.SR.2"
	)
)
rownames(samdf)<-rownames(seqtab)

saveRDS(samdf,file.path(temprds_path,"samdf.RDS"))

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