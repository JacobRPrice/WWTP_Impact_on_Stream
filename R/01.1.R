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

# set paths for raw data
# raw read directory
rawread_path<-file.path(data_path,"rawreads/")
# filtered read directory
filt_path<-file.path(data_path,"filtreads/")
if(!file_test("-d",filt_path)) dir.create(filt_path)

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
# import sequences
########
fns <- sort(list.files(rawread_path, full.names = TRUE))
fnFs <- fns[grepl("R1", fns)]
fnRs <- fns[grepl("R2", fns)]
head(fns)
head(fnFs)
head(fnRs)

#----------------------------------------------------------
########
# trim and filter 
########

########
# visualize quality of sequences
########
# https://www.rdocumentation.org/packages/dada2/versions/1.0.3/topics/plotQualityProfile

# plots of forward reads
plot.qual.f<-NULL
plot.qual.f<-list()
for(i in 1:12) { 
	plot.qual.f[[i]]<-plotQualityProfile(fnFs[i]) + 
	ggtitle("Fwd") + 
	geom_hline(yintercept=30, size=0.5) + 
	geom_hline(yintercept=25, size=0.5) + 
	geom_hline(yintercept=20, size=0.5) +
	geom_vline(xintercept=10, size=0.5, color="red") + 
	geom_vline(xintercept=275, size=0.5, color="green") + 	
	geom_vline(xintercept=70, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=80, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=90, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=100, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=110, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=120, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=130, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=140, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=150, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=160, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=170, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=180, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=190, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=200, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=210, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=220, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=230, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=240, size=0.5, color="gray",linetype="dashed") +
	geom_vline(xintercept=250, size=0.5, color="gray",linetype="dashed") +
	geom_vline(xintercept=260, size=0.5, color="gray",linetype="dashed") +
	geom_vline(xintercept=270, size=0.5, color="gray",linetype="dashed") +
	geom_vline(xintercept=280, size=0.5, color="gray",linetype="dashed") +
	geom_vline(xintercept=290, size=0.5, color="gray",linetype="dashed") +
	xlim(0,300) + ylim(0,40)
}
plot.read.qual.f<-grid.arrange(
	plot.qual.f[[1]],plot.qual.f[[2]],plot.qual.f[[3]],plot.qual.f[[4]],
	plot.qual.f[[5]],plot.qual.f[[6]],plot.qual.f[[7]],plot.qual.f[[8]],
	plot.qual.f[[9]],
	plot.qual.f[[10]],plot.qual.f[[11]],plot.qual.f[[12]],
	nrow=6,ncol=2)
ggsave(file.path(tempfig_path,"read.qual.f.pdf"),plot.read.qual.f,width=2*3,height=6*3,limitsize=FALSE)

# plots of reverse reads
plot.qual.r<-NULL
plot.qual.r<-list()
for(i in 1:12) { 
	plot.qual.r[[i]]<-plotQualityProfile(fnRs[i]) + 
	ggtitle("Rev") + 
	geom_hline(yintercept=30, size=0.5) + 
	geom_hline(yintercept=25, size=0.5) + 
	geom_hline(yintercept=20, size=0.5) +
	geom_vline(xintercept=10, size=0.5, color="red") + 
	geom_vline(xintercept=175, size=0.5, color="green") + 	
	geom_vline(xintercept=70, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=80, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=90, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=100, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=110, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=120, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=130, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=140, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=150, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=160, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=170, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=180, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=190, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=200, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=210, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=220, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=230, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=240, size=0.5, color="gray",linetype="dashed") +
	geom_vline(xintercept=250, size=0.5, color="gray",linetype="dashed") +
	geom_vline(xintercept=260, size=0.5, color="gray",linetype="dashed") +
	geom_vline(xintercept=270, size=0.5, color="gray",linetype="dashed") +
	geom_vline(xintercept=280, size=0.5, color="gray",linetype="dashed") +
	geom_vline(xintercept=290, size=0.5, color="gray",linetype="dashed") +
	xlim(0,300) + ylim(0,40)
}
plot.read.qual.r<-grid.arrange(
	plot.qual.r[[1]],plot.qual.r[[2]],plot.qual.r[[3]],plot.qual.r[[4]],
	plot.qual.r[[5]],plot.qual.r[[6]],plot.qual.r[[7]],plot.qual.r[[8]],
	plot.qual.r[[9]],
	plot.qual.r[[10]],plot.qual.r[[11]],plot.qual.r[[12]],
	nrow=6,ncol=2)
ggsave(file.path(tempfig_path,"read.qual.r.pdf"),plot.read.qual.r,width=2*3,height=6*3,limitsize=FALSE)

########
# carry out trimming and filtering of reads
########
# make list of names with paths for filtered read files
filtFs<-file.path(filt_path,basename(fnFs))
filtRs<-file.path(filt_path,basename(fnRs))
cbind(fnFs,filtFs)
cbind(fnRs,filtRs)

# # carry out trimming and filtering of reads
for(i in seq_along(fnFs)) {
  fastqPairedFilter(
  	c(fnFs[[i]], fnRs[[i]]),
  	c(filtFs[[i]], filtRs[[i]]),
  	rm.phix=TRUE,
  	trimLeft=10, 
  	truncLen=c(275,175),
  	maxN=0, # none allowed for dada2
  	maxEE=2, 
  	truncQ=2,
  	compress=TRUE,
  	verbose=TRUE)
}

#----------------------------------------------------------
########
# infer sequence variants
########

########
# dereplicate forward and reverse reads
########
# https://www.rdocumentation.org/packages/dada2/versions/1.0.3/topics/derepFastq
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
sam.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
cbind(basename(fnFs),basename(fnRs))
cbind(basename(filtFs),basename(filtRs))
cbind(
	basename(filtFs),
	sam.names,
	c(
	"1.USUG.1","1.USUG.2",
	"2.DSUG.1","2.DSUG.2",
	"3.USAmb.1","3.USAmb.2",
	"4.DSAmb.1","4.DSAmb.2",
	"5.BC.1","5.BC.2",
	"6.SR.1","6.SR.2"
	)
)
sam.names<-	c(
	"1.USUG.1","1.USUG.2",
	"2.DSUG.1","2.DSUG.2",
	"3.USAmb.1","3.USAmb.2",
	"4.DSAmb.1","4.DSAmb.2",
	"5.BC.1","5.BC.2",
	"6.SR.1","6.SR.2"
	)
names(derepFs) <- sam.names
names(derepRs) <- sam.names
derepFs[[1]]
derepRs[[1]]
saveRDS(derepFs,file.path(temprds_path,"derepFs.RDS"))
saveRDS(derepRs,file.path(temprds_path,"derepRs.RDS"))

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