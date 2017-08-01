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
pkgs.to.load<-c("phyloseq","ggplot2","gridExtra","vegan","MASS")
analysis_prep(pkgs.to.load)

# set seed
set.seed(20141004)

#----------------------------------------------------------
########
# prepare data
########

########
# qPCR Human-specific bacteroides (B. dorei)
qPCR<-read.csv(file.path(data_path,"qPCR.txt"),sep="\t",header=TRUE)
head(qPCR)
human<-qPCR[,-4]
head(human)
human<-summarySE(human, measurevar="human", groupvars=c("SampleID","site"),na.rm=TRUE)
head(human)

########
# qPCR bacteroides spp.
bacter<-qPCR[,-3]
head(bacter)
bacter<-summarySE(bacter, measurevar="bacter", groupvars=c("SampleID","site"),na.rm=TRUE)
head(bacter)

#----------------------------------------------------------
########
# create and format individual plots
########

########
# qPCR Human-specific bacteroides (B. dorei)
qPCR.human<-
ggplot(human, aes(x=SampleID, y=human, fill=as.factor(site))) +
	geom_bar(position=position_dodge(), stat="identity") +
	geom_errorbar(aes(ymin=human-se, ymax=human+se),
		width=0.2,position=position_dodge(0.9)) +
	theme_bw() +
	ylab("B. dorei (Copy Number/100 mL)") +
	guides(fill=FALSE) +
	theme(axis.text.x=element_text(angle=270,vjust=0.5,hjust=0)) +
	theme(axis.text.y=element_text(angle=90,hjust=0.5)) +
	theme(axis.title.x=element_blank())+ 
	theme(axis.title.y=element_text(size=rel(0.9))) + 
	theme(axis.text.y=element_text(size=rel(0.8))) +
	theme(axis.text.y=element_text(size=rel(0.9))) +
	scale_y_log10()

########
# qPCR bacteroides spp.
qPCR.bacter<-
ggplot(bacter, aes(x=SampleID, y=bacter, fill=as.factor(site))) +
	geom_bar(position=position_dodge(), stat="identity") +
	geom_errorbar(aes(ymin=bacter-se, ymax=bacter+se),
		width=0.2,position=position_dodge(0.9)) +
	theme_bw() +
	ylab("Bacteroides spp. (Copy Number/100 mL)") +
	guides(fill=FALSE) +
	theme(axis.text.x=element_text(angle=270,vjust=0.5,hjust=0)) +
	theme(axis.text.y=element_text(angle=90,hjust=0.5)) +
	theme(axis.title.x=element_blank()) + 
	theme(axis.title.y=element_text(size=rel(0.9))) + 
	theme(axis.text.y=element_text(size=rel(0.8))) +
	theme(axis.text.y=element_text(size=rel(0.9))) +
	scale_y_log10()

#----------------------------------------------------------
########
# make and save plots
########

########
# qPCR - human and bacter 
ggsave(file.path(figs_path,"qPCR-HS&Bacter.eps"),
	grid.arrange(nrow=1,
	qPCR.human + 
	ggtitle("A") +
	theme(plot.title=element_text(hjust=0.5))
	,
	qPCR.bacter + 
	ggtitle("B") +
	theme(plot.title=element_text(hjust=0.5))
	),
	width=190, height=190/2+5,units="mm")

#----------------------------------------------------------
########
# paired t-test for up/down 
########
qPCR<-read.csv(file.path(data_path,"qPCR.txt"),sep="\t",header=TRUE)
human<-qPCR[,-4]
human.up<-human[c(1:6,13:18,37:42,49:54),]
human.down<-human[-c(1:6,13:18,25:36,37:42,49:54,61:72),]
var.test(human.up$human,human.down$human)
t.test(human.up$human,human.down$human,var.equal=FALSE,paired=FALSE)

qPCR<-read.csv(file.path(data_path,"qPCR.txt"),sep="\t",header=TRUE)
bacter<-qPCR[,-3]
bacter.up<-bacter[c(1:6,13:18,37:42,49:54),]
bacter.down<-bacter[-c(1:6,13:18,25:36,37:42,49:54,61:72),]
var.test(bacter.up$bacter,bacter.down$bacter)
t.test(bacter.up$bacter,bacter.down$bacter,var.equal=FALSE,paired=FALSE)

#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
########
# 
########