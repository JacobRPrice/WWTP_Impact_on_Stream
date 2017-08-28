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

########
# relative abundance - Human-specific bacteroides (B. dorei)
psra
psra.human<-subset_taxa(psra, Species %in% c("dorei","dorei/fragilis"))
psra.human

########
# relative abundance - bacteroides spp.
########
psra
psra.bacter<-subset_taxa(psra,Genus=="Bacteroides")
psra.bacter
psra.bacter.glom<-tax_glom(psra.bacter,taxrank="Species",NArm=FALSE)
psra.bacter.glom
psra.bacter.glom.t<-psra.bacter.glom
tax_table(psra.bacter.glom.t)[,7]<-sub("/","/ \n",tax_table(psra.bacter.glom.t)[,7])

#----------------------------------------------------------
########
# plot qPCR and relative abundance together
########

ggsave(
	file.path(figs_path,"qPCR-&-Bact-&-Bdorei.eps"),
	grid.arrange(
		nrow=2,
		ggplot(human, aes(x=SampleID, y=human, fill=as.factor(site))) +
			geom_bar(position=position_dodge(), stat="identity") +
			geom_errorbar(aes(ymin=human-se, ymax=human+se),
										width=0.2,position=position_dodge(0.9)) +
			theme_bw() +
			ylab("B. dorei (Copy Number/100 mL)") +
			guides(fill=FALSE) +
			theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=0)) +
			theme(axis.text.y=element_text(angle=90,hjust=0.5)) +
			theme(axis.title.x=element_blank())+ 
			theme(axis.title.y=element_text(size=rel(0.9))) + 
			theme(axis.text.y=element_text(size=rel(0.9))) +
			scale_y_log10() + 
			ggtitle("A") +
			theme(plot.title=element_text(hjust=0.5))
		,
		ggplot(bacter, aes(x=SampleID, y=bacter, fill=as.factor(site))) +
			geom_bar(position=position_dodge(), stat="identity") +
			geom_errorbar(aes(ymin=bacter-se, ymax=bacter+se),
										width=0.2,position=position_dodge(0.9)) +
			theme_bw() +
			ylab("Bacteroides spp. (Copy Number/100 mL)") +
			guides(fill=FALSE) +
			theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=0)) +
			theme(axis.text.y=element_text(angle=90,hjust=0.5)) +
			theme(axis.title.x=element_blank()) + 
			theme(axis.title.y=element_text(size=rel(0.9))) + 
			theme(axis.text.y=element_text(size=rel(0.9))) +
			scale_y_log10() + 
			ggtitle("B") +
			theme(plot.title=element_text(hjust=0.5))
		,
		plot_bar(psra.human, 
						 x="SampleID", fill="Species") +
			ylab("B. dorei (Rel. Abund.)") +
			theme(legend.position="bottom") +
			theme(axis.title.x=element_blank()) +
			theme(axis.text.y=element_text(angle=90,hjust=0.5)) +
			theme(axis.text.x=element_text(angle=-90,vjust=0.5)) +
			theme(legend.title=element_blank()) + 
			theme(axis.title.y=element_text(size=rel(0.9))) +
			theme(axis.text.y=element_text(size=rel(0.9))) +
			#theme(axis.text.x=element_text(size=rel(0.9))) +
			theme(axis.text.x=element_blank()) +
			guides(fill=guide_legend(nrow=6)) +
			#theme(legend.text=element_text(size=rel(0.75))) +
			theme(legend.text=element_text(size=9)) +
			theme(legend.key.width=unit(0.45,"cm")) +
			theme(legend.key.height=unit(0.46,"cm")) +
			ggtitle("C") +
			theme(plot.title=element_text(hjust=0.5))
		,
		plot_bar(psra.bacter.glom.t, 
						 x="SampleID", fill="Species") +
			ylab("Bacteroides spp. (Rel. Abund.)") +
			theme(legend.position="bottom") +
			theme(axis.title.x=element_blank()) +
			theme(axis.text.y=element_text(angle=90,hjust=0.5)) +
			theme(axis.text.x=element_text(angle=-90,vjust=0.5)) +
			theme(legend.title=element_blank()) + 
			theme(axis.title.y=element_text(size=rel(0.9))) +
			theme(axis.text.y=element_text(size=rel(0.9))) +
			#theme(axis.text.x=element_text(size=rel(0.9))) +
			theme(axis.text.x=element_blank()) +
			guides(fill=guide_legend(nrow=6)) +
			#theme(legend.text=element_text(size=rel(0.75))) +
			theme(legend.text=element_text(size=9)) +
			theme(legend.key.width=unit(0.45,"cm")) +
			theme(legend.key.height=unit(0.46,"cm")) +
			ggtitle("D") +
			theme(plot.title=element_text(hjust=0.5))
	),
	width=190,
	height=240,
	units="mm"
	)

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