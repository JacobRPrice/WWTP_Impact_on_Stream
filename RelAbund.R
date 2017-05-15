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
theme_set(theme_bw())

# set seed
set.seed(20141004)

#----------------------------------------------------------
########
# Phylum Level Bar plot summary of relative abundance 
########
psra
psra.mergeP<-tax_glom(psra,taxrank="Phylum",NArm=FALSE)
psra.mergeP

psra.mergeP.top<-prune_taxa(taxa_sums(psra.mergeP)>.006*12,psra.mergeP)
Phylum.RA<-
plot_bar(
	psra.mergeP.top,
	x="SampleID", 
	fill="Phylum") +
	theme(legend.position="bottom") +
	theme(axis.text.y=element_text(angle=90,hjust=0.5)) +
	theme(axis.title.x=element_blank()) +
	ylab("Top Phyla (Rel. Abund.)")

#----------------------------------------------------------
########
# relative abundance - common gut flora
########

rank_names(psra)
table(tax_table(psra)[,4])
targOrd<-c("Bacteroidales","Bifidobacteriales","Clostridiales")
psra
ps.ind<-subset_taxa(psra, Order %in% targOrd)
ps.ind
Ind_Order.RA<-
plot_bar(
	tax_glom(ps.ind,taxrank="Order",NArm=FALSE),
	x="SampleID",
	fill="Order") +
	theme(legend.position="bottom") +
	theme(axis.text.y=element_text(angle=90,hjust=0.5)) +
	theme(axis.title.x=element_blank()) +
	ylab("Indicator Orders (Rel. Abund.)")

#----------------------------------------------------------
########
# plot relative abundance
########

ggsave(file.path(figs_path,"TopLevelRelAbund.f.eps"),
	plot_bar(psra.mergeP.top, 
		x="SampleID", fill="Phylum") +
		theme(legend.position="bottom") +
		theme(axis.text.y=element_text(angle=90,hjust=0.5)) +
		theme(axis.text.x=element_text(angle=-90,vjust=0.5)) +
		theme(axis.title.x=element_blank()) +
		ylab("Top Phyla (Rel. Abund.)") +
		theme(legend.margin=unit(0,"cm")) + 
		theme(axis.title.y=element_text(size=rel(0.8))) + 
		theme(axis.text.y=element_text(size=rel(0.8))) +
		theme(axis.text.x=element_text(size=rel(0.8),vjust=0.5)) +
		theme(legend.title=element_blank()) + 
		guides(fill=guide_legend(nrow=4)) +
		theme(legend.text=element_text(size=rel(0.9)))
	,
	width=90,
	height=130,
	units="mm"
)

ggsave(file.path(figs_path,"TopLevelRelAbund.i.eps"),
	plot_bar(tax_glom(ps.ind,taxrank="Order",NArm=FALSE), 
		x="SampleID", fill="Order") +
		theme(legend.position="bottom") +
		theme(axis.text.y=element_text(angle=90,hjust=0.5)) +
		theme(axis.text.x=element_text(angle=-90,vjust=0.5)) +
		theme(axis.title.x=element_blank()) +
		ylab("Top Orders (Rel. Abund.)") +
		theme(legend.margin=unit(0,"cm")) + 
		theme(axis.title.y=element_text(size=rel(0.8))) + 
		theme(axis.text.y=element_text(size=rel(0.8))) +
		theme(axis.text.x=element_text(size=rel(0.8),vjust=0.5)) +
		theme(legend.title=element_blank()) + 
		guides(fill=guide_legend(nrow=3)) +
		theme(legend.text=element_text(size=rel(0.9)))
	,
	width=90,
	height=120,
	units="mm"
)

#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
########
# relative abundance - Human-specific bacteroides (B. dorei)
########
psra
psra.human<-subset_taxa(psra, Species %in% c("dorei","dorei/fragilis"))
psra.human

ra.human<-
plot_bar(psra.human,"SampleID",fill="Species") +
	ylab("B. dorei (Rel. Abund.)") +
	theme(legend.position="bottom") +
	theme(axis.text.y=element_text(angle=90,hjust=0.5)) +
	theme(axis.title.x=element_blank()) +
	theme(legend.text=element_text(size=7))

#----------------------------------------------------------
########
# relative abundance - bacteroides spp.
########

psra
psra.bacter<-subset_taxa(psra,Genus=="Bacteroides")
psra.bacter
psra.bacter.glom<-tax_glom(psra.bacter,taxrank="Species",NArm=FALSE)
psra.bacter.glom

ra.bacter<-
plot_bar(psra.bacter.glom,"SampleID",fill="Species") +
	ylab("Bacteroides spp. (Rel. Abund.)") +
	theme(legend.position="bottom") +
	theme(axis.text.y=element_text(angle=90,hjust=0.5)) +
	theme(axis.title.x=element_blank()) +
	theme(legend.text=element_text(size=7))

psra.bacter.glom.t<-psra.bacter.glom
tax_table(psra.bacter.glom.t)[,7]<-sub("/","/ \n",tax_table(psra.bacter.glom.t)[,7])
ra.bacter.t<-
plot_bar(psra.bacter.glom.t,"SampleID",fill="Species") +
	ylab("Bacteroides spp. (Rel. Abund.)") +
	theme(legend.position="bottom") +
	theme(axis.text.y=element_text(angle=90,hjust=0.5)) +
	theme(axis.title.x=element_blank()) +
	theme(legend.text=element_text(size=7))

#----------------------------------------------------------
########
# plot relative abundance of B. dorei and Bacteroides
########

ggsave(file.path(figs_path,"RelAbund-B.dorei_Bacter.eps"),
	grid.arrange(nrow=2,
	ra.human + 
		theme(legend.margin=unit(0,"cm")) +
		theme(legend.title=element_blank()) +
		theme(axis.title.y=element_text(size=rel(0.9))) +
		theme(axis.text.y=element_text(size=rel(0.8))) +
		theme(axis.text.x=element_text(size=rel(0.9))) +
 		theme(legend.position="bottom")  + 
		theme(legend.key.width=unit(0.45,"cm")) +
		theme(legend.key.height=unit(0.46,"cm")) + ggtitle("A")
		,
	ra.bacter.t + 
		theme(legend.margin=unit(0,"cm")) + 
		theme(axis.title.y=element_text(size=rel(0.9))) +
		theme(axis.text.y=element_text(size=rel(0.8))) +
		theme(axis.text.x=element_text(size=rel(0.9))) +
		theme(legend.title=element_blank()) + 
		guides(fill=guide_legend(nrow=6)) + 
  		theme(legend.position="bottom") +
		theme(legend.key.width=unit(0.45,"cm")) +
		theme(legend.key.height=unit(0.46,"cm")) + ggtitle("B")

	)
	, 
	width=90,
	height=200,
	units="mm"
)

#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
########
# 
########