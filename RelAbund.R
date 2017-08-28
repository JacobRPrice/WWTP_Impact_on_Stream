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
psra.mergeP.top

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

#----------------------------------------------------------
########
# plot relative abundance
########

ggsave(file.path(figs_path,"TopLevelRelAbund.f.eps"),
	plot_bar(psra.mergeP.top, 
		x="SampleID", fill="Phylum") +
		theme(legend.position="bottom") +
		theme(axis.text.y=element_text(angle=90,hjust=0.5)) +
		theme(axis.text.x=element_text(angle=90,vjust=0.5)) +
		theme(axis.title.x=element_blank()) +
		ylab("Top Phyla (Rel. Abund.)") +
		#theme(legend.margin=unit(0,"cm")) + 
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
		theme(axis.text.x=element_text(angle=90,vjust=0.5)) +
		theme(axis.title.x=element_blank()) +
		ylab("Top Orders (Rel. Abund.)") +
		#theme(legend.margin=unit(0,"cm")) + 
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
# 
########