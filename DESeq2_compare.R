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
pkgs.to.load<-c("phyloseq","ggplot2","gridExtra")
analysis_prep(pkgs.to.load)
theme_set(theme_bw())

# set seed
set.seed(20141004)

#----------------------------------------------------------
########
# 
########
sigtab<-readRDS(file.path(output_path,"sigtab.RDS"))
sigtab.i<-readRDS(file.path(output_path,"Ind_sigtab.RDS"))

ls.UPvsDown<-rownames(sigtab)
ls.UPvsDown.i<-rownames(sigtab.i)

psra.deseq<-prune_taxa(ls.UPvsDown,psra)
psra.deseq.i<-prune_taxa(ls.UPvsDown.i,psra)

#how many are in common?
length(intersect(taxa_names(psra.deseq),taxa_names(psra.deseq.i)))
# 11 are in common
targOrd<-c("Bacteroidales","Bifidobacteriales","Clostridiales")
ntaxa(psra.deseq)
ntaxa(subset_taxa(psra.deseq, Order %in% targOrd))

ggsave(
	file.path(figs_path,"DESeq.UPvsDown_bar.f.eps"),
	plot_bar(tax_glom(psra.deseq,taxrank="Phylum"),fill="Phylum") +
	ylab("Relative Abundance") +
	theme(axis.title.x=element_blank()) +
	theme(axis.text.y = element_text(angle=90, hjust=0.5)) +
	theme(axis.text.x = element_text(angle=-90, vjust=0.5)) +
	scale_fill_brewer(palette="Paired")
	,
	width=90,height=125,units="mm"
)

ggsave(
	file.path(figs_path,"DESeq.UPvsDown_bar.i.eps"),
	plot_bar(tax_glom(psra.deseq.i,taxrank="Order"),fill="Order") +
	ylab("Relative Abundance") +
	theme(axis.title.x=element_blank()) +
	theme(axis.text.y = element_text(angle=90, hjust=0.5)) +
	theme(axis.text.x = element_text(angle=-90, vjust=0.5)) +
	scale_fill_brewer(palette="Paired")
	,
	width=90,height=125,units="mm"
)


#----------------------------------------------------------
########
# get percentages
########

########
# get total rel abund for full community
sample_sums(psra.deseq)

########
# get total rel abund for indicator community
sample_sums(psra.deseq.i)

#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
########
# 
########