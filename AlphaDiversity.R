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
pkgs.to.load<-c("phyloseq","ggplot2","gridExtra","vegan")
analysis_prep(pkgs.to.load)

# set seed
set.seed(20141004)

#----------------------------------------------------------
########
# rarefaction curves
########
# http://www.fromthebottomoftheheap.net/2015/04/16/drawing-rarefaction-curves-with-custom-colours/

veg.ps<-as.data.frame(otu_table(ps_orig))
dim(veg.ps)
names(veg.ps)[1:3]
rownames(veg.ps)[1:3]
raremax<-min(rowSums(veg.ps))
raremax

col<-colors<-c("red","red","green","green","blue","blue","cyan","cyan","magenta","magenta","orange","orange")
lty<-linetype<-c("solid","dashed","solid","dashed","solid","dashed","solid","dashed","solid","dashed","solid","dashed")

pdf(file.path(figs_path,"rarification.pdf"),width=6,height=6)
par(mar=c(3,3.0,0.5,0.5),
	mgp=c(1.5,0.5,0)
)
out<-rarecurve(veg.ps,
	step=10,
	sample=raremax,
	col=col,
	lty=lty,
	lwd=2,
	label=FALSE)
legend("bottomright",
	legend=rownames(veg.ps),
	col=col,
	lty=lty,
	lwd=2,
	cex=0.8,
	bty="n")
dev.off()

#----------------------------------------------------------
########
# Richness plots
########

AlphaDiversity<-
plot_richness(ps_orig, 
	x="site", 
	color="day", 
	measures=c("observed","Chao1","Shannon")) +
	theme(legend.position="bottom") + 
	theme(axis.text.y=element_text(angle=90,hjust=0.5)) +
  theme(axis.text.x=element_text(angle=90,hjust=0.5)) +
  xlab("Site") + 
  labs(color="Day")
ggsave(file.path(figs_path,"AlphaDiversity.eps"),
	AlphaDiversity,
	width=90,height=120,units="mm")

# gather info for table
sample_sums(ps_orig)
AlphaDiversity[[1]]$variable
AlphaDiversity[[1]]$value
AlphaDiversity[[1]]$se


targOrd<-c("Bacteroidales","Bifidobacteriales","Clostridiales")
ps_orig.ind<-subset_taxa(ps_orig, Order %in% targOrd)
AlphaDiversity.i<-
plot_richness(ps_orig.ind, 
	x="site", 
	color="day", 
	measures=c("observed","Chao1","Shannon")) +
	theme(legend.position="bottom") + 
	theme(axis.text.y=element_text(angle=90,hjust=0.5)) +
  theme(axis.text.x=element_text(angle=90,hjust=0.5)) +
  xlab("Site") + 
  labs(color="Day")
ggsave(file.path(figs_path,"Ind_AlphaDiversity.eps"),
	AlphaDiversity.i,
	width=90,height=120,units="mm")

# gather info for table
sample_sums(ps_orig.ind)
AlphaDiversity.i[[1]]$variable
AlphaDiversity.i[[1]]$value
AlphaDiversity.i[[1]]$se

#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
########
# 
########