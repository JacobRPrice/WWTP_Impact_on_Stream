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
pkgs.to.load<-c("phyloseq","ggplot2","gridExtra","DESeq2","ggrepel")
analysis_prep(pkgs.to.load)

# set seed
set.seed(20141004)

#----------------------------------------------------------
########
# prep data
########
# save unaltered ps 
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

#----------------------------------------------------------
########
# ordinations
########

########
# Principal Components Analysis (PCA)
########
out<-ordinate(ps,method="MDS",distance="euclidean")
evals<-out$values$Eigenvalues
plot.pca<-
plot_ordination(ps, out, color="Location") +  
	coord_fixed(sqrt(evals[2]/evals[1])) +
	theme(legend.position="none") + 
	geom_text_repel(aes(Axis.1, Axis.2, label = SampleID))
ggsave(file.path(figs_path,"PCA.eps"),plot.pca,width=90,height=80,units="mm")

########
# DPCoA 
########
out.dpcoa<-ordinate(ps,method="DPCoA")
#saveRDS(out.dpcoa,file.path(output_path,"out.dpcoa.RDS"))
#out.dpcoa<-readRDS(file.path(output_path,"out.dpcoa.RDS"))
evals<-out.dpcoa$eig
ptmat<-vegan::scores(out.dpcoa,display="sites")
ptdf<-data.frame(labels=rownames(ptmat),ptmat)
ptmap<-aes(xend=Axis1,yend=Axis2,x=0,y=0,shape=NULL,col=NULL,label=labels)
labelmap<-aes(x=10*Axis1,y=10*Axis2,shape=NULL,col=NULL,label=labels)

###
# for plotting seperately
tax_table(ps)[,2]<-sub("_","\n",tax_table(ps)[,2])
p.dpcoa.overlay1<-
plot_ordination(ps,
	out.dpcoa,
	type="taxa",
	color="Phylum") +
	geom_text(labelmap, size=3,fontface="bold", data=ptdf)+
	coord_fixed(sqrt(evals[2]/evals[1])) +
	theme(legend.position="bottom") + 
	theme(legend.title=element_blank()) +
	theme(legend.margin=unit(0,"cm")) +
	theme(legend.text=element_text(size=9)) 

ggsave(file.path(figs_path,"DPCoA_overlay.eps"),p.dpcoa.overlay1,width=190,height=170,units="mm")

#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
########
# INDICATORS
########
#----------------------------------------------------------
########
# ordinations
########

########
# Principal Components Analysis (PCA)
########
out<-ordinate(ps.ind,method="MDS",distance="euclidean")
evals<-out$values$Eigenvalues
plot.pca.ind<-
plot_ordination(ps.ind, out, color="Location") + 
	coord_fixed(sqrt(evals[2]/evals[1])) +
	theme(legend.position="none") + 
	geom_text_repel(aes(Axis.1, Axis.2, label = SampleID))
ggsave(file.path(figs_path,"Ind_PCA.pdf"),plot.pca.ind,width=6.5,height=3.0)

########
# DPCoA 
########
out.dpcoa.i<-ordinate(ps.ind,method="DPCoA")
#saveRDS(out.dpcoa.i,file.path(output_path,"Ind_out.dpcoa.RDS"))
#out.dpcoa.i<-readRDS(file.path(output_path,"Ind_out.dpcoa.RDS"))
evals.i<-out.dpcoa.i$eig
ptmat.i<-vegan::scores(out.dpcoa.i,display="sites")
ptdf.i<-data.frame(labels=rownames(ptmat.i),ptmat.i)
ptmap.i<-aes(xend=Axis1,yend=Axis2,x=0,y=0,shape=NULL,col=NULL,label=labels)
labelmap.i<-aes(x=5*Axis1,y=5*Axis2,shape=NULL,col=NULL,label=labels)

###
# for plotting seperately
#tax_table(ps.ind)[,2]<-sub("_","\n",tax_table(ps.ind)[,2])
p.dpcoa.overlay1.i<-
plot_ordination(ps.ind,
	out.dpcoa.i,
	type="taxa",
	color="Order") +
	geom_text(labelmap.i, size=3,data=ptdf.i)+
	coord_fixed(sqrt(evals.i[2]/evals.i[1])) +
	theme(legend.position="bottom") + 
	theme(legend.title=element_blank()) +
	theme(legend.margin=unit(0,"cm")) +
	theme(legend.text=element_text(size=9)) 

ggsave(file.path(figs_path,"Ind_DPCoA_overlay.eps"),p.dpcoa.overlay1.i,width=190,height=110,units="mm")

#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
########
# 
########