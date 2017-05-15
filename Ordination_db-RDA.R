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
pkgs.to.load<-c("phyloseq","ggplot2","DESeq2","vegan","grid","ggrepel")
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

cap0f<-capscale(
	formula=t(vst) ~ 1, 
	data=sd,
	distance="bray",
	add=TRUE
	)

cap1f<-capscale(
	formula=t(vst) ~ ., 
	data=sd,
	distance="bray",
	add=TRUE
	)

step.capf<-ordistep(cap0f,
	scope=formula(cap1f),
	direction="forward")
step.capf$anova
         # Df     AIC      F Pr(>F)   
# + log_Cl  1 -2.4396 2.2727  0.005 **
# + log_Si  1 -2.8660 2.0168  0.005 **
# + log_Sr  1 -3.5092 1.9713  0.005 **
# + log_Fl  1 -4.3039 1.8357  0.020 * 
plot(step.capf)
plot(step.capf,display=c("wa","cn"))
step.capf
RsquareAdj(step.capf)
# $r.squared
# [1] 0.5769007

# $adj.r.squared
# [1] 0.3351297
anova(step.capf)
         # Df SumOfSqs      F Pr(>F)    
# Model     4  0.45162 2.3861  0.001 ***
# Residual  7  0.33122                
anova(step.capf,by="term")
         # Df SumOfSqs      F Pr(>F)    
# log_Cl    1  0.14497 3.0638  0.001 ***
# log_Si    1  0.11677 2.4679  0.001 ***
# log_Sr    1  0.10302 2.1772  0.002 ** 
# log_Fl    1  0.08686 1.8357  0.026 *  
# Residual  7  0.33122    

ord.f<-ordinate(ps, "CAP",
	formula = ps ~ log_Cl + log_Si + log_Sr + log_Fl)
arrowmat <-vegan::scores(ord.f, display="bp")
arrowdf<-data.frame(labels=rownames(arrowmat),arrowmat)
arrow_map<-aes(xend=0.9*CAP1, yend=0.9*CAP2, x=0, y=0, shape=NULL, color=NULL, label=labels)
label_map<-aes(x = 1.2*0.9*CAP1, y=1.2*0.9*CAP2, shape=NULL, color=NULL, label=labels)
arrowhead<-arrow(length=unit(0.05,"npc"))
evals<-ord.f$CCA$eig

p.cap<-
plot_ordination(ps, ord.f, 
	type = "sample", color="Location"
	) +
	theme(legend.position="none") +
	geom_text(aes(label=SampleID), size=4, vjust=1.5) +
	geom_segment(arrow_map, size =0.5, 
		data = arrowdf, color="gray", arrow = arrowhead) +
	geom_text(label_map, size=4, data = arrowdf)+
	coord_fixed(sqrt(evals[2]/evals[1]))
ggsave(file.path(figs_path,"RDA.f.pdf"),p.cap)

#----------------------------------------------------------
########
# indicators
########
ind<-veganotu(ps.ind)

cap0i<-capscale(
	formula=ind ~ 1, 
	data=sd,
	distance="bray",
	add=TRUE
	)

cap1i<-capscale(
	formula=ind ~ ., 
	data=sd,
	distance="bray",
	add=TRUE
	)

step.capi<-ordistep(cap0i,
	scope=formula(cap1i),
	direction="forward")
step.capi$anova
         # Df    AIC      F Pr(>F)   
# + log_Cl  1 4.0109 3.6952  0.005 **
# + log_Si  1 3.1118 2.4595  0.005 **
# + log_Mg  1 2.8450 1.6633  0.010 **
plot(step.capi)
plot(step.capi,display=c("wa","cn"))
step.capi
RsquareAdj(step.capi)
# $r.squared
# [1] 0.5252421

# $adj.r.squared
# [1] 0.3472079
anova(step.capi)
         # Df SumOfSqs      F Pr(>F)    
# Model     3  0.78544 2.9502  0.001 ***
# Residual  8  0.70994 
anova(step.capi,by="term")
         # Df SumOfSqs      F Pr(>F)    
# log_Cl    1  0.40348 4.5466  0.001 ***
# log_Si    1  0.23435 2.6408  0.014 *  
# log_Mg    1  0.14761 1.6633  0.104    
# Residual  8  0.70994   

########
# plot capscale
ord.i<-ordinate(ps.ind, "CAP",
	formula = ps.ind ~ log_Cl + log_Si + log_Mg)
arrowmat <-vegan::scores(ord.i, display="bp")
arrowdf<-data.frame(labels=rownames(arrowmat),arrowmat)
arrow_map<-aes(xend=1*CAP1, yend=1*CAP2, x=0, y=0, shape=NULL, color=NULL, label=labels)
label_map<-aes(x = 1*1.1*CAP1, y=1*1.1*CAP2, shape=NULL, color=NULL, label=labels)
arrowhead<-arrow(length=unit(0.05,"npc"))
evals<-ord.i$CCA$eig

p.cap.ind<-
plot_ordination(ps.ind, ord.i, 
	type = "sample", color="Location"
	) +
	theme(legend.position="none") +
	geom_text_repel(aes(CAP1,CAP2,label=SampleID,size=4)) +
	geom_segment(arrow_map, size =0.5, 
		data = arrowdf, color="gray", arrow = arrowhead) +
	geom_text(label_map, size=4, data = arrowdf) +
	coord_fixed(sqrt(evals[2]/evals[1])) 
ggsave(file.path(figs_path,"RDA.i.pdf"),p.cap.ind,width=7,height=4.5)

#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
########
# 
########