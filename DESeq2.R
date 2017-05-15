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
pkgs.to.load<-c("phyloseq","ggplot2","DESeq2")
analysis_prep(pkgs.to.load)

# set seed
set.seed(20141004)

#----------------------------------------------------------
########
# Up vs Down
# upstream vs downstream of WWTP
########
# targOrd<-c("Bacteroidales","Bifidobacteriales","Clostridiales")
# ps
# ps.ind<-subset_taxa(ps, Order %in% targOrd)
# ps.ind

sample_names(ps)
keepsamps<-c(
	"1.USUG.1","1.USUG.2",
	"2.DSUG.1","2.DSUG.2",
	"3.USAmb.1","3.USAmb.2",
	"4.DSAmb.1","4.DSAmb.2"
	)
psf<-prune_samples(keepsamps,ps)
psf
psf<-prune_taxa(taxa_sums(psf)>0,psf)
psf

sample_data(psf)$Position
sample_data(psf)$Location

# Test for the effect of Position, controlling for the effect of Location
psfdds<-phyloseq_to_deseq2(psf,~ Location + Position)
psfdds<-DESeq(psfdds, test="Wald", fitType="parametric")
resultsNames(psfdds)
res<-results(psfdds, alpha=0.1)
summary(res)
mcols(res, use.names=TRUE)
mcols(res)$description
# How many otu's had adjusted p-values below 
sum(res$padj<0.05,na.rm=TRUE) 
res.tab<-cbind(as(res, "data.frame"), as(tax_table(psf), "matrix"))
write.csv(as.data.frame(res.tab), file=file.path(output_path,"res.tab.DESeq.UPvsDOWN.txt"))

alpha<-0.05
sigtab<-res[which(res$padj < alpha),]
sigtab<-cbind(as(sigtab, "data.frame"), as(tax_table(psf)[rownames(sigtab),], "matrix"))
dim(sigtab)
colnames(sigtab)
#head(sigtab) 
print.data.frame(sigtab[,c(8,9,10,11,12)],row.names=FALSE) # just taxonomy
saveRDS(sigtab,file=file.path(output_path,"sigtab.RDS"))
write.csv(as.data.frame(sigtab), file=file.path(output_path,"sigtab.DESeq.UPvsDOWN.txt"))

# let's graph the results
#theme_set(theme_bw())
scale_fill_discrete<-function(palname="Set1",...) {
	scale_fill_brewer(palette=palname,...)
}
# Phylum order
x<-tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x<-sort(x,TRUE)
sigtab$Phylum<-factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x<-tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x<-sort(x,TRUE)
sigtab$Genus<-factor(as.character(sigtab$Genus), levels=names(x))
# plot
plot.DESeq.UPvsDOWN<-
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + 
	geom_hline(yintercept=0) +
	theme(axis.text.y = element_text(angle =90, hjust=0.5))
ggsave(file.path(figs_path,"DESeq.UPvsDOWN.eps"),plot.DESeq.UPvsDOWN,width=190,height=120,units="mm")

#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
########
# 
########