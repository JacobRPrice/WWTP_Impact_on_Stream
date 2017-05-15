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
# ps
# ps.ind<-subset_taxa(ps, Order %in% targOrd)
# ps.ind

sample_names(ps)
keepsamps<-c(
	"2.DSUG.1","2.DSUG.2",
	"3.USAmb.1","3.USAmb.2"
	)
psf<-prune_samples(keepsamps,ps)
psf
psf<-prune_taxa(taxa_sums(psf)>0,psf)
psf

sample_data(psf)$Position
sample_data(psf)$Location

# Test for the effect of Position, controlling for the effect of Location
psfdds<-phyloseq_to_deseq2(psf,~ Location)
psfdds<-DESeq(psfdds, test="Wald", fitType="parametric")
resultsNames(psfdds)
res<-results(psfdds, alpha=0.1)
summary(res)
mcols(res, use.names=TRUE)
mcols(res)$description
# How many otu's had adjusted p-values below 
sum(res$padj<0.05,na.rm=TRUE) 
res.tab<-cbind(as(res, "data.frame"), as(tax_table(psf), "matrix"))
write.csv(as.data.frame(res.tab), file=file.path(output_path,"res.tab.DESeq.2vs3.txt"))

alpha<-0.05
sigtab<-res[which(res$padj < alpha),]
sigtab<-cbind(as(sigtab, "data.frame"), as(tax_table(psf)[rownames(sigtab),], "matrix"))
dim(sigtab) # 206 differentially abundant taxa!
colnames(sigtab)
head(sigtab)
print.data.frame(sigtab[,c(8,9,10,11,12)],row.names=FALSE) # just taxonomy
saveRDS(sigtab,file=file.path(output_path,"sigtab.2vs3.RDS"))
write.csv(as.data.frame(sigtab), file=file.path(output_path,"sigtab.DESeq.2vs3.txt"))

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
sigtab.2vs3<-sigtab
# plot
plot.DESeq.2vs3<-
ggplot(sigtab.2vs3, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=4) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + 
	geom_hline(yintercept=0) +
	theme(axis.text.y = element_text(angle =90, hjust=0.5))
ggsave(file.path(figs_path,"DESeq.2vs3.eps"),plot.DESeq.2vs3,width=8,height=5.65)
# ggsave(file.path(figs_path,"DESeq.UPvsDOWN.eps"),plot.DESeq.UPvsDOWN,width=190,height=120,units="mm")


#sigtab.2vs3<-sigtab
#sigtab.2vs3<-readRDS(file.path(output_path,"sigtab.2vs3.RDS"))
sigtab<-readRDS(file.path(output_path,"sigtab.RDS"))
sigtab.i<-readRDS(file.path(output_path,"Ind_sigtab.RDS"))

ls.UPvsDown.2vs3<-rownames(sigtab.2vs3)
ls.UPvsDown<-rownames(sigtab)
ls.UPvsDown.i<-rownames(sigtab.i)

psra.deseq.2vs3<-prune_taxa(ls.UPvsDown.2vs3,psra)
psra.deseq<-prune_taxa(ls.UPvsDown,psra)
psra.deseq.i<-prune_taxa(ls.UPvsDown.i,psra)

dim(sigtab.2vs3) #205 taxa found to be differentially abundant between 2 and 3
names(sigtab.2vs3)
table(sigtab.2vs3$log2FoldChange>0) # almost evenly split between up and down 

table(tax_table(psra.deseq.2vs3)[,2])
 # Actinobacteria Armatimonadetes   Bacteroidetes        Chlorobi   Cyanobacteria   Fibrobacteres 
             # 21               1              55               1              34               1 
     # Firmicutes     Nitrospirae  Planctomycetes  Proteobacteria         SHA-109    Spirochaetae 
              # 2               2               1              77               1               1 
            # TM6 Verrucomicrobia 
              # 2               6 
(21+55+34+77)/205

length(intersect(taxa_names(psra.deseq),taxa_names(psra.deseq.2vs3)))
# 53 of the taxa overlap between differentially abundant up/down vs 2/3
sigtab.2vs3[intersect(taxa_names(psra.deseq),taxa_names(psra.deseq.2vs3)),"log2FoldChange"]
# all of them are negative
# these are the taxa that are 
#	higher in abundance immediately downstream
#	and then decrease in abundance with distance
sigtab[intersect(taxa_names(psra.deseq),taxa_names(psra.deseq.2vs3)),"log2FoldChange"]
# double checking that all of the differentially abundant overlapping taxa are MORE abundant. Which is true in this case. 

# what is the distribution of phyla that are increasing and decreasing 
table(sigtab.2vs3[which(sigtab.2vs3$log2FoldChange>0),"Phylum"])
 # Actinobacteria Armatimonadetes   Bacteroidetes        Chlorobi   Cyanobacteria   Fibrobacteres 
             # 19               1              36               0              20               0 
     # Firmicutes     Nitrospirae  Planctomycetes  Proteobacteria         SHA-109    Spirochaetae 
              # 0               0               0              21               0               0 
            # TM6 Verrucomicrobia 
              # 0               5 

table(sigtab.2vs3[which(sigtab.2vs3$log2FoldChange<0),"Phylum"])
 # Actinobacteria Armatimonadetes   Bacteroidetes        Chlorobi   Cyanobacteria   Fibrobacteres 
              # 2               0              19               1              14               1 
     # Firmicutes     Nitrospirae  Planctomycetes  Proteobacteria         SHA-109    Spirochaetae 
              # 2               2               1              56               1               1 
            # TM6 Verrucomicrobia 
              # 2               1 


psra.intersect<-prune_taxa(intersect(taxa_names(psra.deseq),taxa_names(psra.deseq.2vs3)),psra)
psra.intersect
table(tax_table(psra.intersect)[,6])



length(intersect(taxa_names(psra.deseq.i),taxa_names(psra.deseq.2vs3)))
# 5 taxa are shared 
sigtab.2vs3[intersect(taxa_names(psra.deseq.i),taxa_names(psra.deseq.2vs3)),"log2FoldChange"]
# all of the 2vs3 are negative
sigtab.i[intersect(taxa_names(psra.deseq.i),taxa_names(psra.deseq.2vs3)),"log2FoldChange"]
# and all of the up/down are positive


#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
########
# 
########