#----------------------------------------------------------
#----------------------------------------------------------
# Notes or References:
# 
#----------------------------------------------------------
########
# FUNCTION: read raw data
########
raw_data<-function() {
	readRDS(file.path(data_path,"ps_orig.RDS"))
}

#----------------------------------------------------------
########
# FUNCTION: preprocess data (phyloseq object)
########
preprocessed_data<-function() {
	ps_orig<-raw_data()
	psB<-prune_taxa(taxa_sums(ps_orig)>0,ps_orig)
	ps0<-subset_taxa(psB,!is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
	ps_filtered<-filter_taxa(ps0, function(x) sum(x > 2) >= 2, TRUE)
	ps<-ps_filtered
	psra<-transform_sample_counts(ps,function(x) {x / sum(x)})
	list(ps_orig = ps_orig, ps = ps, psra = psra)
}

#----------------------------------------------------------
########
# FUNCTION: prepare the R environment, load packages, and import phyloseq objects
########
analysis_prep<-function(pkgs) {
	sapply(pkgs, require, character = TRUE)
	attach(preprocessed_data())
}

#----------------------------------------------------------
########
# vegan related functions
########

veganotu <- function(physeq) {
	require("vegan")
	OTU <- otu_table(physeq)
	if (taxa_are_rows(OTU)) {
		OTU <- t(OTU)
	}
	return(as(OTU, "matrix"))
}

vegansd <- function(physeq) {
	require("vegan")
	sd <- sample_data(physeq)
	return(as(sd,"data.frame"))
}

#----------------------------------------------------------
########
# qPCR related functions
# from:
# http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
#Helper functions
########
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}

#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
########
# 
########