#!/bin/bash

# filter raw reads and derep
date
Rscript ./R/01.1.R 

# dada
date
Rscript ./R/01.2-F.R &
Rscript ./R/01.2-R.R &
wait %1 %2

# merge F/R and construct seq table
date
Rscript ./R/01.3.R 

# build phy tree and tax table
# About 24 hrs
date
Rscript ./R/01.4-PhyTree.R &
Rscript ./R/01.4-Tax.R &
wait %1 %2

# prep sample data.frame
date
Rscript ./R/01.4-SampDF.R 

# build phyloseq object
date
Rscript ./R/01.5.R
