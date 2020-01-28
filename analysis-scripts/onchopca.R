library(adegenet)
library(vcfR)
library(dplyr)
library(ggplot2)
library(poppr)
library(spdep)

setwd('/Users/tristanpwdennis/Projects/DNDI/Data/mf-unamplified/ped/sequoia/')

t<-read.vcfR('final.vcf.gz')
obj <- vcfR2genlight(t)
metadata <- read.csv('~/Projects/DNDI/Data/metadata/metadata-mf-males-unamp.csv', header = T)
nm <- obj$ind.names
nm<- as.data.frame(nm)

left_join(nm, metadata,by = c("nm$nm", "metadata$sample_name"))

popinfo<-nm %>% left_join(metadata, by = c("nm"='sample_name'))

strata(obj) <- popinfo                  
setPop(obj) ~Village

