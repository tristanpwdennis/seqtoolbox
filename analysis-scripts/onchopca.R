library(adegenet)
library(vcfR)
library(dplyr)
library(ggplot2)
library(poppr)
library(spdep)

setwd('/Users/tristanpwdennis/Projects/DNDI/Data/mf-unamplified/ped/sequoia/')

t<-read.vcfR('final.vcf.gz')
obj <- vcfR2genlight(t)
metadata <- read.csv('~/Projects/DNDI/Data/metadata/metadata-all.csv', header = T)
nm <- obj$ind.names
nm<- as.data.frame(nm)

left_join(nm, metadata,by = c("nm$nm", "metadata$sample_id"))

popinfo<-nm %>% left_join(metadata, by = c("nm"='sample_name'))

strata(obj) <- popinfo                  
setPop(obj) ~Village
cov <- as.matrix(read.table('~/Projects/DNDI/Data/mf-unamplified/pca/pcangsd-cov.txt'))
pop <- read.csv('~/Projects/DNDI/Data/metadata/metadata-mf-males-unamp.csv')
e<-eigen(cov)

plot(e$vectors[,3:2],col=pop$patient_id)
legend("bottomright",fill=1:30,levels(pop$patient_id))

t<-cbind(metadata, as.data.frame(e$vectors[,1:2]))



s<-t %>% filter(!(sample_name %in% poorlysequencedinds$IND))
ggplot(s, aes(x=V1, y=V2, shape=village, color=parent)) +
  geom_point(aes(size = 5))
