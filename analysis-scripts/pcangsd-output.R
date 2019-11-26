#load externally developed packages
library(adegenet)
library(poppr)
library(vcfR)
library(ggbiplot)
library(dplyr)

#read metadata and select relevant info
pop <- read.csv('~/Projects/DNDI/Data/metadata/metadata-all.csv', header = T)
bamlist <- read.csv('Projects/DNDI/Data/metadata/list.txt', header = T)
bamlist<-filter(bamlist, sample_name != "ii-5")

#read covariance matrix
m <- as.matrix(read.table('~/Desktop/pcangsd-output/pcangsd-test.cov'), header = F)

#compute eigenvalues and eigenvectors
e <- eigen(m)

#plot PCs in base R
dfvectors <- cbind(sample_id = bamlist, as.data.frame(e$vectors))
eigenframe <- left_join(dfvectors, pop, by = "sample_name")
eigenframe <- filter(eigenframe, eigenframe$p_miss < 0.5)

#legend("topright",inset=c(-0.2,0), fill=1:3,levels(as.factor(pop$patient_ID)))

#create pca dataframe from metadata and eigenvectors, make prettier plot in ggplot2

ggplot(eigenframe, aes(x = eigenframe$V2, y = eigenframe$V4))+
  geom_point(aes(col = eigenframe$patient_id), size = 2)

