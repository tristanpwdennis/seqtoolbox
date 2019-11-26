#################################################################################################################################################
#####PACKAGES AND INIT
#################################################################################################################################################

library(dplyr)
library(ggplot2)
library(lmtest)
library(igraph)


sequencing_data <- read.table('~/Projects/DNDI/Data/angsd-output/angsd-output-4/complete-chr/filtered/ngsrelate/newres', header = T)
metadata <- read.csv('~/Projects/DNDI/Data/metadata/metadata-all.csv')

#################################################################################################################################################
#####JOIN METADATA TO SEQUENCING DATA
#################################################################################################################################################

b <- metadata %>% select(sample_name, "doc_b" = mean_doc_prefiltered, 'family_b' = family_id, 'sex_b' = sex)
a <- metadata %>% select(sample_name, "doc_a" = mean_doc_prefiltered, 'family_a' = family_id, 'sex_a' = sex)
f <- left_join(sequencing_data, a, by = c("ida" = "sample_name"))
f <- left_join(f, b, by = c("idb" = "sample_name"))

#use a 'cut' operation to create a factor for each 'division' in our continuous scale (of IBD)
f$founddegree <- cut(f$rab, c(-0.1, 0.007, 0.03, 0.1, 0.125, 0.7, 1), c("unrelated", "third", "second","firstcousin", "fullorhalfsib", "parent-offspring"))

#get total coverage
#and coverage difference
f$total_doc = f$doc_a + f$doc_b
f$diff_doc <- abs(f$doc_a - f$doc_b)

plot(f$cov, f$rab)

t <- filter(f, f$cov > 0.3)



#get total missing
#and missingness difference

f$total_missing = f$missing_a + f$missing_b
f$diff_missing = abs(f$missing_a - f$missing_b)

#################################################################################################################################################
#####ADD UNIFORMITY DATA TO METADATA
#################################################################################################################################################
#READ IN CSVs and store in x
setwd("~/Projects/DNDI/Data/coverage/females_mf/filtered/")
all_files <- list.files(pattern = "*.cov", full.names = TRUE)

x<-NULL
for (file in all_files) {
  the_data <- read.table(file, header = FALSE) %>% subset(V1 == "genome")
  x<-rbind(x,data.frame(IND=gsub("srt.dd.bam.cov|./", "", file),data=the_data))
}

colnames(x) <- c('IND', 'region', 'DP', 'NUM_BASES', 'TOTAL_BASES_IN_GENOME', 'FRAC')

#extract proportion zero covered AND ADD TO SEQUENCING DATA
zerocov <- x %>% filter(DP == 0)

f <- left_join(f, zerocov, by = c("ida" = "IND"))
colnames(f)[colnames(f)=="FRAC"] <- "FRAC_A"
f <- left_join(f, zerocov, by = c("idb" = "IND"))
colnames(f)[colnames(f)=="FRAC"] <- "FRAC_B"
f$total_frac = (f$FRAC_A + f$FRAC_B)

######################################################################################################
#################TEST ONE - TEST SUCCESS OF KINSHIP ESTIMATION - NEED TO TIDY
######################################################################################################

d <- f 

#horrible code to extract known relationships from our matrix
#then add known degree to them
G <- filter(d, grepl("61G", ida) & grepl("61G", idb)) 
G$knowndegree <- paste0("fullorhalfsib", G$degree)
B <- filter(d, grepl("63B", ida) & grepl("63B", idb))
B$knowndegree <- paste0("fullorhalfsib", B$degree)
C <- filter(d, grepl("63C", ida) & grepl("63C", idb))
C$knowndegree <- paste0("fullorhalfsib", C$degree)
ii <- filter(d, grepl("^ii-", ida) & grepl("^ii-", idb))
ii$knowndegree <- paste0("fullorhalfsib", ii$degree)
iii <- filter(d, grepl("iii", ida) & grepl("iii", idb))
iii$knowndegree <- paste0("fullorhalfsib", iii$degree)
relmatrix <- rbind(G, C, ii, iii)

relmatrix %>% select(ida, idb, rab, founddegree, knowndegree)

#do the test
relmatrix$testone <- relmatrix$founddegree == relmatrix$knowndegree

#convert boolean values to numeric '1' for pass and '0' for fail
relmatrix$testone <- as.numeric(relmatrix$testone)

#how many passed?
sum(relmatrix$testone)

list <- relmatrix %>% filter(testone == 1) %>% distinct(ida)
#plot the characteristics of successful identifications

goodmetadata <- left_join(list, metadata, by = c("ida" = "sample_name"))

goodmetadata


######################################################################################################
#################MODELLING UNIFORMITY#################################################################
######################################################################################################

mdata <- relmatrix %>% select(ida, idb, rab, testone,coverage)


plot(mdata$coverage, mdata$testone)

plot(mdata$total_frac, mdata$testone)

m0 <- glm(mdata$testone ~ 1, family = 'binomial')
m1 <- glm(mdata$testone ~ mdata$coverage, family = 'binomial')

1-pchisq(logLik(m1)-logLik(m0), 2)

summary(m1)

exp(-0.8275)


summary(m1)
1-pchisq(243.47, 188)

m2 <- glm(mdata$testone ~mdata$total_frac)




summary(m1)


boxplot(mdata$coverage)
boxplot(mdata$total_frac)











######################################################################################################
#################PLOTS AND MODELS####################################################################
######################################################################################################

#extract data we want
#we want coverage, number of reads, and family
#extract data as two separate dataframes from metadata then perform two separate joins

#data ready to model
mdata <- relmatrix %>% select(total_missing, testone) 
plot(mdata)

plot(metadata$family_id, metadata$p_miss)
plot(metadata$family_id, metadata$F_filtered)
plot(metadata$family_id, metadata$no_miss)




mdata <- na.exclude (mdata)
#two models FOR DOC ON SUCCESS
#null
m0 <- glm(mdata$test~1,family = binomial)
#coverage explaining variation in successful classification
m1 <- glm(mdata$test~1,family = binomial)
m2 <- glm(mdata$test~mdata$total_missing,family = binomial)

lrtest(m1, m2)

summary(m2)

#LRT
logLik(m0)
logLik(m1)
1-pchisq(2*(-139.4155--141.4184),1)
mdata

#two models FOR family ON DOC
#null
mt <- metadata %>% select(family_id, mean_doc, sample_name)
plot(mt, x = mt$family_id, y = mt$mean_doc)
m1 <- glm(mt$family_id ~ mt$mean_doc, family = 'binomial')
m0 <- glm(mt$family_id ~ 1, family = 'binomial')
summary(m0)
summary(m1)
logLik(m0)
logLik(m1)
1-pchisq(2*(-34.37247--37.47765),1)


mt <- relmatrix %>% select(family_a, mean_doc, tests)


#plot family and tests
ggplot(relmatrix, aes(x=family_b, y = test, fill=family_a))+
  geom_bar(stat='identity')

#does family explain variation in successful identifications
m0 <- glm(relmatrix$test ~ 1)
m1 <- glm(relmatrix$test ~ relmatrix$family_a)
summary(m1)


logLik(m0)
logLik(m1)
1-pchisq(2*(-144.4992--148.2308),6)
#no, it is insignificant overall, although there are family specific effects (61G significantly more but overall the model fails to explain this variation)

#Select and plot doc by family
mt <- metadata %>% select(family_id, mean_doc, sample_name)
V <- filter(mt, !grepl("Sample", sample_name))
ggplot(V,aes(x = V$sample_name, y = V$mean_doc, fill = V$family_id, color = V$family_id)) +
  geom_bar(stat = 'identity')
mt <- metadata %>% select(family_id, mean_doc, sample_name)


######################################################################################################
#################PLOTTING KINSHIP NETWORKS############################################################
######################################################################################################
###NEED TO CLEAN THIS UP TO MAKE SURE WE'RE NOT PLOTTING PISH

####
#get clean samples
#
f2<- filter(f, f$coverage > 0.4, f$rab > 0.1)

setwd("~/Desktop/")
write.csv(f2, file = "relmatrix.csv")

nodes <- metadata %>% filter(sample_name %in% f2$ida | sample_name %in% f2$idb)

#get nodes
 #%>% select(sample_name, family_id, parent, patient_id)
#nodes <- nodes %>% filter(!grepl("Sample-5", sample_name))

links <- f2 %>% select(ida, idb, rab, founddegree)
links

net <- graph_from_data_frame(links,vertices = nodes, directed = F)


#check vertex and edge attributeß
vertex_attr(net)
edge_attr(net)

E(net)$color <- as.factor(E(net)$founddegree)
#E(net)$width <- E(net)$ibd*10
V(net)$color <- as.factor((V(net)$patient_id))

net <- delete_edges(net, E(net)[rab < 0.05])

coords <- layout_(net, nicely())

plot(net, vertex.size=5, layout = coords)

sixonegee = c('61G-10', '61G-3',  '61G-5',  '61G-6',  '61G-7','61G-8')
sixthreecee <- c('63C-10', '63C-3',  '63C-4',  '63C-5',  '63C-9')
three <- c('iii-2')
sampleone <- c('Sample-1', 'Sample-2')

c(sixonegee, sixthreecee, three, sampleone)


plot(delete_vertices(net, c(sixonegee, sixthreecee, three, sampleone)), vertex.label=NA)
