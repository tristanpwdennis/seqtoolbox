#tristan dennis jan 2020
#load externally developed packages
library(dplyr)
library(ggplot2)
library(Rmisc)
library(reshape)
library(igraph)
library(ggpedigree)


###########################
#coverage data
##########################
setwd("~/Projects/DNDI/Data/mf-unamplified/coverage/")
metadata <- read.csv("~/Projects/DNDI/Data/metadata/metadata-mf-males-unamp.csv")
#read in all coverage files 'genome' section
#takes filename as a column for identifier
#chops off the file extension. edit the regex in the rbind command if you have different filenames
x<-NULL
for (file in (list.files(pattern = "*cov", full.names = TRUE))) {
  . <- read.table(file, header = FALSE) %>% subset(V1 == "genome")
  x<-rbind(x,data.frame(IND=gsub("_cov|./", "", file),data=.))
}

#add column names
colnames(x) <- c('IND', 'region', 'DP', 'NUM_BASES', 'TOTAL_BASES_IN_GENOME', 'FRAC')


#join metadata to coverage data based on column (you do not need to do this but it is helpful when grouping samples
#and selecting which ones to plot/facet in a plot
x<-left_join(x, metadata, by = c("IND" = "sample_name"))


head(t)

#discretise depth into bins otherwise the plot crashes when coverage is uneven (a lot of small FRAC covered by a lot of bases)
x$DPbin <- cut(x$DP, c(-0.1, 0, 1, 5, 10, 20, 30, 40, 50, 100, 200, 300, 400, 500, 1000, 100000), c("0", "1", "1-5", "5-10", "10-20", "20-30", "30-40", "40-50", "50-100", "100-200", "200-300", "300-400", "400-500", "500-1000", "More than 1000"))

#plot as stacked column plot
#add facet_wrap if you want to split the plot, otherwise remove
#if you remove the facet wrap the plot should work without the metadata

###plot coverage data as fraction/depth
ggplot(x, aes(x = x$IND, y = x$FRAC, fill = as.factor(DPbin))) +
  geom_bar(stat = 'identity') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(~ family_id, scales = "free", space = "free") +
  xlab('Sample name') + ylab('Fraction of genome covered') + labs(fill = 'Depth')

#based on the plots pick a threshold of individuals with relatively poor coverage
#in this case (also make histogram of fractions of coverage==0) let's chuck everything
#below 0.5

poorlysequencedinds <- x %>% filter(DP == 0, FRAC >0.5) %>% select(IND)

##########################
#load kinshipdata
##########################
#set working dir
setwd('/Users/tristanpwdennis/Projects/DNDI/Data/mf-unamplified/forngsrelate/ngsrelate-imputations/test/')

#read in sequencing metadata
metadata <- read.csv('~/Projects/DNDI/Data/metadata/metadata-all.csv')

#get the samples we are interested in
#in this case I have grouped samples by sequencing run (4 = mf, 2 = adult males)
metadata <- metadata %>% filter(sequencing_run == 4 | sequencing_run == 2) %>% select(sample_name, parent, village, nodule_id, sex)

#load all files as a list of dataframes
t <- lapply(list.files(pattern = ".res", full.names = T), function(x) {
  tryCatch(read.delim(x, header = TRUE), error=function(e) NULL)
})


###########################
#munging kinshipdata
##########################

#init df,loop over list of matrices
#calculating ci for each bootstrapped pair comparison
#this will catch every set of samples and we will have to subset using the metadata in a second

df<- data.frame(ida = NA, idb = NA, upper = NA, mean = NA, lower = NA, coverage = NA)
for(i in 1:length(t)) {
  #calculate confidence interval
  confint<-CI(t[[i]]$rab, ci = 0.95)
  #get mean coverage
  #store confint in new df
  df[i,3:5] <- confint
  #ida in new df
  df[i,1] <- paste(t[[i]]$ida[1])
  #idb in new df
  df[i,2] <- paste(t[[i]]$idb[1])
  #coverage too
  df[i,6] <- paste(as.numeric(t[[i]]$coverage[1]))
}


#add metadata to NGSRelate output
df <- left_join(df, metadata, by = c("ida" = "sample_name")) %>% left_join(., metadata, by = c("idb" = "sample_name"))

#add known relationship
#hacky stuff with mutate to get around the fact we have known adults in our dataset
df<-df %>% 
  mutate(., link = ifelse(.$sex.x == 'male' | .$sex.y == 'male', 'malecomp', 'mfcomp')) %>% 
  mutate(., knownrel = ifelse(.$parent.x == .$parent.y & .$link == 'mfcomp', 'fullorhalfsib', 'unknown')) %>% 
  mutate(., withinorbetweenfam = ifelse(.$parent.x == .$parent.y & .$link == 'mfcomp', 'withinfamily', 'betweenfamily'))

#add inferred relationship
df$foundrel <- cut(df$mean, c(-3, 0.0442, 0.0884, 0.177, 0.5, 1), c("unrelated", "third","second", "first", "twin"))

#let's filter out all comparisons between the males as they confuse things
#we can return to them later

df<-filter(df, link == 'mfcomp') 

df
hist(as.numeric(df$coverage))





#get rid of badly covered samples
s<-filter(df, coverage > 0.5)

#kinship data (known and found relationships and demographic info for all samples are now stored in dataframe 't')


s<-filter(s, link == 'mfcomp')

links <- s %>% select(ida, idb, mean, foundrel) 
nodes <- metadata %>% filter(sample_name %in% s$ida | sample_name %in% s$idb)



net <- graph_from_data_frame(links,vertices = nodes, directed = F)
E(net)$color <- as.factor(E(net)$foundrel)
#E(net)$width <- E(net)$ibd*10
V(net)$color <- as.factor((V(net)$parent))
coords <- layout_(net, nicely())
plot(net, vertex.size=5, layout = coords)

#get rid of shit sequence data
#based on poorly covered individuals identified in the coverage section

