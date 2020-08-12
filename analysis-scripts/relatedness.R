#tristan dennis jan 2020
#load externally developed packages
library(dplyr)
library(ggplot2)
library(Rmisc)
library(reshape)
library(igraph)
library(hrbrthemes)
library(viridis)
library(sf)
library(maps)
library("rnaturalearth")
library("rnaturalearthdata")
library(RColorBrewer)
library("ggrepel")
library("ggspatial")


####################################################
#coverage data
####################################################
setwd("~/Projects/DNDI/Data/mf-unamplified/coverage/")
metadata <- read.csv("~/Projects/DNDI/Data/metadata/metadata-all.csv")
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

####################################################
#load kinshipdata
####################################################
#set working dir
setwd('/Users/tristanpwdennis/Projects/DNDI/Data/mf-unamplified/forngsrelate/ngsrelate-imputations/test/')

#read in sequencing metadata
metadata <- read.csv('~/Projects/DNDI/Data/metadata/metadata-all.csv')

#get the samples we are interested in
#in this case I have grouped samples by sequencing run (4 = mf, 2 = adult males)
metadata <- metadata %>% filter(sequencing_run == 4 | sequencing_run == 2) %>% select(sample_name, parent, village, nodule_id, sex, lat, long)

#load all files as a list of dataframes
t <- lapply(list.files(pattern = ".res", full.names = T), function(x) {
  tryCatch(read.delim(x, header = TRUE), error=function(e) NULL)
})


####################################################
#munging kinshipdata
####################################################

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
df$foundrel <- cut(df$mean, c(-3, 0.0442, 0.0884, 0.177, 1), c("unrelated", "third","second", "first"))

#let's filter out all comparisons between the males as they confuse things
#we can return to them later

####################################################
#plotting kinship data
####################################################


hist(as.numeric(df$coverage))

df %>% 
  filter(., link == 'mfcomp') %>% 
  filter(., coverage > 0.8) %>% 
  ggplot(., aes(x = withinorbetweenfam, y = mean)) +
  geom_boxplot(aes(fill = withinorbetweenfam), show.legend = FALSE) +
  theme_ipsum() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  geom_jitter(aes(color=foundrel), size=2, alpha=0.9) +
  xlab("Comparison level") +
  ylab("Mean rxy") + theme_grey(base_size = 22) + 
  ylim(0, 0.6) +
  scale_x_discrete(name ="Comparison Level", breaks=c("betweenfamily","withinfamily"), labels=c("Between Broods", "Within Broods")) +
  labs(colour="Inferred Relationship")
  
  

df %>% 
  filter(., link == 'malecomp') %>% 
  filter(., coverage > 0.8) %>% 
  ggplot(.) +
  geom_boxplot() +
  theme_ipsum() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  geom_jitter(aes(color=foundrel), size=2, alpha=0.9) +
  ylab("Mean coefficient of kinship") + theme_grey(base_size = 22) +
  xlab("Comparison level") + theme_grey(base_size = 22)



#create links and nodes for igraph plot
links <- df %>% 
  select(ida, idb, mean, foundrel, coverage, link, sex.x, sex.y) %>% 
  filter(link == 'malecomp') %>% 
  filter(mean > 0.125) %>% 
  filter(coverage > 0.8)

nodes <- metadata %>% filter(sample_name %in% links$ida | sample_name %in% links$idb)



links
#create graph object
net <- graph_from_data_frame(links,vertices = nodes, directed = F)

#colour links and nodes by whatever you like - in this case relattionship and nodule
E(net)$color <- as.factor(E(net)$foundrel)
V(net)$color <- as.factor((V(net)$nodule_id))
plot(net, vertex.size=5, vertex.label=NA)
####################################################
####################################################
#plotting relatedness over spatial scale
####################################################
world <- ne_countries(scale = "large", returnclass = "sf")
lakes <- ne_download(scale = 10, type = 'lakes', category = 'physical', returnclass = "sf")
river10 <- ne_download(scale = 10, type = 'rivers_lake_centerlines', category = 'physical', returnclass = "sf")
towns <- ne_download(scale = 10, type = 'populated_places', category = 'cultural', returnclass = 'sf')
villages <- st_as_sf((metadata %>% distinct(village, lat, long)), coords = c("long", "lat"), remove = FALSE, crs = 4326, agr = "constant")

ggplot(data = world) +
  geom_sf(fill = "antiquewhite1") +
  geom_sf(data = lakes, colour = '#a6bddb', fill = '#a6bddb') +
  geom_sf(data = river10, colour = '#a6bddb', fill = '#a6bddb') +
  geom_sf(data = villages) +
  scale_fill_viridis_c(trans = "sqrt", alpha = .4) +
  geom_text_repel(data = villages, aes(x = long, y = lat, label = village), fontface = "bold", nudge_x = c(0.5, 0.5, 1, 0.5)) +
  coord_sf(xlim = c(-4, 2), ylim = c(4, 12), expand = FALSE) +
  ggtitle("Sampling Sites", subtitle = "Four Villages in Southeastern Ghana") +
  annotation_scale(location = "bl", width_hint = 0.4) +
  xlab("Longitude") + ylab("Latitude") +
  theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "#a6bddb"))


ggplot(data = world) +
  geom_sf() +
  geom_sf(data = lakes, colour = '#a6bddb', fill = '#a6bddb') +
  geom_sf(data = river10, colour = '#a6bddb', fill = '#a6bddb') +
  geom_sf(data = villages) +
  scale_fill_viridis_c(trans = "sqrt", alpha = .4) +
  geom_text_repel(data = villages, aes(x = long, y = lat, label = village), fontface = "bold") +
  coord_sf(xlim = c(-2, -1.4), ylim = c(5.5, 6.5), expand = FALSE) +
  ggtitle("Sampling Sites", subtitle = "Four Villages in Southeastern Ghana") +
  annotation_scale(location = "bl", width_hint = 0.4) +
  xlab("Longitude") + ylab("Latitude") +
  theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", 
                                        size = 0.5), panel.background = element_rect(fill = "#a6bddb"))




####
?qnbinom()

