#coverage histograms
#adapted with thanks from Aaron Quinlan's Bedtools Protocols
#Tristan Dennis July 2019

library(ggplot2)
library(dplyr)
library(colorbrewer)

setwd("/Users/tristanpwdennis/Projects/DNDI/Data/barcoding_preliminary/coverage/")

#hacky import histograms as df and add headers, individual names
#add as many individuals here as needed
#NOTE: add a loop or something at some point
  
onecov <- read.delim('S1_cov', header = FALSE, sep = "\t")
twocov <- read.delim('S2-cov', header = FALSE, sep = '\t')
threecov <- read.delim('S3-cov', header = FALSE, sep = '\t')
fourcov <- read.delim('S4-cov', header = FALSE, sep = '\t')
fivecov <- read.delim('S5-cov', header = FALSE, sep = '\t')
#add column with individual name
gonecov  <- mutate(onecov, IND = '63B_AF')
gtwocov <- mutate(twocov, IND = '61G_AF')
gthreecov <- mutate(threecov, IND = '61G_MF1')
gfourcov <- mutate(fourcov, IND = '61G_MF2')
gfivecov <- mutate(fivecov, IND = '61G_MF3')

ggonecov = gonecov[gonecov[,1] == 'genome',]
ggtwocov = gtwocov[gtwocov[,1] == 'genome',]
ggthreecov = gthreecov[gthreecov[,1] == 'genome',]
ggfourcov = gfourcov[gfourcov[,1] == 'genome',]
ggfivecov = gfivecov[gfivecov[,1] == 'genome',]


#u can operate on any of the dfs individually
#anything from here is wankiness for ggplot

# plot a quick density function in base R for the genome-wide coverage
#merge data and add headers
cov <- rbind(gonecov, gtwocov, gthreecov, gfourcov, gfivecov)
headers=c("genome", "DP", "NUM_BASES", "TOTAL_BASES", "FRAC", "IND")
colnames(cov) <- headers




#data is already histogram style so need to display as bar chart perhaps
#change subset according to coverage cutoffs as ggplot2 is slow AF

newdata <- subset(cov, DP < 20)
newdata <- subset(newdata, DP > 0)
newdata %>% select(DP, FRAC, IND)

#Quick check
head(newdata)
is.data.frame(newdata)

#helpful to treat DP as numeric thing
newdata$DP <- as.numeric(as.vector(newdata$DP))

#newdata$DP <- as.factor(newdata$DP)

p <- ggplot(newdata, aes(x = DP, y = FRAC)) +
  geom_bar(stat = 'identity', colour = "#3182bd") + facet_grid(IND ~ .) 

  
p











##USEFUL QUICK PLOT FUNCTION FOR EXPLORING
#?plot

#sixtyonegee <- subset(gcov, IND = "61_G_ind")

#plot(sixtyonegee[2:15,2], sixtyonegee[2:15,5], 
#     type='h', col='darkgreen', lwd=10,
#     xlab="Depth", ylab="Fraction of genome at depth")

