library(dplyr)
library(apparent)
library(pegas)
library(adegenet)
library(stringr)
library(sequoia)
library(data.table)
library(tibble)
library(outliers)
library(rcolony)

#read loci and bamlist
bamlist <- read.delim('~/Projects/DNDI/Data/mf-unamplified/forgeno/bamlist', header = FALSE)



#create column header dataframe from text and bamlist then add as header to the loci dataframe
t<-as.data.frame(c('CHROM', 'POS'), header = 'col1')
colnames(t) <- c('col1')
colnames(bamlist) <- c('col1')
loci_header <- rbind(t, bamlist)

beagle <- read.delim('~/Projects/DNDI/Data/mf-unamplified/fordownstream/forcolony.geno', header = F)

colnames(beagle) <- loci_header$col1

#for some reason a column of nas has appeared
#not sure why but it has to die
beagle <- beagle[1:(length(beagle)-1)]

#define a function to add gt delimiters that fussy software want
sub <- function(x) {
  paste(str_sub(x, 1,1), str_sub(x, 1,1), sep = "/")
}

nsub <- function (x) {
  gsub("N/N", "-/-", x)
}

#strip non gt fields
colstokeep <- beagle %>% select(CHROM, POS)
#apply function across geno file
beagle <- beagle %>% 
  select(.,-CHROM, -POS) %>% lapply(., sub) %>% lapply(., nsub)
beagle <- as.data.frame(beagle)
beagle <- cbind(colstokeep, beagle)

#let's reduce the size of this behemoth

PLACEHOLDER <- beagle %>% 
  #REMOVE SHITTY SCAFFOLDS N SEX CHROMS N SYMBIONTS N SHIT
  filter(CHROM == 'OVOC_OM3'|CHROM == 'OVOC_OM1a'|CHROM=='OVOC_OM1b'|CHROM=='OVOC_OM5'|CHROM == 'OVOC_OM4') %>% 
  #REMOVE poorly sequenced individuals
  select(.,-A10_S10, -II_MF_4_S11, -X63C_MF_9_S9, -X63C_MF_8_S8, -X63C_MF_7_S7, -X63C_MF_6_S6, -X63C_MF_11_S22, -X63C_MF_10_S10, -X61G_MF_19_S5)

#tinker with the rowsum until you've lost as much missing data as you care to remove
s <- PLACEHOLDER %>% filter(rowSums(PLACEHOLDER == "-/-")/NCOL(.) < 0.01)

#subsample some data for the apparent analysis        
s <- s %>% 
  sample_n(., 10000) %>% 
  select(., -CHROM, -POS) %>% 
  t(.) %>% 
  as.data.frame(.) %>% 
  rownames_to_column(.)
#remove adult males
adults<-filter(s, !grepl('AF|MF', rowname))
mf<-filter(s, grepl('AF|MF', rowname))
head(s)
#create the column name
off <- data.frame(matrix("Off", nrow = 14, ncol = 1))
colnames(off) <- c('Off')
mffinal <- cbind(s$rowname, off, s[,2:10001])
#run apparent
apparent(InputFile = final, Dyad = F)

##############
#for SEQOUIA
##############

#read in geno file and strip na conotaining column
geno <- read.delim('~/Projects/DNDI/Data/mf-unamplified/forgeno/forgt.geno', strip.white = T)
geno <- geno[1:(length(geno)-1)]
colnames(geno) <- loci_header$col1

#get rid of crappily sequenced individuals
geno <- geno %>% 
  #REMOVE SHITTY SCAFFOLDS N SEX CHROMS N SYMBIONTS N SHIT
  filter(CHROM == 'OVOC_OM3'|CHROM == 'OVOC_OM1a'|CHROM=='OVOC_OM1b'|CHROM=='OVOC_OM5'|CHROM == 'OVOC_OM4') %>% 
  #REMOVE poorly sequenced individuals
  select(.,-A10_S10, -II_MF_4_S11, -'63C_MF_9_S9', -'63C_MF_8_S8', -'63C_MF_7_S7', -'63C_MF_6_S6', -'63C_MF_11_S22', -'63C_MF_10_S10', -'61G_MF_19_S5')

#check distribution of missing data
hist(rowSums(geno == "-1"))

#get rid of missing data
s <- geno %>% filter(rowSums(geno == "-1")/NCOL(.) < 0.01)

#make plink format
#subsample to ease computation times
PLACEHOLDER   <-s%>%  sample_n(., 20000)
PLACEHOLDER$genid <- paste(PLACEHOLDER$CHROM, PLACEHOLDER$POS, sep = '_')
PLACEHOLDER <- select(PLACEHOLDER, -CHROM, -POS)

PLACEHOLDER = setNames(data.frame(t(PLACEHOLDER[,-42])), PLACEHOLDER[,43])
t <- head(PLACEHOLDER,-1) %>% rownames_to_column(.) %>% cbind(PAT = 0, .) %>% cbind(MAT = 0,. ) %>% cbind(SEX = 0,. ) %>% cbind(PHENOTYPE = 0,.) %>% cbind(IID =.$rowname,.) %>% cbind(FID =.$rowname,.) %>% .[ , !(names(.) %in% "rowname")]

t <- filter(t, grepl('AF|MF', FID)) 

geno <- GenoConvert(t)


lhd <- read.csv('~/Projects/DNDI/Data/mf-unamplified/ped/sequoia/lifehistdata.csv')
LifeHistData <- filter(lhd, Sex == 3)
ParOUT <- sequoia(GenoM = geno,  LifeHistData = LifeHistData, MaxSibIter = 0, Err=0.0001, MaxMismatch=50)
rel <- GetMaybeRel(geno, LifeHistData=LifeHistData, SeqList = ParOUT, ParSib = "sib", Err=0.0001, MaxMismatch=50)


