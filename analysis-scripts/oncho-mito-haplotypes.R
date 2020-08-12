library(pegas)
library(ape)
library(adegenet)

input <- c('~/Projects/DNDI/Data/mitochondrial/fasta/10missingness/mito_missingunder10.fna')
metadata <- read.csv('~/Projects/DNDI/Data/metadata/metadata-all.csv')

numetadata <- metadata[metadata$sample_name %in% rownames(d) ,]


d <- ape::read.dna(input, format='fasta')
e <- dist.dna(d)
h <- pegas::haplotype(d)
h <- sort(h, what = "label")
(net <- pegas::haploNet(h))
ind.hap<-with(
  stack(setNames(attr(h, "index"), rownames(h))),
  table(hap=ind, pop=rownames(d)[values])
)
plot(net, size=attr(net, "freq"), scale.ratio=0.2, pie=ind.hap)

legend(-8, 0, colnames(ind.hap), col=rainbow(ncol(ind.hap)), pch=19, ncol=2, cex =0.5, xjust =2, yjust =0.5)

?plot
