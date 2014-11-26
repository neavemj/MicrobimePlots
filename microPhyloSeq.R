
## phyloseq on mothur output for global STP 2.11.14

# install requirements


source("http://bioconductor.org/biocLite.R")
biocLite("phyloseq")

library("devtools")
install_github("phyloseq", "joey711")

install.packages("ggplot2")
install.packages("plyr")
install.packages("vegan")
install.packages("directlabels")
install.packages("knitr")
install.packages("clustsig")

# install latest verion to get bug fixes

library("devtools")
install_github("phyloseq", "joey711")

# load libraries

library("phyloseq")
library("ggplot2")
library("plyr")

# move into directory

setwd('~/microbiome/subprojects/2.16S_MiSeq/2.data/3.mothur/')

# import data into R
# first .shared file or 'otu matrix'

sharedFile = read.table('micro.final.shared')
sharedFile = t(sharedFile)
rownames(sharedFile) = sharedFile[,1]
colnames(sharedFile) = sharedFile[2,]
sharedFile = sharedFile[,2:234]
sharedFile = sharedFile[4:37368,]

class(sharedFile) <- "numeric"

# also import subsampled (7779 seqs) shared file

sharedsubFile = read.table('micro.final.0.03.subsample.shared')
sharedsubFile = t(sharedsubFile)
rownames(sharedsubFile) = sharedsubFile[,1]
colnames(sharedsubFile) = sharedsubFile[2,]
sharedsubFile = sharedsubFile[,2:219]
sharedsubFile = sharedsubFile[4:14762,]

class(sharedsubFile) <- "numeric"

# now my taxonomy file

taxFile = read.table('micro.final.0.03.cons.taxonomy', header=T, sep='\t')
rownames(taxFile) = taxFile[,1]
taxFile = taxFile[,3:9]
taxFile = as.matrix(taxFile)

# import metadata - ie. site, type, etc..

metaFile = read.table('pools.metaData2', header=T, sep='\t')
rownames(metaFile) = metaFile[,1]
metaFile = metaFile[,2:6]

# use phyloseq 

OTU = otu_table(sharedFile, taxa_are_rows = TRUE)
OTUsub = otu_table(sharedsubFile, taxa_are_rows = TRUE)
TAX = tax_table(taxFile)
META = sample_data(metaFile)

physeq = phyloseq(OTU, TAX, META)
physeqSub = phyloseq(OTUsub, TAX, META)

physeq
physeqSub

sample_variables(physeq)
sample_variables(physeqSub)

# do some richness plots on untrimmed data to get acurate estimates

microRich = plot_richness(micro, x = 'species', measures = c('Chao1', 'Shannon', 'Observed'), color = 'site')
microRich + geom_boxplot(data = microRich$data, aes(x = species, y = value, color = NULL), 
                         alpha = 0.1)

plot_richness(microSub, x = 'species', measures = c('Chao1', 'Shannon', 'observed'), color = 'site')

# get rid of any OTUs not present in any samples and get relative abundance

micro <- prune_taxa(taxa_sums(physeq) > 0, physeq)
microSub <- prune_taxa(taxa_sums(physeqSub) > 0, physeqSub)
microSubRel = transform_sample_counts(microSub, function(x) x / sum(x) )
microSubRelFilt = filter_taxa(microSubRel, function(x) mean(x) > 1e-5, TRUE)

# split into species 

microSubRelFiltSpist <- subset_samples(microSubRelFilt, species=='Stylophora pistillata')
microSubRelFiltPdami <- subset_samples(microSubRelFilt, species=='Pocillopora damicornis')
microSubRelFiltPverr <- subset_samples(microSubRelFilt, species=='Pocillopora verrucosa')
microSubRelFiltSea <- subset_samples(microSubRelFilt, species=='seawater')

# some merging to get reasonable graph..

sitesMerged = merge_samples(physeq, "site")
speciesMerged = merge_samples(physeq, "species")

plot_bar(physeq, x="reef", fill="Family")

# try top 50 OTUs in spist

microSubRelFiltSpistTop50names <- names(sort(taxa_sums(microSubRelFiltSpist), TRUE)[1:20])
microSubRelFiltSpistTop50 <- prune_taxa(microSubRelFiltSpistTop50names, microSubRelFiltSpist)
plot_bar(microSubRelFiltSpistTop50, fill="Genus")

# now split into species and seawater etc
# create relative abundance values

top20relAbund = transform_sample_counts(top20otus, function(x) x / sum(x) )
plot_bar(top20relAbund, x="reef", fill="Genus")

top20seawaterTaxa = subset_samples(top20otus, species=='seawater')

# relative abundance of spist otus bar plot

microSubRelFiltSpistFilt = filter_taxa(microSubRelFiltSpist, function(x) mean(x) > 1e-2, TRUE)
plot_bar(microSubRelFiltSpistFilt, fill="Genus", title='Stylophora pistillata')

# relative abundance of pocillopora verrucosa otus bar plot

microSubRelFiltPverrFilt = filter_taxa(microSubRelFiltPverr, function(x) mean(x) > 1e-2, TRUE)
plot_bar(microSubRelFiltPverrFilt, fill="Genus", title='Pocillopora verrucosa')

# relative abundance of pocillopora damicornis otus bar plot

microSubRelFiltPdamiFilt = filter_taxa(microSubRelFiltPdami, function(x) mean(x) > 1e-2, TRUE)
plot_bar(microSubRelFiltPdamiFilt, fill="Genus", title='Pocillopora damicornis')

# relative abundance of seawater otus bar plot

microSubRelFiltSeaFilt = filter_taxa(microSubRelFiltSea, function(x) mean(x) > 1e-2, TRUE)
plot_bar(microSubRelFiltSeaFilt, fill="Genus")

# just spist endozoicomonas OTUs - must put identity in brackets afwaterward (needs exact match)

microSubRelFiltSpistEndo = subset_taxa(microSubRelFiltSpist, Genus=='Endozoicomonas(100)')
plot_bar(microSubRelFiltSpistEndo, fill="Genus", title='Stylophora pistillata')

microSubSpist <- subset_samples(microSub, species=='Stylophora pistillata')
microSubSpistSitesMerged = merge_samples(microSubSpist, "site")
microSubSpistSitesMergedRel = transform_sample_counts(microSubSpistSitesMerged, function(x) x / sum(x) )
microSubSpistSitesMergedRelEndo = subset_taxa(microSubSpistSitesMergedRel, Genus=='Endozoicomonas(100)')
plot_bar(microSubSpistSitesMergedRelEndo, fill="Genus", title='Stylophora pistillata')

# s.pist from aquaria - NOTE: have to do this on non-subsampled data as these get removed in this step

spistAq = prune_samples(c("SPaq", "SPaq-rep2"), physeq)
spistAqrelAbundFilt = filter_taxa(spistAqrelAbund, function(x) mean(x) > 1e-2, TRUE)
plot_bar(spistAqrelAbundFilt, fill="Genus")

# negative controls

neg = prune_samples(c("40-cycles"), physeq)
negrelAbund = transform_sample_counts(neg, function(x) x / sum(x) )
negrelAbundFilt = filter_taxa(negrelAbund, function(x) mean(x) > 1e-2, TRUE)
plot_bar(negrelAbundFilt, fill="Genus")

# compare samples that were replicated

replicatedSamples = prune_samples(c('RP7-rep2', 'RP8-rep2', 'RS5-rep2', 'END74-rep2', 'END59-rep2', 'END62-rep2', '55-2-rep2', '56-2-rep2', '103-rep2', '105-rep2', '250-rep2', '253-rep2', '254-rep2', 'MIC-53-rep2', 'MIC-85-rep2', 'RP7', 'RP8', 'RS5', 'END74', 'END59', 'END62', '55-2', '56-2', '103', '105', '250', '253', '254', 'MIC-53', 'MIC-85'), microSubRelFilt)

replicatedSamplesFilt = filter_taxa(replicatedSamples, function(x) mean(x) > 4e-3, TRUE)
plot_bar(replicatedSamplesFilt, fill="Genus", title='Replicated Samples')

# compare ningaloo samples that were warmed v. PFA preserved

ningalooSamples = prune_samples(c('50', '50-2', '52', '52-2', '53', '53-2'), microSubRelFilt)
ningalooSamplesFilt = filter_taxa(ningalooSamples, function(x) mean(x) > 3e-3, TRUE)
plot_bar(ningalooSamplesFilt, fill="Genus", title='Ningaloo Samples')

# ordination spist

spistOrd <- ordinate(microSubRelFiltSpist, "PCoA", "bray")
spistP1 = plot_ordination(microSubRelFiltSpist, spistOrd, type = "taxa", color = "Phylum", title = "taxa")
print(spistP1)

# ordination by samples

spistP2 = plot_ordination(microSubRelFiltSpist, spistOrd, type = "samples", color = "site")
print(spistP2)
spistP1 + facet_wrap(~Phylum, 4)
spistP2 + geom_polygon(aes(fill = site), alpha=0.5) + geom_point(size = 4) + ggtitle("samples")

# biplot containing both samples and otus..

spistP3 = plot_ordination(microSubRelFiltSpist, spistOrd, type = "biplot", color = "site", shape = "Genus", title = "biplot")

# Some stuff to modify the automatic shape scale
microSubRelFiltSpist.shape.names = get_taxa_unique(microSubRelFiltSpist, "Phylum")
microSubRelFiltSpist.shape <- 15:(15 + length(microSubRelFiltSpist.shape.names) - 1)
names(microSubRelFiltSpist.shape) <- microSubRelFiltSpist.shape.names
microSubRelFiltSpist.shape["samples"] <- 16
spistP3 + scale_shape_manual(values = microSubRelFiltSpist.shape)
show(spistP3)

# ordination poc

p.verrAll = subset_samples(physeq, species=='Pocillopora verrucosa')
p.verrAllrelAbund = transform_sample_counts(p.verrAll, function(x) x / sum(x) )
p.verrAllrelAbundFilt = filter_taxa(p.verrAllrelAbund, function(x) mean(x) > 1e-5, TRUE)
pverrOrd <- ordinate(p.verrAllrelAbundFilt, "NMDS", "bray")
pverrP1 = plot_ordination(p.verrAllrelAbundFilt, pverrOrd, type = "taxa", color = "Genus", title = "taxa")
print(pverrP1)

# ordination by samples

pverrP2 = plot_ordination(p.verrAllrelAbundFilt, pverrOrd, type = "samples", color = "site")
print(pverrP2)
pverrP2 + geom_polygon(aes(fill = site)) + geom_point(size = 4) + ggtitle("samples")


physeqSubOrd <- ordinate(physeqSub, "NMDS", "bray")
physeqP1 = plot_ordination(physeqSub, physeqSubOrd, type = "samples", color = "species")
print(physeqP1)
physeqP1 + geom_polygon(aes(fill = site)) + geom_point(size = 4) + ggtitle("samples")










