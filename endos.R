
# libraries

library("phyloseq")
library("ggplot2")
library('ape')

setwd("./data")

# import data into R
# first .shared file or 'otu matrix'

sharedFile = read.table('micro.endos.shared')
sharedFile = t(sharedFile)
rownames(sharedFile) = sharedFile[,1]
colnames(sharedFile) = sharedFile[2,]
sharedFile = sharedFile[,2:234]
sharedFile = sharedFile[4:824,]
class(sharedFile) <- "numeric"

# now my taxonomy file

taxFile = read.table('micro.endos.taxonomy', header=T, sep='\t')
rownames(taxFile) = taxFile[,1]
taxFile = taxFile[,3:9]
taxFile = as.matrix(taxFile)

# import metadata - ie. site, type, etc..

metaFile = read.table('pools.metaData2', header=T, sep='\t')
rownames(metaFile) = metaFile[,1]
metaFile = metaFile[,2:6]

# import tre file

treeFile = read.tree(file = 'micro.endos.tre')

# use phyloseq 

OTU = otu_table(sharedFile, taxa_are_rows = TRUE)
TAX = tax_table(taxFile)
META = sample_data(metaFile)
TREE = phy_tree(treeFile)

endos = phyloseq(OTU, TAX, META, TREE)

# do some pruning

endosRelAbund = transform_sample_counts(endos, function(x) x / sum(x) )
endosRelAbundFilt = filter_taxa(endosRelAbund, function(x) mean(x) > 1e-2, TRUE)


endosRelAbundFiltSpist <- prune_samples(c("11","49-2","50-2","51-2","52-2","53-2","54-2","57-2","63-2","64-2","65-2","67-2","80-2","75-2","76-2","77-2","101","102","103","105","173","174","175","176","205","206","207","208","209","226","227","228","229","230","231","235","236","237","241","242","243","244","251","253","254","256","257","258","259","260","261","262","MIC-25","MIC-50","MIC-53","MIC-91","MIC-145","MIC-151","MIC-220","MIC-226","MIC-229","MIC-424","MIC-427","MIC-453","RS1","RS2","RS3","RS5","RS6","RS7","RS10","RS12","RS13","RS15","RS16","RS17","RS18","19-rep2","55-2-rep2","56-2-rep2","104-rep2","MIC-85-rep2"), endosRelAbundFilt)

# now plot the tree

plot_tree(endosRelAbundFiltSpist, label.tips = 'taxa_names', color='site', size='abundance')

plot_bar(endosRelAbundFilt) +
  facet_wrap(~site, scales='free')



