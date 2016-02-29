#! /opt/az/local/R/R-3.2.1/installdir/bin/Rscript

######################################################################
# Objective: Program to analyze RNA-seq data at the gene level using DESeq2
# Written by: Austin Dulak
# Date Written: 160218
# Input: Rscript dsloop_v1.R [counts file] [gene tpm file] [targets file] [config file] spotfire-make [output directory]
# Code example: Rscript /gpfs/users/krxw569/bin/dsloop/dsloop_v1.R annotated_combined_merged.counts combined_merged.gene.sf.tpm targets.txt config.txt spotfire-make /gpfs/users/krxw569/bin/dsloop/example/example_from_loop/
######################################################################

################
# A. Start-up
################

print("dsloop_v1.R is running")
print("Checking options")

# Command line arguments into R
options(warn=-1)
args = commandArgs(trailingOnly = TRUE)

# Check if input arguments are correct
if (length(args) < 4) {
	stop('Requires 6 trailing options: Counts file, FPKM file, targets file, config file, make spotfire call, and output dir')
} else {
	args = args
	print("Options are OK")
}

################
# B. Storing input variables
################

counts = as.character(args[1])
norm_dat_file = as.character(args[2])
targets = as.character(args[3])
config = as.character(args[4])
spotfireCall = as.character(args[5])
outputDir = as.character(args[6])

################
# C. Create output directory
################

# Generate directory
if (file.exists(outputDir)){
	outputDir = outputDir
} else {
	dir.create(outputDir)
}
print(paste('Output stored in ', outputDir, sep=""))

################
# 1. Setwd and load dependencies
################

# Load libraries
suppressPackageStartupMessages(library(DESeq2)); library(DESeq2)
suppressPackageStartupMessages(library(stringr)); library(stringr)
suppressPackageStartupMessages(library(gplots)); library(gplots)
suppressPackageStartupMessages(library(reshape2)); library(reshape2)
suppressPackageStartupMessages(library(dplyr)); library(dplyr)

# Load downstream functions
# Correlation plot
myCorPlot = function(dat, mags=c(10,10), ...){
  
  # Correlation
  corDat = cor(dat, method = 'pearson')
  
  # Heatmap2
  heatmap.2(corDat,
            margins = mags, main='Correlation plot',
            Rowv=T, Colv=T, dendrogram='none', scale='none',
            col=colorRampPalette(c('navy', 'white'))(20),
            trace='none',
            key=T, keysize=1, key.title='Scale', key.xlab='Correlation(r)', key.ylab='',
            denscol=NA, cexRow = 0.6, cexCol=0.6)
  
  return(corDat)
}

# PCA plotting and analysis
myPca = function(dat){
	
	# Apply principal components function
	pc = prcomp(dat, scale = FALSE, retx = TRUE)
	
	# Scree plot
	scree = summary(pc)$importance[2,]
	par(las=1, mfrow=c(1,2))
	barplot(scree, names.arg=names(scree), col='cadetblue', border=F, main= 'Scree plot', ylab='Proportion of variance', xlab='Component', width=0.1)

	# Scores plot
	plot(pc$x, xlab=paste('PC1 (', round(scree[1]*100, 1), '%)', sep=''),
		 ylab=paste('PC2 (', round(scree[2]*100, 1), '%)', sep=''),
		 xlim = c(min(pc$x[,1]-80), max(pc$x[,1]+80)),
		 ylim = c(min(pc$x[,2]-80), max(pc$x[,2]+80)), pch='', main = 'Principal components scores')
	text(pc$x, labels= rownames(pc$x), cex=0.7)
	
	# Return PC data
	return(pc)
}

# Clustering of a dataset; ensure genes and/or variables are in column name
myCluster = function(dat, cut = NULL){
 
 # Ward Hierarchical Clustering
 d = dist(dat, method = "euclidean") # distance matrix
 fit = hclust(d, method="ward.D") 
 
 # Plot
 plot(fit, xlab='Clustering: Euclidean distannce = Ward.D', cex=0.6)
 
 if (length(cut)==0){
  cut = NULL # Do nothing
 } else {
  groups = cutree(fit, k=cut) #cut tree into clusters
  rect.hclust(fit, k=cut, border="red")
 }
}

myS2n = function(dat, cls){
 
 # Calculate S2N, mean, and sd (using cls and annotated 1/2 vector)
 s2n = apply(dat, 1, function(x){
  meanZero = mean(as.numeric(x[names(cls[cls==1])]))
  meanOne = mean(as.numeric(x[names(cls[cls==2])]))
  sdZero =  sd(as.numeric(x[names(cls[cls==1])]))       
  sdOne =  sd(as.numeric(x[names(cls[cls==2])]))
  
  # Check and convert means for each group to >=1
  if (meanZero < 1.0){
	tempMeanZero = meanZero + 1
  } else {
	tempMeanZero = meanZero
  }	
  
  if (meanOne < 1.0){
	tempMeanOne = meanOne + 1
  } else {
	tempMeanOne = meanOne
  }	
  
  # Check whether SD > 0.2; if not then covert using abs(mean)*0.2
  if (sdZero < 0.2){
	tempSdZero = abs(tempMeanZero)*0.2
  } else {
	tempSdZero = sdZero
  }
  
  if (sdOne < 0.2){
	tempSdOne = abs(tempMeanOne)*0.2
  } else {
	tempSdOne = sdOne
  }
  
  # Correct for small sd if does not match requirements
  s2n = (meanOne-meanZero) / (tempSdOne+tempSdZero)
  s2n = c(s2n, meanZero, sdZero, meanOne, sdOne)
  return(s2n)
 })
 
 # Process for output
 rownames(s2n) = c('s2n', 'meanZero', 'sdZero', 'meanOne', 'sdOne')
 s2n = t(s2n)
 s2n = s2n[order(s2n[,1]),]
 
 # Return output
 return(s2n)
}

################
# 2. Load data into memory
################

# Targets
targets = read.table(targets, header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=FALSE)
targets$unique_id = apply(targets[,c(3:ncol(targets))], 1, paste, collapse="_") # Generate unique_id
targets$unique_id = gsub(' ', '', targets$unique_id)

# Config
compare = read.table(config, header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=FALSE)

# Counts - generate gene identifier from flat file
rawdata = read.delim(counts, check.names=FALSE, stringsAsFactors=FALSE, row.names=1)
geneNames = data.frame(id=rownames(rawdata), geneSymbol = rawdata$symbol, geneID = paste(rawdata$symbol, rownames(rawdata), sep='_'))
rawdata = rawdata[, -which(colnames(rawdata)=='symbol')]

# FPKM - log2 transform
normdata = read.delim(norm_dat_file,check.names=FALSE, stringsAsFactors=FALSE, row.names=1)
normdata = round(log2((normdata+1)),3)

################
# 3. Convert the CRO sample ids to unique ids generated in the program
################

# Check if CRO names match
cro_match = table(c(colnames(rawdata), colnames(normdata), targets$cro_sample_id))
if (all(cro_match > 2)){
	print("All CRO sample names are correct")
	cro_match = cro_match
} else {
	stop("CRO sample IDs do not match in counts, fpkm, and target files")
}

# Modify colname names
targets = targets[order(targets$cro_sample_id),]
rawdata.sort = sort(colnames(rawdata))
rawdata = rawdata[,rawdata.sort]
normdata.sort = sort(colnames(normdata))
normdata = normdata[,normdata.sort]

# Verify all names are appropriately modified
df = data.frame(targets$cro_sample_id, colnames(rawdata), colnames(normdata))
df$identical = ifelse(df[,1] == df[,2] & df[,2] == df[,3], 'yes', 'no')
if (length(df$identical == 'yes') == nrow(df)){
	colnames(rawdata) = targets$unique_id
	colnames(normdata) = targets$unique_id
	print("All CRO sample names were matched and changed")
} else {
	stop("CRO sample names did not match")
}

###################
# 4. Generate QC plots
###################

# Remove genes with TPM = 0 and with log2 max - min < 1
normdata_qc = normdata[rowSums(normdata) > 0,]
high_variability_count = apply(normdata, 1, function(x){ max(x) - min(x)})
normdata_qc = normdata_qc[names(high_variability_count[high_variability_count >= 1]),]


pdf(paste(outputDir, "/qc_qa_figures.pdf", sep=""), width=11.5, height=11.5)
	# Distribution boxplots
	par(las=2, mar=c(18.1,4.1,4.1,2.1))
	boxplot(normdata_qc,
		ylab='Gene expression (log2TPM)', main = 'Distribution of Gene Expression Values',
		col=topo.colors(ncol(normdata_qc)), cex.axis=0.6)
		
	# Correlation matrix
	par(las=1, mar=c(9.1,4.1,4.1,2.1))
	temp = myCorPlot(normdata_qc, mags=c(18, 18))

	# Cluster plot
	par(mar=c(6.1,4.1,4.1,2.1))
	temp = myCluster(t(normdata_qc))

	# PCA plot
	temp = myPca(t(normdata_qc))
	
	# Raw Counts
	par(mfrow=c(1,1),las=2, mar=c(18.1,4.1,4.1,2.1))
	barplot(colSums(rawdata)/1000000,
		ylab='Library size (10^6)', main = 'Library Sizes',
		col=topo.colors(ncol(rawdata)), cex.names=0.6)
	
invisible(dev.off())

print("Quality report has been generated")

###################
# 5. Start loop for DESeq2
###################

for (i in 1:nrow(compare)){

	# Create directory for testing
	dir.create(paste(outputDir, "/", compare[i,1], "_", compare[i,2], "/", sep=""))
	setwd(paste(outputDir, "/", compare[i,1], "_", compare[i,2], "/", sep=""))
	
	# Print comparison
	print(paste(compare[i,1], "vs.", compare[i,2], "is starting", sep=" "))
	
	##########################
	# Generate classes for deseq2 testing
	##########################
	
	# Subset data for group 1
	identifiers = unlist(strsplit(compare[i,1], split='_'))
	group.names = colnames(rawdata)
	for(j in 1:length(identifiers)){
		names.loop = grep(identifiers[j], group.names)
		group.names = intersect(group.names[names.loop], group.names)
	}
	rawdata.test = rawdata[,group.names]
	normdata.test = normdata[,group.names]
	
	# Store group information
	group.1 = data.frame(sample = group.names, group = compare[i,1])
	
	### Subset data for group 2
	identifiers = unlist(strsplit(compare[i,2], split='_'))
	group.names = colnames(rawdata)
	for(j in 1:length(identifiers)){
		names.loop = grep(identifiers[j], group.names)
		group.names = intersect(group.names[names.loop], group.names)
	}
	rawdata.test = cbind(rawdata.test, rawdata[,group.names])
	normdata.test = cbind(normdata.test, normdata[,group.names])
	
	# Store group information
	group.2 = data.frame(sample = group.names, group = compare[i,2])
	
	### Combine group information
	group = rbind(group.1, group.2)
	
	##########################
	# Generate classes for deseq2 testing
	##########################
	
	# Make individual targets for for testing
	targets2 = data.frame(condition = group$group, row.names=group$sample)
	
	# Generate DESeq2 matrix and analyze
	dds = DESeqDataSetFromMatrix(rawdata.test, DataFrame(targets2), ~ condition)
	dds = DESeq(dds, quiet = T)
	res = results(dds)
	resd = as.data.frame(res)
	resd = resd[complete.cases(resd),]
	
	##########################
	# Generate plots for diagnostic output - volcano and MA
	##########################

	pdf("significance_plots.pdf", width=11.5, height=6.5)	
	# MA plot
	par(mfrow=c(1,2))
	plotMA(res, main=paste0("DESeq2 - ",compare[i,1], "_", compare[i,2]), ylim=c(-3,3))
	
	# Color vector for volcano plot
	colVec = rep('gray80', nrow(resd))
	colVec[resd$padj <= 0.05 & resd$log2FoldChange <= -1] = 'blue'
	colVec[resd$padj <= 0.05 & resd$log2FoldChange >= 1] = 'red'
	plot(resd$log2FoldChange, -log10(resd$padj),
		xlab = 'Log2 fold-change', ylab = '-log10(adjusted p-value)', 
		main = paste0("Volcano - ",compare[i,1], "_", compare[i,2]),
		col=colVec, pch=19)
	invisible(dev.off())

	##########################
	# Generate S2N data for the comparison and write out statistics
	##########################
	
	# Create CLS vector
	cls = as.numeric(group$group)
	names(cls) = group$sample
	
	# Run S2N function
	s2Output = myS2n(dat=normdata.test, cls=cls)
	
	# Annotate column names from s2n
	colnames(s2Output)[2:5] = c(paste('log2TPM',compare[i,1], sep='_'),paste('TPM_sd',compare[i,1], sep='_'),
	paste('log2TPM',compare[i,2], sep='_'), paste('TPM_sd',compare[i,2], sep='_'))
	
	# Merge S2N output with DESeq results
	complete_output = merge(resd, s2Output, by=0)
	
	# Add gene information
	complete_output = merge(geneNames, complete_output, by=1)
	
	# Modify column names
	colnames(complete_output)[1] = "EnsemblGeneID"
	
	# Write table for output
	write.table(complete_output, "deseq2_statistics.txt", na="", col.name=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
	
	# Write preranked files for GSEA input
	dir.create("gsea_input")
	gsea = complete_output[,c('geneSymbol', 's2n')]
	colnames(gsea) = paste("#", colnames(gsea), sep="")
	gsea = gsea[complete.cases(gsea),]
	write.table(gsea, paste(getwd(), "/gsea_input/gsea_s2n_TPM.rnk", sep=""), na="", col.name=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
	
	##########################
	# Generate heatmaps based on S2N
	##########################
	
	# Plot top 50 (25/each direction) by S2N from TPM data
	dir.create("heatmaps")
	
	# Generate ordered dataframe
	myOrd = complete_output[order(complete_output$s2n, decreasing=F),]
	
	# Generate sample list
	mySamp = as.character(group$sample)
	
	# Make deHeatmap
	geneList = c(rev(as.character(head(myOrd$EnsemblGeneID, 25))), as.character(tail(myOrd$EnsemblGeneID, 25)))
	nameList = c(rev(as.character(head(myOrd$geneID, 25))), as.character(tail(myOrd$geneID, 25)))
	
	hm = normdata.test[geneList, mySamp]
	rownames(hm) = nameList
	
	pdf(paste(getwd(),"/heatmaps/heatmap_s2n.pdf", sep=""), height=12, width=12)
	heatmap.2(as.matrix(hm),
		Rowv=FALSE, Colv=FALSE, dendrogram="none",
		margins=c(20,20), cexRow=1, cexCol=1,
		col=colorRampPalette(c("darkblue","white","red"))(50), scale='row',
		key=TRUE, symkey=FALSE,
		density.info="none", trace="none", keysize=1,
		main = 'Heatmap of 50 most \nsignificantly modulated genes')
	invisible(dev.off())

	print(paste(compare[i,1], "vs.", compare[i,2], "edgeR analysis is complete", sep=" "))
}

###################
# 6. Generate tall spotfire file if applicable from command line prompt
###################

if (spotfireCall == 'spotfire-make'){
	print("Generating tall file output for spotfire")
	
	# Generate tall file for plotting in Spotfire - data is log2(tpm + 1)
	spot = normdata
	
	spot = merge(geneNames, spot, by.x=1, by.y=0, all.x=T)
	spot = melt(spot, id.vars = c('id','geneID', 'geneSymbol'))
	colnames(spot)[c(1,4,5)] = c('EnsemblGeneID', 'UniqueSampleID', 'log2(TPM+1)')
	spot = merge(targets, spot, by.x = 'unique_id', by.y = 'UniqueSampleID', all.y = TRUE)
	write.table(spot, paste(outputDir, "/spotfire_tall_file.txt", sep=""), 
		na="", col.name=NA,row.names=TRUE,quote=FALSE,sep="\t")
	print("DESeq2 is complete")
} else {
	stop("DESeq2 is complete")
}