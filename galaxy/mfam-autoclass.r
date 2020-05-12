#!/usr/bin/env Rscript

# ---------- Preparations ----------
# Load libraries
library(MSnbase)
library(xcms)

# Setup R error handling to go to stderr
options(show.error.messages=F, error=function() { cat(geterrmessage(), file=stderr()); q("no",1,F) } )

# Set locales and encoding
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")
loc <- Sys.setlocale(category="LC_ALL", locale="C")
options(encoding="UTF-8")

# Set options
options(stringAsfactors=FALSE, useFancyQuotes=FALSE)



# ---------- Arguments and user variables ----------
# Take in trailing command line arguments
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 26) {
	print("Error! No or not enough arguments given.")
	print("Usage: $0 working_dir mfam_dir metfamily_dir polarity input_msp minimumIntensityOfMaximalMS2peak minimumProportionOfMS2peaks mzDeviationAbsolute_grouping mzDeviationInPPM_grouping doPrecursorDeisotoping mzDeviationAbsolute_precursorDeisotoping mzDeviationInPPM_precursorDeisotoping maximumRtDifference doMs2PeakGroupDeisotoping mzDeviationAbsolute_ms2PeakGroupDeisotoping mzDeviationInPPM_ms2PeakGroupDeisotoping proportionOfMatchingPeaks_ms2PeakGroupDeisotoping mzDeviationAbsolute_mapping minimumNumberOfMS2PeaksPerGroup neutralLossesPrecursorToFragments neutralLossesFragmentsToFragments output_plot_class_abundance output_summary output_validation output_rdata")
	quit(save="no", status=1, runLast=FALSE)
}

# Working directory
setwd(dirname(args[1]))

# MFam directory
mfam_dir <- dirname(args[2])

# MFam directory
metfamily_dir <- dirname(args[3])

# Polarity
polarity <- as.character(args[4])

# MFam Classifier depending on polarity
if (polarity == "positive") {
	classifier_name <- paste0(mfam_dir, "/", "2018-05-29_18_28_56_2018-02-13_09_14_10_LC-MS-MS_pos_21908")
} else {
	classifier_name <- paste0(mfam_dir, "/", "2018-05-29_18_21_05_2018-02-13_09_14_10_LC-MS-MS_neg_11328")
}

# Input MSP containging MS2 spectra
input_msp <- as.character(args[5])

# Input text file containing list of input compound classes that are used to classify spectra
input_classes <- as.character(args[6])

# Parameter: minimumIntensityOfMaximalMS2peak = 10 (100)
minimumIntensityOfMaximalMS2peak <- as.numeric(args[7])

# Parameter: minimumProportionOfMS2peaks = 0.005 (0.05)
minimumProportionOfMS2peaks <- as.numeric(args[8])

# Parameter: mzDeviationAbsolute_grouping = 0.1 (0.01)
mzDeviationAbsolute_grouping <- as.numeric(args[9])

# Parameter: mzDeviationInPPM_grouping = 500 (10)
mzDeviationInPPM_grouping <- as.numeric(args[10])

# Parameter: doPrecursorDeisotoping = TRUE
doPrecursorDeisotoping <- as.logical(args[11])

# Parameter: mzDeviationAbsolute_precursorDeisotoping = 0.001
mzDeviationAbsolute_precursorDeisotoping <- as.numeric(args[12])

# Parameter: mzDeviationInPPM_precursorDeisotoping = 10
mzDeviationInPPM_precursorDeisotoping <- as.numeric(args[13])

# Parameter: maximumRtDifference = 0.02
maximumRtDifference <- as.numeric(args[14])

# Parameter: doMs2PeakGroupDeisotoping = TRUE
doMs2PeakGroupDeisotoping <- as.logical(args[15])

# Parameter: mzDeviationAbsolute_ms2PeakGroupDeisotoping = 0.01
mzDeviationAbsolute_ms2PeakGroupDeisotoping	<- as.numeric(args[16])

# Parameter: mzDeviationInPPM_ms2PeakGroupDeisotoping = 10
mzDeviationInPPM_ms2PeakGroupDeisotoping <- as.numeric(args[17])

# Parameter: proportionOfMatchingPeaks_ms2PeakGroupDeisotoping = 0.9
proportionOfMatchingPeaks_ms2PeakGroupDeisotoping <- as.numeric(args[18])

# Parameter: mzDeviationAbsolute_mapping = 0.01
mzDeviationAbsolute_mapping <- as.numeric(args[19])

# Parameter: minimumNumberOfMS2PeaksPerGroup = 1
minimumNumberOfMS2PeaksPerGroup <- as.numeric(args[20])

# Parameter: neutralLossesPrecursorToFragments = TRUE
neutralLossesPrecursorToFragments <- as.logical(args[21])

# Parameter: neutralLossesFragmentsToFragments = FALSE
neutralLossesFragmentsToFragments <- as.logical(args[22])

# Parameter: output_plot_class_abundance (.pdf)
output_plot_class_abundance <- as.character(args[23])

# Parameter: output_summary (.txt)
output_summary <- as.character(args[24])

# Parameter: output_validation (.tsv)
output_validation <- as.character(args[25])

# Parameter: output_rdata (.rdata)
output_rdata <- as.character(args[26])

# Load MetFamily functions
source(paste0(metfamily_dir, "/", "R_packages.R"))
source(paste0(metfamily_dir, "/", "FragmentMatrixFunctions.R"))
source(paste0(metfamily_dir, "/", "DataProcessing.R"))
source(paste0(metfamily_dir, "/", "Annotation.R"))
source(paste0(metfamily_dir, "/", "Classifiers.R"))



# ---------- Import spectra from MSP ----------
importMs1Ms2data <- function(fileMs1Path, fileMs2Path, minimumIntensityOfMaximalMS2peak, minimumProportionOfMS2peaks, mzDeviationAbsolute_grouping, mzDeviationInPPM_grouping){
	# Box parameters
	parameterSet <- list()
	parameterSet$minimumIntensityOfMaximalMS2peak <- minimumIntensityOfMaximalMS2peak
	parameterSet$minimumProportionOfMS2peaks <- minimumProportionOfMS2peaks
	parameterSet$mzDeviationAbsolute_grouping <- mzDeviationAbsolute_grouping
	parameterSet$mzDeviationInPPM_grouping <- mzDeviationInPPM_grouping
	parameterSet$doPrecursorDeisotoping <- doPrecursorDeisotoping
	parameterSet$mzDeviationAbsolute_precursorDeisotoping <- mzDeviationAbsolute_precursorDeisotoping
	parameterSet$mzDeviationInPPM_precursorDeisotoping <- mzDeviationInPPM_precursorDeisotoping
	parameterSet$maximumRtDifference <- maximumRtDifference
	parameterSet$doMs2PeakGroupDeisotoping <- doMs2PeakGroupDeisotoping
	parameterSet$mzDeviationAbsolute_ms2PeakGroupDeisotoping <- mzDeviationAbsolute_ms2PeakGroupDeisotoping
	parameterSet$mzDeviationInPPM_ms2PeakGroupDeisotoping <- mzDeviationInPPM_ms2PeakGroupDeisotoping
	parameterSet$proportionOfMatchingPeaks_ms2PeakGroupDeisotoping <- proportionOfMatchingPeaks_ms2PeakGroupDeisotoping
	parameterSet$mzDeviationAbsolute_mapping <- mzDeviationAbsolute_mapping
	parameterSet$minimumNumberOfMS2PeaksPerGroup <- minimumNumberOfMS2PeaksPerGroup
	parameterSet$neutralLossesPrecursorToFragments <- neutralLossesPrecursorToFragments
	parameterSet$neutralLossesFragmentsToFragments <- neutralLossesFragmentsToFragments
	
	error <- NULL
	resultObj <- tryCatch(
		{
			convertToProjectFile(
				filePeakMatrix = fileMs1Path, 
				fileSpectra = fileMs2Path, 
				parameterSet = parameterSet, 
				progress = NA
			)
		}, error = function(e) {
			error <<- e
		}
	)
	
	# Errors?
	if(!is.null(error)){
		stop(paste(
			"There occurred an error while processing the input files. Please check the file format and content and try again.", "\n",
			"Occurred error: ", error, sep = ""
		))
		return()
	}
	if(length(resultObj) == 1){
		if(resultObj == "Number of spectra is zero"){
			stop(paste("There are no MS/MS spectra which fulfill the given criteria. Please refine parameter 'Spectrum intensity' and try 'Import MS\u00B9 and MS/MS data' again."))
			return()
		}
	}
	
	lines <- sparseMatrixToString(matrixRows = resultObj$matrixRows, matrixCols = resultObj$matrixCols, matrixVals = resultObj$matrixVals, parameterSet = parameterSet)
	
	# Process project file
	error <- NULL
	dataList <- tryCatch(
		{
			readProjectData(fileLines = lines, progress = FALSE)
		}, error = function(e) {
			error <<- e
		}
	)
	
	if(!is.null(error)){
		stop(paste(
			"There occurred an error while processing the project file. Please check the file format and content and try again.", "\n",
			"Occurred error: ", error, sep = ""
		))
		return()
	}
	
	return(dataList)
}



# ---------- Apply classifier on sample spectra ----------
applyClassifierMs2 <- function(classifierFile, propertiesFile, fileMs1Path = NULL, fileMs2Path, fileClasses, minimumIntensityOfMaximalMS2peak = 2000, minimumProportionOfMS2peaks = 0.05, mzDeviationAbsolute_grouping = 0.01, mzDeviationInPPM_grouping = 10) {
	# Classes of interest
	lines <- readLines(con = fileClasses)
	theseLines <- grepl(pattern = "^Organic compounds", x = lines)
	classesCanonical <- lines[which(theseLines)-1]
	classes <- lines[theseLines]
	classesCanonical <- gsub(x = classesCanonical, pattern = ":", replacement = "")
	classes <- gsub(x = classes, pattern = "/", replacement = "; ")
	
	# Process MS-MS data (msp)
	dataList <- importMs1Ms2data(fileMs1Path, fileMs2Path, minimumIntensityOfMaximalMS2peak, minimumProportionOfMS2peaks, mzDeviationAbsolute_grouping, mzDeviationInPPM_grouping)
	
	propertiesList <- getClassifierProperties(propertiesFile)
	
	# Load and apply classifier
	error <- NULL
	resultObj <- tryCatch(
		{
			doAnnotation(filePath = classifierFile, propertiesList = propertiesList, featureMatrix = dataList$featureMatrix, parameterSet = dataList$importParameterSet, classesWhiteList = classes, progress = FALSE)
		}, error = function(e) {
			error <<- e
		}
	)
	
	if(!is.null(error)){
		stop(paste(
			"There occurred an error while applying the classifiers. Please check the file format and content and try again.", "\n",
			"Occurred error: ", error, sep = ""
		))
		return()
	}
	
	classToSpectra_class	<- resultObj$classToSpectra_class
	properties_class		<- resultObj$properties_class
	mappingSpectraToClassDf <- resultObj$mappingSpectraToClassDf
	
	# Box classifier results to data.frame
	annotationRows <- list()
	for(classIdx in seq_along(classToSpectra_class)){
		
		precursorIndeces <- as.integer(names(classToSpectra_class[[classIdx]]))
		pValues		  <- unname(classToSpectra_class[[classIdx]])
		
		class <- names(classToSpectra_class)[[classIdx]]
		
		for(idx in seq_along(precursorIndeces)){
			precursorIndex	 <- precursorIndeces[[idx]]
			precursorLabel	 <- dataList$precursorLabels[[precursorIndex]]
			mz				 <- dataList$dataFrameInfos[[precursorIndex, "m/z"]]
			rt				 <- dataList$dataFrameInfos[[precursorIndex, "RT"]]
			metaboliteName	 <- dataList$dataFrameInfos[[precursorIndex, "Metabolite name"]]
			
			presentAnnotations <- unlist(dataList$annoArrayOfLists[[precursorIndex]])
			if(length(presentAnnotations) == 0){
				presentAnnotations <- ""
			} else {
				presentAnnotations <- sort(unlist(presentAnnotations))
				presentAnnotations <- paste(presentAnnotations, collapse = "; ")
			}
			
			pValue <- pValues[[idx]]
			
			annotationRows[[length(annotationRows)+1]] <- c(
				"Index" = precursorIndex, 
				"Label" = precursorLabel, 
				"m/z"   = mz, 
				"RT"	= rt, 
				"Metabolite name" = metaboliteName, 
				"Annotation (present)" = presentAnnotations, 
				"Annotation (putative)" = class, 
				"pValue" = pValue
			)
		}
	}
	
	# Box
	head <- c(
		"Index", 
		"Label", 
		"m/z", 
		"RT", 
		"Metabolite name", 
		"Annotation (present)", 
		"Annotation (putative)", 
		"pValue"
	)
	
	if (is.null(unlist(annotationRows))==TRUE) {
		annotationDf <- data.frame(matrix("", ncol=length(head), nrow=5))
		colnames(annotationDf) <- head
	} else {
		annotationDf <- as.data.frame(t(matrix(data = unlist(annotationRows), nrow = length(head))))
		colnames(annotationDf) <- head
	}
	
	return(annotationDf)
}



# #################### Automated in-silico classification ####################
# ---------- Apply MFam classifier ----------
classifiers <- list()

classifiers <- applyClassifierMs2(classifierFile = paste0(classifier_name, ".RData"),
								  propertiesFile = paste0(classifier_name, ".txt"),
								  fileMs1Path = NULL,
								  fileMs2Path = input_msp,
								  fileClasses = input_classes,
								  minimumIntensityOfMaximalMS2peak = minimumIntensityOfMaximalMS2peak,
								  minimumProportionOfMS2peaks = minimumProportionOfMS2peaks,
								  mzDeviationAbsolute_grouping = mzDeviationAbsolute_grouping,
								  mzDeviationInPPM_grouping = mzDeviationInPPM_grouping)

# Diversity of compound classes
div_classes <- data.frame()
obj <- data.frame(mode=i, classes=unique(classifiers[,"Annotation (putative)"]), frequency=0)
for (j in 1:nrow(obj)) obj[j,"frequency"] <- length(which(obj$classes[j] == classifiers[,"Annotation (putative)"]))
div_classes <- rbind(div_classes, obj)

# Plot most abundant classes
pdf(file=output_plot_class_abundance, encoding="ISOLatin1", pointsize=10, width=16, height=10, family="Helvetica")
par(mfrow=c(1,1), mar=c(15,4,4,1), oma=c(0,0,0,0), cex.axis=0.8, cex=0.9)
barplot(div_classes[1:nrow(div_classes),"frequency"], names.arg=gsub('.*\\;','',div_classes[1:nrow(div_classes),"classes"]), las=3, ylab="frequency", main="Most abundant compound classes")
dev.off()



# ---------- Determine how many spectra were classified ----------
spectra_number <- sum(div_classes$frequency)
spectra_classified <- sum(length(unique(classifiers$`Metabolite name`)))

classes_order <- readLines(con=input_classes)
classes_order <- classes_order[which(grepl(pattern="^Organic compounds", x=classes_order))]
classes_order <- gsub(x=classes_order, pattern=":", replacement="")
classes_order <- gsub(x=classes_order, pattern="/", replacement="; ")

classes <- div_classes$classes
classes <- classes[which(grepl(pattern="^Organic compounds", x=classes))]
classes <- gsub(x=classes, pattern=":", replacement="")
classes <- gsub(x=classes, pattern="/", replacement="; ")



# ---------- Write classification summary ----------
output_summary_conn <- file(output_summary)
writeLines(con=output_summary_conn, paste("Number of merged spectra:", spectra_number))
writeLines(con=output_summary_conn, paste("Number of spectra classified:", spectra_classified))
writeLines(con=output_summary_conn, paste("Number of unclassified spectra:", spectra_number - spectra_classified))
writeLines(con=output_summary_conn, paste("Number of classes:", length(classes_order)))
writeLines(con=output_summary_conn, paste("Number of classes with entities:", length(classes)))
writeLines(con=output_summary_conn, paste("Number of classes without entities:", length(classes_order) - length(classes)))
writeLines(con=output_summary_conn, "")
writeLines(con=output_summary_conn, "Classes with entities:")
writeLines(con=output_summary_conn, gsub(x=classes_order[which( (classes_order %in% classes))], pattern=".*; ", replacement=""))
writeLines(con=output_summary_conn, "")
writeLines(con=output_summary_conn, "Classes without entities:")
writeLines(con=output_summary_conn, gsub(x=classes_order[which( ! (classes_order %in% classes))], pattern=".*; ", replacement=""))
close(output_summary_conn)

classes <- classes_order[which(classes_order %in% classes)]



# ---------- Write validation metrics ----------
classifiers_class <- get(load(paste0(classifier_name,".RData")))

# Area under Precision Recall Curve
classifier_auc <- unlist(lapply(X=classes, FUN = function(x) { for (i in 1:length(classifiers_class)) { if (classifiers_class[[i]]$class==x) { y <- classifiers_class[[i]]$AUC; names(y) <- x; return(y) } } } ))

# True Positive Rate for False Positive Rate of 5 Percent
classifier_fpr <- unlist(lapply(X=classes, FUN = function(x) { for (i in 1:length(classifiers_class)) { if (classifiers_class[[i]]$class==x) { y <- classifiers_class[[i]]$TPR_for_FPR_of_5Percent; names(y) <- x; return(y) } } } ))

# Save table with AUC-PR and TPR-FPR rates
classifiers_validation <- data.frame(compound_class=gsub(x=classes, pattern='.*\\;', replacement=''))
classifiers_validation[which(classes %in% names(classifier_pos_auc)), "AUC-PR"] <- classifier_auc
classifiers_validation[which(classes %in% names(classifier_pos_fpr)), "TPR-FPR"] <- classifier_fpr
write.table(x=classifiers_validation, file=output_validation, sep="\t", quote=TRUE, row.names=FALSE, dec=".")



# ---------- Export RData ----------
save.image(file=output_rdata)

#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################



# Count numbers of matched classes in each sample
class_list <- data.frame()
for (i in classes) {
	vec <- as.numeric(c(vec, sum(grepl(pattern=i, x=classifiers[,"Annotation (putative)"]), na.rm=TRUE)))
	class_list <- cbind(class_list, vec)
}
colnames(class_list) <- c("mode", gsub(x=classes, pattern=".*\\; ", replacement=""))
write.table(x=t(class_list), file="classifiers_classes.csv", sep=";", quote=TRUE, row.names=TRUE, dec=".")










# Sunburst plot of classes
source("classifier/sunburst.r")
numberOfSpectra <- as.numeric(unlist(apply(X=class_list[,2:ncol(class_list)], MARGIN=2, FUN = function(x) { sum(x) })))
numberOfSpectra <- numberOfSpectra[which(numberOfSpectra > 0)]
classifierClasses <- classes[which(numberOfSpectra > 0)]
pdf(file="classification_sunburst.pdf", encoding="ISOLatin1", pointsize=8, width=10, height=10, family="Helvetica")
par(mfrow=c(1,1), mar=c(0,0,0,0), oma=c(0,0,0,0), cex.axis=1, cex=1)
sunBurstPlotFromSubstanceClasses(classifierClasses, numberOfSpectra, colorStart=0.0, colorAlpha=0.6)
dev.off()
write.csv(data.frame("Classifier classes"=classifierClasses, "Number of spectra"=numberOfSpectra), file="classification_sunburst.csv", row.names=FALSE)

# Heatmap
heatmap.2(x=as.matrix(class_list[,2:ncol(class_list)]), cexRow=1, cexCol=0.7,
		  Rowv=as.dendrogram(hclust(dist(class_list[,2:ncol(class_list)]))), offsetRow=0,
		  Colv=rev(as.dendrogram(hclust(dist(t(class_list[,2:ncol(class_list)]))))), offsetCol=0,
		  col=colorRampPalette(c('lightgrey','white','darkblue'), alpha=0.1, bias=8)(256),
		  trace="none", margins=c(15,6),
		  key=TRUE, key.title="Color key")

# Treemap
#install.packages('treemap')
library(treemap)

classtree <- data.frame(classname=classes, kingdom=NA, superclass=NA, class=NA, subclass=NA, level5=NA, level6=NA, funcclass=NA, counts=0)
rownames(classtree) <- gsub(x=classes, pattern='.*\\;', replacement='')
for (i in 1:length(classes)) {
	classlist <- strsplit(classes, '\\; ')
	classtree[which(classes==classes[i]), c(2:(1+length(classlist[[i]])))] <- classlist[[i]]
}

class_count <- NULL
for (i in classes) {
	vec <- 0
	for (j in c("pos","neg")) {
		vec <- as.numeric(vec + sum(grepl(pattern=i, x=classifiers[[j]][,"Annotation (putative)"]), na.rm=TRUE))
	}
	class_count <- c(class_count, vec)
}
classtree$counts <- class_count

pdf(file="classification_treemap.pdf", encoding="ISOLatin1", pointsize=8, width=12, height=6, family="Helvetica")
treemap(dtf=classtree, index=c("superclass","class","subclass","level5","level6"),
		vSize="counts", type="index", title="", border.col="white", bg.labels=0, lowerbound.cex.labels=0,
		palette.HCL.options=list(hue_start=0, hue_end=300, luminance=60),
		inflate.labels=FALSE, fontsize.labels=c(16,14,12,10,8), fontcolor.labels="white",
		ymod.labels=c(0.2,0.1,0,-0.1,-0.2),
		force.print.labels=TRUE, overlap.labels=1,
		title.legend="", position.legend="right", fontsize.legend=8)
dev.off()

# Convert data.frame to tree
library(ape)
library(data.tree)
library(plyr)

class_tree <- lapply(strsplit(classifierClasses, "; "), function(x) as.data.frame(t(x)))
class_tree <- rbind.fill(class_tree)
class_tree$pathString <- apply(class_tree, 1, function(x) paste(trimws(na.omit(x)), collapse="/"))
class_tree <- chronos(as.phylo(data.tree::as.Node(class_tree, name="Classes")))

class_nodes <- gsub(x=class_tree$node.label, pattern="_", replacement=" ")
class_tips <- gsub(x=class_tree$tip.label, pattern="_", replacement=" ")
class_nodes_col <- sample(rainbow(length(class_nodes)))
class_nodes_col <- sapply(class_nodes, function(x) { x <- class_nodes_col[grep(x=class_nodes, pattern=x)] } )
class_tips_col <- sample(rainbow(length(class_tips)))
class_tips_col <- sapply(class_tips, function(x) { x <- class_tips_col[grep(x=class_tips, pattern=x)] } )
class_tips_col <- c("#00E8FFFF", "#00E8FFFF", "#00E8FFFF", "#00E8FFFF", "#00E8FFFF", "#FFE800FF", "#FFE800FF", "#FFE800FF")
#, "#5D00FFFF", "#00FFE8FF", "#A200FFFF", "#FF002EFF", "#7400FFFF", "#FF005DFF", "#FFB900FF", "#00FF17FF")

# Normal tree
pdf(file="classification_tree.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=12, family="Helvetica")
plot(class_tree, type="phylogram", use.edge.length=TRUE, root.edge=TRUE, lab4ut="axial", edge.width=0.4,
	 show.node.label=TRUE, node.pos=2,
	 show.tip.label=TRUE, tip.color="black",
	 label.offset=0.01, no.margin=TRUE, cex=1)
dev.off()

# Circular tree
pdf(file="classification_circular_tree.pdf", encoding="ISOLatin1", pointsize=10, width=12, height=12, family="Helvetica")
plot(class_tree, type="radial", use.edge.length=TRUE, root.edge=F, lab4ut="axial", edge.width=1.0,
	 show.node.label=TRUE, node.pos=2, #edge.color=na.omit(as.character(unlist(sapply(funcclass_nodes, function(x) { x <- funcclass_nodes_col[grep(x=funcclass, pattern=x)] } )))),
	 show.tip.label=TRUE, tip.color="black",
	 label.offset=0.01, no.margin=TRUE, cex=1)
dev.off()


