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
# /galaxy/mfam-autoclass.r /data /usr/local/share/MFam /usr/local/share/MetFamily negative Control_NEG_70kNorm_Rep1.msp Control_NEG_70kNorm_Rep1 classes_set.txt 10 0.005 0.01 10 True 0.001 10 0.02 True 0.01 10 0.9 0.01 1 True False output_plot_class_abundance.pdf output_summary.txt output_validation.tsv output_classes_richness.tsv output_rdata.rdata
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 28) {
	print("Error! No or not enough arguments given.")
	print("Usage: $0 working_dir mfam_dir metfamily_dir polarity input_msp sample_name minimumIntensityOfMaximalMS2peak minimumProportionOfMS2peaks mzDeviationAbsolute_grouping mzDeviationInPPM_grouping doPrecursorDeisotoping mzDeviationAbsolute_precursorDeisotoping mzDeviationInPPM_precursorDeisotoping maximumRtDifference doMs2PeakGroupDeisotoping mzDeviationAbsolute_ms2PeakGroupDeisotoping mzDeviationInPPM_ms2PeakGroupDeisotoping proportionOfMatchingPeaks_ms2PeakGroupDeisotoping mzDeviationAbsolute_mapping minimumNumberOfMS2PeaksPerGroup neutralLossesPrecursorToFragments neutralLossesFragmentsToFragments output_plot_class_abundance output_summary output_validation output_classes_richness output_rdata")
	quit(save="no", status=1, runLast=FALSE)
}

# Working directory
setwd(as.character(args[1]))

# MFam directory
mfam_dir <- as.character(args[2])

# MetFamily directory
metfamily_dir <- as.character(args[3])

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

# Sample name
sample_name <- as.character(args[6])

# Input text file containing list of input compound classes that are used to classify spectra
input_classes <- as.character(args[7])

# Parameter: minimumIntensityOfMaximalMS2peak = 10 (100)
minimumIntensityOfMaximalMS2peak <- as.numeric(args[8])

# Parameter: minimumProportionOfMS2peaks = 0.005 (0.05)
minimumProportionOfMS2peaks <- as.numeric(args[9])

# Parameter: mzDeviationAbsolute_grouping = 0.1 (0.01)
mzDeviationAbsolute_grouping <- as.numeric(args[10])

# Parameter: mzDeviationInPPM_grouping = 500 (10)
mzDeviationInPPM_grouping <- as.numeric(args[11])

# Parameter: doPrecursorDeisotoping = TRUE
doPrecursorDeisotoping <- as.logical(args[12])

# Parameter: mzDeviationAbsolute_precursorDeisotoping = 0.001
mzDeviationAbsolute_precursorDeisotoping <- as.numeric(args[13])

# Parameter: mzDeviationInPPM_precursorDeisotoping = 10
mzDeviationInPPM_precursorDeisotoping <- as.numeric(args[14])

# Parameter: maximumRtDifference = 0.02
maximumRtDifference <- as.numeric(args[15])

# Parameter: doMs2PeakGroupDeisotoping = TRUE
doMs2PeakGroupDeisotoping <- as.logical(args[16])

# Parameter: mzDeviationAbsolute_ms2PeakGroupDeisotoping = 0.01
mzDeviationAbsolute_ms2PeakGroupDeisotoping <- as.numeric(args[17])

# Parameter: mzDeviationInPPM_ms2PeakGroupDeisotoping = 10
mzDeviationInPPM_ms2PeakGroupDeisotoping <- as.numeric(args[18])

# Parameter: proportionOfMatchingPeaks_ms2PeakGroupDeisotoping = 0.9
proportionOfMatchingPeaks_ms2PeakGroupDeisotoping <- as.numeric(args[19])

# Parameter: mzDeviationAbsolute_mapping = 0.01
mzDeviationAbsolute_mapping <- as.numeric(args[20])

# Parameter: minimumNumberOfMS2PeaksPerGroup = 1
minimumNumberOfMS2PeaksPerGroup <- as.numeric(args[21])

# Parameter: neutralLossesPrecursorToFragments = TRUE
neutralLossesPrecursorToFragments <- as.logical(args[22])

# Parameter: neutralLossesFragmentsToFragments = FALSE
neutralLossesFragmentsToFragments <- as.logical(args[23])

# Parameter: output_plot_class_abundance (.pdf)
output_plot_class_abundance <- as.character(args[24])

# Parameter: output_summary (.txt)
output_summary <- as.character(args[25])

# Parameter: output_validation (.tsv)
output_validation <- as.character(args[26])

# Parameter: output_classes_richness (.tsv)
output_classes_richness <- as.character(args[27])

# Parameter: output_rdata (.rdata)
output_rdata <- as.character(args[28])

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
obj <- data.frame(classes=unique(classifiers[,"Annotation (putative)"]), frequency=0)
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
cat(file=output_summary, append=FALSE, paste("Number of merged spectra:", spectra_number, "\n"))
cat(file=output_summary, append=TRUE, paste("Number of spectra classified:", spectra_classified, "\n"))
cat(file=output_summary, append=TRUE, paste("Number of unclassified spectra:", spectra_number - spectra_classified, "\n"))
cat(file=output_summary, append=TRUE, paste("Number of classes:", length(classes_order), "\n"))
cat(file=output_summary, append=TRUE, paste("Number of classes with entities:", length(classes), "\n"))
cat(file=output_summary, append=TRUE, paste("Number of classes without entities:", length(classes_order) - length(classes), "\n"))
cat(file=output_summary, append=TRUE, "\n")
cat(file=output_summary, append=TRUE, "Classes with entities: \n")
cat(file=output_summary, append=TRUE, paste(gsub(x=classes_order[which( (classes_order %in% classes))], pattern=".*; ", replacement=""), "\n"))
cat(file=output_summary, append=TRUE, "\n")
cat(file=output_summary, append=TRUE, "Classes without entities: \n")
cat(file=output_summary, append=TRUE, paste(gsub(x=classes_order[which( ! (classes_order %in% classes))], pattern=".*; ", replacement=""), "\n"))

classes <- classes_order[which(classes_order %in% classes)]



# ---------- Write validation metrics ----------
classifiers_class <- get(load(paste0(classifier_name,".RData")))

# Area under Precision Recall Curve
classifier_auc <- unlist(lapply(X=classes, FUN = function(x) { for (i in 1:length(classifiers_class)) { if (classifiers_class[[i]]$class==x) { y <- classifiers_class[[i]]$AUC; names(y) <- x; return(y) } } } ))

# True Positive Rate for False Positive Rate of 5 Percent
classifier_fpr <- unlist(lapply(X=classes, FUN = function(x) { for (i in 1:length(classifiers_class)) { if (classifiers_class[[i]]$class==x) { y <- classifiers_class[[i]]$TPR_for_FPR_of_5Percent; names(y) <- x; return(y) } } } ))

# Save table with AUC-PR and TPR-FPR rates
classifiers_validation <- data.frame(compound_class=gsub(x=classes, pattern='.*\\;', replacement=''))
classifiers_validation[which(classes %in% names(classifier_auc)), "AUC-PR"] <- classifier_auc
classifiers_validation[which(classes %in% names(classifier_fpr)), "TPR-FPR"] <- classifier_fpr
write.table(x=classifiers_validation, file=output_validation, sep="\t", quote=TRUE, row.names=FALSE, dec=".")



# ---------- Export classes richness per sample ----------
# Count numbers of matched classes in each sample
sample_name <- gsub(x=basename(input_msp), pattern='\\..*$', replacement='')
class_list <- NULL
for (i in classes) {
	vec <- as.numeric(sum(grepl(pattern=i, x=classifiers[,"Annotation (putative)"]), na.rm=TRUE))
	class_list <- rbind(class_list, vec)
}
colnames(class_list) <- sample_name
rownames(class_list) <- classes
write.table(x=class_list, file=output_classes_richness, sep="\t", quote=TRUE, row.names=TRUE, dec=".")



# ---------- Export RData ----------
save.image(file=output_rdata)



