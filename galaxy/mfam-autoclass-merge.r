#!/usr/bin/env Rscript

# ---------- Preparations ----------
# Load libraries
#library(MSnbase)
#library(xcms)

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
if (length(args) < 9) {
	print("Error! No or not enough arguments given.")
	print("Usage: $0 working_dir mfam_dir metfamily_dir input_objects output_plot_class_abundance output_summary output_classes_richness output_validation output_rdata")
	quit(save="no", status=1, runLast=FALSE)
}

# Working directory
setwd(as.character(args[1]))

# MFam directory
mfam_dir <- as.character(args[2])

# MetFamily directory
metfamily_dir <- as.character(args[3])

# File with input RData object names
input_objects <- as.character(args[4])

# Parameter: output_plot_class_abundance (.pdf)
output_plot_class_abundance <- as.character(args[5])

# Parameter: output_summary (.txt)
output_summary <- as.character(args[6])

# Parameter: output_classes_richness (.tsv)
output_classes_richness <- as.character(args[7])

# Parameter: output_validation (.tsv)
output_validation <- as.character(args[8])

# Parameter: output_rdata (.rdata)
output_rdata <- as.character(args[9])



# #################### Merge MFam classification objects ####################
# ---------- Load objects ----------
# Read text file with input rdata objects
object_files <- readLines(input_objects)
if (object_files[length(object_files)] == "") object_files <- object_files[-length(object_files)]

# Load RData objects and sample names
rdata_objects <- list()
for (i in 1:length(object_files)) {
	rdata_objects[[i]] <- new.env()
	load(file=object_files[i], envir=rdata_objects[[i]])
}

# Read text file with sample names
sample_names <- list()
for (i in 1:length(object_files)) {
	sample_names[[i]] <- rdata_objects[[i]]$sample_name
}
sample_name <- sample_names

# Polarity
polarity <- list()
classifier_name <- list()
for (i in 1:length(object_files)) {
	polarity[[i]] <- rdata_objects[[i]]$polarity
	
	# MFam Classifier depending on polarity
	if (polarity[[i]] == "positive") {
		classifier_name[[i]] <- paste0(mfam_dir, "/", "2018-05-29_18_28_56_2018-02-13_09_14_10_LC-MS-MS_pos_21908")
	} else {
		classifier_name[[i]] <- paste0(mfam_dir, "/", "2018-05-29_18_21_05_2018-02-13_09_14_10_LC-MS-MS_neg_11328")
	}
}

# Set attribute that this is a merged MFam object
MFam_Merged_RData <- TRUE



# ---------- Plot most abundant classes ----------
# Diversity of compound classes
div_classes <- rdata_objects[[1]]$div_classes
for (i in 2:length(object_files)) {
	div_classes <- merge(div_classes, rdata_objects[[i]]$div_classes, by="classes", all.x=TRUE, all.y=TRUE)
}
colnames(div_classes) <- c("classes", paste0("frequency.",sample_names))
div_classes[is.na(div_classes)] <- 0
div_classes$frequency <- apply(X=div_classes[,c(2:ncol(div_classes))], MARGIN=1, FUN=function(x) { sum(x) })

# Plot most abundant classes
pdf(file=output_plot_class_abundance, encoding="ISOLatin1", pointsize=10, width=16, height=10, family="Helvetica")
par(mfrow=c(1,1), mar=c(15,4,4,1), oma=c(0,0,0,0), cex.axis=0.8, cex=0.9)
barplot(div_classes[1:nrow(div_classes),"frequency"], names.arg=gsub('.*\\;','',div_classes[1:nrow(div_classes),"classes"]), las=3, ylab="frequency", main="Most abundant compound classes")
dev.off()



# ----------  Join classes richness per samples ----------
class_list <- rdata_objects[[1]]$class_list
for (i in 2:length(object_files)) {
	class_list <- merge(class_list, rdata_objects[[i]]$class_list, by="row.names", all.x=TRUE, all.y=TRUE)
}
rownames(class_list) <- class_list$Row.names
class_list$Row.names <- NULL
write.table(x=class_list, file=output_classes_richness, sep="\t", quote=TRUE, row.names=TRUE, dec=".")



# ---------- Write classification summary ----------
# Classes names
classes_order <- rdata_objects[[1]]$classes_order
for (i in 2:length(object_files)) {
	if (length(rdata_objects[[i]]$classes_order) != length(classes_order)) {
		print("Error! Set of reference classes do not match! Please only use one reference classes.txt for merging.")
		quit(status=3)
	}
}
classes <- rownames(class_list)

# Determine how many spectra were classified
spectra_number <- sum(div_classes$frequency)
spectra_classified <- 0
for (i in 1:length(object_files)) spectra_classified <- spectra_classified + sum(length(unique(rdata_objects[[i]]$classifiers$`Metabolite name`)))

# Write classification summary
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
classifiers_validation <- data.frame(compound_class=gsub(x=classes, pattern='.*\\;', replacement=''))

classifier_auc <- list()
classifier_fpr <- list()
for (i in 1:length(unique(polarity))) {
	if (unique(polarity)[i] == "positive") {
		tmp_classifier_name <- paste0(mfam_dir, "/", "2018-05-29_18_28_56_2018-02-13_09_14_10_LC-MS-MS_pos_21908")
	} else {
		tmp_classifier_name <- paste0(mfam_dir, "/", "2018-05-29_18_21_05_2018-02-13_09_14_10_LC-MS-MS_neg_11328")
	}
	tmp_classifiers_class <- get(load(paste0(tmp_classifier_name,".RData")))
	
	# Area under Precision Recall Curve
	classifier_auc[[i]] <- unlist(lapply(X=classes, FUN = function(x) { for (i in 1:length(tmp_classifiers_class)) { if 	(tmp_classifiers_class[[i]]$class==x) { y <- tmp_classifiers_class[[i]]$AUC; names(y) <- x; return(y) } } } ))
	
	# True Positive Rate for False Positive Rate of 5 Percent
	classifier_fpr[[i]] <- unlist(lapply(X=classes, FUN = function(x) { for (i in 1:length(tmp_classifiers_class)) { if (tmp_classifiers_class[[i]]$class==x) { y <- tmp_classifiers_class[[i]]$TPR_for_FPR_of_5Percent; names(y) <- x; return(y) } } } ))
	
	# Write table with AUC-PR and TPR-FPR rates
	if (unique(polarity)[i] == "positive") {
		classifiers_validation[which(classes %in% names(classifier_auc[[i]])), "AUC-PR_positive"] <- classifier_auc[[i]]
		classifiers_validation[which(classes %in% names(classifier_fpr[[i]])), "TPR-FPR_positive"] <- classifier_fpr[[i]]
	} else {
		classifiers_validation[which(classes %in% names(classifier_auc[[i]])), "AUC-PR_negative"] <- classifier_auc[[i]]
		classifiers_validation[which(classes %in% names(classifier_fpr[[i]])), "TPR-FPR_negative"] <- classifier_fpr[[i]]
	}
}

# Save table with AUC-PR and TPR-FPR rates
write.table(x=classifiers_validation, file=output_validation, sep="\t", quote=TRUE, row.names=FALSE, dec=".")



# ---------- Export RData ----------
save.image(file=output_rdata)


