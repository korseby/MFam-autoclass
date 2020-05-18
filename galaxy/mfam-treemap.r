#!/usr/bin/env Rscript

# ---------- Preparations ----------
# Load libraries
library(treemap)

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
if (length(args) < 6) {
	print("Error! No or not enough arguments given.")
	print("Usage: $0 working_dir mfam_dir metfamily_dir input_tsv treemap_method treemap_pdf")
	quit(save="no", status=1, runLast=FALSE)
}

# Working directory
setwd(as.character(args[1]))

# MFam directory
mfam_dir <- as.character(args[2])

# MetFamily directory
metfamily_dir <- as.character(args[3])

# input matrix
input_tsv <- as.character(args[4])

# Parameter: treemap_method (sum | mean | median)
treemap_method <- as.character(args[5])

# Parameter: sunburst_pdf (.pdf)
treemap_pdf <- as.character(args[6])



# #################### MFam Autoclassification treemap plot ####################
# ---------- Load object ----------
class_list <- read.table(file=input_tsv, sep="\t", header=TRUE, stringsAsFactors=TRUE)



# ---------- Treemap plot of classes ----------
# Calculate data matrix
class_list[is.na(class_list)] <- 0
if (sunburst_method == "sum") {
	numberOfSpectra <- as.numeric(unlist(apply(X=class_list, MARGIN=1, FUN = function(x) { sum(x) })))
} else if (sunburst_method == "mean") {
	numberOfSpectra <- as.numeric(unlist(apply(X=class_list, MARGIN=1, FUN = function(x) { mean(x) })))
} else if (sunburst_method == "median") {
	numberOfSpectra <- as.numeric(unlist(apply(X=class_list, MARGIN=1, FUN = function(x) { median(x) })))
} else if (sunburst_method == "diff") {
	numberOfSpectra <- as.numeric(class_list[,1])
	for (i in c(2:ncol(class_list))) {
		numberOfSpectra <- numberOfSpectra - as.numeric(class_list[,i])
	}
	numberOfSpectra <- abs(numberOfSpectra)
}
classifierClasses <- rownames(class_list)

# Treemap data frame
classtree <- data.frame(classname=classifierClasses, kingdom=NA, superclass=NA, class=NA, subclass=NA, level5=NA, level6=NA, counts=0)
rownames(classtree) <- gsub(x=classifierClasses, pattern='.*\\;', replacement='')
for (i in 1:length(classifierClasses)) {
	classlist <- strsplit(classifierClasses, '\\; ')
	classtree[which(classifierClasses==classifierClasses[i]), c(2:(1+length(classlist[[i]])))] <- classlist[[i]]
}
classtree$counts <- numberOfSpectra

# Plot treemap
pdf(file=treemap_pdf, encoding="ISOLatin1", pointsize=8, width=12, height=6, family="Helvetica")
treemap(dtf=classtree, index=c("superclass","class","subclass","level5","level6"),
		vSize="counts", type="index", title="", border.col="white", bg.labels=0, lowerbound.cex.labels=0,
		palette.HCL.options=list(hue_start=0, hue_end=360, luminance=60),
		inflate.labels=FALSE, fontsize.labels=c(16,14,12,10,8), fontcolor.labels="white",
		ymod.labels=c(0.2,0.1,0,-0.1,-0.2),
		force.print.labels=TRUE, overlap.labels=1)#,
		#title.legend="", position.legend="right", fontsize.legend=8)
dev.off()


