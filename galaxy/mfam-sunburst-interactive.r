#!/usr/bin/env Rscript

# ---------- Preparations ----------
# Load libraries
#library(ape)
#library(data.table)
#library(plyr)
library(plotly)
library(htmlwidgets)
library(shiny)
library(sunburstR)

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
if (length(args) < 7) {
	print("Error! No or not enough arguments given.")
	print("Usage: $0 working_dir mfam_dir metfamily_dir input_tsv sunburst_method sunburst_html sunburst_csv")
	quit(save="no", status=1, runLast=FALSE)
}

# Working directory
setwd(as.character(args[1]))

# MFam directory
mfam_dir <- as.character(args[2])

# MetFamily directory
metfamily_dir <- as.character(args[3])

# input data matrix
input_tsv <- as.character(args[4])

# Parameter: sunburst_method (sum | mean | median | diff)
sunburst_method <- as.character(args[5])

# Parameter: sunburst_html (.html)
sunburst_html <- as.character(args[6])

# Parameter: sunburst_tsv (.tsv)
sunburst_csv <- as.character(args[7])



# #################### MFam Autoclassification sunburst plot ####################
# ---------- Load object ----------
class_list <- read.table(file=input_tsv, sep="\t", header=TRUE, stringsAsFactors=TRUE)



# ---------- Sunburst plot of classes ----------
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

# Convert data.frame
class_tree <- data.frame(paths=gsub(x=rownames(class_list), pattern='\\; ', replacement='-'), counts=numberOfSpectra)

# Sunburst plot
sunburst_plot <- sunburst(data=class_tree, count=TRUE)

# Save Sunburst plot
#as_widget()
saveWidget(widget=sunburst_plot, file=sunburst_html, selfcontained=TRUE)

# Save Sunburst table
write.csv(data.frame("Classifier classes"=classifierClasses, "Number of spectra"=numberOfSpectra), file=sunburst_csv, row.names=FALSE)


