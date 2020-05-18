#!/usr/bin/env Rscript

# ---------- Preparations ----------
# Load libraries
library(ape)
library(data.tree)
library(plyr)

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
	print("Usage: $0 working_dir mfam_dir metfamily_dir input_tsv circtree_method circtree_pdf")
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

# Parameter: circtree_method (sum | mean | median)
circtree_method <- as.character(args[5])

# Parameter: circtree_pdf (.pdf)
circtree_pdf <- as.character(args[6])



# #################### MFam Autoclassification circtree plot ####################
# ---------- Load object ----------
class_list <- read.table(file=input_tsv, sep="\t", header=TRUE, stringsAsFactors=TRUE)



# ---------- Circular tree plot of classes ----------
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

# Convert data.frame to tree
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

# Circular tree
pdf(file=circtree_pdf, encoding="ISOLatin1", pointsize=10, width=12, height=12, family="Helvetica")
plot(class_tree, type="radial", use.edge.length=TRUE, root.edge=F, lab4ut="axial", edge.width=1.0,
	 show.node.label=TRUE, node.pos=2, #edge.color=na.omit(as.character(unlist(sapply(funcclass_nodes, function(x) { x <- funcclass_nodes_col[grep(x=funcclass, pattern=x)] } )))),
	 show.tip.label=TRUE, tip.color="black",
	 label.offset=0.01, no.margin=TRUE, cex=1)
dev.off()


