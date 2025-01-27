#!/usr/bin/env Rscript

# ---------- Preparations ----------
# Load libraries
library(circlize)
library(plotrix)

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
	print("Usage: $0 working_dir mfam_dir metfamily_dir input_tsv sunburst_method sunburst_pdf sunburst_csv")
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

# Parameter: sunburst_pdf (.pdf)
sunburst_pdf <- as.character(args[6])

# Parameter: sunburst_tsv (.tsv)
sunburst_csv <- as.character(args[7])



# ---------- sunBurstPlotFromSubstanceClasses ----------
sunBurstPlotFromSubstanceClasses <- function(classifierClasses, numberOfSpectra, colorStart = 0.2, colorAlpha = 0.5){
	level <- unlist(lapply(X = strsplit(x = classifierClasses, split = "; "), FUN = length))
	
	## all class levels
	classesAndSubClasses <- lapply(X = strsplit(x = classifierClasses, split = "; "), FUN = function(x){
		sapply(X = seq_along(x), FUN = function(y){paste(x[1:y], collapse = "; ")})
	})
	classesByLevel <- list()
	labelsByLevel <- list()
	for(levelHere in seq_len(max(level))){
		classesByLevel[[levelHere]] <- sort(unique(unlist(lapply(X = classesAndSubClasses, FUN = function(y){
			if(length(y) < levelHere) return(NULL)
			else return(y[[levelHere]])
		}))))
		labelsByLevel[[levelHere]] <- unlist(lapply(X = strsplit(x = classesByLevel[[levelHere]], split = "; "), FUN = tail, n=1))
	}
	
	## class counts
	countsByLevel <- list()
	for(levelHere in rev(seq_len(max(level)))){
		countsByLevel[[levelHere]] <- unlist(lapply(X = classesByLevel[[levelHere]], FUN = function(class){
			newSpectra <- ifelse(test = class %in% classifierClasses, yes = numberOfSpectra[[which(class == classifierClasses)]], no = 0)
			oldSpectra <- ifelse(test = levelHere < max(level), yes = sum(countsByLevel[[levelHere+1]][grepl(x = classesByLevel[[levelHere+1]], pattern = paste("^", class, sep = ""))]), no = 0)
			return(newSpectra + oldSpectra)
		}))
	}
	rootCount <- sum(countsByLevel[[1]])
	
	## coordinates
	colors <- rainbow(start = colorStart, alpha = colorAlpha, n = 1000)
	startDegreeByLevel <- list()
	spanDegreeByLevel <- list()
	colorByLevel <- list()
	
	for(levelHere in seq_len(max(level))){
		startDegreeByLevel[[levelHere]] <- list()
		spanDegreeByLevel[[levelHere]] <- list()
		colorByLevel[[levelHere]] <- list()
		
		classesToProcess <- classesByLevel[[levelHere]]
		precursorClasses <- NULL
		if(levelHere == 1)  precursorClasses <- ""
		else				precursorClasses <- classesByLevel[[levelHere-1]]
		
		for(precursorClassIdx in seq_along(precursorClasses)){
			precursorClass <- precursorClasses[[precursorClassIdx]]
			classesToProcessHereSelection <- grepl(x = classesToProcess, pattern = precursorClass)
			classesToProcessHere <- classesToProcess[classesToProcessHereSelection]
			startDegree <- ifelse(test = levelHere == 1, yes = 0, no = startDegreeByLevel[[levelHere-1]][[precursorClassIdx]])
			scalingFactor <- ifelse(test = levelHere == 1, yes = 1, no = countsByLevel[[levelHere-1]][[precursorClassIdx]] / sum(countsByLevel[[levelHere]][classesToProcessHereSelection])) ## ambiguous classes
			#startColor  <- ifelse(test = levelHere == 1, yes = 0, no = startDegreeByLevel[[levelHere-1]][[precursorClassIdx]])
			for(classToProcessHere in classesToProcessHere){
				classIdx <- which(classesByLevel[[levelHere]] == classToProcessHere)
				degreeSpan <- 360 * countsByLevel[[levelHere]][[classIdx]] / rootCount * scalingFactor
				startDegreeByLevel[[levelHere]][[classIdx]] <- startDegree
				spanDegreeByLevel [[levelHere]][[classIdx]] <- degreeSpan
				colorByLevel	  [[levelHere]][[classIdx]] <- colors[[(floor(startDegree + degreeSpan / 2) / 360 * length(colors)) + ifelse(test = (floor(startDegree + degreeSpan / 2) / 360 * length(colors))==length(colors), yes = 0, no = 1) ]]
				startDegree <- startDegree + degreeSpan
			}
		}
	}
	
	thereIsNextLevelByLevel <- list()
	for(levelHere in seq_len(max(level))){
		thereIsNextLevelByLevel[[levelHere]] <- list()
		if(levelHere == max(level)){
			thereIsNextLevelByLevel[[levelHere]] <- rep(x = FALSE, times = length(classesByLevel[[levelHere]]))
		} else {
			for(classIdx in seq_along(classesByLevel[[levelHere]]))
				thereIsNextLevelByLevel[[levelHere]][[classIdx]] <- any(grepl(x = classesByLevel[[levelHere+1]], pattern = classesByLevel[[levelHere]][[classIdx]]))
		}
	}
	
	degreeThresholdForDrawing <- 0.5
	#maxLevel <- max(which(unlist(lapply(X = spanDegreeByLevel, FUN = function(x){any(unlist(x) >= degreeThresholdForDrawing)}))))
	
	plotRadius <- 8
	plotCex1 <- 1
	plotCex2 <- 1
	
	plot(1, type="n", xlab="", ylab="", xlim=c(-plotRadius, plotRadius), ylim=c(-plotRadius, plotRadius), axes = FALSE)
	
	## circle segments
	tmp <- sapply(X = seq_along(classesByLevel), FUN = function(levelHere){
		sapply(X = seq_along(classesByLevel[[levelHere]]), FUN = function(classIdx){
			if(spanDegreeByLevel[[levelHere]][[classIdx]] < degreeThresholdForDrawing)  return()
			draw.sector(
				start.degree = startDegreeByLevel[[levelHere]][[classIdx]],
				end.degree = startDegreeByLevel[[levelHere]][[classIdx]] + spanDegreeByLevel[[levelHere]][[classIdx]],
				rou1 = levelHere - 1, rou2 = levelHere, center = c(0,0), clock.wise = FALSE, col = colorByLevel [[levelHere]][[classIdx]], border = "white"
			)
		})
	})
	## segment text
	minimumAngleToShowSegmentText <- 15
	tmp <- sapply(X = seq_along(classesByLevel), FUN = function(levelHere){
		sapply(X = seq_along(classesByLevel[[levelHere]]), FUN = function(classIdx){
			if(spanDegreeByLevel[[levelHere]][[classIdx]] < minimumAngleToShowSegmentText)  return()
			textTokens <- strwrap(x = labelsByLevel[[levelHere]][[classIdx]], width = max(nchar(strsplit(x = "Some text", split = " ")[[1]])))
			#firstOffset <- 1 / (length(textTokens) * 2)
			
			for(idx in seq_along(textTokens)){
				#offset <- firstOffset * (2 * idx - 1)
				middle <- (startDegreeByLevel[[levelHere]][[classIdx]] + spanDegreeByLevel[[levelHere]][[classIdx]]/2)/360 * 2*pi
				isSwitched <- middle > pi/2 & middle < 3 * pi/2
				isSwitched <- middle > pi
				offset <-  ifelse(test = middle > pi, yes = (length(textTokens) - idx + 1) / (length(textTokens) + 1), no = idx / (length(textTokens) + 1))
				#if(isSwitched) textTokens[[idx]] <- rev(textTokens[[idx]])
				#label <- ifelse(test = isSwitched, yes = paste(strsplit(x = labelsByLevel[[levelHere]][[classIdx]], split = " ")[[1]], collapse = " "), no = labelsByLevel[[levelHere]][[classIdx]])
				#middle <- ifelse(test = middle < pi, yes = middle + pi, no = middle)
				arctext(
					x = textTokens[[idx]],
					center = c(0, 0), radius = levelHere - offset - 0.04,
					middle = middle,
					cex = plotCex1, stretch = 1,
					clockwise = !isSwitched
				)
			}
		})
	})
	
	## outer text
	levelMaxHere <- max(level)
	tmp <- sapply(X = seq_along(classesByLevel), FUN = function(levelHere){
		if(levelHere > levelMaxHere)  return()
		sapply(X = seq_along(classesByLevel[[levelHere]]), FUN = function(classIdx){
			if(thereIsNextLevelByLevel[[levelHere]][[classIdx]] | spanDegreeByLevel[[levelHere]][[classIdx]] >= minimumAngleToShowSegmentText)  return()
			#radius <- maxLevel + 0.2
			radius <- levelHere + 0.2
			angle <- (startDegreeByLevel[[levelHere]][[classIdx]] + spanDegreeByLevel[[levelHere]][[classIdx]]/2)/360 * 2*pi
			x <- radius * sin(angle+pi/2)
			y <- radius * cos(angle+pi/2)
			srt <- startDegreeByLevel[[levelHere]][[classIdx]] + spanDegreeByLevel[[levelHere]][[classIdx]]/2
			isSwitched <- srt > 90 & srt < 270
			if(isSwitched) adj <- c(1,0) else adj <- c(0,0.5)
			srt <- ifelse(test = isSwitched, yes = srt + 180, no = srt)
			text(
				x = x, y = -y, labels = labelsByLevel[[levelHere]][[classIdx]], adj = adj,
				srt=srt,
				cex = plotCex2
			)
		})
	})
}



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

# Save Sunburst plot
pdf(file=sunburst_pdf, encoding="ISOLatin1", pointsize=8, width=10, height=10, family="Helvetica")
par(mfrow=c(1,1), mar=c(0,0,0,0), oma=c(0,0,0,0), cex.axis=1, cex=1)
sunBurstPlotFromSubstanceClasses(classifierClasses, numberOfSpectra, colorStart=0.0, colorAlpha=0.6)
dev.off()

# Save Sunburst table
write.csv(data.frame("Classifier classes"=classifierClasses, "Number of spectra"=numberOfSpectra), file=sunburst_csv, row.names=FALSE)


