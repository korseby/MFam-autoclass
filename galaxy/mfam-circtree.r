# ############################## Classification with MetFamily ##############################
#install.packages(c('colourpicker','shinyBS','FactoMineR','slam','cba','squash','plotrix','plotly','circlize'))
source("classifier/ms2_classification.r")
classifiers <- list()

#classifier_name <- "classifier/2019-06-16_MSMS_HR_pos_free__c_13314__s_71354__mFam_MoNA_GNPS_WeizMASS_PlaSMA_Classifier"
classifier_name <- "classifier/2018-05-29_18_28_56_2018-02-13_09_14_10_LC-MS-MS_pos_21908"
classifiers[["pos"]] <- applyClassifierMs2(classifierFile = paste0(classifier_name, ".RData"),
										   propertiesFile = paste0(classifier_name, ".txt"),
										   fileMs1Path = NULL,
										   fileMs2Path = "__ms2_library_spectra_matching_pos.msp",
										   fileClasses = "classifier/classes.txt",
										   minimumIntensityOfMaximalMS2peak = 10, #100
										   minimumProportionOfMS2peaks = 0.005, #0.05
										   mzDeviationAbsolute_grouping = 0.1, #0.01
										   mzDeviationInPPM_grouping = 500 ) #10
classifier_name <- "classifier/2018-05-29_18_21_05_2018-02-13_09_14_10_LC-MS-MS_neg_11328"
classifiers[["neg"]] <- applyClassifierMs2(classifierFile = paste0(classifier_name, ".RData"),
										   propertiesFile = paste0(classifier_name, ".txt"),
										   fileMs1Path = NULL,
										   fileMs2Path = "__ms2_library_spectra_matching_neg.msp",
										   #fileMs2Path = paste0("__neg.msp"),
										   fileClasses = "classifier/classes.txt",
										   minimumIntensityOfMaximalMS2peak = 10, #100
										   minimumProportionOfMS2peaks = 0.005, #0.05
										   mzDeviationAbsolute_grouping = 0.1, #0.01
										   mzDeviationInPPM_grouping = 500 ) #10

# Diversity of compound classes
div_classes <- data.frame()
for (i in c("pos","neg")) {
	obj <- data.frame(mode=i, classes=unique(classifiers[[i]][,"Annotation (putative)"]), frequency=0)
	for (j in 1:nrow(obj)) obj[j,"frequency"] <- length(which(obj$classes[j] == classifiers[[i]][,"Annotation (putative)"]))
	div_classes <- rbind(div_classes, obj)
}

# Plot most abundant classes
pdf(file="classification_classcounts.pdf", encoding="ISOLatin1", pointsize=10, width=10, height=10, family="Helvetica")
par(mfrow=c(1,1), mar=c(15,4,4,1), oma=c(0,0,0,0), cex.axis=0.8, cex=0.9)
barplot(div_classes[1:nrow(div_classes),"frequency"], names.arg=gsub('.*\\;','',div_classes[1:nrow(div_classes),"classes"]), las=3, ylab="frequency", main="Most abundant compound classes")
dev.off()

# Determine how many spectra were classified
spectra_number <- 27 + 41 # neg + pos

spectra_classified <- sum(length(unique(classifiers[["pos"]]$`Metabolite name`)) +
						  length(unique(classifiers[["neg"]]$`Metabolite name`)))

classes_order <- readLines(con="classifier/classes.txt")
classes_order <- classes_order[which(grepl(pattern="^Organic compounds", x=classes_order))]
classes_order <- gsub(x=classes_order, pattern=":", replacement="")
classes_order <- gsub(x=classes_order, pattern="/", replacement="; ")

classes <- div_classes$classes
classes <- classes[which(grepl(pattern="^Organic compounds", x=classes))]
classes <- gsub(x=classes, pattern=":", replacement="")
classes <- gsub(x=classes, pattern="/", replacement="; ")

print(paste("Number of merged spectra:", spectra_number))
print(paste("Number of spectra classified:", spectra_classified))
print(paste("Number of unclassified spectra:", spectra_number - spectra_classified))
print(paste("Number of classes:", length(classes_order)))
print(paste("Number of classes with entities:", length(classes)))
print(paste("Number of classes without entities:", length(classes_order) - length(classes)))
print("Classes with entities:")
print(gsub(x=classes_order[which( (classes_order %in% classes))], pattern=".*; ", replacement=""))
print("Classes without entities:")
print(gsub(x=classes_order[which( ! (classes_order %in% classes))], pattern=".*; ", replacement=""))

classes <- classes_order[which(classes_order %in% classes)]

# Count numbers of matched classes in each sample
class_list <- data.frame(mode=c("pos","neg"))
for (i in classes) {
	vec <- NULL
	for (j in c("pos","neg")) {
		vec <- as.numeric(c(vec, sum(grepl(pattern=i, x=classifiers[[j]][,"Annotation (putative)"]), na.rm=TRUE)))
	}
	class_list <- cbind(class_list, vec)
}
colnames(class_list) <- c("mode", gsub(x=classes, pattern=".*\\; ", replacement=""))
rownames(class_list) <- c("pos","neg")
write.table(x=t(class_list), file="classifiers_classes.csv", sep=";", quote=TRUE, row.names=TRUE, dec=".")

# Classifiers
classifiers_pos_class <- get(load("classifier/2018-05-29_18_28_56_2018-02-13_09_14_10_LC-MS-MS_pos_21908.RData"))
classifiers_neg_class <- get(load("classifier/2018-05-29_18_21_05_2018-02-13_09_14_10_LC-MS-MS_neg_11328.RData"))

# Area under Precision Recall Curve
classifier_pos_auc <- unlist(lapply(X=classes, FUN = function(x) { for (i in 1:length(classifiers_pos_class)) { if (classifiers_pos_class[[i]]$class==x) { y <- classifiers_pos_class[[i]]$AUC; names(y) <- x; return(y) } } } ))
classifier_neg_auc <- unlist(lapply(X=classes, FUN = function(x) { for (i in 1:length(classifiers_neg_class)) { if (classifiers_neg_class[[i]]$class==x) { y <- classifiers_neg_class[[i]]$AUC; names(y) <- x; return(y) } } } ))

# True Positive Rate for False Positive Rate of 5 Percent
classifier_pos_fpr <- unlist(lapply(X=classes, FUN = function(x) { for (i in 1:length(classifiers_pos_class)) { if (classifiers_pos_class[[i]]$class==x) { y <- classifiers_pos_class[[i]]$TPR_for_FPR_of_5Percent; names(y) <- x; return(y) } } } ))
classifier_neg_fpr <- unlist(lapply(X=classes, FUN = function(x) { for (i in 1:length(classifiers_neg_class)) { if (classifiers_neg_class[[i]]$class==x) { y <- classifiers_neg_class[[i]]$TPR_for_FPR_of_5Percent; names(y) <- x; return(y) } } } ))

# Save table with AUC-PR and TPR-FPR rates
classifiers_validation <- data.frame(compound_class=gsub(x=classes, pattern='.*\\;', replacement=''))
classifiers_validation[which(classes %in% names(classifier_pos_auc)), "AUC-PR_pos"] <- classifier_pos_auc
classifiers_validation[which(classes %in% names(classifier_neg_auc)), "AUC-PR_neg"] <- classifier_neg_auc
classifiers_validation[which(classes %in% names(classifier_pos_fpr)), "TPR-FPR_pos"] <- classifier_pos_fpr
classifiers_validation[which(classes %in% names(classifier_neg_fpr)), "TPR-FPR_neg"] <- classifier_neg_fpr
write.table(x=classifiers_validation, file="classifiers_validation.csv", sep=";", quote=TRUE, row.names=FALSE, dec=".")

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

