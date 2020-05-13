# ############################## Classification with MetFamily ##############################



print("---------------------------------------")
print(classifiers_validation)
print("---------------------------------------")
quit()

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

