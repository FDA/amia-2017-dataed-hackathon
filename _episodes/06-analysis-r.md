---
title: "Recreating the Analysis with R"
teaching: 40
exercises: 20
questions:
- "How can R be used for hierarchical clustering and visualization of gene expression data?"
- "Can Supporting Figure 6 be recreated in R?"
objectives:
- "Learn to use R for hierarchical clustering and visualizing gene expression data."
- "Attempt to recreate Supporting Figure 6."
keypoints:
- "We can use R to write an uncentered correlation function, but the similarity calculations are not identical to those in the original paper."
- "We can use R to create dendrograms that show the clustering of arrays and genes, but the clusters are not identical to those in the original paper."
- "Our replicate heat map is not identical to the original heat map because of differences in clustering and color coding."
---

In this lesson we will explore how to cluster and visualize the gene expression data from Supporting Figure 6 with R. R has packages for hierarchical clustering and production of heat maps with dendrograms.

As discussed in [Lesson 5]({{ page.root }}/05-analysis-cluster/), the Cluster program was originally used by the authors to perform hierarchical clustering, and TreeView was used to visualize a heat map with dendrograms showing the clustering results. The original data file that was inputted in Cluster is not available, so we again choose the output files from Cluster as the starting point of this lesson. Therefore we will not add steps to transform, center, or normalize the data in R. We will read the SupplFigure6.cdt file exactly as provided, create a similarity metric, perform hierarchical clustering, and visualize the heat map and clustering results.

# Reading the Data

We can begin by loading the additional packages `dendextend` ([see also](https://cran.r-project.org/web/packages/dendextend/vignettes/introduction.html)), `dendextendRcpp`, and `NMF` that will be used in this lesson (they should installed following the instructions on the [Setup]({{ page.root }}/setup/) page first).

~~~
library(dendextend)
library(dendextendRcpp)
library(NMF)
~~~
{: .source}

As in [Lesson 4]({{ page.root }}/04-figure-hive-plot/), we will save the original SupplFigure6.cdt file in the current working directory specified in R, which can be determined with the function `getwd`. Then we can read in this file with the function `read.delim`. We select only the relevant rows and columns and store the gene expression data in a matrix format for use in subsequent steps. The weights for individual genes and arrays are also stored for future use. (Recall from [Lesson 5]({{ page.root }}/05-analysis-cluster/) that EWEIGHT is 1 for all arrays in the data file, but GWEIGHT takes one of four values for genes: 0.3, 0.35, 0.5, 1). These weights will influence the calculation of similarity between genes and arrays, respectively.

~~~
getwd()
orig_data <- read.delim("Suppl_Figure6.cdt", 
header = FALSE, 
as.is = TRUE)

expressions <- orig_data[4:nrow(orig_data), 5:ncol(orig_data)]

mat <- matrix(as.numeric(unlist(expressions)), nrow = nrow(expressions),
ncol = ncol(expressions))

dimnames(mat) <- list(orig_data[4:nrow(orig_data), 3],
orig_data[1, 5:ncol(orig_data)])

gweight <- as.numeric(orig_data[4:nrow(orig_data), 4])

eweight <- as.numeric(orig_data[3, 5:ncol(orig_data)])
~~~
{: .source}

# Creating a Similarity Metric

Next, we write a function `sim` to calculate similarity between pairs of genes and pairs of arrays. Sørlie et al. selected a similarity metric and an average linkage function in Cluster.  Although the details of the similarity metric are not discussed in the paper, they refer to the work of [Eisen et al](http://www.pnas.org/content/95/25/14863.full).  As in [Lesson 5]({{ page.root }}/05-analysis-cluster/), we assume that the authors selected uncentered correlation as the similarity metric for both genes and arrays. To write an uncentered correlation function in R, we refer to Section 3.4 of the Cluster 3.0 [user manual](http://bonsai.hgc.jp/~mdehoon/software/cluster/manual/index.html), which explains that missing entries are pairwise deleted during calculation of similarity. We choose to write an if-else statement to get a pair of complete entries. We use the weighted covariance function `cov.wt` in R to incorporate weights.

~~~
sim<-function(x,y,w){

	if(sum(is.na(x+y))>0){
		complete <-cbind(x,y,w)[-which(is.na(x+y)),]
	}
	else{
		complete <-cbind(x,y,w)
	}
	
	return(cov.wt(complete[,1:2], wt= complete[,3],center=FALSE,cor=TRUE)$cor[1,2])
}
~~~
{: .source}

Now that we are equipped with a function to calculate similarity, we can create distance matrices for arrays and genes. Starting with arrays, we create a square matrix for the arrays and store the similarity of every pair. We create a *distance* matrix by defining distance = 1 - similarity for each pair of arrays. Similarity can take values ranging from -1 to 1. Therefore, pairs of items with high similarity (close to 1), will have a distance close to zero, and pairs with low similarity (close to -1) will have a distance close to 2. Note that `gweight` is used for weighting the indices of the array-pair (not `eweight`).

~~~
a.sim.mat<-matrix(ncol=ncol(mat), nrow=ncol(mat))
dimnames(a.sim.mat)<-list(colnames(mat),colnames(mat))

for(i in 1:ncol(mat)){
	for(j in 1:ncol(mat)){
		if(is.na(a.sim.mat[i,j])){
			a.sim.mat[i,j]<-sim(mat[,i], mat[,j], gweight)
			a.sim.mat[j,i]<-a.sim.mat[i,j]
		}
	}
} 

a.dist<-1-as.dist(a.sim.mat)
~~~
{: .source}

Similarly, we create a square matrix for the genes and store the similarity of every pair. We again define distance = 1 - similarity for each pair of genes. Note that `eweight` is used for weighting the indices of the gene-pair.

~~~
g.sim.mat<-matrix(ncol=nrow(mat), nrow=nrow(mat))
dimnames(g.sim.mat)<-list(rownames(mat),rownames(mat))


for(i in 1:nrow(mat)){
	for(j in 1:nrow(mat)){
		if(is.na(g.sim.mat[i,j])){
			g.sim.mat[i,j]<-sim(mat[i,], mat[j,], eweight)
			g.sim.mat[j,i]<-g.sim.mat[i,j]
		}
	}
} 

g.dist<-1-as.dist(g.sim.mat)
~~~
{: .source}

We can inspect our similarity calculations to see if our process is identical to the one originally used in Cluster. For arrays, we can check our `a.sim.mat` matrix for the similarity of *NormBreast2* (corresponding to ARRY120X in the original SupplFigure6.atr file) with *NormBreast3* (corresponding to ARRY121X in the original SupplFigure6.atr file).

~~~
a.sim.mat[which(rownames(a.sim.mat)=="NormBreast2"), which(rownames(a.sim.mat)=="NormBreast3")]
~~~
{: .source}

We observe that in our replicate in R, the similarity between these arrays is 0.7835408, which is identical (after rounding) to the result of our re-clustering in [Lesson 5]({{ page.root }}/05-analysis-cluster/). However, as discussed in [Lesson 5]({{ page.root }}/05-analysis-cluster/) the similarity between ARRY120X and ARRY121X is 0.783569336 in the original calculation. Therefore we have already shown that our similarity calculations are not identical to the originals.

For genes, we can compare *116219 GATA3 GATA binding protein 3 Hs.169946 H72474* (corresponding to GENE25X in the original .gtr file) with *101778 GATA3 GATA binding protein 3 Hs.169946 R31441* (corresponding to GENE7X).

~~~
g.sim.mat[which(rownames(g.sim.mat)=="116219 GATA3 GATA binding protein 3 Hs.169946 H72474 "), which(rownames(g.sim.mat)=="101778 GATA3 GATA binding protein 3 Hs.169946 R31441 ")]
~~~
{: .source}

In our replicate, the similarity is 0.9766245, which is identical after rounding to the result of our re-clustering in [Lesson 5]({{ page.root }}/05-analysis-cluster/). As discussed in [Lesson 5]({{ page.root }}/05-analysis-cluster/), the similarity between GENE7X and GENE25X is 0.97668457 in the original calculation.  We conclude that we have not identically recreated the original similarity calculations. It is possible that we have recreated our own steps from [Lesson 5]({{ page.root }}/05-analysis-cluster/) that were performed in Cluster 3.0; however, the similarity matrices are not available from Cluster 3.0 for performing a complete comparison. This lesson does not aim to recreate the similarity matrices from Cluster 3.0.

# Performing Hierarchical Clustering

We can proceed to clustering after calculation of the distance matrices. We can cluster arrays and genes independently by using the function `hclust` in R with an average linkage. We re-order the resulting clustering, where possible, to mimic the original order that was reported by the authors. Recall from [Lesson 2]({{ page.root }}/02-hierarchical/) that we can re-order the two children of a node without changing the clustering. Finally, we create dendrograms for the arrays and genes for use in the next step, visualization.

~~~
a.hclust<-hclust(a.dist, method="average")

a.dend.plot<-as.dendrogram(rotate(a.hclust, order=match(1:122,a.hclust$order)))

g.hclust<-hclust(g.dist, method="average")

g.dend.plot <- as.dendrogram(rotate(g.hclust, order=match(1:552,g.hclust$order)))
~~~
{: .source}

# Visualizing the Results

Before plotting the heat map, we can plot the array and gene dendrograms, respectively. We want to use these plots to inspect whether each dendrogram matches the original one. For arrays, we will again inspect the positions of Norway 83-BE (corresponding to ARRY103X in the original SupplFigure6.atr file) and Norway FU02-BE (corresponding to ARRY9X). We will plot and color these arrays in red for clarity. The array dendrogram will be saved as *Lesson6_arrDend.png* in the current working directory.

~~~
col.arr <-rep("black", 122)
col.arr[which(labels(a.dend.plot)=="Norway 83-BE")]<-"red"
col.arr[which(labels(a.dend.plot)=="Norway FU02-BE")]<-"red"

png("Lesson6_arrDend.png", width=1920, height=960)
a.dend.plot %>% color_branches(col = col.arr) %>% color_labels(col = col.arr) %>% plot(main="Array Cluster Dendrogram")
dev.off()
~~~
{: .source}

<a href="{{ page.root }}/fig/r-array-dend-large.png"><img src="{{ page.root }}/fig/r-array-dend.png" alt="The array dendrogram reclustered by R" /></a>

As in [Lesson 5]({{ page.root }}/05-analysis-cluster/), we find that Norway 83-BE and Norway FU02-BE are clustered as children under the same node in the replicate dendrogram, whereas they are children under different nodes in the original dendrogram (we can see this by inspecting the original SuppleFigure6.atr file).

For genes, we will inspect the positions of *108924  \*\*Homo sapiens cDNA FLJ34425 fis, clone HHDPC2008297 Hs.120638 N81017* (corresponding to  GENE31X in the original SupplFigure6.gtr file) and *120444 N33 Putative prostate cancer tumor suppressor Hs.71119 H13424* (corresponding to GENE12X). We will plot and color these genes in red for clarity. The gene dendrogram will be saved as *Lesson6_geneDend.png* in the current working directory.

~~~
col.genes <-rep("black", 552)
col.genes[which(labels(g.dend.plot)=="108924  **Homo sapiens cDNA FLJ34425 fis, clone HHDPC2008297 Hs.120638 N81017 ")]<-"red"
col.genes[which(labels(g.dend.plot)=="120444 N33 Putative prostate cancer tumor suppressor Hs.71119 H13424 ")]<-"red"

g.dend.plot <- set(g.dend.plot,"labels", substr(labels(g.dend.plot), start=1, stop=10))
png("Lesson6_geneDend.png", width=3600, height=1000)
g.dend.plot %>% assign_values_to_leaves_nodePar(c(rep(0.7, 552)), "lab.cex") %>% color_branches(col = col.genes)%>% color_labels(col = col.genes) %>% hang.dendrogram() %>% plot(main="Gene Cluster Dendrogram")
dev.off()
~~~
{: .source}

<a href="{{ page.root }}/fig/r-gene-dend-large.png"><img src="{{ page.root }}/fig/r-gene-dend.png" alt="The gene dendrogram reclustered by R" /></a>

In the replicate figure, we observe that these genes are children under different nodes, whereas in the original figure they are children of the same node (we can see this by inspecting the original SupplFigure6.gtr file). They are also ordered far apart from each other in the replicate figure, whereas they are ordered next to each other in the original figure. With these counterexamples, we conclude we have not recreated the identical dendrograms presented by Sørlie et al.

Next, we plot a heat map using the function `aheatmap` in the package `NMF`. Our replicate dendrograms have been used to order the rows and columns of the heat map accordingly.

~~~
png("Lesson6_heatmap.png", width=960, height=2400)
aheatmap(mat, Rowv=g.dend.plot, Colv=a.dend.plot, col=c("green", "black", "red"))
dev.off()
~~~
{: .source}

<a href="{{ page.root }}/fig/r-heat-map-large.png"><img src="{{ page.root }}/fig/r-heat-map.png" alt="The heat map generated by R" /></a>

Due to the differences between the replicate and original dendrograms, it follows that the replicate heat map is not ordered identically to the original figure. We also observe that the color coding in the replicate heat map is not identical to the original figure. For example, it is unclear how gene expression values that exceed 5.6 overexpression or 5.6 underexpression were represented with the original color scale. The replicate heat map uses a color scale with a greater range of values. In this lesson, we have not aimed to direct the values of the heat map to achieve the same color coding. Additional exploration is encouraged outside of this lesson.

# Examination Questions

> ## 1
> Write another similarity metric in R and compare the results of this function with the similarity calculations appearing in the original SupplFigure6.atr and SuppleFigure6.gtr files. For example, does a centered correlation function better replicate the original calculations?
{: .challenge}

> ## 2
> Explore whether an adjustment of the color coding in the heat map achieves results that are closer to the original figure. For example, restrict the color scale to a 5.6 magnitude difference in each direction. How does the new figure compare to the original one?
{: .challenge}
