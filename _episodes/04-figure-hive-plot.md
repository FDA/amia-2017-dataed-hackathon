---
title: "Recreating Supporting Figure 6 with Hive Plots"
teaching: 40
exercises: 20
questions:
- "What is a hive plot?"
- "How can a hive plot be created using R?"
- "How can gene expression data be represented with a hive plot?"
objectives:
- "Explain how hive plots represent network data."
- "Demonstrate creation of a hive plot for gene expression data."
keypoints:
- "Hive plots are intended to represent large network data."
- "The expression profiles for a large number of genes can be visualized using a hive plot."
---

In this lesson, we will take a detour away from investigating reproducibility to explore an alternative visualization approach for microarray analyses.  In place of a heat map showing gene expression, we will represent the microarray results as a network.  Networks attempt to visualize complex relationships between entities (displayed as nodes) by connecting them with edges.

In this lesson, we will use the R language to create a network based on the expression values in Supporting Figure 6 and visualize this network as a hive plot.  You should have installed R following the instructions on the [Setup]({{ page.root }}/setup/) page, which will allow you to run all of the commands in this lesson.

# Introduction to Hive Plots

Large and complex networks are not easily represented by traditional networks visualizations, as they can become cluttered and difficult to interpret.  [Hive Plots](http://www.hiveplot.net/) are a relatively new method for displaying large network data, where nodes are placed along meaningful axes, providing a quantitative aspect to the visualization that is often lacking in other network representations.

## The HiveR Package

We will be using the [HiveR package](http://academic.depauw.edu/~hanson/HiveR/HiveR.html) in R to create the hive plots in this lesson.  For detailed information about this package, refer to the [HiveR Reference Manual](https://cran.r-project.org/web/packages/HiveR/HiveR.pdf).

> ## Other Languages
> Hive plot implementations exist in other commonly used languages as well, including Java and Python.
{: .callout}

# Converting the Figure to a Hive Plot

The original data for the heat map in Supporting Figure 6 can be used to construct a hive plot that demonstrates many of the same key points.  We will walk through the steps needed to convert the original data into a suitable format and then use the HiveR package to draw the hive plot.

## Planning the Hive Plot

Before starting to work on the visualization, it is worth considering the nature of the data and deciding on a useful hive plot representation.  In this case, we have a list of arrays that have each been evaluated for a large number of genes to find their expression levels relative to a baseline.  The important associations are thus between arrays (or groups of arrays) and genes.  The association also takes two forms: overexpression and underexpression.

In a network, associations are represented as the edges that connect two nodes.  For hive plots, the nodes are on various axes, and the edges fill the spaces in between.  Here, we will represent the arrays on one axis and create two copies of the gene axis.  This will let us show both the overexpressed and underexpressed associations between arrays and genes using multiple axes.

## Reading the Data

The first step is reading the SupplFigure6.cdt data file (the contents of this file are discussed in [Lesson 5]({{ page.root }}/05-analysis-cluster/)), but first we have to make sure that R knows where it is.  You can check the current working directory in R with the `getwd` command.

~~~
getwd()
~~~
{: .source}

~~~
[1] "C:/Users/User.Name/Documents"
~~~
{: .output}

By default, R will be running in your Documents folder on Windows, so let's put the data files there for now.  You'll need at least the SupplFigure6.cdt file (instructions are on the [Setup]({{ page.root }}/setup/) page) for this lesson.

> ## Change Working Directory
> Note that you can change the working directory in R with the `setwd` command if you'd like to have R read the files from a different location.
{: .callout}

Now we're able to read in the file and begin to set up our R data structure, which will be a matrix containing the numeric expression values.

~~~
orig_data <- read.delim("SupplFigure6.cdt",
                        header=FALSE,
                        as.is=TRUE)

expressions <- orig_data[4:nrow(orig_data), 5:ncol(orig_data)]

mat <- matrix(as.numeric(unlist(expressions)),
              nrow=nrow(expressions), ncol=ncol(expressions))

dimnames(mat) <- list(orig_data[4:nrow(orig_data),3],
                      orig_data[1, 5:ncol(orig_data)])
~~~
{: .source}

Then, we will duplicate each row of the matrix, so that genes are represented twice and can be placed on two separate axes.  We'll differentiate the names of the duplicated rows by putting the letters "dup" before the gene name.

~~~
mat.thresh <- mat[rep(seq_len(nrow(mat)), each = 2), ]
indices <- seq(2, nrow(mat.thresh), by=2)
rownames(mat.thresh)[indices] <- lapply(rownames(mat.thresh)[indices], function(x) paste0("dup", x))
~~~
{: .source}

Next, we will retain only the strong associations between genes and arrays, i.e. those with a measured overexpression or underexpression greater than some threshold.  We set this threshold (`express.thresh`) at 4.0 at this time because it leads to visualizations that are not very cluttered, but you are encouraged to experiment with other values.  We will also separate the gene listings for positive and negative values at this time.  We retain the large positive values in the data frame for the rows based on the original gene names and retain the large negative values for the duplicated rows.  We accomplish this by indexing the odd and even rows separately.

~~~
express.thresh <- 4.0
mat.thresh[mat.thresh <= express.thresh & mat.thresh >= -express.thresh] <- NA

indices <- which(mat.thresh < 0, arr.ind = TRUE)
mat.thresh[indices[seq(1, nrow(indices), by=2), ] ] <- NA

indices <- which(mat.thresh > 0, arr.ind = TRUE)
mat.thresh[indices[seq(2, nrow(indices), by=2), ] ] <- NA
~~~
{: .source}

For this data to represent a network, we need to transform the matrix so it looks like an adjacency matrix.  An adjacency matrix (***A***) lists all the nodes in a network as both column names and row names, creating a square matrix.  The value *a*<sub>*ij*</sub> within the matrix represents the weight of the edge between the nodes corresponding to row *i* and column *j*.  Fortunately, we don't have to add any new values to create the proper adjacency matrix, but we do need to add new columns (representing the gene nodes) and new rows (representing the array nodes).

~~~
# First add all the gene names from the rows as new columns (on the right)
extra <- matrix(ncol = nrow(mat.thresh), nrow = nrow(mat.thresh))
colnames(extra) <- rownames(mat.thresh)
mat.sq <- cbind(mat.thresh, extra)

# Now add all the array names from the columns as new rows (at the top)
extra <- matrix(ncol = ncol(mat.sq), nrow = ncol(mat.thresh))
rownames(extra) <- colnames(mat.thresh)
mat.sq <- rbind(extra, mat.sq)
~~~
{: .source}

Finally, we replace values in the matrix that cannot be read as edge weights by HiveR.  This means replacing `NA` values with zeros and taking the absolute value of the matrix to remove negative edge weights.  Since we have already prepared our second gene axis to only contain negative values, we will still be able to differentiate overexpressed and underexpressed genes.

~~~
mat.sq[is.na(mat.sq)] <- 0
mat.sq <- abs(mat.sq)
~~~
{: .source}

## Creating the Hive Plot

With the adjacency matrix ready, we can invoke the function `adj2HPD` from the HiveR package to create a hive plot data (HPD) object, which we will soon be able to use to create the hive plot visualization.

~~~
require(HiveR)
f6hive <- adj2HPD(M = mat.sq,
                 axis.cols = c("white"),
                 type = "2D",
                 desc = "the array and gene expression data")
~~~
{: .source}

If we examine the structure of the data object with the function `str`, we see that nodes and edges are both represented with several attributes.

~~~
str(f6hive)
~~~
{: .source}

~~~
List of 5
 $ nodes    :'data.frame':	1226 obs. of  6 variables:
  ..$ id    : int [1:1226] 1 2 3 4 5 6 7 8 9 10 ...
  ..$ lab   : chr [1:1226] "Norway 51-BE" "Stanford 16" "Norway 39-BE" "Norway 17-BE" ...
  ..$ axis  : int [1:1226] 1 1 1 1 1 1 1 1 1 1 ...
  ..$ radius: num [1:1226] 1 1 1 1 1 1 1 1 1 1 ...
  ..$ size  : num [1:1226] 1 1 1 1 1 1 1 1 1 1 ...
  ..$ color : chr [1:1226] "black" "black" "black" "black" ...
 $ edges    :'data.frame':	410 obs. of  4 variables:
  ..$ id1   : int [1:410] 129 133 133 135 135 135 136 139 139 141 ...
  ..$ id2   : int [1:410] 79 82 90 60 73 88 45 60 86 60 ...
  ..$ weight: num [1:410] 4.57 4.03 4.08 4.7 4.01 ...
  ..$ color : chr [1:410] "gray" "gray" "gray" "gray" ...
 $ desc     : chr "the array and gene expression data"
 $ axis.cols: chr "white"
 $ type     : chr "2D"
 - attr(*, "class")= chr "HivePlotData"
~~~
{: .output}

To create even an initial version of this hive plot, the most important attributes to set are `nodes$axis` and `nodes$radius`.  Setting the axes that nodes will appear on will be easy because we have already planned our 3 axes.  Here, we assign the axes by taking advantage of the order in which the node names appear in the matrix, with the array nodes listed first and the gene nodes alternating between the original and duplicate names.  The array nodes will be assigned to the first axis, the original gene nodes to the second, and the duplicate gene nodes to the third.

~~~
f6hive$nodes$axis <- c(rep(as.integer(1), ncol(mat)),
                       rep(c(as.integer(2), as.integer(3)), nrow(mat)))
~~~
{: .source}

The `nodes$radius` variable defines how far out along the axis each node should appear.  Instead of setting this manually, we will use the built-in ranking feature, which spaces nodes according to the order in which they appeared in the input data.

~~~
f6hive <- manipAxis(f6hive, "rank")
~~~
{: .source}

With these variables set, it is now possible to produce a basic hive plot showing over and underexpression links between arrays and genes.

~~~
plotHive(f6hive, bkgnd = "white")
~~~
{: .source}

<img src="{{ page.root }}/fig/hive-plot-1.png" alt="A poorly visualized hive plot showing strong connections between genes and arrays from Supporting Figure 6" />

The array axis is at the top, and you can see that edges connect certain arrays on this axis to gene locations along the two gene axes in the lower left and lower right.  There is, however, clearly room for improvement in the visuals of this plot, so let's make a few adjustments.  We will reweight the edges to make them thinner and easier to distinguish, as well as reducing the size of the nodes, so they don't appear as a thick, solid line.  We will also stretch the array axis, so it provides better detail.

~~~
f6hive$edges$weight <- f6hive$edges$weight / max(f6hive$edges$weight)
f6hive$nodes$size <- 0.4
f6hive <- manipAxis(f6hive, "scale", action = c(3.75, 1.0, 1.0))
plotHive(f6hive, bkgnd = "white")
~~~
{: .source}

<img src="{{ page.root }}/fig/hive-plot-2.png" alt="A more easily readable version of the hive plot showing strong connections between genes and arrays from Supporting Figure 6" />

This is instantly a much more readable plot, and more information about the connections can be discerned.  By placing the nodes along the axis in the same order as they appear in Supporting Figure 6, we have preserved some of the layout information from the heat map for this plot.  Each axis in the hive plot is ordered out from the center, so the arrays near the left side of the heat map will be in the lower part of the array axis (the top axis) on the hive plot.  Similarly, genes listed near the top of the heat map will appear close to the center of the hive plot.  They are mirrored on both of the gene axes (lower left and lower right).  With this hive plot, we can see some strong relationships between arrays and genes that were grouped near each other in the heat map.

But the biggest question right now is, "which side is which?"  It's not yet clear which edges represent overexpressed genes and which underexpressed.  A quick comparison of the values in the data file against Supporting Figure 6 reveals that positive values correspond to red coloring in the figure, so let's set our positive edges to red and our negative edges to green.  Recall that we retained positive values in the odd rows with the original gene names and negative values in the even rows with the duplicated gene names.

~~~
f6hive$edges$color <- ifelse(f6hive$edges$id1 %% 2 == 0, "green", "red")
plotHive(f6hive, bkgnd = "white")
~~~
{: .source}

<img src="{{ page.root }}/fig/hive-plot-3.png" alt="A version of the hive plot showing strong connections between genes and arrays from Supporting Figure 6 with colored edges representing the type of connection" />

Now we can identify many of these connections as the same strong regions from Supporting Figure 6.  For example, the red curves that connect the top part of the array axis to genes near the midpoint of the gene axis represent the red region on the right side of the heat map about halfway down (cutout E in Figure 1 in the paper).

# More Representations

There is more we can do to make this hive plot clearly display the same information as the heat map.  For example, we can use the subtypes that the authors provide on the array axis to show how the types correlate with different genes.  You can obtain the names of the arrays that were classified into each subtype from [Supporting Table 3](http://www.pnas.org/content/suppl/2003/06/16/0932692100.DC1/2692Table3.pdf) from the paper.  We have correlated these array names with R colors in [this CSV file]({{ page.root }}/data/subtype_colors.csv), which you can download and place into the R working directory.  Now R can assign colors to specific nodes on the array axis.

~~~
f6hive$nodes$color <- "gray60"
subclass.colors <- read.csv("subtype colors.csv",
                            header=TRUE, stringsAsFactors = FALSE)
replace <- subclass.colors$color[match(f6hive$nodes$lab,
                                       subclass.colors$array)]
f6hive$nodes$color[!is.na(replace)] <- replace[!is.na(replace)]
~~~
{: .source}

<img src="{{ page.root }}/fig/hive-plot-4.png" alt="A version of the hive plot showing strong connections between genes and subtype colored arrays from Supporting Figure 6 with colored edges representing the type of connection" />

Let's also bring some of that color distinction to the edges, to make it clear how the array subtypes are connected to various genes.  This line will replace each edge's color with the color used for the node stored in its `id2` position.  We have also post-processed the final image to add labels.

~~~
f6hive$edges$color <- f6hive$nodes$color[match(f6hive$edges$id2, f6hive$nodes$id)]
plotHive(f6hive, bkgnd = "white")
~~~
{: .source}

<img src="{{ page.root }}/fig/hive-plot-labeled.png" alt="A version of the hive plot showing subtype colored connections between genes and arrays from Supporting Figure 6" />

# Further Hive Plot Topics

Hive plots can be useful for a variety of other bioinformatics topics as they specialize in highlighting connections in relatively large networks.  Refer to the main [hive plot website](http://www.hiveplot.net/) for many more examples and inspiration for further use.

# Examination Questions

> ## 1
> Recreate the plots with the expression threshold for visible connections set to 2.0, 3.0, and 5.0.  How do these plots differ from the plots shown in the lesson (which used a threshold of 4.0)?  What are the challenges in using a low threshold?  What are the challenges in using a high threshold?  Which plot do you think is the most informative?
{: .challenge}

> ## 2
> Think of at least one other way to represent the data from Supporting Figure 6 as a hive plot by defining the axes or the edge connections differently.  Don't be constrained to 3 axes; hive plots can have anywhere from 2 to 6 axes!  Think about the kinds of strong relationships that could show up in your visualization and exactly what they would look like.  If you have time, go ahead and make the plot to see if it looks like what you expect.
{: .challenge}
