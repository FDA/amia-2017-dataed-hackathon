---
title: "Recreating the Analysis with Cluster and TreeView"
teaching: 40
exercises: 20
questions:
- "What tools were used for analysis in the original paper?"
- "Can certain analysis steps be reproduced?"
objectives:
- "Learn how to analyze the original data in Cluster 3.0 by defining the same parameters described in the publication."
- "Create the new visualization using the new clustering results and observe the differences between the new and the original heat map/dendrograms."
keypoints:
- "The reclustering of the original data with Cluster 3.0 does not support the recreation of Supporting Figure 6."
---

We successfully recreated a visualization of the final clustered data in [Lesson 3]({{ page.root }}/03-figure-treeview/) with TreeView 3.0.  However, we would also like to examine the clustering procedure performed by Sørlie et al. to see if this analysis step can also be reproduced.  In this lesson, we will use the same tools as the authors to perform the hierarchical clustering and compare the results

# Analysis Tools

The primary paper by Sørlie et al. cites previous analysis published in "[Gene expression patterns of breast carcinomas distinguish tumor subclasses with clinical implications](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC58566/)."  The "Microarray analysis" paragraph in the methodology section of that paper describes the analysis of a subset of the same data used in the current paper.  Based on the description in the two papers, it is concluded that the average-linkage hierarchical clustering was applied in both cases.  The Cluster program supported the clustering process, and the results were displayed using the TreeView tool described in [Lesson 4]({{ page.root }}/04-figure-hive-plot/).  Follow the instructions on the [Setup]({{ page.root }}/setup/) page to obtain both of these tools.

Cluster 3.0 contains algorithms for the analysis of gene expression data. As described in the Cluster 3.0 online manual, the tool contains routines that support hierarchical, k-means and k-medians clustering as well as the construction of 2D self-organizing maps.

> ## Cluster 3.0
> Cluster 3.0 is available for Microsoft Windows, Mac OS X and Linux/Unix platforms, and routines are available in Python and Perl.
{: .callout}

The default user interface for Cluster 3.0 is shown below.

<img src="{{ page.root }}/fig/cluster-interface.png" alt="The initial user interface for Cluster 3.0 before importing any data" />

# Data Pre-Processing

As mentioned on the [Setup]({{ page.root }}/setup/) page, the authors have provided the following three data files, which will be discussed in this lesson:

 * SupplFigure6.cdt
 * SupplFigure6.gtr
 * SupplFigure6.atr

These files have the characteristic extensions (.cdt, .gtr, and .atr) of output from the Cluster software.

Cluster 3.0 accepts a formatted .txt file with gene expressions as its input. After loading the gene expression data in Cluster 3.0, the user has the option to perform pre-clustering steps (such as filtering and adjustment of the data), followed by the selection of clustering parameters and the actual clustering step. Cluster 3.0 independently clusters genes and microarrays and outputs .cdt, .gtr, and .atr data files containing information about each clustering. The .cdt file contains the data on which the clustering was performed (following any optional pre-clustering steps by the user). The .gtr and .atr files contain the gene tree and array tree, respectively, that were created by Cluster 3.0. The [user manual](http://bonsai.hgc.jp/~mdehoon/software/cluster/manual/index.html) for Cluster 3.0 discusses the formatting of the input file, clustering parameters, and output files.

The original .txt input file that the authors imported in Cluster is not available, so we have chosen the .cdt output file of Cluster as the starting point of this lesson. As per the Cluster 3.0 manual, the .cdt file contains the values used to perform the clustering (i.e. the values that existed immediately prior to the clustering but after any data transformation, centering, and normalization by the user). The authors have discussed that pre-clustering steps included log<sub>2</sub> normalization and median-centering of rows. This lesson does not aim to explore any other filtering or adjustments that may have been made by the authors prior to clustering. We expect that re-processing the data as it appeared immediately prior to clustering with an identical selection of parameters will produce an identical hierarchical clustering result.

The .cdt output file follows most, but not all, of the specifications required by Cluster 3.0 for an input file.

Prior to using any of the tool functionalities to perform the clustering process, you will need to pre-process the data distributed with the Sørlie et al. paper and generate a readable format for Cluster 3.0.  Here you will use tab-delimited text files in a particular format as shown below (found under Help>File format).

<img src="{{ page.root }}/fig/cluster-format.png" alt="The window in Cluster 3.0 that describes the expected format for input data" />

Rows represent genes and columns represent samples or observations.  The file must contain at least the gene names (in red along the left side), the names of the observations (also in red along the top), and the DATA (shaded in green).  The other columns and rows in black are optional.  DATA entries may be blank.

We will use a spreadsheet application to edit the original data file (SupplFigure6.cdt) from the paper to match this format.  Check the [Setup]({{ page.root }}/setup/) page for instructions on obtaining this data file, as well as the two accompanying files that will be needed later in this lesson.  Open the data file in the spreadsheet application.  Once open, you’ll see that the file looks very similar to the format required for Cluster, but that there are some extra entries that should ideally be removed.  We have highlighted the rows and columns that **should be deleted** below.  These areas contain information associated with either the results of the clustering performed during the initial analysis or other optional information.

<img src="{{ page.root }}/fig/cluster-input-highlight.png" alt="Thee first few rows and columns in the data file with entries that should be deleted highlighted in yellow" />

We particularly delete the first column and second row that include the gene and array codes that correlate with the clustering information found in the other two data files for the paper (SupplFigure6.gtr and SupplFigure6.atr).  Recall that these files supported the construction of the dendrograms shown in Supporting Figure 6 as described in [Lesson 2]({{ page.root }}/02-hierarchical/).  The second column is not necessary either and should be removed. The re-formatted file should look like this:

<img src="{{ page.root }}/fig/cluster-input.png" alt="The first few rows and columns of the pre-processed data file that is ready to be used in Cluster 3.0" />

We will retain the GWEIGHT and EWEIGHT information, even though they are optional for Cluster 3.0, because the original paper includes the description, "Replicate clones representing one single UniGene cluster were weighted before clustering," in the Supporting Materials and Methods.  It is therefore important to include the GWEIGHT column that includes four distinct values for all genes: 0.3, 0.35, 0.5, and 1.  Although the EWEIGHT column contains the same value for all arrays ("1") and does not appear to contribute to the overall clustering, it is retained in the re-formatted file.

We save the re-formatted file as a tab-delimited text file (name it as Suppl_Figure6.txt) and place it into a new subfolder under the folder that was created for [Lesson 3]({{ page.root }}/03-figure-treeview/) ("Cluster_Results\new_clustering").

# Reclustering and Comparison

We open the tab-delimited text file containing the data in Cluster 3.0 by selecting the *File>Open data file* command.  The program checks the number of columns in each row of the file to see if there are any discrepancies, which would be reported in an error message.  If each row appears correct, the data will be loaded into the main interface as shown.

<img src="{{ page.root }}/fig/cluster-interface-filled.png" alt="The user interface for Cluster 3.0 after preparing to perform the reclustering" />

The path for the loaded file is shown at the top.  This lesson does not aim to make any selections under the tabs "Filter Data" or "Adjust Data" because the imported .txt file contains data post-filtering and post-adjustment by the authors. Therefore, we can proceed directly to the "Hierarchical" tab to select the parameters for clustering. To run the clustering process, check the Cluster boxes in the Gene and Array frames, retain the default "Job name" as Suppl_Figure6, and click on the Average Linkage button.  (Note, it is not specified if the authors used centered or uncentered similarities. In this lesson, we use “Correlation (centered)” for both genes and arrays.) Cluster 3.0 will generate three output files and place them in the same folder: *Suppl_Figure6.cdt*, *Suppl_Figure6.gtr*, and *Suppl_Figure6.atr*.

The *Suppl_Figure6.cdt* output file includes the data used in the hierarchical clustering process as well as the gene and array IDs generated by the algorithm to link this with the next two files that contain the appropriate information to build the dendrograms.  Here's a sample of the data file; the codes from the original data file are shown in red.

<img src="{{ page.root }}/fig/recluster-compare-1.png" alt="A portion of the main output file of Cluster 3.0 from the reclustering with some original data in red for comparison" />

The *Suppl_Figure6.gtr* and *Suppl_Figure6.atr* output files (shown below as A and B, respectively) contain gene tree and array tree information generated by the clustering algorithm.  Again, we are showing the first few lines from the paper's original data files in red for comparison.

<img src="{{ page.root }}/fig/recluster-compare-2.png" alt="Portions of both of the tree output files of Cluster 3.0 from the reclustering with some original data in red for comparison" />

These comparisons show that there is a different ordering to the genes and arrays between the original clustering data files and the newly clustered data.  However, the above figures do not illustrate the actual differences between the dendrograms as arrays or genes could still belong to the same clusters. To better explore the differences and similarities between the dendrograms, we import the new clustering results into TreeView 3.0 following the process described in [Lesson 3]({{ page.root }}/03-figure-treeview/).

Let's compare the generated heat map and dendrograms for the new clustering with the versions for the original data.  This can be accomplished by opening two separate instances of TreeView and importing the different datasets into them.

<a href="{{ page.root }}/fig/recluster-compare-3-large.png"><img src="{{ page.root }}/fig/recluster-compare-3.png" alt="A side-by-side comparison of the new clustering with the original clustering as shown in TreeView 3.0" /></a>

The comparison between the two clearly shows that not only the order of arrays and genes in the dendrograms but also the clustering itself differs. For example, the "Norway 83-BE" array is clustered with the "Norway FU02-BE" array based on the new results; however, it forms its own cluster in the original clustering presented in the paper. An overview of the remaining results indicates that even though there are some similarities, there are also noticeable differences between the old and the new clustering. We can therefore argue that **it is not possible to fully reproduce Supporting Figure 6 based on the reclustering of the original dataset distributed with the paper**.

# Examination Questions

> ## 1
> Open two instances of TreeView 3.0, one with the three original data files and one with the three new data files.  Choose two arrays and evaluate whether they are clustered differently in each instance of TreeView. Repeat the same for two genes. Keep the TreeView instances for the next exercise.
{: .challenge}

> ## 2
> Delete the GWEIGHT column and the EWEIGHT row from the original SupplFigure6.cdt file distributed with the paper and repeat the clustering process. Open a new, third, instance of TreeView 3.0 and compare it with the two previous instances from exercise 1. What are the differences you observe, especially between the second and third instances of the visualization (the new clustering with and without weights)?
{: .challenge}
