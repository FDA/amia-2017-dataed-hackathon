---
layout: lesson
root: .
---

The "Reproducibility of Microarray and Gene Expression Analysis" course focuses on the reproduction of a figure included in the paper "[Repeated observation of breast tumor subtypes in independent gene expression data sets](http://www.pnas.org/content/100/14/8418.full)" by Sørlie et al. [Supporting Figure 6](http://www.pnas.org/content/suppl/2003/06/16/0932692100.DC1/2692Fig6.pdf) shows a heat map of the gene expression patterns across 122 breast tissue samples and two dendrograms representing the gene and array clustering.  The reproducibility of published results is significant in Biomedical Informatics and Data Science.  Information visualizations have become a key element in scientific presentations and it is therefore important to be able to recreate them by reusing the original data and following the same steps.

The education material supports the primary objective of the course and has been organized accordingly in six main lessons.  It particularly provides full guidance to any interested scientist with minimum to medium knowledge on the topic, as illustrated in the User Stories presented in the next section.  The course starts with an initial [Setup]({{ page.root }}/setup/) phase for data retrieval and software installation.  [Lesson 1]({{ page.root }}/01-expression/) and [Lesson 2]({{ page.root }}/02-hierarchical/) discuss key concepts in microarray analysis and hierarchical clustering that are necessary for understanding the heat map and dendrogram visualization illustrated in Supporting Figure 6.  [Lesson 3]({{ page.root  }}/03-figure-treeview/) steps through the reconstruction of the figure using the most recent version of the [TreeView tool](https://bitbucket.org/TreeView3Dev/treeview3) described in the paper's methodology to visualize the clustered [data](http://genome-www.stanford.edu/breast_cancer/robustness/data.shtml) distributed by the authors.  [Lesson 4]({{ page.root  }}/04-figure-hive-plot/) describes the use of hive plots for visualization of the clustered data.  Hive plots present an interesting alternative to the heat map-dendrogram combination of the original figure and can be created in R using the original data and code provided in the lesson.  The last two lessons take a step back and start from the semi-processed version of the data.  [Lesson 5]({{ page.root }}/05-analysis-cluster/) describes the reformatting of the semi-processed data for processing in the Cluster 3.0 tool (an earlier version of this tool was used for the analysis presented in the paper), the actual clustering process, and the visualization of results in TreeView 3.0.  The new figure and the clustering results are compared with the original findings and any differences or similarities are further discussed.  [Lesson 6]({{ page.root }}/06-analysis-r/) follows the same structure; however, the clustering and visualization are conducted in R.  The complete details at each step and the corresponding R code are provided in this lesson.

The Software Carpentry templates have been used to structure lessons into coherent blocks with estimates for the time required.

> ## Prerequisites
>
> Extensive background knowledge on gene expression or data analysis is not required,
> though participants should be familiar with basic scientific principles regarding
> data quality, research procedures, and the publication process.  Prior experience
> using the R language would be helpful in following certain lessons, but fully operable
> code is provided at all times, which will allow anyone to successfully reach the results.
>
{: .prereq}
