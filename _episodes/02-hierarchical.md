---
title: "Introduction to Hierarchical Clustering"
teaching: 10
exercises: 10
questions:
- "What is the basic principal behind hierarchical clustering?"
- "How can we interpret the branch levels?"
- "How are the cancer subtypes defined and illustrated in Supporting Figure 6?"
objectives:
- "Understand the basic principles behind hierarchical clustering."
- "Learn how to interpret a dendrogram."
keypoints:
- "Hierarchical clustering groups observations on a continuous scale and allows a dissimilarity threshold to determine the number of clusters."
---

Understanding hierarchical clustering is essential to interpreting not just Supporting Figure 6, but also the conclusions reached in the Sørlie et al. paper.  Looking at Supporting Figure 6 (or the sample heat map from [Lesson 1]({{ page.root }}/01-expression/)) you can see two dendrograms, one at the very top and one to the left of the microarray heat map.  This lesson will teach you the basic principles behind hierarchical clustering and how this method was used to create the dendrograms in Supporting Figure 6.  At the end of this lesson, you will understand how this technique allowed the authors to divide both genes and tumor samples into groups.

# Principles of Hierarchical Clustering

Hierarchical clustering groups observations and displays them in dendrograms.  A sample dendrogram is shown here.  Clustering is accomplished by calculating a dissimilarity measure between all observations and then combining the observations that are least dissimilar.

<img src="{{ page.root }}/fig/hier-cluster-sample.png" alt="A sample dendrogram showing hierarchical clustering of hypothetical data" />

At the bottom of the dendrogram, each observation is shown as its own cluster.  Moving up the dendrogram, the observations that are least dissimilar to each other are clustered, which reduces the number of observations for future clustering (starting with *n* observations, clustering two observations leads to *n*−1 observations for future clustering).  In this example dendrogram, the first two observations to merge are the ones labeled 4 and 9.  This process of merging the least dissimilar clusters repeats until all observations merge into a single cluster at the top of the dendrogram.

The height at which observations merge is the dissimilarity measure used to assess differences between clusters.  The lower in the dendrogram two observations merge, the more similar (i.e. less dissimilar) the two observations (or clusters of observations) are to each other.  The higher in the dendrogram two observations merge, the more dissimilar the observations are.  The dissimilarity metric value that defines the height can be seen in the scale bar to the left of the dendrogram in the sample above.  It is important to note however that not all dendrograms contain such a scale bar.

As hierarchical clustering represents data hierarchically, the number of possible clusters can range from one to *n* depending on the number of observations (*n*) and the dissimilarity value chosen.  To determine the number of clusters you will use for further analysis, pick an appropriate dissimilarity value and make a horizontal cut across the dendrogram.  The number of clusters beneath the height of this cut is the final number of clusters.

## Important Considerations

A few important points should be considered before using hierarchical clustering to draw conclusions from your data.  First, the dissimilarity metric should be carefully chosen to reflect the underlying relationship within the data.  For example, clustering using Euclidean distance may yield different results than clustering using correlation-based distance; choosing a dissimilarity measure that reflects relationships in your data is essential.  Second, conclusions about the similarity of two observations should only be drawn from the vertical axis (which measures dissimilarity).  The horizontal distance between two observations is not informative, and horizontal proximity should be never be confused with similarity.  For example, in the sample dendrogram above, observation 6 is equally as similar to observation 5 as to observation 9 (dissimilarity = 2.0).  Finally, it is important to note that the assumption of hierarchical structure within the data may be incorrect.  If the "true" clusters are not nested, so that the best division into three groups does not result from taking the best division into two groups and splitting up one of the groups, then the resulting clusters will lack validity.  Determining whether or not your data can be represented hierarchically should be considered before performing hierarchical clustering

# Hierarchical Clustering in Supporting Figure 6

Now that you know more about the principles of hierarchical clustering, we will delve into how this clustering method is illustrated in Supporting Figure 6.  There are two dendrograms shown in Supporting Figure 6, one at the very top and one to the left of the microarray heat map.  The top dendrogram clusters tumor samples (shown along the top of the microarray heat map).  Branches of this dendrogram are colored purple, blue, pink, red, green, gray, or black.  The left dendrogram clusters genes (shown to the right of the microarray heat map), and all branches are colored black.

Within the paper, Sørlie et al. mention that both genes and tumor samples were clustered using correlation-based dissimilarity metrics.  Therefore, the pattern with which gene specific microarray expression varied across all tumor samples determined the gene cluster while tumor samples were clustered according to distinct variations in gene expression patterns.

<img src="{{ page.root }}/fig/supp-fig-6-dendrogram.png" alt="The dendrogram showing array clustering at the top of Supporting Figure 6" />

Clustering the tumor samples led to five tumor subtypes, which have also been discussed in the authors' [previous work](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC58566/); this course does not aim to explore the selection of these five subtypes.  Parent branches for the five cancer subtypes are colored in black.  The other branches in this dendrogram are colored according to their correlation with the average expression value of each gene within a tumor subtype.  This was accomplished by creating a subtype profile consisting of the average expression for each of the 551 genes within each subtype.  Branches were colored if the tumor sample was strongly correlated with subtype profile while tumor samples with low correlations are shown in gray.

It is clear from Supporting Figure 6 that hierarchical clustering played a major role in the definition of cancer subtypes and in clustering genes.  As this clustering method forms the backbone of the conclusions reached later in this paper, examining the details of the methodology is critical to reproducing both Supporting Figure 6 and the work of Sørlie et al.  We will expand on these specifics in later sections, but the information you learned in this section should supply you with the background knowledge needed for later lessons.

# Examination Questions

> ## 1
> Refer to the figure below.  How many cluster would exist if you made a horizontal cut at similarity = 0.75?  Which observations would fall into each cluster?
> <img src="{{ page.root }}/fig/hier-cluster-sample.png" alt="A sample dendrogram showing hierarchical clustering of hypothetical data" />
{: .challenge}

> ## 2
> Refer to the array dendrogram for Supporting Figure 6 that is recreated in the lesson above.  Which of the five colored cancer subtypes has the most homogeneous tumors?
{: .challenge}
