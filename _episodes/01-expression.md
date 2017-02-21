---
title: "Introduction to Gene Expression"
teaching: 10
exercises: 10
questions:
- "What is a microarray?"
- "What biological concepts underlie microarray analysis?"
- "How are the results of microarray analysis visualized?"
objectives:
- "Understand the basic principles behind microarray analysis."
- "Learn how to interpret heat maps."
- "Interpret the heat map in Supporting Figure 6."
keypoints:
- "Microarray heat maps can be used to compare and identify gene expression patterns."
---

Understanding the principles of microarray analysis is essential for interpreting Supporting Figure 6. To produce this figure, Sørlie et al. used the results of DNA microarray analysis to both produce the heat map and provide input for the hierarchical clustering represented by the tree-like structures (or dendrograms) shown in the figure. This lesson will teach you the underlying concepts of microarray analysis and how these results are visualized in a heat map. At the end of this lesson you will understand how microarray analysis provides quantification of gene expression as well as how to interpret a heat map. You will learn about hierarchical clustering and dendrograms in the [next lesson]({{ page.root }}/02-hierarchical/).

# Microarray Analysis

Microarray analysis is an important tool in clinical medicine that allows us to measure the expression levels of thousands of genes and study the effects of diseases and treatments on gene expression. The microarray itself is an array on which single strand DNA, RNA, or polynucleotide probes are fixed in a specific manner. Since each cell in the array corresponds to a known gene, the microarray serves as a genetic reference standard to analyze the gene expression pattern of unknown samples.

DNA Microarray analysis is performed by first extracting mRNA from the unknown tissue and using reverse transcriptase to synthesize its complementary DNA (cDNA). Fluorescent dye is introduced into the cDNA as it is being synthesized and this fluorescent cDNA is then run over the microarray. If the unknown tissue contains genes corresponding to a probe in the microarray, the fluorescent cDNA will bind to that probe. After allowing time for binding and washing the microarray to remove unbound cDNA, the microarray is illuminated by a laser that causes the bound fluorescent cDNA molecules to fluoresce proportional to the bound quantity. As the strength of the signal from each cell is dependent on the quantity of cDNA, the signal collected from each cell can be used to estimate the expression level of a gene. In practice, to standardize this quantification when using multiple microarrays, experimental samples are compared with a known cell line on the same microarray, using different frequency fluorescent dyes to distinguish the known and unknown sample. More information on microarray analysis is described by [Tarca, Romero, and Draghici](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2435252/), by [Babu](http://www.mrc-lmb.cam.ac.uk/genomes/madanm/microarray/), and by [Schena](http://www.wiley.com/WileyCDA/WileyTitle/productCd-0471414433.html), as well as in this [short video clip](https://www.youtube.com/watch?v=VNsThMNjKhM).

# Visualizing the Microarrays as a Heat Map

The results of microarray analysis are often visualized as a heat map consisting of a rectangular array of blocks with varying colors and intensities. Here we present a sample heat map visualizing the hypothetical expression of genes 1-10 (rows) among arrays 1-10 (columns).

<img src="{{ page.root }}/fig/sample-heat-map.png" alt="A sample heat map showing hypothetical data" />

The color of each cell represents the expression level of a gene for an array. Shades of red are used to represent overexpressed or upregulated genes while shades of green are used to represent under-expressed or downregulated genes. Cells colored in black indicate there was no change in expression. A gradient color bar on the upper right hand side of the above figure shows the how color intensity is used to quantify differences in gene expression with brighter reds indicating more upregulated gene expression and brighter green indicating more downregulated gene expression.

This course focuses on the heat map shown in Supporting Figure 6 by Sørlie et al. which illustrates the various gene expression patterns in 122 tissue samples (columns) for over 500 genes. The same color conventions are used so that intense red reflects overexpressed genes while intense green reflects underexpressed genes. Black color indicates little deviation from the median gene expression across all samples, and gray represents missing data. The gradient intensity bar is shown in the top right and ranges from 1 to 5.6 times the median intensity in both directions. Also shown in Supporting Figure 6 are tree-like visualizations called dendrograms which will be discussed in the [next lesson]({{ page.root }}/02-hierarchical/). Now that you understand how to interpret a heat map and how microarray analysis provided quantification of gene expression levels, you are ready to further explore how this data was used for clustering.

# Examination Questions

Analyze the gene expressions in the microarray heat map shown above and answer the following questions:

> ## 1
> Which genes are most overexpressed in array 4 and array 1?
{: .challenge}

> ## 2
> Which genes are most underexpressed in array 4 and array 8?
{: .challenge}
