---
title: "Macrophage Gene Signatures"
author: "Giuliano Netto Flores Cruz"
date: "October 9, 2018"
bibliography: library.bib
csl: abnt-ufrgs-initials.csl
output: 
  word_document: 
    keep_md: yes
---





## Macrophage Genomic Integrative Analysis

  Gene expression integrative analysis is a computationally-expensive and rather complex task. First, one must integrate different technologies (_i.e_ microarray and RNA sequencing) from several distinct platforms [@Walsh2015]. Second, although greater number of samples yields higher statistical power, the potential confounding factors must be taken into account. When it comes to lung injury, a wide range of animal models and human conditions have already been tested, and their respective data sets should be treated with care. For instance, the widely used bleomycin-triggered IPF model shows enrichment of traditionally M1-associated genes at very early stages [@Bauer2015]. Fungal infection models, on the other hand, show divergent genomic markers with potentially protective roles associated to genes from the M2 spectrum [@Bhatia2011; @Margalit2015]. Both cases, though, may lead to pulmonary fibrosis through macrophage activity [@Iwasaki2016; @Gieseck2018; @Wynn2016]. Other challenges include the adequacy of sample sizes, pre-processing techniques, statistical analysis, validation, and experimental design, as well as the lack of a comprehensive framework for the execution of gene expression meta-analyses [@Ramasamy2008]. Of note, the term "meta-analysis" refers to when the researcher analyzes each data set separately and draws conclusions based on the final statistical results combined, whereas "cross-platform normalization" is used to describe the integration of raw data from multiple sources for combined downstream analysis [@Walsh2015]. Here, we use "integrative analysis" to denote both terms interchangeably as not all data sets analyzed are suitable for merging.
  
  The applicability of integrative analysis for elucidating reproducible macrophage gene signatures and even predicting clinical outcome based on the enrichment of those signatures has been previously tested [@Becker2015a]. Using data from human-derived macrophages challenged with two sets of activation stimuli, namely "classical" (IFN-g + LPS; TNF-a) and "alternative" (IL-4; IL-13), the authors were able to establish prognostic values in diverse clinical settings such as viral infections and asthma. Noteworthy, however, is the gene signatures being relatively limited by the M1 versus M2 paradigm, which hinders interpretation at the microenvironment level. After all, how to understand the heterogeneity within both groups of cells and, furthermore, how to address potentially overlapping gene signatures from macrophage subsets that were overlooked (or that are yet to be described)? How comprehensive should an integrated analysis be to assure wide reproducibility of gene expression patterns?
  Recently, an elegant work integrated several data sets from wide variety of human diseases as well as several mouse strains within the context of LPS exposure[@Buscher2017]. Surprising was not the high level of variance across strains, but the ability to nevertheless infer the degree of polarization of macrophages from a vast sort of patients .
  
  
  
  


## Macrophage Gene Signatures

  Macrophages have been demonstrated to develop highly complex activation profiles in a diverse set of microenvironments [@Ginhoux2016]. Specifically, these cells are key players in conditions related to idiopathic pulmonary fibrosis (IPF) [@Bauer2015; @Wynn2011; @Venosa2016a]. Although many genetic markers are known to play important roles within the macrophage biology context, defining a robust set of gene signatures for the currently known phenotypic subsets remains a challenging task [@Martinez2014]. Recent cytometric and genomic approaches have revealed major limitations regarding the classic M1 versus M2-polarization model, which is no longer accepted as suitable to explain the biological dynamics of macrophage response [@Martinez2014].
  
  In the pursuit of standardization towards reproducible research, back in 2014 a group of specialists suggested nomenclatures and experimental guidelines for the macrophage activation profiles well-established by then [@Murray2014]. However, the recent abundance of genomic data has challenged the classical protein-level techniques used to sort macrophage subsets and therefore novel approaches emerged. In the same year that Murray's paper was published, a multi-center work attributed a much higher heterogeneity to macrophages through machine learning algorithms applied to single-cell transcriptome analysis [@Xue2014]. 
  
  Assessing the transcriptomes from almost 300 _in-vitro_ stimulated human macrophages, Xue and colleagues used weighted gene co-expression network analysis (WGCNA) to identify 49 coexpression modules, each of which ranging from less than 30 to over 800 distinct genes of size [@Xue2014]. Based on Pearson correlation, WGCNA defines gene clusters, known as _transcriptional modules_, which present specific co-expression patterns across each treatment condition [@Langfelder2008]. As an example, these modules can then be used to visualize the comprehensiveness of the M1 versus M2 model. As noted by the authors, stimuli not M1- or M2-associated showed prominent patterns consistent with a rather dynamic spectrum model of cell activation.
  
  In order to achieve high sensivity for potential macrophage phenotypes, in this study the 49 transcriptional modules produced by Xue and colleagues were used as gene sets for further analysis. The Figure 1 shows the distribution of number of genes across the different modules.
  
![](Defining_Macrophage_Gene_signatures_files/figure-docx/Figure-1.png)<!-- -->


  
  
  
  
  
  
  
  
  
  