---
title: "Development"
bibliography: config/library.bib
csl: config/abnt-ufrgs-initials.csl
output: 
  word_document:
    reference_docx: config/template.docx
    keep_md: true
---




# IPF animal models at the genomic level

  The bleomycin-induced IPF animal model is widely used to understand lung fibrosis pathology, regardless of its limited capability of mimicking the actual human disease [@Moeller2008]. Bauer and colleagues studied this question when comparing microarray data from one hundred lung samples from IPF patients with rat lungs sampled several time points after bleomycin exposure [@Bauer2015]. Although they were able to identify disease-relevant translational gene markers, the point of highest rat-human gene expression commonality was at day 7 after rat lung aggression. The authors suggest that these gene signatures can be used to identify IPF patients and to stratify these according to disease severity. Here, we reanalyze their data in order to further understand time course patterns in gene expression and their relation with cellular pathological activity.
  
  Using the arrayQualityMetrics R package, we were able to identify five outliers based on overall expression data and these were removed from further analysis - although the original paper indicated 17 outliers [@Bauer2015]. As morphological and cytometric analyses indicate that bleomycin model shows time-related pathological events [@Venosa2016; @Izbicki2002], we cut the original data into 5 supposedly divergent phases: namely, "Healthy" for untreated samples, "Injury" for rats killed at early exposure time points (3 and 7 days), "Early Fibrosis" (day 14), "Late Fibrosis" (days 21 through 28), and "Healing" (days 42 through 56). This generally arbitrary classification successfully showed descriptive gene expression patterns in principal component analysis - Figure 4. 
  
![Figure 4 - Principal Component Analysis of gene expression from Bauer and colleagues data (2016).](Development_files/figure-docx/pca_bauer2015.png)
  
  The first three principal components separate control and bleomycin samples. Most importantly, the samples at early time points - _i.e._ when the recent aggression induces major inflammatory responses - fall well separated from later times as well as from control samples. Notably, the injury-labeled samples fall further beyond others, followed by early fibrosis, late fibrosis, and finally healing-labeled and control samples - almost mimicking the actual time course experimental design and suggesting the impact of measurement times on IPF animal model assessment. Even though the authors indicate that day 7 (injury phase) is the point with maacross time must not be neglected, especially regarding assessment of IPF candidates and early-diagnosis procedures. 
  
  Once disease is installed, one may expect reproducible gene expression patterns, even though this understanding is hindered by the idiopathic characteristic of the condition. However, those patients with developing histopathological characteristics that are yet to be diagnosed as typical IPF may not reflect such genomic patterns. Furthermore, it has been reported that gene signatures differ significantly across IPF patients with progressive and stable conditions [@Boon2009]. Thus, it is important to note the importance of longituginal studies regarding genomic signatures as these may prove themelves helpful when predicting disease onset, progression, and stabilization.

  Regarding macrophage biology, several approaches are possible to assess their dynamics in animal models. As previously noted, Venosa and colleagues were able to describe macrophage activity in an animal model of IPF induced by nitrogen mustard [@Venosa2016]. Using data from cytometric, qRT-PCR, and other nonmolecular assays, the authors demonstrated the inflammatory profile of infiltrating cells at early time points, while anti-inflammatory and healing profiles where dominant at later times. The proposed kynetics related well with gene expression patterns, although high-throughput technologies were not used. 
  
  Here, the first macrophage characterizaition addresses the M1 versus M2 paradigm. As proposed by Buscher and colleagues, the Polarization Factor Ratio (PFR) is intended to describe the degree of macrophage polarization towards M1 or M2 spectra [@Buscher2017]. As a simple model, it derives from the expression levels of M1- and M2-markers. Using Bauer's data, two-way mixed design anova with multilevel modeling did not reveal any time-dependent patterns (p > 0.05). However, it did reveal significant differences across treatment groups (p < 0.005), which was further confirmed by t-student (p < 0.01) and Exact Wilcoxon-Mann-Whitney tests (p < 0.01). Figure 5 shows that time courses for both groups fail to trend any direction significantly - perhaps due to noisy data points. However, the courves do not overlap completely, and the boxplots illustrate the distribution differences. In the original paper, Buscher and colleagues demonstrated higher prediction power for the PFR built over IL12b and Arginase 1 expression levels - when comparing to the same score constructed with inducible nitric oxide synthase (iNOS - NOS2 gene, also M1-related) instead of the cited interleukin. Here, the effect size comparison seemed te shifted, and the PFR (iNOS/Arg1) showed better separation between bleomycin- and PBS-treated animals. Taken together, these data initialy indicate that PFR is capable to reflect a slight overall macrophage polarization towards an M1 spectrum in the given pulmonary fibrosis animal model. As a two-gene model, however, such a conclusion is clearly an oversimplification of macrophage and IPF biology.
  
  ![Figure 5 - PFR with IPF animal model data from Bauer and colleagues. ** p = 0.007983; *** p = 4.335e-07 (Exact Wilcoxon-Mann-Whitney Test).](Development_files/figure-docx/pfr_bauer2015.png)
  
  In order to further characterize the time-course differences in overall gene expression of bleomycin- and PBS-treated rats, we performed differential expression analysis on Bauer's data using a two-step statistical method which is especially designed for time course data and is implemented in the Bioconductor package, maSigPro [@Conesa2006]. First, the procedure fits a global model for all genes in a given dataset. Then, it applies step-wise regression as a means of variable selection so that it can detect significant differences across study groups and consistent expression profiles across time points. Although the dataset tested contained eight time points for each group, here we relied on a cubic regression model, Higher polynomial degrees have yielded high noisy fitting and possibly high rates of type I error (data not shown), which is somehow expected when working with overly complex polynomials [@Conesa2006]. One could argue the use of splines, but these are not available in the maSigPro prackage. To avoid underfitting, the time points 42 and 56 were excluded from this analysis. 
  
  MaSigPro also conducts cluster analysis with several strategies to identify similar expression profiles across time. Using the algorithm from mclust R package (Normal Mixture Modeling for Model-Based Clustering, Classification, and Density Estimation - available on CRAN), it can group the time courses into an optimal k number of clusters based on finite normal mixture modeling [@Scrucca2016]. However, hierarchical clustering showed similar results (with k = 9) and these were taken for further analysis. Figure 6 shows the nine clusters produced by maSigPro and their time course profiles. The dashed lines represent the fitted models, while solid lines show the true median expression values. 
  
  ![Figure 6 - Hierarchical clustering reveals time course expression profiles of differentially expressed genes in IPF animal model.](Development_files/figure-docx/hclust_bauer2015.png)
  
  
  Notably, there are fairly similar clusters ( _e.g._ clusters 6 and 7). Here, however, we are particularly interested in those which show overexpression either at early or later time points, as these may be representative of eventual macrophage polarization patterns. For instance, one may speculate cluster 1 to be filled with genes related to the M1 spectrum, while cluster 4 seems to follow a transitory course and, finally, cluster 9 may represent an M2-polarized enrironment. Cleary, these are limited speculations once overall expression patterns greatly overlook description of macrophage dynamics. Therefore, deeper characterization required assessment of gene profiles on a case-by-case basis.
  
  Based on recent literature, we sought to find immune-associated genes that have been previously reported as macrophage- or IPF-related. Figure 7 shows the time course profiles for a selected set of significantly differentiated chemokines. 
  
  
  

