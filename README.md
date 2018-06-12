## Repository structure


This repository contains the following files: 

`main.r`: Main function (in-line code) that follows the analysis pipenile developed to download, normalize, merge and run ordinal regression on the gene expression data for 21 healthy samples and 30 Systemic Lupus Erythematous samples. This function will create a `RESULT` folder to store all the results from running the analysis. It also writes any generated plot to a `PLOT` folder if exists (it will create one otherwise).

`RMA.r`: Function to perform RMA normalization separately on the 2 different platforms contained in the dataset and to merge the normalized results into a single expression matrix.

`PCAplots.r`: Series of wrapper functions to inspect the data and look for batch effect problems.

`QuantileNorm.r`: Function to quantile normalize the expression matrix.

`Combat.r`: Function to mitigate the bach effect introduced by the use of different platforms.

`OrdinalRegression.r`: Function to perform ordinal regression on the normalised and batch corrected gene expression matrix. It returns a list of significant genes (default significance thresold is 1e-06) ordered according to the slope of their beta coefficients. High positive beta means that the gene progressively upregulates with disease progression whereas small negative beta means that the gene progressively downregulated with disease progression. 


`MergedDatasetsInfo.csv`: Clinical information, publically available from [Nephroseq](https://www.nephroseq.org/resource/login.html)




### Identification of genes associated to severity of kidney damage in Systemic Lupus Erythematous

The search for disease-associated genes and biomarkers relies on the discovery of statistical links between gene expression and disease phenotype. In most methods, clinical metrics are treated as binary data (e.g., disease vs. control). However, in many cases, even the most basic clinical data provide a richer description of the disease process <sup>1</sup>. Systematic integration of ordinal clinical metrics (such as Tumor stage, neurodegeneration levels or more in general measurements of disease progression) with gene expression data may lead to identifying a subset of the genes that play a critical role in disease progression. Once experimentally validated, these genes could be important candidates for therapeutic targets.<sup>1</sup>
To develop an approach that can take advantage of information on the severity of the disease, we analyzed gene expression data from glomeruli tissues in Systemic Lupus Erythematous (SLE) and Living Donors (LD). Patients who suffer of SLE can be classified in five categories that, identified using the estimated filtration rate (eGFR) numbers, depict the severity and pattern of decline of the kidney filtering function.


### Data 



The dataset used, GEO accession number **GSE104948**, contains glomerular transcriptome from European Renal cDNA Bank subjects and living donors. It is composed of 196 samples. Out of them all, we collected 30 samples of Systemic Lupus Erythematous (SLE) and 21 Living Donors (LD). 



|                | Platform HG-U133A |Platform HG-U133_Plus_2 | 
| -----------    | ----------------- | ---------------------- |
|  **Condition** |  
|      LD        |         3         |           18           |
|      SLE       |         30        |            0           |


Clinical metrics corresponding to disease phenotype for these samples are publicly available. Specifically, we categorized the patients based on the reported severity of kidney demage.
The gene expression compendium was normalised using Robust Multi Array (RMA) normalisation.
After RMA normalization, as the number of genes differ due to the two platform used, we merge the genes from HG-U133A (13768) and HG-U133_Plus_2 (22048). The final results yield a gene expressium compendium of 13768 genes and 51 samples.
After merging, we inspect PCAs looking for possible clusters which denote a batch effect problem.


<p align="center">
    Samples grouped by condition
    <img src="https://github.com/cecca46/Ordinal-regression-in-Lupus-/blob/master/PLOT/PCAbyGroup.png" width="700" height="700">
</p>




<p align="center">
    Samples grouped by platform
    <img src="https://github.com/cecca46/Ordinal-regression-in-Lupus-/blob/master/PLOT/PCAbyPlatform.png" width="700" height="700">
</p>


Clearly from the latter plots, there is strong batch effect driving the samples from the two different platforms far apart. The 3 healthy samples shown in the pictures, wrongly cluster with SLE samples as they were analysed using the same chip.
We then tried to reduce platform specific effect using Combat, after quantile normalization on the expression matrix, while adjusting for the sample condition (SLE-LD).

<p align="center">
    Samples grouped by condition after Combat
     <img src="https://github.com/cecca46/Ordinal-regression-in-Lupus-/blob/master/PLOT/PCAbyGroupCombat%20.png" width="700" height="700">
</p>

<p align="center">
    Samples grouped by platform after Combat
    <img src="https://github.com/cecca46/Ordinal-regression-in-Lupus-/blob/master/PLOT/PCAbyPlatformCombat.png" width="700" height="700">
</p>



In the latter plots, the 3 living donor samples from platform HG-U133A cluster with the living donors from platform HG-U133_Plus_2. From these plots seem that the differences among the samples is now driven by biological differences, namely the SLE vs LD condition, and not by technical variation. 

### Ordinal analysis

Using these samples, we sought to identify and rank genes that are linked strongly to the severity of SLE progression. For this purpose, we adapted a statistical model known as ordinal regression to integrate real-valued expression data with the stage of disease. In this regression model, the values of the response variable have an ordinal relationship (stage 1 < stage 2 < stage 3 < stage 4 < stage 5). The ordinal regression model makes no assumptions about the relative quantitative value of the scale; this is essential as clinical phenotypes have a qualitative and not quantitative nature (the distance between one class and another can be unknown) and a severity score of four does not represent twice the severity of a score of two. To increase the statistical power of the analysis, we grouped the samples according to their proximity of SLE progression in the way shown in the table:


|     Ordinal Category	   |      Conditions   |       Sample size      | 
| -------------------------| ----------------- | ---------------------- |
|     Ordinal Category 1   |        Control    |            21          |
|     Ordinal Category 2   |    Stage 1 or 2   |            18          |
|     Ordinal Category 3   |    Stage 3A or 3B |            7           |
|     Ordinal Category 4   |    Stage 4 or 5   |            5           |


Using this approach, we identified 506 genes (termed phenotype-associated genes, PAGs) whose transcriptional dysregulation was significantly associated with renal pathological severity (Wald test, p-value of all the fitted parameters <1e-6). The magnitude of the β parameter, which is fitted for each gene, correlates to the rate by which the expression of a gene is altered with the increase of severity in disease progression. We identified 279 consistently upregulated PAGs, and 227 consistently downregulated PAGs due to progressive kidney decline function.
Using the significance scores calculated using Ordinal regression, we performed gene set enrichment using the tool Piano on the whole gene list. 


<p align="center">
    Piano Distinct Directional Down 
    <img src="https://github.com/cecca46/Ordinal-regression-in-Lupus-/blob/master/PLOT/DistinctDOWN.png" width="700" height="700">
</p>

<p align="center">
    Piano Distinct Directional Up 
    <img src="https://github.com/cecca46/Ordinal-regression-in-Lupus-/blob/master/PLOT/DistinctUP.png" width="700" height="700">
</p>

<p align="center">
    Piano Mixed Directional Down 
    <img src="https://github.com/cecca46/Ordinal-regression-in-Lupus-/blob/master/PLOT/MixedDOWN.png" width="700" height="700">
</p>

<p align="center">
    Piano Mixed Directional Up
    <img src="https://github.com/cecca46/Ordinal-regression-in-Lupus-/blob/master/PLOT/MixedUP.png" width="700" height="700">
</p>


<p align="center">
    Piano excluded directionality (non directional)
    <img src="https://github.com/cecca46/Ordinal-regression-in-Lupus-/blob/master/PLOT/NonDirectional.png" width="700" height="700">
</p>


The following Venn Diagrams show the number of common enriched processes among the the non directional, directional up, mixed up and non directional, directional down, mixed down respectively:

<p align="center">
     Non directional, directional up, mixed up intersection
    <img src="https://github.com/cecca46/Ordinal-regression-in-Lupus-/blob/master/PLOT/UPIntersection.png" width="500" height="500">
</p>

**Upregulated common biological processes**: *immune system process, innate immune response,neutrophil degranulation,protein binding,extracellular region,apoptotic process,inflammatory response,leukocyte migration,regulation of immune response,defense response to virus,viral process,neutrophil chemotaxis,immune response,regulation of megakaryocyte differentiation,interferon-gamma-mediated signaling pathway,lipopolysaccharide-mediated signaling pathway,negative regulation of viral genome replication,type I interferon signaling pathway,regulation of inflammatory response,signaling pattern recognition receptor activity,CENP-A containing nucleosome assembly,response to virus,negative regulation of gene expression, epigenetic,phagocytosis, engulfment,B cell receptor signaling pathway,toll-like receptor TLR1:TLR2 signaling pathway,cellular response to triacyl bacterial lipopeptide,positive regulation of interleukin-8 secretion,chromatin silencing at rDNA,DNA replication-dependent nucleosome assembly,telomere organization,antigen processing and presentation of exogenous peptide antigen via MHC class I,2'-5'-oligoadenylate synthetase activity,myeloid dendritic cell activation involved in immune response,regulation of B cell apoptotic process*


<p align="center">
    Non directional, directional down, mixed down intersection
    <img src="https://github.com/cecca46/Ordinal-regression-in-Lupus-/blob/master/PLOT/DOWNIntersection.png" width="700" height="500">
</p>

**Downregulated common biological processes**: *homophilic cell adhesion via plasma membrane adhesion molecules, regulation of cellular ketone metabolic process*


As the slope of the fitted lines indicates the rate by which the expression of a gene changes with progressive SLE, we hypothesized that genes with larger slopes would be likely to have important roles in the disease. Thus, we ranked the genes according to the magnitude of the β parameter.

|       Rank             |     Gene name	        |    β parameter    |       pvalue           | 
| ---------------------  | -------------------------| ----------------- | ---------------------- |
|       1                |     TCF3                 |    11.15402638    |            3.37E-17    |
|       2                |     CENPI                |    8.914415045    |            1.73E-16    |
|       3                |     GUCY2C               |    8.344231823    |            1.11E-15    |
|       4                |     RPL27                |    7.771460644    |            1.88E-15    |
|       5                |     MVP                  |    7.453617767    |            2.48E-18    |
|       6                |     ZDHHC3               |    7.306487705    |            1.87E-17    |
|       7                |     ZDHHC18              |    6.839532265    |            3.99E-18    |
|       8                |     TINF2                |    6.760683612    |            8.26E-17    |
|       9                |     ABHD2                |    6.617354198    |            2.23E-16    |
|       10               |     RBL1                 |    6.527797283    |            7.63E-16    |
|       11               |     RNASEH2A             |    6.22245656     |            2.84E-19    |
|       12               |     TIMELESS             |    6.063382403    |            2.04E-17    |

### xCell enrichment 

xCell is a webtool that performs cell type enrichment analysis from gene expression data for 64 immune and stroma cell types. It is a gene signatures-based method learned from thousands of pure cell types from various sources.<sup>2</sup> We provided as input to xCell our merged gene expression matrix composed by 51 samples and 13768 genes.

<p align="center">
    <img src="https://github.com/cecca46/Ordinal-regression-in-Lupus-/blob/master/PLOT/xCelHeatMap.png" width="700" height="700">
</p>

From the xCell output, among others general knwon trends can be observed for Macrophages and Monocytes; they are consistently upregulated in SLE condition with respect to healty controls. To assess the significance of these results, we correlated the enrichment scores for Macrophages and Monocytes to eGFR numbers. 
As the activity of SLE is measured by eGFR we assumed a negative correlation between eGFR levels and enrichment scores obtained for Macrophages and Monocytes using xCell. The lower eGFR, the higher the enrichment score for Macrophages. We tested this hyphotesis using spearman correlation between eGFR and enrichment scores.

Correlation between eGFR number and Macrophages yields a rho coefficient of -0.545 with a pvalue of 3.519e-05 while correlation between eGFR number and Monocytes yields a rho coefficient of -0.451 with a pvalue of 0.0008873.

Interestingly, the strongest correlation coefficent was found to be between eGFR number and Myocytes (rho = 0.564, pvalue = 1.618e-05) suggesting that enrichment in Myocytes corresponds to a better SLE prognosis. 

### License

Distributed under the GNU GPLv3 License. See accompanying file LICENSE.txt or copy at http://www.gnu.org/licenses/gpl-3.0.html.


### Reference

1. [Identifying therapeutic targets by combining transcriptional data with ordinal clinical measurements, Pirhaji et al. 2017](https://www.nature.com/articles/s41467-017-00353-6#Sec10)

2. [xCell: digitally portraying the tissue cellular heterogeneity landscape, Dvir Aran et al. 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5688663/)

