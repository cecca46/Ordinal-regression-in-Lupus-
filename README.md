## Repository Structure

This repository contains the following files: 

`main.r`: Main function (in-line code) that follows the analysis pipenile developed to download, normalize, merge and run ordinal regression on the gene expression data for 21 healthy samples and 30 Systemic Lupus Erythematous samples. This function will create a RESULT folder to store all the results (plots included) from running the analysis.

`RMA.r`: Function to perform RMA normalization separately on the 2 different platforms contained in the dataset and to merge the normalized results into a single expression matrix.

`PCAplots.r`: Series of wrapper functions to inspect the data and look for batch effect problems.

`QuantileNorm.r`: Function to quantile normalize the expression matrix.

`Combat.r`: Function to mitigate the bach effect introduced by the use of different platforms.

`OrdinalRegression.r`: Function to perform ordinal regression on the normalised and batch corrected gene expression matrix. It returns a list of significant genes (default significance thresold is 1e-06) ordered according to the slope of their beta coefficients. High positive beta means that the gene progressively upregulates with disease progression whereas small negative beta means that the gene progressively downregulated with disease progression. 

`MergedDatasetsInfo.csv`: Clinical information, publically available from [Nephroseq](https://www.nephroseq.org/resource/login.html)




### Identification of genes associated to severity of kidney damage in Systemic Lupus Erythematous

The search for disease-associated genes and biomarkers relies on the discovery of statistical links between gene expression and disease phenotype. In most methods, clinical metrics are treated as binary data (e.g., disease vs. control). However, in many cases, even the most basic clinical data provide a richer description of the disease process. Systematic integration of ordinal clinical metrics (such as Tumor stage, neurodegeneration levels or more in general measurements of disease progression) with gene expression data may lead to identifying a subset of the genes that play a critical role in disease progression. Once experimentally validated, these genes could be important candidates for therapeutic targets.
To develop an approach that can take advantage of information on the severity of the disease, we analyzed gene expression data from glomeruli tissues in Systemic Lupus Erythematous (SLE) and Living Donors (LD). Patients who suffer of SLE can be classified in five categories that, identified using the estimated filtration rate (eGFR) numbers, depict the severity and pattern of decline of the kidney filtering function.


### Data 

The dataset used, GEO accession number **GSE104948**, contains glomerular transcriptome from European Renal cDNA Bank subjects and living donors. It is composed of 196 samples. Out of them all, we collected 30 samples of Systemic Lupus Erythematous (SLE) and 21 Living Donors (LD).  Clinical metrics corresponding to disease phenotype for these samples are publicly available. Specifically, we categorized the patients based on the reported severity of kidney demage.
The gene expression compendium was normalised using Robust Multi Array (RMA) normalisation.
After RMA normalization, as the number of genes differ due to the two platform used, we merge the genes from HG-U133A (13768) and HG-U133_Plus_2 (22048). The final results yield a gene expressium compendium of 13768 genes and 51 samples.
After merging, we inspect PCAs looking for possible clusters which denote a batch effect problem.

![](https://github.com/cecca46/Ordinal-regression-in-Lupus-/blob/master/PLOT/PCAbyGroup.pdf)
