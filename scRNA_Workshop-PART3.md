---
title: "Single Cell RNAseq Part 3"
author: "Bioinformatics Core"
output:
    html_document:
      keep_md: TRUE
---
## Load libraries

```r
library(Seurat)
library(ggplot2)
```

## Load the Seurat object

```r
load(file="pre_sample_corrected.RData")
experiment.aggregate
```

```
## An object of class Seurat 
## 12811 features across 2681 samples within 1 assay 
## Active assay: RNA (12811 features)
```

```r
experiment.test <- experiment.aggregate
set.seed(12345)
rand.genes <- sample(1:nrow(experiment.test), 500,replace = F)
mat <- as.matrix(GetAssayData(experiment.test, slot="data"))
mat[1:5,1:5]
```

```
##        ACTCTAATGTGGGTATG-UCD_Adj_VitE AGGCTGGTCAATCACAC-UCD_Adj_VitE
## Xkr4                        0.0000000                       0.000000
## Sox17                       0.0000000                       0.000000
## Mrpl15                      0.0000000                       0.000000
## Lypla1                      0.0000000                       1.096085
## Tcea1                       0.7540035                       1.096085
##        ATGACTAGCACATGACT-UCD_Adj_VitE AAGCGTCGTCTCTAAGG-UCD_Adj_VitE
## Xkr4                                0                       0.000000
## Sox17                               0                       0.000000
## Mrpl15                              0                       0.000000
## Lypla1                              0                       2.098709
## Tcea1                               0                       0.000000
##        ACATCGGGTCCATGCTC-UCD_Adj_VitE
## Xkr4                                0
## Sox17                               0
## Mrpl15                              0
## Lypla1                              0
## Tcea1                               0
```

```r
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
experiment.test@meta.data %>% count(batchid, orig.ident) # All 808 cells of ADJ_VITE group are in Batch2
```

```
## # A tibble: 3 x 3
##   batchid orig.ident        n
##   <chr>   <fct>         <int>
## 1 Batch1  UCD_Supp_VitE   947
## 2 Batch1  UCD_VitE_Def    926
## 3 Batch2  UCD_Adj_VitE    808
```

```r
# adding extra signal to a random set of 500 genes in Batch2
mat[rand.genes,experiment.test$batchid=="Batch2"] <- mat[rand.genes,experiment.test$batchid=="Batch2"] + 0.18
experiment.test = SetAssayData(experiment.test, slot="data", new.data= mat )
```



## Exploring Batch effects 3 ways, none, Seurat [vars.to.regress] and COMBAT

First lets view the data without any corrections

## PCA in prep for tSNE

ScaleData - Scales and centers genes in the dataset. 

```r
?ScaleData
experiment.test.noc <- ScaleData(object = experiment.test)
```

```
## Centering and scaling data matrix
```

### Run PCA

```r
experiment.test.noc <- RunPCA(object = experiment.test.noc, verbose = F)
#DimPlot(object = experiment.test.noc, group.by = "batchid", reduction = "pca")
#DimPlot(object = experiment.test.noc, group.by = "batchid", dims = c(2,3), reduction = "pca")
```

PCA Elbow plot to determine how many principal components to use in downstream analyses.  Components after the "elbow" in the plot generally explain little additional variability in the data.


```r
ElbowPlot(experiment.test.noc)
```

![](scRNA_Workshop-PART3_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

We use 10 components in downstream analyses. Using more components more closely approximates the full data set but increases run time.

### TSNE Plot

```r
pcs.use <- 10
experiment.test.noc <- RunTSNE(object = experiment.test.noc, dims = 1:pcs.use)
DimPlot(object = experiment.test.noc,  group.by = "batchid")
```

<img src="scRNA_Workshop-PART3_files/figure-html/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />




## Correct for sample to sample differences (seurat)

Use vars.to.regress to correct for the sample to sample differences and percent mitochondria

```r
experiment.test.regress <- ScaleData(object = experiment.test, 
                    vars.to.regress = c("batchid"), model.use = "linear")
```

```
## Regressing out batchid
```

```
## Centering and scaling data matrix
```

```r
experiment.test.regress <- RunPCA(object =experiment.test.regress, verbose = F)
#DimPlot(object = experiment.test.regress, group.by = "batchid", reduction.use = "pca")
```

### Corrected TSNE Plot

```r
experiment.test.regress <- RunTSNE(object = experiment.test.regress, dims = 1:pcs.use)
DimPlot(object = experiment.test.regress, group.by = "batchid", reduction = "tsne")
```

<img src="scRNA_Workshop-PART3_files/figure-html/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />





## COMBAT corrected
https://academic.oup.com/biostatistics/article-lookup/doi/10.1093/biostatistics/kxj037

```r
library(sva)
?ComBat
m = as.matrix(GetAssayData(experiment.test))
#m[1:5,1:5]
```

Why are they performing batch correction across the "orig.ident" and not "batchid" that was used in Seurat ScaleData vars to regress?

```r
com = ComBat(dat=m, batch=as.numeric(as.factor(experiment.test$orig.ident)), prior.plots=FALSE, par.prior=TRUE)
```

```
## Found3batches
```

```
## Adjusting for0covariate(s) or covariate level(s)
```

```
## Standardizing Data across genes
```

```
## Fitting L/S model and finding priors
```

```
## Finding parametric adjustments
```

```
## Adjusting the Data
```

```r
experiment.test.combat <- experiment.test
experiment.test.combat <- SetAssayData(experiment.test.combat, new.data = as.matrix(com))
experiment.test.combat = ScaleData(experiment.test.combat)
```

```
## Centering and scaling data matrix
```

Principal components on ComBat adjusted data

```r
experiment.test.combat <- RunPCA(object = experiment.test.combat, verbose =F)
#DimPlot(object = experiment.test.combat, group.by = "batchid", reduction = "pca")
```

### TSNE plot on ComBat adjusted data

```r
experiment.test.combat <- RunTSNE(object = experiment.test.combat, dims = 1:pcs.use)
DimPlot(object = experiment.test.combat, group.by = "batchid", reduction = "tsne")
```

![TSNE plot, ComBat adjusted ](scRNA_Workshop-PART3_files/figure-html/unnamed-chunk-12-1.png)


## Paul ComBat correction on Batchid not Orig.Ident

```r
com_paul = ComBat(dat=m, batch=as.numeric(as.factor(experiment.test$batchid)), prior.plots=FALSE, par.prior=TRUE)
```

```
## Found2batches
```

```
## Adjusting for0covariate(s) or covariate level(s)
```

```
## Standardizing Data across genes
```

```
## Fitting L/S model and finding priors
```

```
## Finding parametric adjustments
```

```
## Adjusting the Data
```

```r
mycombattest <- experiment.test
mycombattest <- SetAssayData(mycombattest, new.data = as.matrix(com_paul))
mycombattest = ScaleData(mycombattest)
```

```
## Centering and scaling data matrix
```

```r
mycombattest <- RunPCA(object = mycombattest, verbose =F)
mycombattest <- RunTSNE(object = mycombattest, dims = 1:pcs.use)
DimPlot(object = mycombattest, group.by = "batchid", reduction = "tsne")
```

![](scRNA_Workshop-PART3_files/figure-html/unnamed-chunk-13-1.png)<!-- -->



#### Question(s)

1. Try a couple of PCA cutoffs (low and high) and compare the TSNE plots from the different methods.  Do they look meaningfully different?


```r
# data
#DimPlot(RunTSNE(experiment.test.noc, dims = 1:5), group.by = "batchid", reduction = "tsne")
DimPlot(RunTSNE(experiment.test.noc, dims = 1:10), group.by = "batchid", reduction = "tsne")  + ggtitle("No correction, pca 10")
```

![](scRNA_Workshop-PART3_files/figure-html/unnamed-chunk-14-1.png)<!-- -->

```r
DimPlot(RunTSNE(experiment.test.noc, dims = 1:20), group.by = "batchid", reduction = "tsne")  + ggtitle("No correction, pca 20")
```

![](scRNA_Workshop-PART3_files/figure-html/unnamed-chunk-14-2.png)<!-- -->

```r
#DimPlot(RunTSNE(experiment.test.regress, dims = 1:5), group.by = "batchid", reduction = "tsne") + ggtitle("Seurat regress, pca 5")
DimPlot(RunTSNE(experiment.test.regress, dims = 1:10), group.by = "batchid", reduction = "tsne") + ggtitle("Seurat regress, pca 10")
```

![](scRNA_Workshop-PART3_files/figure-html/unnamed-chunk-14-3.png)<!-- -->

```r
DimPlot(RunTSNE(experiment.test.regress, dims = 1:20), group.by = "batchid", reduction = "tsne") + ggtitle("Seurat regress, pca 20")
```

![](scRNA_Workshop-PART3_files/figure-html/unnamed-chunk-14-4.png)<!-- -->



## Get the next Rmd file

```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2019-single-cell-RNA-sequencing-Workshop-UCD_UCSF/master/scrnaseq_analysis/scRNA_Workshop-PART4.Rmd", "scRNA_Workshop-PART4.Rmd")
```

## Session Information

```r
sessionInfo()
```

```
## R version 3.6.2 (2019-12-12)
## Platform: x86_64-apple-darwin15.6.0 (64-bit)
## Running under: macOS Catalina 10.15.3
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] sva_3.34.0          BiocParallel_1.20.1 genefilter_1.68.0  
## [4] mgcv_1.8-31         nlme_3.1-143        dplyr_0.8.3        
## [7] ggplot2_3.2.1       Seurat_3.1.2       
## 
## loaded via a namespace (and not attached):
##   [1] backports_1.1.5      sn_1.5-4             plyr_1.8.5          
##   [4] igraph_1.2.4.2       lazyeval_0.2.2       splines_3.6.2       
##   [7] listenv_0.8.0        TH.data_1.0-10       digest_0.6.23       
##  [10] htmltools_0.4.0      gdata_2.18.0         fansi_0.4.1         
##  [13] magrittr_1.5         memoise_1.1.0        cluster_2.1.0       
##  [16] ROCR_1.0-7           limma_3.42.0         globals_0.12.5      
##  [19] annotate_1.64.0      RcppParallel_4.4.4   matrixStats_0.55.0  
##  [22] R.utils_2.9.2        sandwich_2.5-1       colorspace_1.4-1    
##  [25] blob_1.2.0           ggrepel_0.8.1        xfun_0.12           
##  [28] crayon_1.3.4         RCurl_1.98-1.1       jsonlite_1.6        
##  [31] zeallot_0.1.0        survival_3.1-8       zoo_1.8-7           
##  [34] ape_5.3              glue_1.3.1           gtable_0.3.0        
##  [37] leiden_0.3.1         future.apply_1.4.0   BiocGenerics_0.32.0 
##  [40] scales_1.1.0         mvtnorm_1.0-12       DBI_1.1.0           
##  [43] bibtex_0.4.2.2       Rcpp_1.0.3           metap_1.2           
##  [46] plotrix_3.7-7        viridisLite_0.3.0    xtable_1.8-4        
##  [49] reticulate_1.13      bit_1.1-15.1         rsvd_1.0.2          
##  [52] SDMTools_1.1-221.2   stats4_3.6.2         tsne_0.1-3          
##  [55] htmlwidgets_1.5.1    httr_1.4.1           gplots_3.0.1.2      
##  [58] RColorBrewer_1.1-2   TFisher_0.2.0        ica_1.0-2           
##  [61] pkgconfig_2.0.3      XML_3.99-0.3         R.methodsS3_1.7.1   
##  [64] farver_2.0.1         uwot_0.1.5           utf8_1.1.4          
##  [67] tidyselect_0.2.5     labeling_0.3         rlang_0.4.2         
##  [70] reshape2_1.4.3       AnnotationDbi_1.48.0 munsell_0.5.0       
##  [73] tools_3.6.2          cli_2.0.1            RSQLite_2.2.0       
##  [76] ggridges_0.5.2       evaluate_0.14        stringr_1.4.0       
##  [79] yaml_2.2.0           npsurv_0.4-0         knitr_1.27          
##  [82] bit64_0.9-7          fitdistrplus_1.0-14  caTools_1.17.1.4    
##  [85] purrr_0.3.3          RANN_2.6.1           pbapply_1.4-2       
##  [88] future_1.15.1        R.oo_1.23.0          compiler_3.6.2      
##  [91] plotly_4.9.1         png_0.1-7            lsei_1.2-0          
##  [94] tibble_2.1.3         stringi_1.4.5        highr_0.8           
##  [97] lattice_0.20-38      Matrix_1.2-18        multtest_2.42.0     
## [100] vctrs_0.2.1          mutoss_0.1-12        pillar_1.4.3        
## [103] lifecycle_0.1.0      Rdpack_0.11-1        lmtest_0.9-37       
## [106] RcppAnnoy_0.0.14     data.table_1.12.8    cowplot_1.0.0       
## [109] bitops_1.0-6         irlba_2.3.3          gbRd_0.4-11         
## [112] R6_2.4.1             KernSmooth_2.23-16   gridExtra_2.3       
## [115] IRanges_2.20.2       codetools_0.2-16     MASS_7.3-51.5       
## [118] gtools_3.8.1         assertthat_0.2.1     withr_2.1.2         
## [121] sctransform_0.2.1    mnormt_1.5-5         multcomp_1.4-12     
## [124] S4Vectors_0.24.3     parallel_3.6.2       grid_3.6.2          
## [127] tidyr_1.0.0          rmarkdown_2.1        Rtsne_0.15          
## [130] numDeriv_2016.8-1.1  Biobase_2.46.0
```
