
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MyRPackage

<!-- badges: start -->
<!-- badges: end -->

R functions used at my work.

## Installation

You can install the development version of MyRPackage from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("fentouxungui/MyRPackage")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(MyRPackage)
```

### scRNAseq

#### Predict Cluster location from bulk RNA-seq

**method 1: Region top Genes in binary mode**

``` r
data(scRNA)
data(bulkRNA)
score.list <- scRNAseq_Score_Region(scRNA, bulkRNA)
scRNAseq_Score_Region_evaluate(score.list)
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

``` r
# correlation of each parameter combination
# scRNAseq_Score_Region_evaluate2(score.list)
```

``` r
scRNAseq_Score_Region_plot(score.list)
#> Using UMI Cutoff: 20; Genes Used: 10
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

``` r
scRNAseq_Score_Region_plot(score.list, 100, 100)
```

<img src="man/figures/README-unnamed-chunk-3-2.png" width="100%" />

**method 2: Expression correlation**

``` r
score.matrix <- scRNAseq_Score_Region2(scRNA, bulkRNA, Method = "spearman")
pheatmap::pheatmap(score.matrix)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

``` r
score.matrix <- scRNAseq_Score_Region2(scRNA, bulkRNA, Method = "spearman", Genes.Selection = "Top")
pheatmap::pheatmap(score.matrix)
```

<img src="man/figures/README-unnamed-chunk-4-2.png" width="100%" />

**compare results from two methods**

``` r
head(scRNAseq_Score_Compare(score.list,score.matrix),20)
#>   UMI-10-Genes-400  UMI-100-Genes-400 UMI-1000-Genes-400 UMI-1500-Genes-300 
#>          0.7387318          0.7387318          0.7387318          0.7123555 
#>   UMI-10-Genes-300  UMI-100-Genes-300 UMI-1000-Genes-300   UMI-10-Genes-500 
#>          0.7095885          0.7095885          0.7095885          0.7006579 
#>  UMI-100-Genes-500 UMI-1000-Genes-500    UMI-50-Genes-50   UMI-500-Genes-50 
#>          0.7006579          0.7006579          0.6971685          0.6971685 
#> UMI-1500-Genes-200    UMI-50-Genes-30   UMI-500-Genes-30   UMI-10-Genes-200 
#>          0.6921613          0.6844487          0.6844487          0.6827991 
#>  UMI-100-Genes-200 UMI-1000-Genes-200   UMI-30-Genes-200  UMI-200-Genes-300 
#>          0.6827991          0.6827991          0.6770857          0.6748740
```
