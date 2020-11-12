
# SCISSOR

## Shape Changes In Selecting Sample Outliers in RNA-seq

## Overview

`SCISSOR` (shape changes in selecting sample outliers in RNA-seq) aims for unsupervised screening of a range of structural alterations in RNA-seq data. `SCISSOR` considers a novel shape property of aligned short read data through a base-level pileup file. This intact and uncompressed view of RNA-seq profile enables the unbiased discovery of structural alterations by looking for anomalous shapes in expression. This approach holds promise for identifying otherwise obscured genetic aberrations. As a result, `SCISSOR` identifies known as well as novel aberrations including abnormal splicing, intra-/intergenic deletions, small indels, alternative transcription start/termination. 

## Documentation

The documentation is available at [https://hyochoi.github.io/SCISSOR](https://hyochoi.github.io/SCISSOR)

## System Requirements

This package has been tested on the following systems:

* macOS: Mojave (10.14.6)  
* Linux: CentOS (7.5)

## R package dependencies

`SCISSOR` requires the following packages:

* `BiocManager`
* `Rsamtools`   
* `GenomicRanges`  
* `refGenome`   
* `RColorBrewer`   
* `wesanderson`   
* `nloptr`  
* `zoo`


## Installation Guide

1. install `devtools`:

```r
if ("devtools" %in% rownames(installed.packages()) == FALSE) {install.packages("devtools")}
library(devtools)
```

2. install `SCISSOR`:

```r
install_github("hyochoi/SCISSOR")
library(SCISSOR)
```

Installation takes 5-10 mins. 

**Note:**  

If you see errors like the following: 
```r
‘rlang’ 0.3.4 is already loaded, but >= 0.4.0 is required
```
try this:
```r
devtools::install_github("r-lib/rlang", build_vignettes = TRUE)
```




