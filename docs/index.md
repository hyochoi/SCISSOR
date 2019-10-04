# SCISSOR 

`SCISSOR` is an R package containing tools for statistical analysis and visualization of base-level RNA-seq data. 

## Overview

`SCISSOR` (shape changes in selecting sample outliers in RNA-seq) aims for unsupervised screening of a range of structural alterations in RNA-seq data. `SCISSOR` considers a novel shape property of aligned short read data through a base-level pileup file. This intact and uncompressed view of RNA-seq profile enables the unbiased discovery of structural alterations by looking for anomalous shapes in expression. This approach holds promise for identifying otherwise obscured genetic aberrations. As a result, `SCISSOR` identifies known as well as novel aberrations including abnormal splicing, intra-/intergenic deletions, small indels, alternative transcription start/termination. 

## Statistical model

With the goal of detecting samples exhibiting anomalous shapes, `SCISSOR` models base-level read counts using a high-dimensional latent variable framework that is naturally integrated to its normalization, abnormal feature extraction, and quantification. A latent variable is used to model an underlying abnormal trajectory, i.e. an outlier direction in a high-dimensional space, that is interrogated for outliers. An outlier case with shape changes then can be a data point that is strongly involved in one or multiple abnormal trajectories, which enables modeling complex structural variation. 

## Statistical method

`SCISSOR` extracts a latent space associated with abnormal sequencing coverage and quantifies the level of abnormality in a robust way for determining the cases with shape changes. As the type of structure of interest is outlying/abnormal, it uses a projection pursuit approach to measure how outlying a sample is in the most extreme one-dimensional direction. At each gene under consideration, the resulting statistic is an outlyingness score for each sample with larger values indicating more severe deviation from other samples in the dataset. For each outlier, `SCISSOR` produces the most outlying direction as a single best trajectory that describes abnormalities of the corresponding outlier, which can be used to recover the latent space of outlying outlier directions. 

## Free software

`SCISSOR` is free software and available on [Github](https://github.com/hyochoi/SCISSOR).


## Documentation

* [Getting Started](installation.md)  
* [Tutorial](tutorial.md)   
* [Demo](example.md)


