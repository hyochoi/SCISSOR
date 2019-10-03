
# SCISSOR (Shape Changes In Selecting Sample Outliers in RNA-seq)

## Description

High-throughput sequencing protocols such as RNA-seq have made it possible to interrogate the structure and sequence variation of mRNA in addition to quantitative measures of transcript expression made available by more established microarray and other molecular techniques. While many computational tools have been proposed for identifying mRNA variation by considering splice junctions or exon-level expression, the promise of RNA-seq remains largely unrealized.  In the current report, we propose a novel framework and techniques for unbiased and robust aberrant RNA structure discovery using short read sequencing data. Shape changes in selecting sample outliers in RNA-seq (SCISSOR), is a series of procedures for transformation and normalization of base level RNA sequencing data in a transcript independent manner, as well as a statistical framework for its analysis. The resulting high dimensional object is amenable to both supervised and unsupervised analysis of structural alterations across RNA-seq cohorts in a manner that independently recaptures known variants (such as splice site mutations in tumor suppressor genes) as well as novel variants difficult to identify by any existing methodology such as recurrent alternate start sites and recurrent complex deletions in 3’ UTR.

## System Requirements

* Hardware requirements

* Software requirements


## Installation Guide

1. install devtools:

```r
if ("devtools" %in% rownames(installed.packages()) == FALSE) {install.packages("devtools")}
library(devtools)
```

2. install SCISSOR:

```r
install_github("hyochoi/SCISSOR")
library(SCISSOR)
```

Installation takes ---- mins. 

