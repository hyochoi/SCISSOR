
# Tutorial 

This tutorial is for a single gene analysis with `SCISSOR`. For more information, you can find R documentation for each function using `?function`. E.g.:

```r
?read_BAM
```

## Get gene annotation

SCISSOR needs gene annotation (genomic ranges for exons) formatted as "chr1:1-100,200-300:+". You can use `build_gaf` to obtain gene annotation from GTF file. Make sure that you specify a gene symbol, e.g. TP53, as an input in `build_gaf`.

```r
Gene = "TP53"
regions = build_gaf(Gene=Gene,GTF.file="./Homo_sapiens.GRCh37.87.gtf")
```

Suppose that GTF file have information about our hypothetical gene, "TOY". (See [Demo](example.md)) Then, the `build_gaf` gives the exons for the gene "TOY", and the output will be something like:

```r
> regions
[1] "chrQ:7571719-7572198,7574858-7575157,7598088-7598437:-"
```

## Get coverage from BAM files

SCISSOR takes base-level pileup for a single gene as an input. If you want to get the pileup data from BAM files, you can use `read_BAM`. Suppose that your BAM files are located under the directory, ~/bamDir/. You can read the part of the BAM files for particular regions of interest into R. Then, the resulting data object **pileupData** is a matrix where samples are in columns and genomic coordinates are in rows. If you have IDs for your samples, you can specify them for the argument `caseIDs` with the same order as the `BAMfiles`. 

```r
BAMfiles = list.files(path="~/bamDir/")
BAMfilesPath = as.character(sapply(BAMfiles,function(x) paste(getwd(),x,sep="/")))
pileupData = read_BAM(BAMfiles=bamfilesPath,caseIDs=case.barcodes,
                      symbol=Gene,regions=regions,outputType=outputType)
```

`outputType` is to set the type of intronic region that will be included in the pileup output, with choices "whole_intron", "part_intron" (default), or "only_exon". If you want SCISSOR to look for changes in intronic regions, it is good to set `outputType`="part_intron". Gene models are modified to include a portion but not all of intronic regions to facilitate the common alterations that involve intron-exon boundaries. The omission of large portions of introns is reasonable because they complicate the visualization of RNA pileups and add little to biologic signal. To determine which parts of introns to be included in the model, a basic rule is that the total lengths of bases for all exons and all introns at a gene to be approximately equal for the current SCISSOR application. This helps to make variations of expression at exonic regions and intronic regions comparable. More details on the rule for determining the intronic regions are in Methods of the paper (cite). `outputType`="only_exon" may be useful when you are only interested in changes in exons (e.g. exon skipping, alternative exon, deletion, etc.). We do not recommend `outputType`="whole_intron" for the statistical analysis of SCISSOR because it is very likely to add a large amount of noise rather than signals. However, `outputType`="whole_intron" can be useful when you want to visualize coverage with whole intronic regions. 


## Get genomic ranges for SCISSOR analysis

An important step of SCISSOR is to get genomic ranges for the gene of interest using `get_Ranges`. Wait, we already obtained `regions` previously... Why do we need this step? This step is needed to specify what regions are included in the analysis because our data object (pileup) might include part of the introns, whole exons, or only exons, which can be different from `regions`. Determined by the argument `outputType`, we get our new annotation to be used in the downstream analysis. 

```r
geneRanges = get_Ranges(Gene=Gene,regions=regions,outputType="part_intron")
```

Using this new annotation `geneRanges` as an input in other core functions, we let them know our genomic ranges of interest.

## Plot coverage

Let's plot coverage using `plot_pileup` for the TOY gene. Here, we randomly chose a subset of samples (randomSamples). 

```r
plot_pileup(Pileup=pileupData,Ranges=geneRanges,cases=randomSamples,
            main="Raw coverage")
```

![toyrawpileup](images/TOY_raw.png)

Let's plot the log-transformed coverage ($ \log (pileup + c) $) where $c$ is a pseudo count that is added before the log-transformation. To print labels for raw read count instead of log-transformed read count on the y-axis, you can specify the pseudo-count ($c$) for the argument `logcount`. 

```r
plot_pileup(Pileup=log10(pileupData+1),Ranges=geneRanges,cases=randomSamples,
            main="Log-transformed coverage",logcount=1)
```

![toyrawpileup](images/TOY_log.png)


## Run SCISSOR

`Scissor` is all-in-one function performing transformation, normalization, and the statistical analysis. 

We have base-level pileup data (as the object, **pileupData**) and genomic ranges (as the object, **geneRanges**) from the previous steps. `Scissor` takes these as inputs with other options to  identify various types of structural changes such as abnormal splicing (exon skipping and intron retention), alternative transcription start or termination, small deletions, and etc. You can use `Scissor` as the following simple command:

```r
ScissorOutput=Scissor(pileupData=pileupData,Ranges=geneRanges)
```

`Scissor` performs:

* logarithmic transformation by automatically choosing the log shift parameter
* base-level normalization   
* global shape change detection by exploring all possible low-dimensional space   
* local shape change detection by exploring residual space   



`Scissor` provides:

* shape changes identified
* outlyingness scores (global and local)  
* cutoff values (global and local)    
* most outlyingness directions for the identified shape changes   

For more information, see the R documentaion: `?Scissor` or `help(Scissor)`. 






