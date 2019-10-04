## Install from Github

First, install `devtools`:

```r
if ("devtools" %in% rownames(installed.packages()) == FALSE) {install.packages("devtools")}
library(devtools)
```

Next, install `SCISSOR`:

```r
install_github("hyochoi/SCISSOR")
library(SCISSOR)
```

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

