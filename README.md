# locuscomparer

## 1. Installation
The package locuscomparer is an R package for visualization of GWAS-eQTL colocalization events. 


Use the following commands to install LocusCompareR:

```
install.packages("devtools")
library(devtools)
install_github("boxiangliu/locuscomparer")
library(locuscomparer)
```

In addition, you will need to install [`plink`](https://www.cog-genomics.org/plink2). You will also need to download [1000 Genomes genotypes](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) for LD calculation. You can run this [file](https://raw.githubusercontent.com/boxiangliu/locuscomparer/master/src/download_1000g.sh) to download the 1000 Genomes genotypes. Note that they are very large. 


## 2. Input 
The input format to locuscompare is a two-column tab-delimited text file. Here is an example file:

```
rsid	pval
rs62156064	0.564395
rs7562234	0.399642
rs11677377	0.34308
rs35076156	0.625237
```

## 3. Example

TODO: give an example and a plot. 
