# locuscomparer

## 1. Installation
The package locuscomparer is an R package for visualization of GWAS-eQTL colocalization events. 


Use the following commands to install LocusCompareR:

```
install.packages("devtools")
library(devtools)
install_github("boxiangliu/locuscomparer")
```

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

To illustrate the use of locuscompare, we use the Nikpay et al. (2015) and the coronary artery eQTL dataset from GTEx. First download these example [GWAS](https://raw.githubusercontent.com/boxiangliu/locuscomparer/master/data/PHACTR1_gwas.tsv) and [eQTL](https://raw.githubusercontent.com/boxiangliu/locuscomparer/master/data/Coronary_Artery_PHACTR1_eqtl.tsv) datasets. 

Then run the following commands: 
```
library(locuscomparer)
gwas_fn = system.file('extdata','gwas.tsv', package = 'locuscomparer')
eqtl_fn = system.file('extdata','eqtl.tsv', package = 'locuscomparer')
main(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'GWAS', title2 = 'eQTL')
```

The file `ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz` can be downloaded from 1000 Genomes through [FTP](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/). 

You should see a figure like this. 

![locuscompare](https://raw.githubusercontent.com/boxiangliu/locuscomparer/master/fig/locuscompare.png)
