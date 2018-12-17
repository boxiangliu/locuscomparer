# locuscomparer

## 1. Installation
LocusCompareR is an R package for visualization of GWAS-eQTL colocalization events. 


Use the following commands to install LocusCompareR. If you don't have `devtools`, uncomment the first line to install it. 

```
# install.packages("devtools")
devtools::install_github("boxiangliu/locuscomparer")
```

## 2. Example

To illustrate the use of locuscompare, we use the GWAS dataset from Nikpay et al. (2015) and the coronary artery eQTL dataset from GTEx v7 at the *PHACTR1* locus: 

```
library(locuscomparer)
gwas_fn = system.file('extdata','gwas.tsv', package = 'locuscomparer')
eqtl_fn = system.file('extdata','eqtl.tsv', package = 'locuscomparer')
main(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title1 = 'CAD GWAS', title2 = 'Coronary Artery eQTL')
```

The output from the `main` function is a figure like the following:

![](https://raw.githubusercontent.com/boxiangliu/locuscomparer/master/fig/locuscompare.png)

The labeled SNP is the lead SNP (in this case for both studies), and other SNPs are colored according to their LD $r^2$ with the lead SNP.

## 3. Using your own dataset:

The input to `locuscompare::main()` is a two-column tab-delimited text file with two columns: 

1. rsid 
2. pval 

Here is an example file:

```
rsid	pval
rs62156064	0.564395
rs7562234	0.399642
rs11677377	0.34308
rs35076156	0.625237
```

You can download the example files below:  [GWAS](https://raw.githubusercontent.com/boxiangliu/locuscomparer/master/inst/extdata/gwas.tsv) and [eQTL](https://raw.githubusercontent.com/boxiangliu/locuscomparer/master/inst/extdata/eqtl.tsv) datasets. 

Then run the following commands: 
```
library(locuscomparer)
gwas_fn = 'path/to/gwas.tsv'
eqtl_fn = 'path/to/eqtl.tsv'
main(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'GWAS', title2 = 'eQTL')
```

