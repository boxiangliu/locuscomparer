context('Test main')

in_fn1 = system.file('extdata','gwas.tsv', package = 'locuscomparer')
in_fn2 = system.file('extdata','eqtl.tsv', package = 'locuscomparer')
marker_col1 = 'rsid'
marker_col2 = 'rsid'
pval_col1 = pval_col2 = 'pval'
snp = NULL
population = 'EUR'
title1 = 'GWAS'; title2 = 'eQTL'
legend = combine = TRUE
lz_ylab_linebreak = FALSE
legend_position = 'bottomright'

p = main(in_fn1 = in_fn1, in_fn2 = in_fn2)
p

