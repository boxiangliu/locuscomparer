context('Test read metal')

in_fn1 = system.file('extdata','gwas.tsv', package = 'locuscomparer')
marker_col1 = 'rsid'
pval_col1 = 'pval'
d1 = read_metal(in_fn1, marker_col1, pval_col1)

test_that('headers are correct',{
    expect_equal(colnames(d1), c('rsid','pval','logp'))
})
