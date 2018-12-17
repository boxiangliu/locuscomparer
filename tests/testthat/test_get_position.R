context('test get_position')

in_fn1 = system.file('extdata','gwas.tsv', package = 'locuscomparer')
in_fn2 = system.file('extdata','eqtl.tsv', package = 'locuscomparer')
marker_col1 = marker_col2 = 'rsid'
pval_col1 = pval_col2 = 'pval'
snp = NULL
population = 'EUR'

d1 = read_metal(in_fn1, marker_col1, pval_col1)
d2 = read_metal(in_fn2, marker_col2, pval_col2)
merged = merge(d1, d2, by = "rsid", suffixes = c("1", "2"), all = FALSE)

res = get_position(merged)
test_that('get_position returns the right position',{
    expect_equal(res[1,'chr'], '6')
    expect_equal(res[1,'pos'], 12366743)
})
