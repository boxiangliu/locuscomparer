context('get_lead_snp lead selection')

test_that('min(pval1+pval2) for clean p-values (backward compatible)', {
    m = data.frame(rsid = c('a','b','c'),
                   pval1 = c(0.1, 0.001, 0.5), pval2 = c(0.2, 0.01, 0.5),
                   logp1 = c(1, 3, 0.3), logp2 = c(0.7, 2, 0.3),
                   stringsAsFactors = FALSE)
    expect_equal(get_lead_snp(m), 'b')
})

test_that('falls back to max(logp1+logp2) when a p-value is 0', {
    m = data.frame(rsid = c('a','b'),
                   pval1 = c(0, 0.3), pval2 = c(0.99, 0.3),
                   logp1 = c(307.6, 0.523), logp2 = c(0.0044, 0.523),
                   stringsAsFactors = FALSE)
    # min(pval): b (0.6) < a (0.99) -> 'b';  max(logp): a (307.6) > b (1.05) -> 'a'
    expect_equal(get_lead_snp(m), 'a')
})

test_that('falls back to max(logp1+logp2) when p-values are NA (z/logp input)', {
    m = data.frame(rsid = c('a','b'),
                   pval1 = NA_real_, pval2 = NA_real_,
                   logp1 = c(10, 2), logp2 = c(1, 8),
                   stringsAsFactors = FALSE)
    expect_equal(get_lead_snp(m), 'a')  # a = 11 > b = 10
})

test_that('explicit snp is returned and validated', {
    m = data.frame(rsid = c('a','b'), pval1 = c(0.1, 0.2), pval2 = c(0.1, 0.2),
                   logp1 = c(1, 0.7), logp2 = c(1, 0.7), stringsAsFactors = FALSE)
    expect_equal(get_lead_snp(m, 'b'), 'b')
    expect_error(get_lead_snp(m, 'zzz'), 'not found')
})
