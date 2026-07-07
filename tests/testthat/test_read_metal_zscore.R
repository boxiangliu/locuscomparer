context('read_metal z-score input')

test_that('zscore_col yields the correct -log10 p (z = 1.96 ~ p 0.05)', {
    res = read_metal(data.frame(rsid = 'rs1', z = 1.959964, stringsAsFactors = FALSE), zscore_col = 'z')
    expect_equal(res$logp, -log10(0.05), tolerance = 1e-4)
    expect_true(is.na(res$pval))
})

test_that('z-score is sign-agnostic (two-sided)', {
    res = read_metal(data.frame(rsid = c('a','b'), z = c(5, -5), stringsAsFactors = FALSE), zscore_col = 'z')
    expect_equal(res$logp[1], res$logp[2])
})

test_that('large z stays finite (no Inf underflow)', {
    res = read_metal(data.frame(rsid = 'a', z = 50, stringsAsFactors = FALSE), zscore_col = 'z')
    expect_true(is.finite(res$logp) && res$logp > 500)
})

test_that('cannot supply both zscore_col and logp_col', {
    expect_error(
        read_metal(data.frame(rsid = 'a', z = 1, lp = 1, stringsAsFactors = FALSE),
                   zscore_col = 'z', logp_col = 'lp'),
        'only one')
})
