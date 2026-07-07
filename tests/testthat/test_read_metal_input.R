context('read_metal flexible input')

test_that('pval == 0 is floored at -log10(xmin) with a warning', {
    d = data.frame(rsid = c('rs1','rs2'), pval = c(1e-10, 0), stringsAsFactors = FALSE)
    expect_warning(res <- read_metal(d), 'underflow')
    expect_true(all(is.finite(res$logp)))
    expect_equal(res$logp[2], -log10(.Machine$double.xmin))
})

test_that('normal p-values are unchanged and silent', {
    d = data.frame(rsid = c('rs1','rs2'), pval = c(0.5, 1e-8), stringsAsFactors = FALSE)
    expect_silent(res <- read_metal(d))
    expect_equal(res$logp, -log10(c(0.5, 1e-8)))
})

test_that('logp_col is used directly and pval becomes NA', {
    d = data.frame(rsid = c('rs1','rs2'), neglogp = c(5, 320), stringsAsFactors = FALSE)
    res = read_metal(d, logp_col = 'neglogp')
    expect_equal(res$logp, c(5, 320))
    expect_true(all(is.na(res$pval)))
})

test_that('a missing named column errors clearly', {
    expect_error(read_metal(data.frame(rsid = 'a', pval = 0.1), logp_col = 'nope'), 'not found')
})
