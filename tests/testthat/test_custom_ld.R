context('locuscompare custom LD')

test_that('supplying ld without snp errors', {
    ld = data.frame(SNP_A = 'rs1', SNP_B = 'rs2', R2 = 0.5, stringsAsFactors = FALSE)
    expect_error(locuscompare('a.tsv', 'b.tsv', ld = ld), "specify the lead 'snp'")
})

test_that('malformed ld (missing required columns) errors', {
    bad = data.frame(A = 'rs1', B = 'rs2', r2 = 0.5, stringsAsFactors = FALSE)
    expect_error(locuscompare('a.tsv', 'b.tsv', snp = 'rs1', ld = bad), 'SNP_A, SNP_B, R2')
})

test_that('ld whose SNP_A never matches the lead snp warns', {
    ld = data.frame(SNP_A = 'rsX', SNP_B = 'rsY', R2 = 0.5, stringsAsFactors = FALSE)
    expect_warning(
        tryCatch(locuscompare('a.tsv', 'b.tsv', snp = 'rs1', ld = ld), error = function(e) NULL),
        'low-LD')
})
