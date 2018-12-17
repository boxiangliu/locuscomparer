context('Test assign color')

SNP_B = paste0('rs', 0:10)
R2 = seq(1, 0, -0.1)
ld = data.frame(SNP_A = 'rs1', SNP_B, R2)
rsid = paste0('rs',0:11)

res = assign_color(rsid, 'rs1', ld)

test_that('assign_color',{
    expect_equal(unname(res['rs1']),'purple')
    expect_equal(unname(res['rs11']),'blue4')
    expect_equal(unname(res['rs0']),'red')
})
