context('Test locuszoom')

in_fn1 = system.file('extdata','gwas.tsv', package = 'locuscomparer')
in_fn2 = system.file('extdata','eqtl.tsv', package = 'locuscomparer')
marker_col1 = marker_col2 = 'rsid'
pval_col1 = pval_col2 = 'pval'
snp = NULL
chr='6'
population = 'EUR'


d1 = read_metal(in_fn1, marker_col = 'rsid', pval_col = 'pval')
d2 = read_metal(in_fn2, marker_col = marker_col2, pval_col = pval_col2)
merged = merge(d1, d2, by = "rsid", suffixes = c("1", "2"), all = FALSE)
merged = get_position(merged)
chr = unique(merged$chr)

snp = get_lead_snp(merged, snp)
ld = retrieve_LD(chr, snp, population)
color = assign_color(merged$rsid, snp, ld)

shape = ifelse(merged$rsid == snp, 23, 21)
names(shape) = merged$rsid

size = ifelse(merged$rsid == snp, 3, 2)
names(size) = merged$rsid

merged$label = ifelse(merged$rsid == snp, merged$rsid, '')
metal = merged[, c('rsid', 'logp1', 'chr', 'pos', 'label')]
colnames(metal)[which(colnames(metal) == 'logp1')] = 'logp'
p = make_locuszoom(metal,title,chr,color,shape,size,ylab_linebreak=FALSE)
p
