context('Test combined plot')

in_fn1 = system.file('extdata','gwas.tsv', package = 'locuscomparer')
in_fn2 = system.file('extdata','eqtl.tsv', package = 'locuscomparer')
marker_col1 = marker_col2 = 'rsid'
pval_col1 = pval_col2 = 'pval'
chr='6'
snp = NULL
population = 'EUR'

d1 = read_metal(in_fn1, marker_col1, pval_col1)
d2 = read_metal(in_fn2, marker_col2, pval_col2)
merged = merge(d1, d2, by = "rsid", suffixes = c("1", "2"), all = FALSE)
merged = get_position(merged)
snp = get_lead_snp(merged, snp)

ld = retrieve_LD(chr, snp, population)
color = assign_color(merged$rsid, snp, ld)

shape = ifelse(merged$rsid == snp, 23, 21)
names(shape) = merged$rsid

size = ifelse(merged$rsid == snp, 3, 2)
names(size) = merged$rsid

merged$label = ifelse(merged$rsid == snp, merged$rsid, '')

p = make_combined_plot(merged, 'GWAS', 'eQTL', ld, chr, snp = NULL,
                       combine = TRUE, legend = TRUE,
                       legend_position = 'bottomright', lz_ylab_linebreak=FALSE)
p
