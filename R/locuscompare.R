# Make locuscatter plots
# Boxiang Liu
# 2017-12-07
#' @import cowplot
#' @importFrom data.table data.table

#' @export
read_metal=function(in_fn,marker_col='rsid',pval_col='pval'){
    if (is.character(in_fn)){
        if (grepl('.gz',in_fn)){
            d = data.table::fread(sprintf('gunzip -c %s',in_fn))
        } else {
            d = data.table::fread(in_fn)
        }

        data.table::setnames(d,c(marker_col,pval_col),c('rsid','pval'))

    } else if (is.data.frame(in_fn)){
        d = in_fn
    } else {
        stop('The argument "in_fn" must be a string or a data.frame')
    }

    data.table::setDT(d)
    d = d[,list(rsid,pval,logp=-log10(pval))]
    return(d)
}

#' @export
retrieve_LD = function(chr,snp,population){
    data(config)
    on.exit(rm(config))

    conn = RMySQL::dbConnect(RMySQL::MySQL(),"locuscompare",config$b,config$c,config$a)
    on.exit(RMySQL::dbDisconnect(conn))

    res1 = DBI::dbGetQuery(
        conn = conn,
        statement = sprintf(
            "select SNP_A, SNP_B, R2
            from tkg_p3v5a_ld_chr%s_%s
            where SNP_A = '%s';",
            chr,
            population,
            snp
        )
    )

    res2 = DBI::dbGetQuery(
        conn = conn,
        statement = sprintf(
            "select SNP_B as SNP_A, SNP_A as SNP_B, R2
            from tkg_p3v5a_ld_chr%s_%s
            where SNP_B = '%s';",
            chr,
            population,
            snp
        )
    )

    res = rbind(res1,res2)
    data.table::setDT(res)
    return(res)
}

#' @export
assign_color=function(rsid,snp,ld){
    all_snps=unique(ld$SNP_A)
    color_dt=ld[SNP_A==snp,list(rsid=SNP_B,color=cut(R2,breaks=c(0,0.2,0.4,0.6,0.8,1),
                                                     labels=c('blue4','skyblue','darkgreen','orange','red'),

                                                     include.lowest=TRUE))]
    color_dt=rbind(color_dt,data.table(rsid=all_snps[!all_snps%in%color_dt$rsid],color='blue4'))
    color_dt[rsid==snp,color:='purple']
    color_dt=rbind(color_dt,data.table(rsid=rsid[!rsid%in%all_snps],color='grey'))
    color=as.character(color_dt$color)
    names(color)=color_dt$rsid
    return(color)
}

#' @export
make_combined_plot = function (merged, title1, title2, ld, snp = NULL, combine = TRUE,
                               legend = TRUE, legend_position = c('bottomright','topright','topleft'),lz_ylab_linebreak=FALSE) {
    snp = get_lead_snp(merged, snp)
    print(sprintf("INFO - %s", snp))
    color = assign_color(merged$rsid, snp, ld)
    shape = ifelse(merged$rsid == snp, 23, 21)
    names(shape) = merged$rsid
    size = ifelse(merged$rsid == snp, 3, 2)
    names(size) = merged$rsid
    merged[, `:=`(label, ifelse(rsid == snp, rsid, ""))]
    p1 = make_locuscatter(merged, title1, title2, ld, color,
                          shape, size, legend, legend_position)
    p2 = make_locuszoom(merged[, list(rsid, logp = logp1, label)],
                        title1, ld, color, shape, size, lz_ylab_linebreak)
    p3 = make_locuszoom(merged[, list(rsid, logp = logp2, label)],
                        title2, ld, color, shape, size, lz_ylab_linebreak)
    if (combine) {
        p2 = p2 + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
        p4 = plot_grid(p2, p3, align = "v", nrow = 2, rel_heights=c(0.8,1))
        p5 = plot_grid(p1, p4)
        return(p5)
    }
    else {
        return(list(locuscompare = p1, locuszoom1 = p2, locuszoom2 = p3))
    }
}

#' @export
make_locuscatter = function (merged, title1, title2, ld, color, shape, size, legend = TRUE, legend_position = c('bottomright','topright','topleft')) {
    p = ggplot(merged, aes(logp1, logp2)) + geom_point(aes(fill = rsid,
                                                           size = rsid, shape = rsid), alpha = 0.8) + geom_point(data = merged[label !=
                                                                                                                                   ""], aes(logp1, logp2, fill = rsid, size = rsid, shape = rsid)) +
        xlab(bquote(.(title1) ~ -log[10] * '(P)')) + ylab(bquote(.(title2) ~
                                                                 -log[10] * '(P)')) + scale_fill_manual(values = color, guide = "none") +
        scale_shape_manual(values = shape, guide = "none") +
        scale_size_manual(values = size, guide = "none") + geom_text_repel(aes(label = label))
    if (legend == TRUE) {
        legend_position = match.arg(legend_position)
        if (legend_position == 'bottomright'){
            legend_box = data.frame(x = 0.8, y = seq(0.4, 0.2, -0.05))
        } else if (legend_position == 'topright'){
            legend_box = data.frame(x = 0.8, y = seq(0.8, 0.6, -0.05))
        } else {
            legend_box = data.frame(x = 0.2, y = seq(0.8, 0.6, -0.05))
        }
        p = ggdraw(p) + geom_rect(data = legend_box, aes(xmin = x,
                                                         xmax = x + 0.05, ymin = y, ymax = y + 0.05), color = "black",
                                  fill = rev(c("blue4", "skyblue", "darkgreen", "orange",
                                               "red"))) + draw_label("0.8", x = legend_box$x[1] +
                                                                         0.05, y = legend_box$y[1], hjust = -0.3, size = 10) +
            draw_label("0.6", x = legend_box$x[2] + 0.05, y = legend_box$y[2],
                       hjust = -0.3, size = 10) + draw_label("0.4",
                                                             x = legend_box$x[3] + 0.05, y = legend_box$y[3],
                                                             hjust = -0.3, size = 10) + draw_label("0.2", x = legend_box$x[4] +
                                                                                                       0.05, y = legend_box$y[4], hjust = -0.3, size = 10) +
            draw_label(parse(text = "r^2"), x = legend_box$x[1] +
                           0.05, y = legend_box$y[1], vjust = -2, size = 10)
    }
    else {
        NULL
    }
    return(p)
}

#' @export
make_locuszoom=function(metal,title,ld,color,shape,size,ylab_linebreak=FALSE){
    data=merge(metal,unique(ld[,list(chr=CHR_A,pos=BP_A,rsid=SNP_A)]),by='rsid')
    chr=unique(data$chr)
    p = ggplot(data,aes(x=pos,logp))+
        geom_point(aes(fill=rsid,size=rsid,shape=rsid),alpha=0.8)+
        geom_point(data=data[label!=''],aes(x=pos,logp,fill=rsid,size=rsid,shape=rsid))+
        scale_fill_manual(values=color,guide='none')+
        scale_shape_manual(values=shape,guide='none')+
        scale_size_manual(values=size,guide='none')+
        scale_x_continuous(labels=function(x){sprintf('%.1f',x/1e6)})+
        geom_text_repel(aes(label=label))+
        xlab(paste0('chr',chr,' (Mb)'))+
        ylab(bquote(.(title)~-log[10]*'(P)'))+
        theme(plot.margin=unit(c(0.5, 1, 0.5, 0.5), "lines"))
    if (ylab_linebreak==TRUE){
        p = p + ylab(bquote(atop(.(title),-log[10]*'(P)')))
    }
    return(p)
}


#' @export
get_lead_snp = function(merged, snp = NULL){
    if (is.null(snp)) {
        snp = merged[which.min(pval1 + pval2), rsid]
    }
    else {
        if (!snp %in% merged$rsid) {
            stop(sprintf("%s not found in the intersection of %s and %s", snp, in_fn1, in_fn2))
        }
    }
    return(snp)
}

in_fn1 = system.file('extdata','gwas.tsv', package = 'locuscomparer')
in_fn2 = system.file('extdata','eqtl.tsv', package = 'locuscomparer')
marker_col1 = marker_col2 = 'rsid'
pval_col1 = pval_col2 = 'pval'
chr='6'
snp = NULL
population = 'EUR'

#' @export
main = function (in_fn1, marker_col1 = "rsid", pval_col1 = "pval", title1 = "eQTL",
                 in_fn2, marker_col2 = "rsid", pval_col2 = "pval", title2 = "GWAS",
                 snp = NULL, population = "EUR", vcf_fn, combine = TRUE, legend = TRUE,
                 legend_position = c('bottomright','topright','topleft'),
                 lz_ylab_linebreak = FALSE) {
    d1 = read_metal(in_fn1, marker_col1, pval_col1)
    d2 = read_metal(in_fn2, marker_col2, pval_col2)
    merged = merge(d1, d2, by = "rsid", suffixes = c("1", "2"), all = FALSE)
    snp = get_lead_snp(merged, snp)
    ld = retrieve_LD(chr, snp, population)
    p = make_combined_plot(merged, title1, title2, ld, snp, combine,
                           legend, legend_position, lz_ylab_linebreak)
    return(p)
}

