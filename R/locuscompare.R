#' Read association summary statistics from file and append column.
#' The file must contain 2 columns: markers, i.e SNPs, and p-value.
#' The marker column should contain SNP rsIDs.
#'
#' @param in_fn (string) Path to the input file.
#' @param marker_col (string, optional) Name of the marker column. Default: 'rsid'.
#' @param pval_col (string, optional) Name of the p-value column. Default: 'pval'.
#' @examples
#' in_fn = system.file('extdata', 'gwas.tsv', package = 'locuscomparer')
#' d1 = read_metal(in_fn, marker_col = 'rsid', pval_col = 'pval')
#' @export
read_metal=function(in_fn,marker_col='rsid',pval_col='pval'){
    # message('Reading ', in_fn)

    if (is.character(in_fn)){

        d = read.table(in_fn, header = TRUE, stringsAsFactors = FALSE)
        colnames(d)[which(colnames(d) == marker_col)] = 'rsid'
        colnames(d)[which(colnames(d) == pval_col)] = 'pval'

    } else if (is.data.frame(in_fn)){

        d = in_fn

    } else {

        stop('The argument "in_fn" must be a string or a data.frame')

    }

    d$logp = -log10(d$pval)
    return(d[,c('rsid','pval','logp')])
}

#' Append two columns, chromosome (chr) and position (pos), to the input data.frame.
#'
#' @param x (data.frame) Input data.frame.
#' @param genome (string, optional) Genome assembly, either 'hg19' or 'hg38'. Default: 'hg19'.
#' @examples
#' in_fn = system.file('extdata', 'gwas.tsv', package = 'locuscomparer')
#' d1 = read_metal(in_fn, marker_col = 'rsid', pval_col = 'pval')
#' genome = "hg19"
#' get_position(d1, genome)
#' @export
get_position=function(x, genome = c('hg19','hg38')){

    data(config)
    on.exit(rm(config))

    conn = RMySQL::dbConnect(RMySQL::MySQL(),"locuscompare",config$b,config$c,config$a)
    on.exit(RMySQL::dbDisconnect(conn))

    stopifnot('rsid' %in% colnames(x))

    genome = match.arg(genome)

    cmd = sprintf("select rsid, chr, pos from tkg_p3v5a_%s where rsid in ('%s')",genome,paste0(x$rsid,collapse="','"))
    res = DBI::dbGetQuery(conn = conn, statement = cmd)
    y=merge(x,res,by='rsid')
    return(y)
}


#' Retrive SNP pairwise LD from database.
#' SNP pairwise lD are calculated based on 1000 Genomes Project Phase 3 version 5.
#' For storage-efficiency, the output will only include SNPs with r2 > 0.2 with the
#' input SNP.
#' @param chr (string) Chromosome name. e.g. '22'. Notice that the name should not contain 'chr'.
#' @param snp (string) SNP rsID.
#' @param population (string) One of the 5 popuations from 1000 Genomes: 'AFR', 'AMR', 'EAS', 'EUR', and 'SAS'.
#' @examples
#' retrieve_LD('6', 'rs9349379', 'AFR')
#'
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
    return(res)
}

#' Get the lead SNP from the list of SNPs in input data.frame
#' The lead SNP is defined as the SNP with the lowest sum of p-values
#' from the two studies.
#' @param merged (data.frame) Input data.frame, which is a result by merging two
#' association studies.
#' @param snp (string, optional) SNP rsID. If NULL, the function will select the
#' lead SNP based on the sum of p-values from the two studies. If an rsID is supplied,
#' the function will simply return the rsID.
#' @examples
#' # Select the lead SNP
#' in_fn_1 = system.file('extdata', 'gwas.tsv', package = 'locuscomparer')
#' d1 = read_metal(in_fn_1, marker_col = 'rsid', pval_col = 'pval')
#' in_fn_2 = system.file('extdata', 'eqtl.tsv', package = 'locuscomparer')
#' d2 = read_metal(in_fn_2, marker_col = 'rsid', pval_col = 'pval')
#' merged = merge(d1, d2, by = "rsid", suffixes = c("1", "2"), all = FALSE)
#' get_lead_snp(merged)
#' @export
get_lead_snp = function(merged, snp = NULL){
    if (is.null(snp)) {
        snp = merged[which.min(merged$pval1 + merged$pval2), 'rsid']
    }
    else {
        if (!snp %in% merged$rsid) {
            stop(sprintf("%s not found in the intersection of in_fn1 and in_fn2.", snp))
        }
    }
    return(as.character(snp))
}

#' Assign color to each SNP according to LD.
#' @param rsid (character vector) A vector of rsIDs on which to assign color.
#' @param snp (string) rsID for lead SNP. This SNP will be colored purple.
#' Other SNPs will be assigned color based on their LD with the lead SNP.
#' @param ld (data.frame) The output from `retrieve_LD()`.
#' @examples
#' # the data.frame merged comes from the example for `get_lead_snp()`.
#' # the data.frame ld comes from the example for `retrieve_LD()`.
#' in_fn_1 = system.file('extdata', 'gwas.tsv', package = 'locuscomparer')
#' d1 = read_metal(in_fn_1, marker_col = 'rsid', pval_col = 'pval')
#' in_fn_2 = system.file('extdata', 'eqtl.tsv', package = 'locuscomparer')
#' d2 = read_metal(in_fn_2, marker_col = 'rsid', pval_col = 'pval')
#' merged = merge(d1, d2, by = "rsid", suffixes = c("1", "2"), all = FALSE)
#' ld = retrieve_LD('6', 'rs9349379', 'AFR')
#' color = assign_color(rsid = merged$rsid, snp = 'rs9349379', ld)
#' @export
assign_color=function(rsid,snp,ld){
    ld = ld[ld$SNP_A==snp,]
    ld$color = as.character(cut(ld$R2,breaks=c(0,0.2,0.4,0.6,0.8,1), labels=c('blue4','skyblue','darkgreen','orange','red'), include.lowest=TRUE))

    color = data.frame(rsid, stringsAsFactors = FALSE)
    color = merge(color, ld[, c('SNP_B', 'color')], by.x = 'rsid', by.y = 'SNP_B', all.x = TRUE)
    color[is.na(color$color),'color'] = 'blue4'
    if (snp %in% color$rsid){
        color[rsid == snp,'color'] = 'purple'
    } else {
        color = rbind(color, data.frame(rsid = snp, color = 'purple'))
    }

    res = color$color
    names(res) = color$rsid

    return(res)
}


#' Add a column of SNP labels to input data.frame
#' @param merged (data.frame) Input data.frame, which is a result by merging two
#' association studies. See the example under `get_lead_snp()` for generation of
#' such data.frame.
#' @param snp (character vector) A vector of SNP rsIDs. If only labeling one SNP,
#' this can also be a single string.
#' @examples
#' # The data.frame merged comes from the example for `get_lead_snp()`.
#' in_fn_1 = system.file('extdata', 'gwas.tsv', package = 'locuscomparer')
#' d1 = read_metal(in_fn_1, marker_col = 'rsid', pval_col = 'pval')
#' in_fn_2 = system.file('extdata', 'eqtl.tsv', package = 'locuscomparer')
#' d2 = read_metal(in_fn_2, marker_col = 'rsid', pval_col = 'pval')
#' merged = merge(d1, d2, by = "rsid", suffixes = c("1", "2"), all = FALSE)
#' merged = add_label(merged, 'rs9349379')
#' @export
add_label = function(merged, snp){
    merged$label = ifelse(merged$rsid %in% snp, merged$rsid, '')
    return(merged)
}


#' Make a scatter plot (called the LocusCompare plot).
#' Each axis of the LocusCompare plot represent the -log10(p-value) from
#' an association study. Each point thus represent a SNP. By default, the lead SNP
#' is a purple diamond, whereas the other SNPs are colored according to
#' their LD with the lead SNP.
#' @import ggplot2
#' @import cowplot
#' @param merged (data.frame) An input data.frame which has the following
#' columns: rsid, pval1 (p-value for study 1), logp1 (p-value for study 2),
#' logp1 (log p-value for study 1), logp2 (log p-value for study 2), chr, pos.
#' See the example for `get_lead_snp()` on how to generate this data.frame.
#' @param title1 (string) The title for the x-axis.
#' @param title2 (string) The title for the y-axis.
#' @param color (data.frame) The output from `assign_color()`.
#' @param shape (data.frame) Specification of the shape of each SNP. See example blow on how to generate this data.frame.
#' @param size (data.frame) Specification of the size of each SNP. See example below on how to generate this data.frame.
#' @param legend (boolean) Whether to include the legend.
#' @param legend_position (string, optional) Either 'bottomright','topright', or 'topleft'. Default: 'bottomright'.
#' @examples
#' # The data.frame `merged` comes from the example of `add_label()`.
#' # The data.frame `color` comes from the example of `assign_color()`.
#' in_fn_1 = system.file('extdata', 'gwas.tsv', package = 'locuscomparer')
#' d1 = read_metal(in_fn_1, marker_col = 'rsid', pval_col = 'pval')
#' in_fn_2 = system.file('extdata', 'eqtl.tsv', package = 'locuscomparer')
#' d2 = read_metal(in_fn_2, marker_col = 'rsid', pval_col = 'pval')
#' merged = merge(d1, d2, by = "rsid", suffixes = c("1", "2"), all = FALSE)
#' snp = 'rs9349379'
#' shape = ifelse(merged$rsid == snp, 23, 21)
#' names(shape) = merged$rsid
#' size = ifelse(merged$rsid == snp, 3, 2)
#' names(size) = merged$rsid
#' ld = retrieve_LD('6', 'rs9349379', 'AFR')
#' color = assign_color(rsid = merged$rsid, snp = 'rs9349379', ld)
#' merged = add_label(merged, snp)
#' make_scatterplot(merged, title1 = 'GWAS', title2 = 'eQTL', color, shape, size)
#' @export
make_scatterplot = function (merged, title1, title2, color, shape, size, legend = TRUE, legend_position = c('bottomright','topright','topleft')) {

    p = ggplot(merged, aes(logp1, logp2)) +
        geom_point(aes(fill = rsid, size = rsid, shape = rsid), alpha = 0.8) +
        geom_point(data = merged[merged$label != "",],
                   aes(logp1, logp2, fill = rsid, size = rsid, shape = rsid)) +
        xlab(bquote(.(title1) ~ -log[10] * '(P)')) +
        ylab(bquote(.(title2) ~ -log[10] * '(P)')) +
        scale_fill_manual(values = color, guide = "none") +
        scale_shape_manual(values = shape, guide = "none") +
        scale_size_manual(values = size, guide = "none") +
        ggrepel::geom_text_repel(aes(label = label))+
        theme_classic()

    if (legend == TRUE) {
        legend_position = match.arg(legend_position)
        if (legend_position == 'bottomright'){
            legend_box = data.frame(x = 0.8, y = seq(0.4, 0.2, -0.05))
        } else if (legend_position == 'topright'){
            legend_box = data.frame(x = 0.8, y = seq(0.8, 0.6, -0.05))
        } else {
            legend_box = data.frame(x = 0.2, y = seq(0.8, 0.6, -0.05))
        }

        p = ggdraw(p) +
            geom_rect(data = legend_box,
                      aes(xmin = x, xmax = x + 0.05, ymin = y, ymax = y + 0.05),
                      color = "black",
                      fill = rev(c("blue4", "skyblue", "darkgreen", "orange", "red"))) +
            draw_label("0.8", x = legend_box$x[1] + 0.05, y = legend_box$y[1], hjust = -0.3, size = 10) +
            draw_label("0.6", x = legend_box$x[2] + 0.05, y = legend_box$y[2], hjust = -0.3, size = 10) +
            draw_label("0.4", x = legend_box$x[3] + 0.05, y = legend_box$y[3], hjust = -0.3, size = 10) +
            draw_label("0.2", x = legend_box$x[4] + 0.05, y = legend_box$y[4], hjust = -0.3, size = 10) +
            draw_label(parse(text = "r^2"), x = legend_box$x[1] + 0.05, y = legend_box$y[1], vjust = -2, size = 10)
    }

    return(p)
}

#' Make a locuszoom plot.
#' Details see http://locuszoom.org/.
#' @import ggplot2
#' @import cowplot
#' @param metal (data.frame) input file with two column, rsID and p-value.
#' See `read_metal()` for more details.
#' @param title (string) y-axis title.
#' @param chr (string) chromosome.
#' @param color (data.frame) The output from `assign_color()`.
#' @param shape (data.frame) Specification of the shape of each SNP. See example blow on how to generate this data.frame.
#' @param size (data.frame) Specification of the size of each SNP. See example below on how to generate this data.frame.
#' @param ylab_linebreak (boolean, optional) Whether to break the line of y-axis. If FALSE, the y-axis title and '-log10(p-value)'
#' will be on the same line. Default: FALSE.
#' @examples
#' \dontrun{  # clearly the function works but example setup is tricky
#' in_fn_1 = system.file('extdata', 'gwas.tsv', package = 'locuscomparer')
#' d1 = read_metal(in_fn_1, marker_col = 'rsid', pval_col = 'pval')
#' in_fn_2 = system.file('extdata', 'eqtl.tsv', package = 'locuscomparer')
#' d2 = read_metal(in_fn_2, marker_col = 'rsid', pval_col = 'pval')
#' merged = merge(d1, d2, by = "rsid", suffixes = c("1", "2"), all = FALSE)
#' snp = 'rs9349379'
#' merged = add_label(merged, snp)
#' d1 = add_label(d1, snp)
#' shape = ifelse(merged$rsid == snp, 23, 21)
#' names(shape) = merged$rsid
#' size = ifelse(merged$rsid == snp, 3, 2)
#' names(size) = merged$rsid
#' chr = '6'
#' d1 = get_position(d1)
#' ld = retrieve_LD(chr, snp, "AFR")
#' color = assign_color(rsid = merged$rsid, snp = 'rs9349379', ld)
#' make_locuszoom(d1, title = 'GWAS', chr, color, shape, size)
#' }
#' @export
make_locuszoom=function(metal,title,chr,color,shape,size,ylab_linebreak=FALSE){

    p = ggplot(metal,aes(x=pos,logp))+
        geom_point(aes(fill=rsid,size=rsid,shape=rsid),alpha=0.8)+
        geom_point(data=metal[metal$label!='',],aes(x=pos,logp,fill=rsid,size=rsid,shape=rsid))+
        scale_fill_manual(values=color,guide='none')+
        scale_shape_manual(values=shape,guide='none')+
        scale_size_manual(values=size,guide='none')+
        scale_x_continuous(labels=function(x){sprintf('%.1f',x/1e6)})+
        ggrepel::geom_text_repel(aes(label=label))+
        xlab(paste0('chr',chr,' (Mb)'))+
        ylab(bquote(.(title)~-log[10]*'(P)'))+
        theme_classic()+
        theme(plot.margin=unit(c(0.5, 1, 0.5, 0.5), "lines"))

    if (ylab_linebreak==TRUE){
        p = p + ylab(bquote(atop(.(title),-log[10]*'(P)')))
    }
    return(p)
}


#' Generated a combined plot with two locuszoom plots and a locuscompare
#' plot. Each locuszoom plot represent an association study.
#' @param merged (data.frame) An input data.frame which has the following
#' columns: rsid, pval1 (p-value for study 1), logp1 (p-value for study 2),
#' logp1 (log p-value for study 1), logp2 (log p-value for study 2), chr, pos.
#' See the example for `get_lead_snp()` on how to generate this data.frame.
#' @param title1 (string) The title for the x-axis.
#' @param title2 (string) The title for the y-axis.
#' @param ld (data.frame) The output from `retrieve_LD()`.
#' @param chr (string) Chromosome name. e.g. '22'. Notice that the name should not contain 'chr'.
#' @param snp (string, optional) SNP rsID. If NULL, the function will select the lead SNP. Default: NULL.
#' @param combine (boolean, optional) Should the three plots be combined into one plot? If FALSE, a list of
#' three plots will be returned. Default: TRUE.
#' @param legend (boolean, optional) Should the legend be shown? Default: TRUE.
#' @param legend_position (string, optional) Either 'bottomright','topright', or 'topleft'. Default: 'bottomright'.
#' @param lz_ylab_linebreak (boolean, optional) Whether to break the line of y-axis of the locuszoom plot.
#' @param genome character(1) for get_position
#' If FALSE, the y-axis title and '-log10(p-value)'. will be on the same line. Default: FALSE.
#' @examples
#' in_fn_1 = system.file('extdata', 'gwas.tsv', package = 'locuscomparer')
#' d1 = read_metal(in_fn_1, marker_col = 'rsid', pval_col = 'pval')
#' in_fn_2 = system.file('extdata', 'eqtl.tsv', package = 'locuscomparer')
#' d2 = read_metal(in_fn_2, marker_col = 'rsid', pval_col = 'pval')
#' merged = merge(d1, d2, by = "rsid", suffixes = c("1", "2"), all = FALSE)
#' snp = 'rs9349379'
#' merged = add_label(merged, snp)
#' chr = '6' 
#' ld = retrieve_LD(chr, snp, "AFR")
#' make_combined_plot(merged, 'GWAS', 'eQTL', ld, chr)
#' @export
make_combined_plot = function (merged, title1, title2, ld, chr, snp = NULL,
                               combine = TRUE, legend = TRUE,
                               legend_position = c('bottomright','topright','topleft'),
                               lz_ylab_linebreak=FALSE, genome="hg19") {
    snp = get_lead_snp(merged, snp)
    # print(sprintf("INFO - %s", snp))
    color = assign_color(merged$rsid, snp, ld)

    shape = ifelse(merged$rsid == snp, 23, 21)
    names(shape) = merged$rsid

    size = ifelse(merged$rsid == snp, 3, 2)
    names(size) = merged$rsid

    merged = add_label(merged, snp)
    if (!("pos" %in% names(merged))) merged = get_position(merged, genome)

    p1 = make_scatterplot(merged, title1, title2, color,
                          shape, size, legend, legend_position)

    metal1 = merged[,c('rsid', 'logp1', 'chr', 'pos', 'label')]
    colnames(metal1)[which(colnames(metal1) == 'logp1')] = 'logp'
    p2 = make_locuszoom(metal1, title1, chr, color, shape, size, lz_ylab_linebreak)

    metal2 = merged[,c('rsid', 'logp2', 'chr', 'pos', 'label')]
    colnames(metal2)[which(colnames(metal2) == 'logp2')] = 'logp'
    p3 = make_locuszoom(metal2, title2, chr, color, shape, size, lz_ylab_linebreak)

    if (combine) {
        p2 = p2 + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
        p4 = cowplot::plot_grid(p2, p3, align = "v", nrow = 2, rel_heights=c(0.8,1))
        p5 = cowplot::plot_grid(p1, p4)
        return(p5)
    }
    else {
        return(list(locuscompare = p1, locuszoom1 = p2, locuszoom2 = p3))
    }
}

#' Make a locuscompare plot.
#' @param in_fn1 (string) Path to the input file for study 1.
#' @param in_fn2 (string) Path to the input file for study 2.
#' @param marker_col1 (string, optional) Name of the marker column. Default: 'rsid'.
#' @param pval_col1 (string, optional) Name of the p-value column. Default: 'pval'.
#' @param title1 (string) The title for the x-axis.
#' @param marker_col2 (string, optional) Name of the marker column. Default: 'rsid'.
#' @param pval_col2 (string, optional) Name of the p-value column. Default: 'pval'.
#' @param title2 (string) The title for the y-axis.
#' @param snp (string, optional) SNP rsID. If NULL, the function will select the lead SNP. Default: NULL.
#' @param population (string, optional) One of the 5 popuations from 1000 Genomes: 'AFR', 'AMR', 'EAS', 'EUR', and 'SAS'. Default: 'EUR'.
#' @param combine (boolean, optional) Should the three plots be combined into one plot? If FALSE, a list of
#' three plots will be returned. Default: TRUE.
#' @param legend (boolean, optional) Should the legend be shown? Default: TRUE.
#' @param legend_position (string, optional) Either 'bottomright','topright', or 'topleft'. Default: 'bottomright'.
#' @param lz_ylab_linebreak (boolean, optional) Whether to break the line of y-axis of the locuszoom plot.
#' @param genome (string, optional) Genome assembly, either 'hg19' or 'hg38'. Default: 'hg19'.
#' @examples
#' in_fn1 = system.file('extdata','gwas.tsv', package = 'locuscomparer')
#' in_fn2 = system.file('extdata','eqtl.tsv', package = 'locuscomparer')
#' locuscompare(in_fn1 = in_fn1, in_fn2 = in_fn2)
#' @export
locuscompare = function(in_fn1, in_fn2, marker_col1 = "rsid", pval_col1 = "pval",
                 title1 = "eQTL",marker_col2 = "rsid", pval_col2 = "pval", title2 = "GWAS",
                 snp = NULL, population = "EUR", combine = TRUE, legend = TRUE,
                 legend_position = c('bottomright','topright','topleft'),
                 lz_ylab_linebreak = FALSE, genome = c('hg19','hg38')) {
    d1 = read_metal(in_fn1, marker_col1, pval_col1)
    d2 = read_metal(in_fn2, marker_col2, pval_col2)

    merged = merge(d1, d2, by = "rsid", suffixes = c("1", "2"), all = FALSE)
    genome = match.arg(genome)
    merged = get_position(merged, genome)

    chr = unique(merged$chr)
    if (length(chr) != 1) stop('There must be one and only one chromosome.')

    snp = get_lead_snp(merged, snp)
    ld = retrieve_LD(chr, snp, population)
    p = make_combined_plot(merged, title1, title2, ld, chr, snp, combine,
                           legend, legend_position, lz_ylab_linebreak)
    return(p)
}

