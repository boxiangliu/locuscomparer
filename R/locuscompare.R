# Make locuscatter plots
# Boxiang Liu
# 2017-12-07
library(R.utils)
library(data.table)
library(cowplot)
library(ggrepel)
library(stringr)

read_metal=function(in_fn,marker_col='rsid',pval_col='pval'){
    if (is.character(in_fn)){
        if (grepl('.gz',in_fn)){
            d=fread(sprintf('gunzip -c %s',in_fn))
        } else {
            d=fread(in_fn)
        }

        setnames(d,c(marker_col,pval_col),c('rsid','pval'))

    } else if (is.data.frame(in_fn)){
        d=in_fn
    } else {
        stop('in_fn must be a string or a data.frame')
    }
    setDT(d)
    d=d[,list(rsid,pval,logp=-log10(pval))]
    return(d)
}

read_full_metal=function(in_fn,marker_col='rsid',pval_col='pval',a1_col,a2_col,effect_col,se_col){
    if (is.character(in_fn)){
        if (grepl('.gz',in_fn)){
            d=fread(sprintf('gunzip -c %s',in_fn))
        } else {
            d=fread(in_fn)
        }

        setnames(d,c(marker_col,a1_col,a2_col,effect_col,se_col,pval_col),c('rsid','a1','a2','effect','se','pval'))

    } else if (is.data.frame(in_fn)){
        d=in_fn
    } else {
        stop('in_fn must be a string or a data.frame')
    }
    setDT(d)
    d=d[,list(rsid,a1,a2,effect,se,pval,logp=-log10(pval))]
    return(d)
}
extract_population=function(population,
                            out_file,
                            panel='/srv/persistent/bliu2/shared/1000genomes/phase3v5a/integrated_call_samples_v3.20130502.ALL.panel'){
    panel=fread(panel)
    x=panel[super_pop==population,list(sample,sample)]
    fwrite(x,out_file,sep='\t',col.names=FALSE)
}


subset_vcf=function(vcf_in,rsid,population,vcf_out_prefix){
    pop_fn=system.file('extdata/population/',sprintf('%s.txt',population),package = 'locuscomparer')
    out_dir=tempdir()
    rsid_fn=sprintf('%s/rsid.txt',out_dir)

    write.table(rsid,rsid_fn,sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)

    command=sprintf('plink --vcf %s --keep-allele-order --keep %s --extract %s --recode vcf-iid --out %s',vcf_in,pop_fn,rsid_fn,vcf_out_prefix)
    print(command)
    system(command)
    vcf_out = sprintf('%s.vcf',vcf_out_prefix)
    return(vcf_out)
}

get_position=function(vcf_in,x){
    y = read.table(
        file = vcf_in,
        header = FALSE,
        stringsAsFactors = FALSE
        )[,1:3]
    colnames(y) = c('chr','pos','rsid')

    stopifnot('rsid' %in% colnames(x))
    x = merge(x,y,by='rsid')
    return(x)
}

get_rsid=function(vcf_in,x){
    y = read.table(
        file = vcf_in,
        header = FALSE,
        stringsAsFactors = FALSE
        )[,1:3]
    colnames(y) = c('chr','pos','rsid')

    stopifnot(all(c('chr','pos') %in% colnames(x)))
    x = merge(x,y,by=c('chr','pos'))
    return(x)
}

calc_LD=function(rsid,vcf_in){
    out_fn_prefix=tempfile()
    out_fn=sprintf('%s.ld',out_fn_prefix)

    command=sprintf('plink --vcf %s --keep-allele-order --r2 --ld-window 9999999 --ld-window-kb 9999999 --out %s',vcf_in,out_fn_prefix)
    print(command)
    system(command)

    ld=fread(out_fn)
    ld2=ld[,list(CHR_A=CHR_B,BP_A=BP_B,SNP_A=SNP_B,CHR_B=CHR_A,BP_B=BP_A,SNP_B=SNP_A,R2)]
    ld=rbind(ld,ld2)
    return(ld)
}


assign_color=function(rsid,snp,ld){
    all_snps=unique(ld$SNP_A)

    color_dt=ld[SNP_A==snp,list(rsid=SNP_B,color=cut(R2,breaks=c(0,0.2,0.4,0.6,0.8,1),
                                                     labels=c('blue4','skyblue','darkgreen','orange','red'),
                                                     include.lowest=TRUE))]
    color_dt=rbind(color_dt,data.table(rsid=all_snps[!all_snps%in%color_dt$rsid],color='blue4'))
    color_dt=rbind(color_dt,data.table(rsid=rsid[!rsid%in%all_snps],color='grey'))
    color_dt[rsid==snp,color:='purple']
    color=as.character(color_dt$color)
    names(color)=color_dt$rsid
    return(color)
}


make_combined_plot=function(merged,title1,title2,ld,snp=NULL,combine=TRUE,legend=TRUE){
    if (is.null(snp)){
        snp=merged[which.min(pval1+pval2),rsid]
    } else {
        if(!snp%in%merged$rsid){
            stop(sprintf('%s not found in %s',snp,in_fn1))
        }
    }
    print(sprintf('INFO - %s',snp))

    color=assign_color(merged$rsid,snp,ld)
    shape=ifelse(merged$rsid==snp,23,21)
    names(shape)=merged$rsid
    size=ifelse(merged$rsid==snp,3,2)
    names(size)=merged$rsid
    merged[,label:=ifelse(rsid==snp,rsid,'')]

    p1=make_locuscatter(merged,title1,title2,ld,color,shape,size,legend)
    p2=make_locuszoom(merged[,list(rsid,logp=logp1,label)],title1,ld,color,shape,size)
    p3=make_locuszoom(merged[,list(rsid,logp=logp2,label)],title2,ld,color,shape,size)

    if (combine){
        p2=p2+theme(axis.text.x=element_blank(),axis.title.x=element_blank())
        p4=plot_grid(p2,p3,align='v',nrow=2)
        p5=plot_grid(p1,p4)
        return(p5)
    } else {
        return(list(locuscompare=p1,locuszoom1=p2,locuszoom2=p3))
    }
}


make_locuscatter=function(merged,title1,title2,ld,color,shape,size,legend=TRUE){
    p=ggplot(merged,aes(logp1,logp2))+
        geom_point(aes(fill=rsid,size=rsid,shape=rsid),alpha=0.8)+
        geom_point(data=merged[label!=''],aes(logp1,logp2,fill=rsid,size=rsid,shape=rsid))+
        xlab(bquote(.(title1)~-log[10](P)))+ylab(bquote(.(title1)~-log[10](P)))+
        scale_fill_manual(values=color,guide='none')+
        scale_shape_manual(values=shape,guide='none')+
        scale_size_manual(values=size,guide='none')+
        geom_text_repel(aes(label=label))

    if (legend==TRUE){
        legend_box=data.frame(x=0.8,y=seq(0.4,0.2,-0.05))
        p=ggdraw(p)+geom_rect(data=legend_box,aes(xmin=x,xmax=x+0.05,ymin=y,ymax=y+0.05),color='black',fill=rev(c('blue4','skyblue','darkgreen','orange','red')))+
            draw_label('0.8',x=legend_box$x[1]+0.05,y=legend_box$y[1],hjust=-0.3,size=10)+
            draw_label('0.6',x=legend_box$x[2]+0.05,y=legend_box$y[2],hjust=-0.3,size=10)+
            draw_label('0.4',x=legend_box$x[3]+0.05,y=legend_box$y[3],hjust=-0.3,size=10)+
            draw_label('0.2',x=legend_box$x[4]+0.05,y=legend_box$y[4],hjust=-0.3,size=10)+
            draw_label(parse(text='r^2'),x=legend_box$x[1]+0.05,y=legend_box$y[1],vjust=-2.0,size=10)
    } else {
        NULL
    }
    return(p)
}

make_locuszoom=function(metal,title,ld,color,shape,size){
    data=merge(metal,unique(ld[,list(chr=CHR_A,pos=BP_A,rsid=SNP_A)]),by='rsid')
    chr=unique(data$chr)
    ggplot(data,aes(x=pos,logp))+
        geom_point(aes(fill=rsid,size=rsid,shape=rsid),alpha=0.8)+
        geom_point(data=data[label!=''],aes(x=pos,logp,fill=rsid,size=rsid,shape=rsid))+
        scale_fill_manual(values=color,guide='none')+
        scale_shape_manual(values=shape,guide='none')+
        scale_size_manual(values=size,guide='none')+
        scale_x_continuous(labels=function(x){sprintf('%.1f',x/1e6)})+
        geom_text_repel(aes(label=label))+
        xlab(paste0('chr',chr,' (Mb)'))+
        ylab(bquote(.(title)~-log[10](P)))+
        theme(plot.margin=unit(c(0.5, 1, 0.5, 0.5), "lines"))
}

main=function(in_fn1,marker_col1='rsid',pval_col1='pval',title1='eQTL',
              in_fn2,marker_col2='rsid',pval_col2='pval',title2='GWAS',
              snp=NULL,population='EUR',vcf_fn,combine=TRUE,legend=TRUE){

    d1=read_metal(in_fn1,marker_col1,pval_col1)
    d2=read_metal(in_fn2,marker_col2,pval_col2)

    merged=merge(d1,d2,by='rsid',suffixes=c('1','2'),all=FALSE)

    sub_vcf_prefix = tempfile()
    sub_vcf_fn = subset_vcf(vcf_in = vcf_fn, rsid = merged$rsid, population = population, vcf_out_prefix = sub_vcf_prefix)
    merged = get_position(vcf_in = sub_vcf_fn,x = merged)
    ld = calc_LD(rsid = merged$rsid, vcf_in = sub_vcf_fn)

    p=make_combined_plot(merged,title1,title2,ld,snp,combine,legend)
    return(p)
}

