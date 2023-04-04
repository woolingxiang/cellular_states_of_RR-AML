#' Install R packages.
#'
#' Install R packages. Can be specified for samples or clusters to avoid confusion.
#'
#' @param pcg A vector of names of R package to be installed.
#'
#' @return A vector of information.

install_packages = function(pcg) {
	if(length(grep('^BiocManager$',installed.packages()))==0){
	  install.packages("BiocManager",repos=c(
	    'http://mirrors.ustc.edu.cn/CRAN/',
	    'http://mirror.lzu.edu.cn/CRAN/',
	    'http://mirrors.tuna.tsinghua.edu.cn/CRAN/'
	  ))
	}
	new = pcg[!(pcg %in% installed.packages()[, "Package"])]
  	if (length(new)) BiocManager::install(new)
  	sapply(pcg, library,character.only = T)
}

#' Determine the color scheme.
#'
#' Determine the color scheme. Can be specified for samples or clusters to avoid confusion.
#'
#' @param type Type of scheme ("samples" or "clusters").
#'
#' @return A vector of colors.
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggsci pal_d3 pal_igv
get_color_scheme = function(type = "clusters") {
  if (type == "samples") {
    color_scheme = c(brewer.pal(5, "Set1"), brewer.pal(8, "Dark2"), pal_igv("default")(51))
  }
  if (type == "clusters") {
    color_scheme = c(pal_d3("category10")(10), pal_d3("category20")(20), pal_d3("category20b")(20), pal_d3("category20c")(20),pal_igv("default")(51))
  }
  return(color_scheme)
}


#' Determine the point size for tSNE plots (smaller for larger datasets).
#'
#' @param num_cells Number of cells (points on the plot).
#'
#' @return Numeric point size.
get_dr_pt_size = function(num_cells) {

  pt_size = 1.8
  if (num_cells > 1000) pt_size = 1.2
  if (num_cells > 5000) pt_size = 1.0
  if (num_cells > 10000) pt_size = 0.8
  if (num_cells > 25000) pt_size = 0.6
  if (num_cells > 50000) pt_size = 0.4
  if (num_cells > 100000) pt_size = 0.2
  return(pt_size)

}


# plot colored by specified variable
plot_scatter_group = function(metadata_tbl, x_var = "UMAP_1", y_var = "UMAP_2", aspect_ratio = 1, color_var, color_scheme) {

  ggplot(metadata_tbl, aes(x = !!sym(x_var), y = !!sym(y_var))) +
    geom_point(aes(color = !!sym(color_var)), size = get_dr_pt_size(metadata_tbl)) + theme_few() + 
    theme(
      aspect.ratio = aspect_ratio,
      axis.ticks = element_blank(),
      axis.text = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    scale_color_manual(values = color_scheme)

}


# plot split by specified variable
plot_scatter_split = function(metadata_tbl, x_var = "UMAP_1", y_var = "UMAP_2", aspect_ratio = 1, rows_var = NULL, cols_var = NULL, color_var, color_scheme) {

  gp =
    ggplot(metadata_tbl, aes(x = !!sym(x_var), y = !!sym(y_var))) +
    geom_point(aes(color = !!sym(color_var)), size = get_dr_pt_size(metadata_tbl)) + theme_few() + 
    theme(
      aspect.ratio = aspect_ratio,
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      strip.background = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    scale_color_manual(values = color_scheme)

  if (is.null(rows_var)) {
    gp + facet_grid(cols = vars(!!sym(cols_var)))
  } else if (is.null(cols_var)) {
    gp + facet_grid(rows = vars(!!sym(rows_var)))
  } else {
    gp + facet_grid(rows = vars(!!sym(rows_var)), cols = vars(!!sym(cols_var)))
  }
}


# density plot split by specified variable
# calculate density normalized to 1, independently for each facet variable
plot_density_split = function(metadata_tbl, x_var, y_var, split_var, num_bins,low='white',high='darkred') {

  # ran into some issues with merging split geom_hex
  ggplot(metadata_tbl, aes(x = !!sym(x_var), y = !!sym(y_var))) + theme_few() + 
    # geom_hex(aes(fill = stat(ndensity)), bins = num_bins) +
    stat_bin_2d(aes(fill = stat(ndensity)), bins = num_bins) +
    theme(
      aspect.ratio = 1,
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      strip.background = element_blank()
    ) +
    scale_fill_gradient2(low = low, high = high) +
    facet_wrap(vars(!!sym(split_var)))

}

# get table for density plot, split by stage
get_density_diff_table = function(metadata_tbl, x_var, y_var, split_var, num_bins) {

  # generate a density plot split by stage
  density_plot = plot_density_split(metadata_tbl = metadata_tbl, x_var = x_var, y_var = y_var, split_var = split_var, num_bins = num_bins)

  # produce an object that can be rendered
  density_plot_tbl = ggplot_build(density_plot)

  # panel labels
  panels_tbl =
    tibble(
      PANEL = density_plot_tbl$layout$layout$PANEL,
      stage = density_plot_tbl$layout$layout[[split_var]]
    )

  # merge panel contents and panel names
  density_tbl = density_plot_tbl$data[[1]]
  density_tbl = density_tbl %>% full_join(panels_tbl, by = "PANEL")

  return(density_tbl)

}

# density plot split by specified variable
# split normalization (adding norm_split_var) may not work
plot_density_diff = function(metadata_tbl, x_var = "UMAP_1", y_var = "UMAP_2", split_var, num_bins, 
                             group_pos, group_neg, interpolate = FALSE,low = "#053061", mid = "gray80", high = "#E41A1C") {

  density_tbl = get_density_diff_table(metadata_tbl = metadata_tbl, x_var = x_var, y_var = y_var, split_var = split_var, num_bins = num_bins)

  min_density = quantile(density_tbl$density, 0)

  density_pos_tbl =
    density_tbl %>%
    filter(stage == group_pos) %>%
    select(x, y, cells_pos = count, density_pos = density)
  density_neg_tbl =
    density_tbl %>%
    filter(stage == group_neg) %>%
    select(x, y, cells_neg = count, density_neg = density)

  density_split_tbl = full_join(density_pos_tbl, density_neg_tbl, by = c("x", "y"))
  density_split_tbl[is.na(density_split_tbl)] = min_density
  density_split_tbl = density_split_tbl %>% mutate(density_diff = density_pos - density_neg)
  density_split_tbl = density_split_tbl %>% mutate(density_ratio = log(density_pos/density_neg))

  min_density_diff = density_split_tbl %>% pull(density_diff) %>% quantile(0.01)
  max_density_diff = density_split_tbl %>% pull(density_diff) %>% quantile(0.99)
  min_density_ratio = density_split_tbl %>% pull(density_ratio) %>% quantile(0.01)
  max_density_ratio = density_split_tbl %>% pull(density_ratio) %>% quantile(0.99)

  density_split_tbl =
    density_split_tbl %>%
    mutate(
      cells = cells_pos + cells_neg,
      log_density = log(density_pos + density_neg),
      density_ratio = if_else(density_ratio < min_density_ratio, min_density_ratio, density_ratio),
      density_ratio = if_else(density_ratio > max_density_ratio, max_density_ratio, density_ratio)
    ) %>%
    filter(cells > 0)

  ggplot(density_split_tbl, aes(x = x, y = y)) +
    # geom_tile(aes(fill = density_ratio)) +
    geom_tile(aes(fill = density_ratio), interpolate = interpolate) + theme_few() + 
    theme(
      aspect.ratio = 1,
      axis.ticks = element_blank(),
      axis.text = element_blank()
    ) +
    labs(title = glue("{group_pos} vs {group_neg}"), x = x_var, y = y_var) +
    scale_fill_gradient2(low = low, mid = mid, high = high)
}


# density plot split by specified variable
# split normalization (adding norm_split_var) may not work
get_density_split_tbl = function(metadata_tbl, x_var = "UMAP_1", y_var = "UMAP_2", split_var, num_bins, group_pos, group_neg) {

  density_tbl = get_density_diff_table(metadata_tbl = metadata_tbl, x_var = x_var, y_var = y_var, split_var = split_var, num_bins = num_bins)

  min_density = quantile(density_tbl$density, 0)

  density_pos_tbl =
    density_tbl %>%
    filter(stage == group_pos) %>%
    select(x, y, cells_pos = count, density_pos = density)
  density_neg_tbl =
    density_tbl %>%
    filter(stage == group_neg) %>%
    select(x, y, cells_neg = count, density_neg = density)

  density_split_tbl = full_join(density_pos_tbl, density_neg_tbl, by = c("x", "y"))
  density_split_tbl[is.na(density_split_tbl)] = min_density
  density_split_tbl = density_split_tbl %>% mutate(density_diff = density_pos - density_neg)
  density_split_tbl = density_split_tbl %>% mutate(density_ratio = log(density_pos/density_neg))

  min_density_diff = density_split_tbl %>% pull(density_diff) %>% quantile(0.01)
  max_density_diff = density_split_tbl %>% pull(density_diff) %>% quantile(0.99)
  min_density_ratio = density_split_tbl %>% pull(density_ratio) %>% quantile(0.01)
  max_density_ratio = density_split_tbl %>% pull(density_ratio) %>% quantile(0.99)

  density_split_tbl =
    density_split_tbl %>%
    mutate(
      cells = cells_pos + cells_neg,
      log_density = log(density_pos + density_neg),
      density_ratio = if_else(density_ratio < min_density_ratio, min_density_ratio, density_ratio),
      density_ratio = if_else(density_ratio > max_density_ratio, max_density_ratio, density_ratio)
    ) %>%
    filter(cells > 0)
    
    return(density_split_tbl)
}


# violin plot split by specified group.
# default group is orig.ident
plot_violin = function(metadata_tbl, color_scheme, y_var, x_var = "orig.ident") {

  violin_plot = ggplot(metadata_tbl, aes(x = !!sym(x_var), y = !!sym(y_var), fill = !!sym(x_var))) +
    geom_violin() +
    xlab(x_var) +
    ylab(y_var) +
    scale_fill_manual(values = color_scheme,
                      name = x_var) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  return(violin_plot)
}


# fitted curve plot of marker
markerCurveFit = function(rds,marker,assays='RNA',mycols=NULL,order=NULL,time.dat){
	library(pROC); library(ROCR); library(tidyverse); library(broom)

	marker.dat = t(as.data.frame(rds[[assays]]@data[intersect(marker,rownames(rds[[assays]]@data)),]))
	cm = intersect(rownames(time.dat),rownames(marker.dat))
	tmp = cbind(time.dat[cm,],marker.dat[cm,])

	co = sapply(order[1:(length(order)-1)],function(x){
		roc = roc(response = ifelse(tmp$type==x, "1", "0"), predictor = as.numeric(tmp$Time))
		return(pROC::coords(roc, "best", best.method=NULL, best.weights=c(1, 0.5), transpose=FALSE)$threshold)})
	pp_list = lapply(marker,function(i){
		expr_tmp = data.frame(Gene = tmp[,i], Time= tmp$Time,type=tmp$type)
		model = loess(Gene ~ Time, data = expr_tmp, span = 0.75)
		pred = predict(model, newdata = expr_tmp$Time, se=T)

		ci = pred$se.fit * qt(0.95 / 2 + 0.5, pred$df)
		expr_tmp = data.frame(expr_tmp, fit = pred$fit, ymin = pred$fit - ci, ymax = pred$fit + ci, se = pred$se.fit)

		expr_tmp$Type = cut(expr_tmp$Time, breaks=c(min(expr_tmp$Time), as.numeric(co), max(expr_tmp$Time)), labels = order)
		p = ggplot(data=expr_tmp) + 
		  geom_line(aes(x=Time, y=fit, color=Type)) +
		  geom_ribbon(aes(x=Time, ymin = ymin, ymax = ymax, fill=Type), alpha=0.2)+
		  scale_color_manual(values=mycols)+
		  scale_fill_manual(values=mycols)+
		  theme(axis.text.x = element_text(angle = 315, hjust = 0, vjust = 0.5),
		  		panel.grid.major =element_blank(), panel.grid.minor = element_blank(), 
		  		panel.background = element_blank(),axis.line = element_line(colour = "black"), 
		  		axis.text=element_text(colour="black"), legend.position="none")+
		  labs(x="Time", y=i)
		  return(p)
	})
	return(pp_list)
}


# Sankey plot 
Sankey = function(proportion,x,y,stratum,alluvium,fill,label,cols,alpha=0.5){
	library(ggalluvial)
	p = ggplot(proportion,
       	aes(x = proportion[,x], y = proportion[,y], stratum = proportion[,stratum], 
       		alluvium = proportion[,alluvium], fill = proportion[,fill], label = proportion[,label])) +
  	   	scale_x_discrete(expand = c(.1, .1)) + scale_fill_manual(values = cols)+ylab(y)+ xlab(x)+
  	   	geom_flow() + geom_stratum(alpha = alpha) + geom_text(stat = "stratum", size = 2.5) +
  	   	theme(panel.border=element_rect(color='black', fill=NA), panel.grid.major =element_blank(), 
       		 panel.grid.minor = element_blank(), panel.background = element_blank(), 
       		 axis.line = element_line(colour = "black"),axis.text=element_text(colour="black")) 
  	return(p)
}

# boxplot of markers
markerBoxplot = function(rds,marker,assays='RNA',mycols=NULL,order=NULL,var,width=12,height=12){
        marker.dat = t(as.data.frame(as.matrix(rds[[assays]]@data[intersect(marker,rownames(rds[[assays]]@data)),])))
        tmp = data.frame(State = as.character(rds[[var]][,1]),CellName = colnames(rds), row.names=colnames(rds))
        cm = intersect(rownames(tmp),rownames(marker.dat))
        tmp = cbind(tmp[cm,],marker.dat[cm,])

        if(is.null(order)) order = as.character(unique(tmp$State))
        else tmp$State = factor(tmp$State,levels=order)

        pdf(paste0(assays,'_marker_boxplot.pdf'),useDingbats=F,width=width,height=height)
        par(mfrow=c(2,2))

        for(i in marker){
                boxplot(tmp[,i]~State, data = tmp,ylab=i,main=i,xlab='',notch=TRUE,outline=FALSE,col=mycols)
                lines(1:length(order),tapply(tmp[,i],tmp$State,median),col='black')}
        dev.off()
}


# boxplot of markers split of group
markerBoxplotGroup = function(rds,marker,assays='RNA',mycols=NULL,order=NULL,fill,var,xvar,sample){
	marker.dat = t(as.data.frame(rds[[assays]]@data[intersect(marker,rownames(rds[[assays]]@data)),]))
	tmp = data.frame(State = as.character(rds[[var]][,1]),CellName = colnames(rds), source = as.character(rds[[sample]][,1]), row.names=colnames(rds))
	cm = intersect(rownames(tmp),rownames(marker.dat))
	tmp = cbind(tmp[cm,],marker.dat[cm,])

	if(is.null(order)) order = as.character(unique(tmp$State))
	else tmp$State = factor(tmp$State,levels=order)

	p.list = lapply(marker,function(x){
		p = ggplot(tmp, aes(x = tmp[,xvar], y = tmp[,x])) + 
			geom_boxplot(aes(fill = tmp[,fill]),position=position_dodge(0.8),width=0.6) + 
			scale_fill_manual(values = mycols) + ylab(x) + 
				theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
				  	  panel.background=element_blank(),axis.line=element_line(colour="black"))
		return(p)})
	return(p.list)
}

### barplot of markers
markerBarplot = function(rds,marker,assays='RNA',mycols=NULL,order=NULL,var,cutoff){
	marker.dat = t(as.data.frame(rds[[assays]]@data[intersect(marker,rownames(rds[[assays]]@data)),]))
	tmp = data.frame(State = as.character(rds[[var]][,1]),CellName = colnames(rds), row.names=colnames(rds))
	cm = intersect(rownames(tmp),rownames(marker.dat))
	tmp = cbind(tmp[cm,],marker.dat[cm,])

	if(is.null(order)) order = as.character(unique(tmp$State))
	else tmp$State = factor(tmp$State,levels=order)

	p.list = lapply(marker,function(x){
		por = data.frame(tmp,Gene = ifelse(as.numeric(as.character(tmp[,x]))>cutoff,1,0)) %>% 
						group_by(State) %>% summarise(all = n(),expr_count = sum(Gene))
		por$percentage = por$expr_count/por$all
		p = ggplot(por, aes(x=State, y=percentage)) + geom_bar(stat="identity", position="dodge", aes(fill=State)) + 
			scale_fill_manual(values=mycols) + ylab(x) + 
			theme(panel.border=element_rect(color='black', fill=NA), panel.grid.major =element_blank(), 
       		 panel.grid.minor = element_blank(), panel.background = element_blank(), 
       		 axis.line = element_line(colour = "black"),axis.text=element_text(colour="black")) 
		return(p)
	})
	return(p.list)
}


### survival curve
survival.plot = function(dat,marker,time,status,group.stand='median',xlim=50,high.pct=0.5,low.pct=0.5,noIntermediate=TRUE){
  library(survival); library(survminer)
  if(group.stand=='median') dat$Group = ifelse(dat[,marker]>=median(dat[,marker]),'High','Low')
  if(group.stand=='mean') dat$Group = ifelse(dat[,marker]>=mean(dat[,marker]),'High','Low')
  if(group.stand=='quantile'){
    dat = dat[order(dat[,marker],decreasing = T),]
    high = ceiling(high.pct*nrow(dat))
    if(high.pct+low.pct==1) low = nrow(dat) - high
    else low = ceiling(low.pct*nrow(dat))
    dat$Group = c(rep('High',times = high),rep('Intermediate',times = (nrow(dat)-high-low)),rep('Low',times = low))}
  if(noIntermediate) dat = dat[which(dat$Group %in% c('High','Low')),]
  dat$Time = as.numeric(as.character(dat[,time])) ; dat$Status = as.numeric(as.character(dat[,status]))
  fit = survfit(Surv(Time,Status) ~ Group,data = dat)
  if(length(unique(dat$Group))==2 & length(grep('High',dat$Group))>0 & length(grep('Low',dat$Group))>0){
    p = ggsurvplot(fit, data = dat, 
                   risk.table = TRUE,pval = TRUE,
                   risk.table.height=0.25,xlim=c(0,xlim),
                   palette = c('red','blue'),
                   legend.labs = c('High','Low'),
                   ggtheme =  theme(panel.grid.major=element_blank(),
                                    panel.grid.minor=element_blank(),
                                    panel.background=element_blank(),
                                    axis.line=element_line(colour="black"))) + ggtitle(marker)
  }else{
    p = ggsurvplot(fit, data = dat, 
                   risk.table = TRUE,pval = TRUE,
                   risk.table.height=0.25,xlim=c(0,xlim),
                   ggtheme =  theme(panel.grid.major=element_blank(),
                                    panel.grid.minor=element_blank(),
                                    panel.background=element_blank(),
                                    axis.line=element_line(colour="black"))) + ggtitle(marker)
  }
  
  return(p)
}


### relapse curve
relapse.plot = function(dat,marker,time,status,group.stand='median',xlim=50,high.pct=0.5,low.pct=0.5,noIntermediate=TRUE){
  library(survival); library(survminer)
  if(group.stand=='median') dat$Group = ifelse(dat[,marker]>=median(dat[,marker]),'High','Low')
  if(group.stand=='mean') dat$Group = ifelse(dat[,marker]>=mean(dat[,marker]),'High','Low')
  if(group.stand=='quantile'){
    dat = dat[order(dat[,marker],decreasing = T),]
    high = ceiling(high.pct*nrow(dat))
    if(high.pct+low.pct==1) low = nrow(dat) - high
    else low = ceiling(low.pct*nrow(dat))
    dat$Group = c(rep('High',times = high),rep('Intermediate',times = (nrow(dat)-high-low)),rep('Low',times = low))}
  if(noIntermediate) dat = dat[which(dat$Group %in% c('High','Low')),]
  dat$Time = as.numeric(as.character(dat[,time])) ; dat$Status = as.numeric(as.character(dat[,status]))
  fit = survfit(Surv(Time,Status) ~ Group,data = dat)
  if(length(unique(dat$Group))==2 & length(grep('High',dat$Group))>0 & length(grep('Low',dat$Group))>0){
    p = ggsurvplot(fit, data = dat, surv.median.line = "hv",fun = "event",
                   risk.table = TRUE,pval = TRUE,
                   risk.table.height=0.25,xlim=c(0,xlim),
                   palette = c('red','blue'),
                   legend.labs = c('High','Low'),
                   ggtheme =  theme(panel.grid.major=element_blank(),
                                    panel.grid.minor=element_blank(),
                                    panel.background=element_blank(),
                                    axis.line=element_line(colour="black"))) + ggtitle(marker)
  }else{
    p = ggsurvplot(fit, data = dat, surv.median.line = "hv",fun = "event",
                   risk.table = TRUE,pval = TRUE,
                   risk.table.height=0.25,xlim=c(0,xlim),
                   ggtheme =  theme(panel.grid.major=element_blank(),
                                    panel.grid.minor=element_blank(),
                                    panel.background=element_blank(),
                                    axis.line=element_line(colour="black"))) + ggtitle(marker)
  }
  
  return(p)
}


### Monocle trajectory analysis
newMonocleAnalysis = function(seurat_object,assays='integrated',savePath,cols, ordering_genes=NULL){

    #Extract data, phenotype data, and feature data from the SeuratObject
    data = as(as.matrix(seurat_object[[assays]]@data), 'sparseMatrix')
    pd = new('AnnotatedDataFrame', data = seurat_object@meta.data)
    fData = data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    fd = new('AnnotatedDataFrame', data = fData)

    #Construct monocle cds
    monocle = newCellDataSet(data,phenoData = pd,featureData = fd,expressionFamily = uninormal())

    #Run ordering algorithm
    if(!is.null(ordering_genes)) ordering_genes = ordering_genes
    if(is.null(ordering_genes)) ordering_genes = seurat_object[[assays]]@var.features

    monocle = setOrderingFilter(monocle, ordering_genes)

    # reduce dimension 
    monocle = reduceDimension(monocle,norm_method="none", reduction_method="DDRTree",
                            max_components=3,pseudo_expr = 0,verbose=TRUE)
    # order cells
    monocle = orderCells(monocle)

    pdf(savePath,useDingbats=F,width=10,height=10)
    print(plot_cell_trajectory(monocle,  color_by = "seurat_clusters",theta = 10,show_branch_points = FALSE,
                         show_tree = TRUE,cell_size = 1.5) + 
    		scale_color_manual(breaks = c("X", "Y", "Z"), values=cols) + theme(legend.position = "right"))
    dev.off()

    return(monocle)
}

### calculate entropy
Entropy <- function(d){
  res <- 0
  for(i in 1:length(d)){
    if(d[i]!=0) res <- res + d[i]*log(d[i])}
  return (-res)
}



### cellcycle score: generate the cellcycle signatures list
ccScore.generate = function(dat.gene,index1,index2){
	rs = list()
	for(i in unique(dat.gene[,index1])){
		rs[[i]] = as.character(dat.gene[which(dat.gene[,index1]==i),index2])}
	return(rs)
}

### cellcycle score：calculate cellcycle scores 
ccScore.analyze = function(dat.expr,dat.gene.list,k=10,bin=25,controlSize=100){

	gene.mean = data.frame(Gene = rownames(dat.expr), GeneMean = apply(dat.expr,1,mean),row.names=rownames(dat.expr))
	gene.mean = gene.mean[order(gene.mean$GeneMean,decreasing=T),]

	message('Splicing bin...')
	binlength = floor(nrow(gene.mean)/bin)
	tmp = c()
	for(j in 1:bin){
		if(j != bin){ tmp = c(tmp,rep(paste0('BIN',j),times=binlength))}
		else { tmp = c(tmp,rep(paste0('BIN',j),times=(nrow(dat.expr)-(bin-1)*binlength)))}}
	gene.mean$BIN = tmp

	rs = c()
	for(mod in names(dat.gene.list)){

		dat.gene = dat.gene.list[[mod]]
		set.seed(12345)
		cellcycleScore = c() # row-cell;col-cellcyclescoreGene
		dat.gene2 = c()
		for(gene in dat.gene){
			if(length(intersect(gene,rownames(dat.expr)))>0){
				dat.gene2 = c(dat.gene2,gene)
				index = gene.mean$BIN[which(gene.mean$Gene==gene)]
				rangeGene = gene.mean$Gene[which(gene.mean$BIN==index)]
				scoreCell = c()
				for(kfold in 1:k){
					controlGene = sample(setdiff(rangeGene,gene),size=controlSize,replace=F)
					scoreCellTmp = apply(dat.expr[c(gene,controlGene),],2,function(x){return(mean(x[1])-mean(x[2:length(x)]))})
					scoreCell = cbind(scoreCell,scoreCellTmp)}	
				cellcycleScore = cbind(cellcycleScore,apply(scoreCell,1,mean))}}
		colnames(cellcycleScore) = dat.gene2
		print(length(dat.gene2))
		if(is.null(rs)) rs = data.frame(ccScore=apply(cellcycleScore,1,mean),row.names=rownames(cellcycleScore))
		else rs = cbind(rs,data.frame(ccScore=apply(cellcycleScore,1,mean),row.names=rownames(rs)))
	}
	colnames(rs) = paste0(names(dat.gene.list),'_raw')
	rs.norm = as.data.frame(apply(rs,1,function(x){(x-min(x))/(max(x)-min(x))}))
	if(dim(rs.norm)[1]!=dim(rs)[1] | dim(rs.norm)[2]!=dim(rs)[2]) rs.norm = t(rs.norm)
	colnames(rs.norm) = gsub('_raw','_norm',colnames(rs))
	
	return(cbind(rs,rs.norm))
}


### estimate the stemness for each cell
CytoTRACE <- function(mat, batch = NULL, enableFast = TRUE,ncores = 1,subsamplesize = 1000){

  range01 <- function(x){(x-min(x))/(max(x)-min(x))}

  #inputs
  a1 <- mat
  a2 <- batch
  if(ncol(mat) < 3000){
    enableFast = FALSE
    message("The number of cells in your dataset is less than 3,000. Fast mode has been disabled.")
    } else {
  message("The number of cells in your dataset exceeds 3,000. CytoTRACE will now be run in fast mode (see documentation). You can multi-thread this run using the 'ncores' flag. To disable fast mode, please indicate 'enableFast = FALSE'.")
  }
  #Checkpoint: NAs and poor quality genes
  pqgenes <- is.na(rowSums(mat>0)) | apply(mat, 1, var) == 0
  num_pqgenes <- length(which(pqgenes == TRUE))
  mat <- mat[!pqgenes,]
  if(num_pqgenes>0){
    warning(paste(num_pqgenes, "genes have zero expression in the matrix and were filtered"))
  }

  #Subsample routine
  if(enableFast == FALSE){
    size <- ncol(mat)
  } else if (enableFast == TRUE & subsamplesize < ncol(mat)){
    size <- subsamplesize
  } else if (enableFast == TRUE & subsamplesize >= ncol(mat)){
    stop("Please choose a subsample size less than the number of cells in dataset.")
  }

  chunk <- round(ncol(mat)/size)
  subsamples <- split(1:ncol(mat), sample(factor(1:ncol(mat) %% chunk)))
  message(paste("CytoTRACE will be run on", chunk, "sub-sample(s) of approximately",
                round(mean(unlist(lapply(subsamples, length)))), "cells each using", min(chunk, ncores),"/", ncores, "core(s)"))

  message(paste("Pre-processing data and generating similarity matrix..."))
  batches <- parallel::mclapply(subsamples, mc.cores = min(chunk, ncores), function(subsample){
    #Checkpoint: log2-normalization
    mat <- mat[,subsample]
    batch <- batch[subsample]

    if(max(mat)<50){
      mat <- 2^mat - 1
    }

    #Checkpoint: ERCC standards
    if(length(grep("ERCC-", rownames(mat)))>0){
      mat <- mat[-grep("ERCC-", rownames(mat)),]
    }

    #Checkpoint: Sequencing depth normalization
    mat <- t(t(mat)/apply(mat, 2, sum))*1000000

    #Checkpoint: NAs and poor quality cells
    pqcells <- is.na(apply(mat>0, 2, sum)) | apply(mat>0, 2, sum) <= 10
    num_pqcells <- length(which(pqcells == TRUE))
    mat <- mat[,!pqcells]

    #Checkpoint: log2-normalize
    mat <- log(mat+1,2)
    mat <- data.matrix(mat)

    #Calculate pre-batch corrected gene counts
    counts <- apply(mat>0, 2, sum)
    #Checkpoint: Batch correction
    if(ncol(a1) == length(a2)){
      #filter poor quality cells from batch vector
      batch <- batch[!pqcells]

      #Run Combat
      suppressMessages(mat <- sva::ComBat(mat, batch, c()))
      mat <- data.matrix(mat)

      #Replace negative values after batch correction with zeroes for compatibility with downstream steps
      mat[which(mat<0)] <- 0
    }
    #Rescale each single cell with gene counts to convert relative transcript abundances to absolute RNA content prior to cell lysis (credit: Census, Qiu et al., 2017)
    census_normalize <- function(mat, counts) {
      xnl <- 2^data.matrix(mat) - 1
      rs <- apply(xnl, 2, sum)
      rnorm <- t(t(xnl) * counts/rs)
      A <- log(rnorm+1,2)
      return(A)
    }

    mat2 <- census_normalize(mat, counts)
    #Function to identify the most variable genes
    mvg <- function(matn) {
      A <- matn
      n_expr <- rowSums(A > 0);
      A_filt <- A[n_expr >= 0.05 * ncol(A),];
      vars <- apply(A_filt, 1, var);
      means <- apply(A_filt, 1, mean);
      disp <- vars / means;
      last_disp <- tail(sort(disp), 1000)[1];
      A_filt <- A_filt[disp >= last_disp,];

      return(A_filt)
    }

    #Filter out cells not expressing any of the 1000 most variable genes
    mat2.mvg <- mvg(mat2)
    rm1 <- colSums(mat2.mvg) == 0
    mat2 <- mat2[, !rm1]
    counts <- counts[!rm1]

    #Calculate similarity matrix
    similarity_matrix_cleaned <- function(similarity_matrix){
      D <- similarity_matrix
      cutoff <- mean(as.vector(D))
      diag(D) <- 0;
      D[which(D < 0)] <- 0;
      D[which(D <= cutoff)] <- 0;
      Ds <- D
      D <- D / rowSums(D);
      D[which(rowSums(Ds)==0),] <- 0
      return(D)
    }
    D <- similarity_matrix_cleaned(HiClimR::fastCor(mvg(mat2)))

    return(list(mat2 = mat2,counts = counts, D = D))
  }
  )
  #Prepare for downstream steps
  mat2 <- do.call(cbind, lapply(batches, function(x) x$mat2))
  counts <- do.call(c, lapply(batches, function(x) x$counts))
  filter <- colnames(a1)[-which(colnames(a1) %in% colnames(mat2))]
  if(length(filter)>0){
    warning(paste(length(filter), "poor quality cells were filtered based on low or no expression. See 'filteredCells' in returned object for names of filtered cells."))
  }
  #Calculate gene counts signature (GCS) or the genes most correlated with gene counts
  message("Calculating gene counts signature...")
  ds2 <- sapply(1:nrow(mat2), function(x) ccaPP::corPearson(mat2[x,],counts))
  names(ds2) <- rownames(mat2)
  gcs <- apply(mat2[which(rownames(mat2) %in% names(rev(sort(ds2))[1:200])),],2,mean)

  samplesize <- unlist(lapply(lapply(batches, function(x) x$counts), length))
  gcs2 <- split(gcs, as.numeric(rep(names(samplesize), samplesize)))
  D2 <- lapply(batches, function(x) x$D)

  #Regress gene counts signature (GCS) onto similarity matrix
  regressed <- function(similarity_matrix_cleaned, score){
    out <- nnls::nnls(similarity_matrix_cleaned,score)
    score_regressed <- similarity_matrix_cleaned %*% out$x
    return(score_regressed)
  }

  #Apply diffusion to regressed GCS using similarity matrix
  diffused <- function(similarity_matrix_cleaned, score, ALPHA = 0.9){
    vals <- score
    v_prev <- rep(vals);
    v_curr <- rep(vals);

    for(i in 1:10000) {
      v_prev <- rep(v_curr);
      v_curr <- ALPHA * (similarity_matrix_cleaned %*% v_curr) + (1 - ALPHA) * vals;

      diff <- mean(abs(v_curr - v_prev));
      if(diff <= 1e-6) {
        break;
      }
    }
    return(v_curr)
  }

  message("Smoothing values with NNLS regression and diffusion...")
  cytotrace <- parallel::mclapply(1:length(D2), mc.cores = ncores, function(i) {
    gcs_regressed <- regressed(D2[[i]], gcs2[[i]])
    gcs_diffused <- diffused(D2[[i]], gcs_regressed)
    cytotrace <- rank(gcs_diffused)
  }
  )

  cytotrace <- cytotrace_ranked <- unlist(cytotrace)
  cytotrace <- range01(cytotrace)

  #Calculate genes associated with CytoTRACE
  cytogenes <- sapply(1:nrow(mat2),
                          function(x) ccaPP::corPearson(mat2[x,], cytotrace))
  names(cytogenes) <- rownames(mat2)
  message("Calculating genes associated with CytoTRACE...")

  #Final steps
  names(cytotrace) <- names(gcs) <- names(counts) <- colnames(mat2)
  cytotrace <- cytotrace[colnames(a1)]; cytotrace_ranked <- cytotrace_ranked[colnames(a1)]; gcs <- gcs[colnames(a1)]; counts <- counts[colnames(a1)]

  mat2 <- t(data.frame(t(mat2))[colnames(a1),])
  names(cytotrace) <- names(cytotrace_ranked) <- names(gcs) <- names(counts) <- colnames(mat2) <- colnames(a1)

  message("Done")
  return(list(CytoTRACE = cytotrace, CytoTRACErank = cytotrace_ranked, cytoGenes = sort(cytogenes, decreasing = T), GCS = gcs, gcsGenes = sort(ds2, decreasing = T),
              Counts = counts, filteredCells = filter, exprMatrix = mat2))
}



plot_Circ_heatmap <- function(mat, cluster, distmethod = 'euclidean') {
    # ref: https://zhuanlan.zhihu.com/p/136138642
    #@mat: row(sample or groups) X col(pair info or genes)
    #@cluster: set cluster number in hclust function

    library(dendextend)
    library("circlize")
    library(RColorBrewer)

    mat=scale(mat, center = TRUE, scale = TRUE)
    dend <-as.dendrogram(hclust(dist(t(mat),method = distmethod)))
    n=1
    dend <-dend %>% set("branches_k_color", k = n) 
    #par(mar=c(7.5,3,1,0))
    #plot(dend)
    # 聚类后的样本信息
    mat2 = mat[, order.dendrogram(dend)]
    lable1=row.names(mat2)
    lable2=colnames(mat2)
    nr = nrow(mat2)
    nc = ncol(mat2)
    col_fun = colorRamp2(c(-1.5, 0, 1.5), c("skyblue", "white", "red"))
    col_mat = col_fun(mat2)
    par(mar=c(0,0,0,0))
    circos.clear()
    circos.par(
        canvas.xlim =c(-2,2),
        canvas.ylim = c(-2,2),
        cell.padding = c(0,0,0,0), 
        gap.degree =90
    )
    factors = "a"
    circos.initialize(factors, xlim = c(0, ncol(mat2)))
    circos.track(
        ylim = c(0, nr),bg.border = NA,track.height = 0.05*nr, 
        panel.fun = function(x, y) {
            for(i in 1:nr) {
                circos.rect(xleft = 1:nc - 1, ybottom = rep(nr - i, nc),
                    xright = 1:nc, ytop = rep(nr - i + 1, nc),
                    border = "black",
                    col = col_mat[i,]
                )
                circos.text(x = nc,
                    y = 6.4 -i,
                    labels = lable1[i],
                    facing = "downward", niceFacing = TRUE,
                    cex = 0.6,
                    adj = c(-0.2, 0))
            }
        }
    )
    for(i in 1:nc){
        circos.text(x = i-0.04,
        y = 11,
        labels = lable2[i],
        facing = "clockwise", niceFacing = TRUE,
        cex = 0.5,adj = c(0, 0))
    }
    #添加树
    max_height <-max(attr(dend, "height"))
    circos.track(ylim = c(0, max_height),bg.border = NA,track.height = 0.3, 
        panel.fun = function(x, y){
        circos.dendrogram(dend = dend,
        max_height = max_height)
    })
    circos.clear()
    # 添加图例
    library(ComplexHeatmap)
    lgd <- Legend(at = c(-2,-1, 0, 1, 2), col_fun = col_fun, 
                title_position = "topcenter",title = "Z-score")
    draw(lgd, x = unit(0.7, "npc"), y = unit(0.7, "npc"))
}

# single-cell Analysis

library(Seurat)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(ssgsea.GBM.classification)
library(monocle)
library(infercnv)
library(reshape2)
# library(SCANER)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggthemes);library(ggpubr); library(SingleR)
library(SingleCellExperiment);library(tidyverse);library(reticulate);library(tidyverse)
library(ggplot2);library(openxlsx);library(ggsci);library(ROGUE)
library(Seurat);library(ROCR);library(cluster);library(parallel)
library(rlang); library(ggplot2); library(ggthemes); library(glue);library(ggalluvial)
library(fmsb);library(harmony)
source('~/R/ssgseaMOD.r')
source('~/R/ccScore.R')

suppressMessages(library(ggplot2))
suppressMessages(library(tidyverse))