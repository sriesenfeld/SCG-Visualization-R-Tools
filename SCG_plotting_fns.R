## S.J. Riesenfeld
## Code for visualization and analysis of single-cell genomic data.
## Available here: https://github.com/sriesenfeld/SCG-Visualization-R-Tools
## Under the BSD 3-Clause "New" or "Revised" License

library("ggpubr") # ggarrange() wrapper for cowplot::plot_grid()
library("cowplot") # ggdraw(), plot_grid() 
library("gridExtra") # marrangeGrob()
library("scales") #viridis_pal()
library("reshape2") #melt()
library("plyr")
library("dplyr")
# library("Seurat") # Not required but typical

###### Includes code for top level functions
# barp.cnts(), plotViolins(), image_matrix(), plot.embed.fac(), plot.embed.genes(), plot.embed.features()
###### and the utility functions
# get.colors.wcap(), get.viridis.col(), multiPanelPlot(), and high.loadings()
######

# barp.cnts(): plots data as a ggplot2 stacked bar plot of either counts or proportions.
## ARGUMENTS
# groupvars: vector of up to 3 names of factors in data, excluding color.by,
## or 4 names including color.by, for splitting the groups
# color.by: the name of a factor in data for the main grouping variable
# colors.v: optional named list of colors corresponding to the levels in the color.by factor 
# do.prop: Boolean; if TRUE, plot proportions, otherwise counts
# title: string giving plot title
# bar.width: numeric value passed as width to geom_bar()
# show.legend: Boolean; if FALSE, the color legend is removed
# n.col: number of panel columns desired (ignored if length(groupvars)>2
# n.col.max: maximum number of panel columns desired (ignored if length(groupvars)>2) 
barp.cnts <- function(data, groupvars, color.by=groupvars[1], 
                      colors.v=NULL, do.prop=F,
                      title=NULL,
                      bar.width=0.7, show.legend=TRUE, n.col=NULL, n.col.max=10) {
  if (length(groupvars)>3) { warning(paste0("Only the first 3 variables in groupvars will be used"))}
  gvars=setdiff(groupvars, color.by)
  groupvars=c(color.by, gvars)
  if (is.null(title)) { 
    title=paste0(color.by, ifelse(do.prop, " proportions", " counts"), 
                 " within ", gvars[1],
                 ifelse(length(gvars)>1, 
                        paste0(", split by ", paste0(
                          gvars[2:min(length(gvars),3)], collapse=", ")), "")) 
  }
  if (! all (groupvars %in% colnames(data))) {
    stop(paste0("Cannot find ", paste0(groupvars[!groupvars%in%colnames(data)], collapse=", "), 
                "in column names of data"))
  }
  if (! all (sapply(groupvars, function(x) { is.factor(data[,x])}))) {
    stop("Expecting all variables in groupvars and color.by to name factors in data")
  }
  p=ggplot(data, aes_string(x=gvars[1], fill=color.by))
  if (do.prop) {
     p=p + geom_bar(position=position_fill(), width=bar.width) + ylab("proportion")
  } else {
    p=p + geom_bar(position=position_stack(), width=bar.width) 
  }
  if (!is.null(colors.v)) { p= p+  scale_fill_manual(values=colors.v) }
  if (length(gvars)==2) {
      n.levels=nlevels(data[,gvars[2]])
      if (is.null(n.col)) { n.col=ifelse(n.levels<n.col.max, n.levels, 
                                         min(n.col.max, ceiling(sqrt(n.levels)))) }
      p=p+facet_wrap(as.formula(paste0("~", gvars[2])), 
                     ncol=n.col) 
  } else if (length(gvars)>2) {
    p=p+facet_grid(as.formula(paste0(gvars[3], "~", gvars[2])))
  }
  p= p+  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) + ggtitle(title)
  if (!show.legend) { p = p+theme(legend.position="none")}
  return(p)
}

# plotViolins(): function for plotting numeric data as violins split by groups
## ARGUMENTS:
# measurevars: vector of names of rows in num.data corresponding to numeric variables to be plotted 
# groupvars: vector of names of columns in metadata corresponding to grouping factors
# ob: optional Seurat object
# num.data: numeric matrix, whose column names match the row names of metadata
# metadata: data.frame of features of the column names of num.data
# measurevar.name: name of the type of variable being plotted, e.g., "gene"
# measure.name: name of the type of numeric data being plotted, e.g., "expression"
# color.by: name of factor in metadata by which the violins should be colored 
# colors.v: named vector of colors corresponding to the levels of color.by factor
# color.def: default color of violins if color.by is NULL
# quantiles: vector quantiles to be plotted on top of violins, passed as draw_quantiles to
##  geom_violin(); can be finicky, default of the median usually works
# adjust: numeric value passed as adjust to geom_violin()
# scale: passed as scale to geom_violin()
# facet.scales: passed as scales to facet_wrap()
# plot.means: Boolean; if TRUE, plot mean values on top of violins
# means.shape, means.size, means.fill, means.color: passed as 
##  shape, size, fill, and color, respectively to geom_point() for plotting mean values
# n.col: desired number of panel columns (ignored if there is more than one facet)
# n.col.max: max desired number of panel columns (ignored if there is more than one facet)
# strip.text.size, strip.text.face: font size and face for panel titles
# axis.title.size: font size for axis titles
# axis.x.text.size, axis.y.text.size: font sizes for x and y axis text, respectively
# rotate.axis.text.x: Boolean, if TRUE, rotate x axis text
# panel.spacing: if not NULL, passed to theme() as panel.spacing=unit(panel.spacing, "points")
# size.geom: passed as size to ggplot()
# cc.ylim: if not NULL, ylimites of coordinates are adjusted by coord_cartesian(ylim=cc.ylim) 
# no.legend: Boolean, if TRUE, remove color legend
# base_theme: ggplot2 base theme to use
# no.x.axis.tick.labels: Boolean, if TRUE, remove x axis tick labels
plotViolins<- function(measurevars, groupvars, ob=NULL,
                       num.data=GetAssayData(ob), metadata=ob[[]],
                       measurevar.name="variable", measure.name="value",
                       color.by=NULL, colors.v=NULL, color.def="gray50", 
                       quantiles=c(0.5), adjust=1.25, scale="width", facet.scales="free_y", 
                       plot.means=T, #plot.means.df=NULL, 
                       means.shape=23, means.size=1.5, 
                       means.fill="white", means.color="black",
                       n.col=NULL, n.col.max=10, # only applies if there is at most 1 facet
                       strip.text.size=10, strip.text.face="italic", 
                       axis.title.size=8, 
                       axis.x.text.size=8, axis.y.text.size=8, 
                       rotate.axis.text.x=T, #rotate.axis.text.y=F, 
                       panel.spacing=NULL, size.geom=0.5, 
                       cc.ylim=NULL, 
                       no.legend=F, base_theme=theme_bw(),
                       no.x.axis.tick.labels=((!no.legend) && (!is.null(color.by)) 
                                              && (color.by %in% groupvars))) {
  if (!all(measurevars %in% rownames(num.data))) { 
    warning(paste0("Cannot find ", paste0(measurevars[!measurevars %in% rownames(num.data)], collapse=", "), 
                   " in row names of num.data, so removing it")) 
    measurevars=measurevars[measurevars%in% rownames(num.data)]
  }
  if (length(measurevars)==0) { stop("Need measurevars to have length at least 1")}
  if (!all(colnames(num.data)==rownames(metadata))) { 
    stop("Column names of num.data should match row names of metadata")}
  if (!all(groupvars %in% colnames(metadata))) { 
    stop(paste0("Cannot find ", paste0(groupvars[!groupvars %in% colnames(metadata)], collapse=", "),
                " in column names of metadata" )) }
  if (!all(sapply(metadata[,groupvars,drop=FALSE], is.factor))) { 
    stop(paste0("Expecting ", groupvars[!sapply(metadata[,groupvars,drop=FALSE], is.factor)], 
                " to name factors in metadata")) }
  if (!is.null(color.by)) { 
    if (!color.by %in% colnames(metadata)) { 
      stop(paste0("Cannot find ", color.by ," in column names of metadata" )) }
    if (!is.factor(metadata[,color.by])) { stop("Expecting ", color.by , " to name a factor in metadata") }
  }
  xvar=measurevar.name; groupvars.fct=NULL
  if (length(groupvars)>0) { 
    xvar=groupvars[1] 
    if (length(measurevars)>1) {
      groupvars.fct=measurevar.name
      if (length(groupvars)>=2) { 
        if (length(groupvars)>2) { warning("If measurevars has length >1, then only the first two entries of groupvars are used") }
        groupvars.fct=c(groupvars[2], measurevar.name)
      } 
    } else if (length(groupvars)>1) {
      if (length(groupvars) > 3) { warning("Only the first three entries of groupvars are used") }
      groupvars.fct=groupvars[2:min(3,length(groupvars))]
    }
  } 
  regroup.vars=setdiff(unique(c(xvar, groupvars.fct, color.by)), measurevar.name)
  data=melt(data.frame(t(num.data[measurevars,,drop=FALSE]),
                       metadata[,regroup.vars,drop=FALSE]), id.vars=regroup.vars)
  colnames(data)[colnames(data)=="variable"]=measurevar.name
  colnames(data)[colnames(data)=="value"]=measure.name
  #if (is.null(plot.means.df)) {
  plot.means.df <- summarySE(data, measurevar=measure.name, groupvars=c(regroup.vars, measurevar.name))
  #}  
  if (is.null(color.by)) {
    p = ggplot(data, aes_string(x=xvar, y=measure.name), fill=color.by, size=size.geom)
  } else {
    p = ggplot(data, aes_string(x=xvar, y=measure.name, fill=color.by), size=size.geom) 
  }
  p=p+geom_violin(adjust=adjust, scale=scale, draw_quantiles=quantiles)
  if (!is.null(colors.v)) { p = p + scale_fill_manual(values=colors.v) }
  p=p+base_theme
  if (length(groupvars.fct)==1) { 
    n.levels=nlevels(data[,groupvars.fct])
    if (is.null(n.col)) { n.col=ifelse(n.levels<n.col.max, n.levels, 
                                       min(n.col.max, ceiling(sqrt(n.levels)))) }
    p=p+facet_wrap(groupvars.fct[1], ncol=n.col, scales=facet.scales) 
  } else if (length(groupvars.fct)>1) { 
    p=p+facet_grid(paste0(groupvars.fct[2], "~", groupvars.fct[1]),  scales=facet.scales) 
  }
  if (plot.means) {
    p= p+ geom_point(data=plot.means.df, stat="identity", color=means.color, fill=means.fill, shape=means.shape, size=means.size) 
  }
  if (rotate.axis.text.x){ p=p+ theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) }
  p=p+theme(axis.text.x= element_text(size=axis.x.text.size), 
            axis.text.y= element_text(size=axis.y.text.size),
            axis.title=element_text(size=axis.title.size),
            strip.text = element_text(size = strip.text.size, face=strip.text.face))
  if (no.x.axis.tick.labels) {
    p=p+ theme(axis.text.x = element_blank())
  }
  if (!is.null(panel.spacing)) { 
    p=p+ theme(panel.spacing=unit(panel.spacing, "points"))
  }
  if (no.legend) { p=p+theme(legend.position="none") }
  if (!is.null(cc.ylim)) { 
    p=p+coord_cartesian(ylim=cc.ylim) } ## does a zoom without throwing out any data
  return(p)
}


# image_matrix(): Function for plotting a matrix as an image
## ARGUMENTS:
# mm: numeric or factor matrix
# colors.v: color values
# no.axis.text.x; Boolean, if TRUE, remove x axis text
# no.axis.text.y: Boolean, if TRUE< remove y axis text
# rotate.axis.text.x: Boolean, if TRUE rotate x axis text
# axis.pos.x: string specifying the x axis position
# categ: Boolean, if TRUE, treat data in mm as discrete
# xlab, ylab: strings giving x and y axis labels, respectively
# title: string giving plot title
image_matrix <- function(mm, colors.v=NULL,
                         no.axis.text.x=FALSE, no.axis.text.y=TRUE,
                         rotate.axis.text.x=TRUE,
                         axis.pos.x=c("top", "bottom"),
                         categ=TRUE, xlab="", ylab="", title="") {
  axis.pos.x=match.arg(axis.pos.x)
  mm.df=do.call(rbind, lapply(1:ncol(mm), function(i) {
    return(data.frame(
      "cols"=rep(i,nrow(mm)),
      "rows"=rev(1:nrow(mm)),
      "value"=mm[,i])) })) 
  print(dim(mm.df))
  mm.df$cols=factor(mm.df$cols, labels=colnames(mm))
  mm.df$rows=factor(mm.df$rows, labels=rownames(mm))
  if (categ) { 
    mm.df$value=factor(mm.df$value) 
  }
  p=ggplot(mm.df, aes(x=cols,y=rows,fill=value))
  if (!is.null(colors.v)) {
    if (categ) { 
      p=p+ scale_fill_manual(values=colors.v)
    } else {
      p=p+ scale_fill_gradientn(colors = colors.v)
    }
  }
  p=p+scale_x_discrete(name=xlab, expand=expansion(0), position=axis.pos.x) +
    scale_y_discrete(name=ylab, expand=expansion(0))
  p=p+geom_raster() + theme_bw()  
  if (no.axis.text.x) {
    p=p+ theme(axis.text.x = element_blank()) }
  if (no.axis.text.y) {
    p=p+ theme(axis.text.y = element_blank()) }
  if (rotate.axis.text.x){ 
    hjust=ifelse(axis.pos.x=="top", 0, 1)
    vjust=ifelse(axis.pos.x=="top", 0, 1)
    p=p+ theme(axis.text.x = 
                 element_text(angle=60, vjust=vjust, hjust=hjust)) }
  p=p+ggtitle(title)
  return(p)
}

# plot.embed.fac(): plot a factor (i.e., categorical variable) on a 2d embedding
##  (optionally) as a grid of plots highlighting individual factor levels. 
##  Need at minimum fac.v and embed.xy to be specified.
## ARGUMENTS:
# ob: Seurat object (optional; if NULL, then fac.v and embed.xy must be given)
# fac.name: name of factor in the metadata of the Seurat object 
##  (if NULL, then fac.v must be provided)
# fac.v: factor data to be plotted, same length as number of rows of embed.xy 
# embed.xy: matrix or data.frame of embedding coordinates (named columns 1 and 2)
# do.facet.gray:  Boolean; if TRUE, then plot each level of the factor separately 
##  against grayed out background; if FALSE, plot factor levels together; 
# facet.levels: levels of the factor to include in the plot
# colorsv: vector of colors of same length as facet.levels, used only if 
##  do.facet.gray==FALSE
# color.gray: color for the non-highlighted factor level, used only if 
##  do.facet.gray==TRUE
# color.highlight: color for the highlighted factor level, used only if 
##  do.facet.gray==TRUE
# return.plot.list: Boolean; if TRUE, then return list of individual factor
##  level ggplots rather than one that combines those into a single plot
# do.save: Boolean; save the combined plot into file given by filename argument
# filename: name of file, used only if do.save==TRUE, passed to multiPanelPlot
# do.label: Boolean; if TRUE, factors are labeled on the plot
# label.size: Passed to geom_text if factors are labeled on the plot
# alpha, size, shape: Passed to geom_point
# leg.title: Passed as title to guide_legend 
# no.legend: Boolean; if TRUE, include a legend
# leg.wid: Passed as legend.key.width to theme
# leg.ncol: Passed as ncol to guide_legend
# leg.override.size: Passed as size in override.aes to guide_legend
# leg.pos: Passed as legend.position to theme
# leg.dir: Passed as legend.direction to theme
# text.rel.size: Relative size of axis and legend title and text
# panel.title.rel.size: Relative size of panel title text
# x.lab, y.lab: x axis label, y axis label
# leg.labels: legend labels
# xlim, ylim: if not NULL, passed to xlim and ylim respectively
# base_theme_min: Boolean; if TRUE, use theme_light as base theme rather than
##  default theme
# common.legend, title, subtitle, title.size, subtitle.size, base_height, 
##  base_asp, rel.height.title, align, n.col.max, n.row.max, n.col, and n.row: 
##  Passed to multiPanelPlot()
plot.embed.fac = function(ob=NULL, fac.name=NULL, fac.v=ob[[fac.name]][,1],  
                         embed.xy=Embeddings(ob, "umap"), do.facet.gray=TRUE, 
                         facet.levels=levels(fac.v), colorsv=NULL, 
                         color.gray="gray85", color.highlight="orangered1", 
                         return.plot.list=TRUE, 
                         do.save=TRUE, filename="./fac.multipanel.pdf",
                         do.label=F, label.size=5, 
                         alpha=0.9, size=0.5, shape=20, 
                         leg.title=fac.name, no.legend=do.facet.gray, 
                         leg.wid=25, leg.ncol=NULL,
                         leg.override.size=3, leg.pos="top", 
                         leg.dir=ifelse(leg.pos%in%c("left","right"), "vertical","horizontal"),
                         text.rel.size=1.0, panel.title.rel.size=1.2, 
                         x.lab=colnames(embed.xy)[1], y.lab=colnames(embed.xy)[2],
                         leg.labels=levels(fac.v), xlim=NULL, ylim=NULL,
                         base_theme_min=TRUE,
                         common.legend=FALSE, 
                         title="", subtitle="", title.size=14, subtitle.size=9, 
                         base_height=3.75, base_asp=ifelse(no.legend, 1.0, 1.618),
                         rel.height.title=0.1,
                         align="hv", #pmarg=unit(c(1,1,1,1), "pt"), 
                         n.col.max=10, n.row.max=15, 
                         n.col=ifelse(length(facet.levels)<n.col.max, length(facet.levels), 
                                      min(n.col.max, ceiling(sqrt(length(facet.levels))))), 
                         n.row=min(n.row.max, ceiling(length(facet.levels)/n.col))
                         ) {
  if (!is.factor(fac.v)) {
    stop("Expecting a factor")
  }
  if (no.legend) { leg.pos="none"}
  facet.other=F; o=order(fac.v)
  data.plot=data.frame(embed.xy, "gp"=fac.v)[o,]; #print(head(data.plot)) #print(levels(data.plot$gp))
  colnames(data.plot)[c(1,2)]=c("x","y")
  plots.l=NULL
  if (do.facet.gray) {
    if (is.null(colorsv) || (length(colorsv)==0)) {
      colorsv=rep(color.highlight, nlevels(fac.v)); names(colorsv)=levels(fac.v)
    } else if (length(colorsv) != length(facet.levels)) {
      stop(paste0("If colorsv vector is specified, its length must be equal to the length of facet.levels"))
    } 
    plots.l=lapply(facet.levels, function(c) {
      plot.df=data.plot; plot.df$gp=factor(plot.df$gp, levels=c(setdiff(levels(plot.df$gp), c),c))
      plot.df=plot.df[order(plot.df$gp),]
      if (!is.null(colorsv)) { colorsv[names(colorsv)!=c]=color.gray }
      p=ggplot(plot.df, aes(x=x, y=y)) + geom_point(aes(color=gp), alpha=alpha, size=size, shape=shape) + 
        ggtitle(paste0(c))
      if (!is.null(colorsv)) {
        p=p+scale_color_manual(values=colorsv)
      } else {
        p=p+scale_color_discrete(labels=leg.labels)
      }
      return(p) })
  } else {
    p=ggplot(data.plot, aes(x=x, y=y)) + geom_point(aes(color=gp), alpha=alpha, size=size, shape=shape)
    if (do.label) {
      data.plot %>% dplyr::group_by(gp) %>% dplyr::summarize(x = median(x), y = median(y)) -> centers
      p=p+geom_text(data=centers, aes(x=x,y=y,label=gp), size = label.size, color="black")
    }
    if (!is.null(colorsv)) {
      p=p+scale_color_manual(values=colorsv)#
    } else {
      p=p+scale_color_discrete(labels=leg.labels)
    }
    n.col=1; n.row=1
    plots.l=list(p)
  }
  plots.l=lapply(plots.l, function(p) {
    if (base_theme_min) {p = p+theme_light()}
    if (!is.null(xlim)) { p=p+ xlim(xlim)}
    if (!is.null(ylim)) { p=p+ ylim(ylim)}
    p = p + xlab(x.lab) + ylab(y.lab)
    if (leg.pos!="none") {
      p=p+guides(color=guide_legend(title=leg.title, 
                                    ncol=leg.ncol, override.aes = list(alpha = 1, size=leg.override.size)))
    }    
    if (!is.null(leg.wid)) {
      p <-p + theme(legend.key.width=unit(leg.wid, "point"))
    }
    if (!do.save && return.plot.list) {
      p=p+ggtitle(title)
    }
    p = p + theme(aspect.ratio=1.0, legend.position = leg.pos, legend.direction=leg.dir,
                  axis.title=element_text(size=rel(text.rel.size)), 
                  axis.text=element_text(size=rel(text.rel.size)),
                  legend.title=element_text(size=rel(text.rel.size)),
                  legend.text=element_text(size=rel(text.rel.size)), 
                  plot.title=element_text(size=rel(panel.title.rel.size))
    )
    return(p) 
  })
  if (do.save || (!return.plot.list)) {
    p=multiPanelPlot(plots.l,filename=filename, do.save=do.save, 
                     return.plot=!return.plot.list, common.legend=common.legend,
                     title=title, subtitle=subtitle, title.size=title.size, 
                     subtitle.size=subtitle.size,
                     base_height=base_height, base_asp=base_asp,
                     rel.height.title=rel.height.title, 
                     align=align, n.col.max=n.col.max, n.row.max=n.row.max, 
                     n.col=n.col, n.row=n.row)
    if (!return.plot.list) {
      return(p)
    }
  } 
  return(plots.l)
}

# plot.embed.genes(): create plots of expression of genes on 2d embedding
##  Need at minimum genes.plot, data, and embed.xy to be specified.
## ARGUMENTS:
# genes: vector of names of genes, which should correspond to names of rows 
##  in data
# ob: Seurat object (optional; if NULL, then data and embed.xy must be given)
# embed.xy: matrix or data.frame of embedding coordinates (named columns 1 and 2)
# data.type: the assay type in ob, used only if ob is specified and data is not
# data: matrix of gene expression
# ord: an option that specifies the order in which points are plotted, i.e., 
##  if ord=="increasing", the points plotted first have the smallest
##  expression values and the points plotted last (on top) have the biggest
# log.first: Boolean; if TRUE, log2(data+1) plotted rather than data
# mid.frac, mid.val, cap.probs: Passed to get.colors.wcap
# return.plot.list: Boolean; if TRUE, then return list of individual factor
##  level ggplots rather than one that combines those into a single plot
# do.save: Boolean; save the combined plot into file given by filename argument
# filename: name of file, used only if do.save==TRUE, passed to multiPanelPlot
# alpha, size, shape: Passed to geom_point
# leg.title: Passed as title to guide_legend 
# leg.margin: Passed as legend.margin to theme
# leg.pos: Passed as legend.position to theme
# leg.dir: Passed as direction to to guide_colorbar in scale_color_gradientn
# leg.raster: Boolean, passed as raster to guide_colorbar in scale_color_gradientn
# leg.nbin: Passed as nbin to guide_colorbar in scale_color_gradientn
# leg.bar.wid, leg.bar.ht: Passed as relative values for barwidth and barheight 
##  to guide_colorbar in scale_color_gradientn
# leg.ticks.linewidth: Passed as ticks.linewidth to guide_colorbar in scale_color_gradientn
# leg.override.size: Passed as size in override.aes to guide_legend
# text.rel.size: Relative size of axis and legend title and text
# panel.title.rel.size: Relative size of panel title text
# panel.title.plot.face: Face of font of panel titles
# x.lab, y.lab: x axis label, y axis label
# no.axis.text: Boolean; if TRUE, all axis text is blanked out
# xlim, ylim: if not NULL, passed to xlim and ylim respectively
# base_theme_min: Boolean; if TRUE, use theme_light as base theme rather than
##  default theme
# common.legend, title, subtitle, title.size, subtitle.size, base_height, 
##  base_asp, rel.height.title, align, n.col.max, n.row.max, n.col, and n.row: 
##  Passed to multiPanelPlot()
# ...: additional arguments passed to get.colors.wcap()
plot.embed.genes = function(genes, ob=NULL, embed.xy=Embeddings(ob, "umap"), 
                           data.type="SCT",
                           data=GetAssayData(ob,assay=data.type, slot="data"), 
                           ord=c("increasing", "decreasing", "none"),
                           log.first=FALSE, 
                           mid.frac=0.5, mid.val=NULL, cap.probs=c(0.01, 0.99),
                           return.plot.list=TRUE,
                           do.save=TRUE, filename="./genes.multipanel.pdf", 
                           alpha=0.9, size=0.5, shape=16, 
                           leg.title="", leg.margin=margin(1,1,1,1, unit="pt"),
                           leg.pos="right", leg.dir=ifelse(leg.pos%in%c("left","right"), "vertical","horizontal"), 
                           leg.raster=TRUE, leg.nbin=30,
                           leg.bar.wid=ifelse(leg.dir=="vertical", rel(0.3), rel(2.5)), 
                           leg.bar.ht=ifelse(leg.dir=="vertical", rel(2.5), rel(0.3)),
                           leg.ticks.linewidth=0.5,
                           text.rel.size=0.6, panel.title.rel.size=1.2,
                           panel.title.plot.face="bold.italic", 
                           xlab=NULL, ylab=NULL, no.axis.text=FALSE,
                           xlim=NULL, ylim=NULL,
                           base_theme_min=TRUE,
                           common.legend=FALSE, title="", subtitle="",  
                           title.size=14, subtitle.size=9, 
                           base_height=3.5, base_asp=1.1,
                           rel.height.title=0.1, 
                           align="hv", pmarg=unit(c(1,1,1,1), "pt"), 
                           n.col.max=10, n.row.max=15, 
                           ...) {
  ord=match.arg(ord )
  genes.plot=genes[genes %in% rownames(data)]
  if (length(genes)!=length(genes.plot)) {
    warning(paste0("Some genes were not found in row names: ", 
                   paste0(setdiff(genes, genes.plot), collapse=", ")))
  }
  if (length(genes.plot)==0) {
    stop("No genes to plot")
  }
  data=data[genes.plot,,drop=FALSE]
  if (log.first) {
    data=log2(data+1)
  }
  if (ncol(data)!=nrow(embed.xy)) { 
    stop("Expecting number of columns of data to correspond to number of rows of embed.xy")}
  plots.l=lapply(genes.plot, function(g) { #print(g);
    mfr=mid.frac; mva=mid.val
    if ((!is.null(mid.frac)) && (length(mid.frac)>1)) { mfr=mid.frac[g] }
    if ((!is.null(mid.val)) && (length(mid.val)>1)) { mva=mid.val[g] }
    colrs.cap=get.colors.wcap(data[g,,drop=F], mid.frac=mfr,mid.val=mva, 
                              cap.probs=cap.probs,...)
    if (is.null(colrs.cap)) {warning("Cannot create plot for ", g, "; skipping"); return(NULL)}
    colrs=colrs.cap[["colrs"]]
    plot.mat=colrs.cap[["dat"]]; 
    dat=data.frame(cbind(embed.xy, t(plot.mat)))
    if (ord %in% c("increasing", "decreasing")) {
      # print(paste0("Putting values in ", ord, " order"))
      o=order(plot.mat[1,], decreasing=(ord=="decreasing")); 
      dat = dat[o,] 
    }     
    plot.df=melt(dat, "id.vars"=colnames(embed.xy))
    p=ggplot(plot.df, aes_string(x=colnames(embed.xy)[1], y=colnames(embed.xy)[2])) + 
      geom_point(aes(color=value), alpha=alpha, size=size, shape=shape)
    if (base_theme_min) { p=p+theme_light() }
    p= p + ggtitle(g) +
      scale_color_gradientn(colours=colrs, 
                            guide=guide_colorbar(title=leg.title, override.aes=list(alpha=1), direction=leg.dir, ticks.linewidth=leg.ticks.linewidth,
                                                 barwidth=rel(leg.bar.wid), barheight=rel(leg.bar.ht), raster=leg.raster, nbin=leg.nbin)) 
    p=p+theme(legend.position=leg.pos, aspect.ratio=1, axis.text=element_text(size=rel(text.rel.size)),
              axis.title=element_text(size=rel(text.rel.size)), legend.title=element_text(size=rel(text.rel.size)),
              legend.text=element_text(size=rel(text.rel.size)), 
              plot.title=element_text(size=rel(panel.title.rel.size), face=panel.title.plot.face),
              plot.margin=pmarg, legend.margin=leg.margin)
    if (!is.null(xlim)) { p=p+ xlim(xlim)}
    if (!is.null(ylim)) { p=p+ ylim(ylim)}
    if (!is.null(xlab)) { p=p+ xlab(xlab)}
    if (!is.null(ylab)) { p=p+ ylab(ylab)}
    if (no.axis.text) { p=p+theme(axis.text=element_blank()) }  
    return(p)
  })
  
  if (do.save || !return.plot.list) { 
    p=multiPanelPlot(plots.l,filename=filename, do.save=do.save, return.plot=!return.plot.list, common.legend=common.legend,
                     title=title, subtitle=subtitle, title.size=title.size, subtitle.size=subtitle.size,
                     base_height=base_height, base_asp=base_asp, rel.height.title=rel.height.title, align=align, 
                     n.col.max=n.col.max, n.row.max=n.row.max)
    if (!return.plot.list) {
      return(p)
    }
  }
  if (return.plot.list) {
    return(plots.l)
  } 
}

# plot.embed.features(): wrapper function for plot.embed.genes() that makes it
##  easy to plot meta.data features with the same functionality
plot.embed.features = function(features, ob=NULL, 
                            data=t(ob[[]]), # get metadata
                            filename="./features.multipanel.pdf", 
                            ...) {
  plot.embed.genes(genes=features, ob=ob, data=data, filename=filename, ...)
}

# get.colors.wcap(): Utility function to get a palette that is centered 
##  appropriately and that reflects potentially capped values in the data.
## ARGUMENTS:
# data: vector or matrix of values to be plotted
# mid.val: numeric value in data that should represent the center of the palette;
##  optional; if NULL, then ctr or mid.frac determine mid.val
# ctr: Boolean; if TRUE, then mid.val, if NULL, is set to 0 
# mid.frac: a fraction of the range of data (after capping) used to determine 
##  mid.val, only if mid.val is NULL and ctr==FALSE.
# cap: numeric vector of length 2, minimum and maximum thresholds below and above 
##  which values in data should be capped.
# cap.probs: quantiles used to determine cap, only if cap is NULL
# add.full: Boolean; if TRUE, use capped values to determine the color range and
##  palette, but then add back in the full range so it shows in the legend;
##  values above cap[2] will all have the same color, and values below cap[1] will 
##  also have the same color.
# eps.div: Jitter parameter; in case the range is degenerate, i.e., consists of  
##  just one value, that value is divided by eps.div, and the range is expanded 
##  equally in both directions by the result 
# na.col: Color used for NA values
get.colors.wcap <- function(data, mid.val=NULL, ctr=F, mid.frac=0.5, 
                            cap=NULL, cap.probs=c(0, 1), add.full=FALSE, eps.div=10, na.col="gray50",
                            col.fn=c("get.viridis.col"),...) {
  col.fn=get(match.arg(col.fn))
  full.range=range(data, na.rm=TRUE)
  orig.data=data # save orig input
  num.v=as.vector(data)
  if (is.null(cap)) {
    cap=signif(quantile(num.v, cap.probs, na.rm=T),2)
  }
  if (is.infinite(cap[1])) {
    if (cap[1]>0) {stop("Lower cap cannot be Inf")
    } else { cap[1]=min(num.v[!is.infinite(num.v)]); warning(print(paste0("Set lower cap to min non-infinite value: ", cap[1]))) }
  }
  if (is.infinite(cap[2])) {
    if (cap[2]<0) {stop("Upper cap cannot be -Inf")
    } else { cap[2]=max(num.v[!is.infinite(num.v)]); warning(print(paste0("Set upper cap to max non-infinite value: ", cap[2]))) }
  }
  if (cap[2]==cap[1]) {
    warning(print(paste0("Caps are equal: ", cap[1], "; using full range instead!: ", 
                         paste0(signif(range(num.v),3), collapse=", "))))
    cap=range(num.v); 
    eps=abs(cap[2])/(2*eps.div); 
    if (cap[1]==cap[2]) {
      warning(print(paste0("Caps are still equal; expanding range by ", signif(eps)))) #return(NULL) 
      new.cap=cap[2]+eps
      if (cap[2]<0) {cap[2]=min(new.cap,0)}
      new.cap=cap[1]-eps
      if (cap[1]>0) {cap[1]=max(new.cap,0)}
    }
  }
  data=squish(data, cap, only.finite=F)
  range.val=range(data, na.rm=TRUE)
  if (range.val[2]==range.val[1]) {range.val=cap}
  if (is.null(mid.val)) {
    if (ctr) { mid.val=0 
    } else {
      mid.val=signif((mid.frac*(range.val[2]-range.val[1]) + range.val[1]),2)
    }
  } 

  if (is.na(mid.val)) { 
    colrs=na.col; names(colrs)="NA"
  } else {
    colrs=do.call(col.fn, c(list(range.val=range.val, mid.val=mid.val),list(...))); 
  }
  if (!is.na(mid.val) && add.full && ((full.range[1]!=range.val[1]) || (full.range[2]!=range.val[2]))) {
    frac.full=(range.val[2]-range.val[1])/(full.range[2]-full.range[1]); 
    print(paste0("limited fraction of full range: ", signif( frac.full, 4)))
    if (!is.null(names(colrs))) {
      names(colrs)=as.numeric(names(colrs))*frac.full;
    }
    frac.init=((range.val[1]-full.range[1])/(full.range[2]-full.range[1]))
    step.unit=frac.full/length(colrs) #; print(paste0("step.unit: ", signif(step.unit, 4)))
    rep.init=round(frac.init/step.unit)
    rep.end=round((1-frac.init-frac.full)/step.unit)
    if ((rep.init>0) && (full.range[1]!=range.val[1])) {
      print("Including minimum of full range of values in names of colorscale")
      colrs=c(rep(colrs[1], rep.init), colrs)
      if (!is.null(names(colrs))) {
        names(colrs)=as.numeric(names(colrs))+frac.init
        names(colrs)[1:rep.init]=seq(0,by=step.unit, length.out=rep.init)
      }
    }
    if ((rep.end>0) && (full.range[2]!=range.val[2])) {
      print("Including maximum of full range of values in names of colorscale")
      colrs=c(colrs, rep(colrs[[length(colrs)]], rep.end))
      if (!is.null(names(colrs))) {
        names(colrs)[(length(colrs)-rep.end+1):(length(colrs))]=
          c(frac.init+frac.full+seq(step.unit, by=step.unit, length.out=(rep.end-1)), 1.0)
      }
    }
  }
  ret.l=list("colrs"=colrs, "dat"=data, "mid"=mid.val)
  if (add.full) { ret.l[["dat"]]=orig.data }
  return(ret.l)
}

# get.viridis.col(): Utility function for getting a palette based on the 
## viridis color package, where the right number of steps is selected to 
## place the center color at the specified mid value.
## ARGUMENTS:
# range.val: length 2 numeric vector, range of values to be covered
# mid.val: value in the range where the palette should be centered
# n.steps.final: number of steps in the final palette
# pal.fn: function for creating the palette, defaults to viridis_pal()  
# direction: 1 or -1, specifies color order; passed to pal.fn
# begin: hue in [0,1] where color map should begin; passed to pal.fn
# end: hue in [0,1] where color map should end; passed to pal.fn
# pal.opt: color map to use, passed as option to pal.fn
get.viridis.col <- function(range.val=NULL, mid.val=NULL, n.steps.final=11, 
                            pal.fn=viridis_pal,
                            direction=1, begin=0, end=1, pal.mid=0.5, pal.opt=
                              c("viridis", "magma", "inferno", "plasma", "cividis")) {
  frac.names=TRUE # needed for plotly, doesn't seem to hurt
  pal.opt=match.arg(pal.opt)
  if (is.null(range.val) && !is.null(mid.val)) {
    stop("In get.viridis.col(): If mid.val is specified, range.val should also be specified")
  }
  frac.low=0.5; frac.high=0.5
  n.steps.final=min(n.steps.final, max(n.steps.final, 3), 255)
  if (!is.null(range.val)) {
    if (is.null(mid.val)) {
      mid.val=range.val[1]+(range.val[2]-range.val[1])/2.0
    }
    ## fraction of actual values in the low vs. high ranges
    frac.low=(mid.val-range.val[1])/(range.val[2]-range.val[1]);  
    frac.high=1.0-frac.low;
  }
  pal.mid=(end-begin)*pal.mid; 
  n.steps.low=min(max(2, ceiling(n.steps.final*frac.low)), n.steps.final-1)
  n.steps.high=n.steps.final-n.steps.low+1
  cols.low=pal.fn(begin=begin, end=pal.mid, direction=direction, option=pal.opt)(n.steps.low)
  cols.high=pal.fn(begin=pal.mid, end=end, direction=direction, option=pal.opt)(n.steps.high) 
  if (direction==1) {
    cols=c(cols.low, cols.high[2:length(cols.high)])
  } else {
    cols=c(cols.high, cols.low[2:length(cols.low)])
  }
  if (!is.null(range.val)) {
    if (frac.names) { 
      names(cols)=c( 
        seq(from=0, to=frac.low, length.out=n.steps.low), # plotly needs normalized values as names
        seq(from=frac.low, to=1.0, length.out=n.steps.high)[2:n.steps.high]) # plotly needs normalized values as names
    } 
  }
  #print(cols)
  return(cols)
}


# multiPanelPlot(): a utility function for combining ggplots into a single plot
## ARGUMENTS:
# plots.l: list of ggplot plots
# filename: where plot is saved, if do.save==TRUE
# do.save: Boolean, if TRUE, save plot to file named by filename argument
# return.plot: Boolean; if TRUE, return the plot that is created; else return nothing
# common.legend: Boolean; if TRUE, create a common legend for all the plots
# title: title for combined plot
# subtitle: subtitle for combined plot
# title.size: font size for title
# subtitle.size: font size for subtitle
# base_height: height for individual panels, passed as base_height to save_plot
# base_asp: the aspect ratio (width/height) of one subplot; 
##  default of 1.618 (the golden ratio) works well for figures with a legend 
# rel.height.title: relative height for title space in combined plot
# align: Passed as align argument to ggarrange
# n.col.max: Maximum number of columns desired in 1-page grid 
# n.row.max: Maximum number of rows desired in 1-page grid;
##  if length(plots.l) > n.col.max*n.row.max, then a pdf of multiple pages 
##    will be created, as many as required to include all plots in plots.l
multiPanelPlot <- function(plots.l, filename="./multipanel.pdf", do.save=TRUE, return.plot=FALSE,
                           common.legend=FALSE, 
                           title="", subtitle="", title.size=11, subtitle.size=7, 
                           base_height=3.71, base_asp=1.618,
                           rel.height.title=0.1,
                           align="hv", #pmarg=unit(c(1,1,1,1), "pt"), 
                           n.col.max=10, n.row.max=10, 
                           n.col=ifelse(length(plots.l)<n.col.max, length(plots.l), 
                                        min(n.col.max, ceiling(sqrt(length(plots.l))))), 
                           n.row=min(n.row.max, ceiling(length(plots.l)/n.col))) {
  n.pp=n.row*n.col;#   print(n.pp)
  #print(length(plots.l))
  multip.flag=FALSE
  if (length(plots.l)>n.pp) { 
    multip.flag=TRUE
  }
  ggarrange.args=list(ncol=n.col, common.legend=common.legend, align=align)
  title <- ggdraw() + draw_label(title, fontface='bold', size=title.size)
  plots.pp=do.call(ggarrange, c(ggarrange.args, list(plotlist=plots.l, nrow=n.row)))
  if (!multip.flag) { plots.pp=list(plots.pp) }
  plots.pp=lapply(plots.pp, function(p) { add_sub(p, label=subtitle, size=subtitle.size, fontface="bold") })
  plots.pp=lapply(plots.pp, function(p) {plot_grid(title, p, ncol=1, rel_heights=c(rel.height.title, 1))})
  p=do.call(marrangeGrob, c(list(plots.pp), list(nrow = 1, ncol = 1), top=""))
  if (do.save) {
    print(paste0("Saving plot to file ", filename))
    if(multip.flag) { print(paste0("File of multiple pages required"))}
    save_plot(filename=filename, plot=p, base_height=base_height, base_asp=base_asp,
              ncol=n.col, nrow=n.row)
  } 
  if (return.plot) { return(p) 
  } else {return()}
}

# summarySE(): Utility function for computing statistics on the data by groups
## This function needs to be updated as it still uses plyr functions instead of dplyr.
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=TRUE,
                      conf.interval=.95, .drop=TRUE) {
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  colnames(datac)[which(colnames(datac)=="mean")]=measurevar
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  return(datac)
}

# high.loadings(): Utility function for getting the variables associated with each PC.
high.loadings <- function(pca.loadings, n.vars=min(10, nrow(pca.loadings)),
                          n.pcs=min(10, ncol(pca.loadings)), 
                          which.pcs=1:n.pcs, var.thr=1.0, max.n.vars=30, do.print=TRUE,
                          balanced=TRUE) {
  use.pcs=colnames(pca.loadings)[which.pcs]
  
  genes_l=lapply(use.pcs, function(pc_x) {
    if (!is.null(n.vars)) {
      n=n.vars
    } else {
      y=cumsum(sort(pca.loadings[,pc_x]^2, decreasing=TRUE))
      n=min(match(FALSE, y<var.thr)-1, max.n.vars, na.rm=TRUE )
    }
    if (is.na(n) || (n==0)) { 
      genes=NULL 
    } else {
      if (balanced && !is.null(n.vars)) {
        n.top=ceiling(n/2); n.bot=floor(n/2)
        ord.pc=order(pca.loadings[,pc_x], decreasing=TRUE)
        genes=unique(c(rownames(pca.loadings)[ord.pc][1:n.top],
                       rownames(pca.loadings)[ord.pc][(length(ord.pc)-n.bot+1):length(ord.pc)]))
      } else {
        ord.pc=order(abs(pca.loadings[,pc_x]), decreasing=TRUE)
        genes=rownames(pca.loadings)[ord.pc][1:n]
      }
    }
    if (do.print) { 
      print(pc_x)
      if (is.null(n.vars)) { print(paste0("Number of vars: ", n)) }
      print(genes)
      cat ("\n")
    }
    genes
  })
  names(genes_l)=use.pcs
  return(genes_l)
}

