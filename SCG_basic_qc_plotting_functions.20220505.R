## S.J. Riesenfeld
# Code for visualizing basic QC of single-cell genomics data
library("ggplot2")
library("dplyr")

# plot.histn(): Create a histogram (or density) plot of a data feature.
## Arguments:
## n.per.cell: numeric data vector (e.g., total number of nUMI per cell);
## cutoffs: optional, a numeric vector of length 2 (min/max thresholds to visualize);
## do.log: Boolean, if TRUE, log2(n.per.cell+1) (and log2(cutoffs+1), if provided) 
##   will be plotted;
## sample.name: optional string, used for plot title;
## label: optional string to specify a name of the data in n.per.cell;
## bins: number of bins for the histogram;
## label.cutoffs: Boolean, if TRUE and if !is.null(cutoffs), then numeric text
##   annotations are provided on the plot for the cutoffs;
## do.density: Boolean, if TRUE, make a density plot instead of a histogram;
## return.plot: Boolean, if TRUE, return the ggplot object rather than printing it;
## highlight.color: Color used to highlight and annotate the cutoffs.
plot.histn <- function(n.per.cell, cutoffs=NULL, do.log=T, sample.name="", label="", bins=100, label.cutoffs=T,
                       do.density=F, return.plot=FALSE, highlight.color="orangered1") {
  if (do.log) {
    if (!is.null(cutoffs)) { cutoffs=log2(cutoffs+1) }
    n.per.cell=log2(n.per.cell+1)
    xlab=paste0("log2 #", label, " per cell")
    tit=paste0(sample.name, ": Hist. #", label, " Per Cell (Log Scale)")
  } else {
    tit=paste0(sample.name, ": Hist. #", label, " Per Cell")
    xlab=paste0("#", label, " per cell")
  }
  pct.label=""
  if (!is.null(cutoffs)) { 
    cutoffs.plot=signif(cutoffs, 3)
    n.pass=sum((n.per.cell>cutoffs[1]) & (n.per.cell<cutoffs[2]))
    pct.label=paste0(round(n.pass/length(n.per.cell)*100, 1),
                   "% (", n.pass, " of ", length(n.per.cell) , ") cells meet #", label, " cutoffs (", cutoffs.plot[1], ", ", cutoffs.plot[2], ")")
    print(pct.label)
  }
  plot.df=data.frame("n.per.cell"=n.per.cell); 
  p=ggplot(plot.df, aes(n.per.cell)) + xlab(xlab) + 
    ggtitle(tit, subtitle=pct.label) + theme_light()
  if (do.density) {
    p=p+geom_density()
  } else {
    p=p+ geom_histogram(bins=bins)
  }
  if (!is.null(cutoffs)) {
    p=p+ geom_vline(xintercept=cutoffs.plot, col=highlight.color, linetype="dashed", alpha=0.8) 
    if (label.cutoffs) {
      p=p+geom_text(data=data.frame("cutoffs"=cutoffs.plot), aes(x=cutoffs.plot, label=cutoffs.plot), y=-1, 
                    col=highlight.color, alpha=0.8)
    }
  }
  if (return.plot) {
    return(p)
  } else {
    print(p)
  }
}

# plot.ft1.v.ft2(): Plot two data/feature vectors, ft1 and ft2, against each other; 
# compute which cells are outliers with respect to each marginal distribution and
# also with respect to a loess fit of the joint distribution for the non-outlier cells.
## Arguments:
## ft1: numeric data vector (e.g., total number of UMIs per cell);
## ft2: numeric data vector (e.g., total number of genes detected per cell);
## ft1.name: optional, string to use as the name of the data in ft1 in the plot;
## ft2.name: analogous to ft1.name ;
## ft1.cutoffs: optional, a numeric vector of length 2 giving min/max thresholds 
##    for the data in ft1;
## ft2.cutoffs: analogous to ft1.cutoffs
## do.log: Boolean, if TRUE, log2(ft1+1) and log2(ft2+1) will be plotted, as well 
##    as log2(ft1.cutoffs+1) and log2(ft2.cutoffs+1), if provided; 
## sample.name: optional, string used in plot title;
## return.plot: Boolean, if TRUE, changes the returned value. See "Returns" below.
## do.loess: Boolean, if TRUE, then the non-outlier cells wrt cutoffs are used to 
##    make a loess approximation of ft2~ft1;
## loess.resid.probs: a numeric vector of length two, min and max quantiles for 
##    loess residuals; cells below the min or above the max are outliers;
## pt.size: point size for the scatter plot;
## pt.alpha: alpha value for the points in the scatter plot;
## linewdith: line size for the threshold demarcations;
## pass.color: color of points that are within all QC thresholds
## highlight.color: color of points that will be excluded by at least one QC measure
## contour.alpha: alpha value for the 2d density contour lines;
## contour.col: color for the 2d density contour lines;
## contour.linewidth: line size for the 2d density contours;
## ... arguments passed on to stat_density_2d()
# Returns: If return.plot is FALSE, returns a vector of indices of 
### cells that meet all the QC measures. If return.plot is TRUE, returns a 
### list with an entry named "cells2keep", which is the vector of indices of 
### cells that meet all the QC measures, and an entry named "p", which is the 
### created ggplot object.
plot.ft1.v.ft2 <- function(ft1, ft2, ft.names=c("ft1", "ft2"),
                           ft.cutoffs=list(NULL, NULL),
                           do.log=T, sample.name="", return.plot=F, 
                           do.loess=T, loess.resid.probs=c(0.005, 0.995),
                           pt.size=0.5, pt.alpha=0.5, linewidth=0.5, 
                           contour.alpha=0.75, contour.col="deepskyblue", 
                           contour.linewidth=0.25,
                           pass.color="gray30", highlight.color="orangered1",
                           ...) {
  performing.qc=FALSE
  if ((FALSE %in% sapply(ft.cutoffs, is.null)) 
      || (loess.resid.probs[1]>0) || (loess.resid.probs[2]<1)) {
    performing.qc=TRUE
  }
  legend.title="QC"
  cutoff.line.color=highlight.color
  pf.colors=c("Pass"=pass.color, "Fail"=highlight.color)
  log.tag=""
  if (do.log) {
    ft1=log2(ft1+1); ft2=log2(ft2+1)
    ft.cutoffs=lapply(ft.cutoffs, function(x) {if (!is.null(x)) {log2(x+1)} else {x}})
    log.tag="log2 "
  }
  #print(summary(ft1)); print(summary(ft2))
  #print(ft.cutoffs)
  ft.pass.l=lapply(1:2, function(i) { if (!is.null(ft.cutoffs[[i]])) {
    ft=list(ft1,ft2)[[i]]; thr.v=ft.cutoffs[[i]]
    ((ft>thr.v[1]) & (ft<thr.v[2])) } else {NULL} })
 
  cells2keep= (ft.pass.l[[1]] & ft.pass.l[[2]]) # a Boolean vector
  if (sum(cells2keep)==0){
    warning("All cells failed QC")
  } else if (do.loess) {
    ft2.loess=loess(ft2[cells2keep]~ft1[cells2keep]) # fit only for the QC'd cells
    ft2.resid.limits=quantile(ft2.loess$residuals, prob=loess.resid.probs)
    loess.pass=((ft2.loess$resid>ft2.resid.limits[1]) & (ft2.loess$resid<ft2.resid.limits[2])) 
    cells2keep[cells2keep] = (cells2keep[cells2keep] & loess.pass)
  }
  n.pass=sum(cells2keep)
  pct.label=NULL
  if (performing.qc) {
    pct.label=paste0(round(n.pass/length(cells2keep)*100,1), "% (",
    n.pass, " of ", length(cells2keep), ") cells jointly meet input QC cutoffs")
    print(pct.label)
  }
  qc.pass=factor(cells2keep, levels=c("FALSE", "TRUE"), labels=c("Fail", "Pass"))
  qc.df=data.frame(qc.pass, ft1, ft2); colnames(qc.df)=c("QC.pass", ft.names)
  p=ggplot(qc.df, aes_string(x=ft.names[1], y=ft.names[2])) 
  if (!performing.qc) {
    p = p+ geom_point(alpha=pt.alpha, shape=16, size=pt.size) + 
      ggtitle(sample.name) + theme_light() +
      xlab(paste0(log.tag, ft.names[1], " per cell")) + ylab(paste0(log.tag, ft.names[2], " per cell"))
  } else {
    p=p + geom_point(aes(color=QC.pass), alpha=pt.alpha, shape=16, size=pt.size) +
      scale_color_manual(values=pf.colors, guide=guide_legend(title=legend.title, override.aes=list(alpha=1, size=1))) +
      geom_vline(xintercept=ft.cutoffs[[1]], color=cutoff.line.color, linetype="dashed", size=linewidth) +
      geom_hline(yintercept=ft.cutoffs[[2]], color=cutoff.line.color, linetype="dashed", size=linewidth) + 
      ggtitle(paste0(sample.name), subtitle=pct.label) 
    if (!is.null(ft.cutoffs[[1]])) {
      p=p+ annotate("text", x=signif(ft.cutoffs[[1]],3), label=signif(ft.cutoffs[[1]],3),
                    y=min(ft2)-0.2, col=pf.colors["Fail"])
    }
    if (!is.null(ft.cutoffs[[2]])) {
      p=p+annotate("text", y=signif(ft.cutoffs[[2]],3), label=signif(ft.cutoffs[[2]],3),
                   x=min(ft1)-0.3, col=pf.colors["Fail"])
    }
  }
  p = p + stat_density_2d(color=contour.col, alpha=contour.alpha, 
                          size=contour.linewidth, ...)
  ret.v=which(cells2keep) # indices of cells passing QC
  if (return.plot) {
    return(list("cells2keep"=ret.v, "p"=p))
  } else {
    print(p)
    return(ret.v)
  }
}
