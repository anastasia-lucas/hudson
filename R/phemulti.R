#' phemulti
#'
#' Create mulitple Manhattan plots for PheWAS
#' Dependencies: ggplot2, gridExtra, grid, RColorBrewer
#' Suggested: ggrepel
#' @param phemulti list of data frames, if not plato or plink format, must contain PHE, SNP, CHR, POS, pvalue, columns, optional Shape
#' @param phegroup optional grouping file for phenotypes, must contain PHE and Group columns
#' @param tline optional pvalue threshold to draw red line at in top plot
#' @param bline optional pvalue threshold to draw red line at in bottom plot
#' @param log10 plot -log10() of pvalue column, boolean
#' @param yaxis label for y-axis in the format c("top", "bottom"), automatically set if log10=TRUE
#' @param opacity opacity of points, from 0-1, useful for dense plots
#' @param annotate_snp vector of RSIDs to annotate
#' @param annotate_p list of pvalue thresholds to annotate in the order of c(p_top, p_bottom)
#' @param plot_titles optional vector of strings for plot titles
#' @param chrcolor1 first alternating color for chromosome
#' @param chrcolor2 second alternating color for chromosome
#' @param highlight_snp vector of snps to highlight
#' @param highlight_p list of pvalue thresholds to annotate in the order of c(p_top, p_bottom)
#' @param highlighter color to highlight
#' @param groupcolors named vector of colors where names correspond to data in 'PHE' or 'Group' column
#' @param freey allow y-axes to scale with the data
#' @param background variegated or white
#' @param chrblocks boolean, turns on x-axis chromosome marker blocks
#' @param file file name of saved image
#' @param hgt height of plot in inches
#' @param hgtratio height ratio of plots, equal to top plot proportion
#' @param wi width of plot in inches
#' @param res resolution of plot in pixels per inch
#' @param ymax y-axis max value, if log10 is TRUE then log10 value must be specified
#' @param ymin y-axis min value, if log10 is TRUE then log10 value must be specified
#' @param legend_pos location of legend, values allowed  = ("bottom", "right","top","left")
#' @return gtable (grob) object. Draw with \code{\link[grid]{grid.draw}}.
#' @import ggplot2
#' @importFrom gridExtra grid arrangeGrob grid.arrange
#' @export
#' @examples
#' data(phewas.t)
#' data(phewas.b)
#' phemulti(list(phewas.t, phewas.b), 
#' plot_titles = c("PheWAS Example: Data 1", "PheWAS Example: Data 2"))
phemulti <- function(datalist, phegroup, tline, bline, log10=TRUE, yaxis, opacity=1, annotate_snp, annotate_p, highlight_snp, highlight_p, highlighter="red", plot_titles=NULL, chrcolor1="#AAAAAA", chrcolor2="#4D4D4D", groupcolors, freey=FALSE, background="variegated", chrblocks=TRUE, file=NULL, hgtratio=0.5, hgt=7, wi=12, res=300, ymax = NULL, ymin = NULL, legend_pos = "bottom") {
  
  #Get all chromosomes and phenotypes from all datasets
  chrs <- c()
  phes <- c()
  for (data in datalist) {
    chrs <- unique(c(chrs, unique(data$CHR)))
    phes <- unique(c(phes, unique(data$PHE)))
  }
  phe_dummy_data <- data.frame(phes)
  colnames(phe_dummy_data) = c("PHE")
  phe_dummy_data$SNP <- 1:nrow(phe_dummy_data)
  phe_dummy_data$CHR <- 1
  phe_dummy_data$POS <- 1
  phe_dummy_data$pvalue <- 1
  
  #Set up colors
  nchrcolors <- nlevels(factor(chrs))
  if(!missing(groupcolors)){
    dcols <- c(rep(x=c(chrcolor1, chrcolor2), length.out=nchrcolors, each=1), "#FFFFFF", "#EBEBEB")
    names(dcols) <-c(levels(factor(chrs)), "shade_ffffff", "shade_ebebeb")
    topcols <- c(dcols, groupcolors)
  } else {
    #Top Colors
    ngroupcolors <- nlevels(factor(phes))
    if(ngroupcolors>15){
      if (!requireNamespace(c("RColorBrewer"), quietly = TRUE)==TRUE) {
        stop("Please install RColorBrewer to add color attributes for more than 15 colors.", call. = FALSE)
      } else {
        getPalette = grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))
        topcols <- c(rep(x=c(chrcolor1, chrcolor2), length.out=nchrcolors, each=1), getPalette(ngroupcolors), "#FFFFFF", "#EBEBEB")
      }
    } else {
      pal <- pal <- c("#009292", "#920000", "#490092", "#db6d00", "#24ff24", 
                      "#ffff6d", "#000000", "#006ddb", "#004949","#924900", 
                      "#ff6db6", "#6db6ff","#b66dff", "#ffb6db","#b6dbff")
      topcols <- c(rep(x=c(chrcolor1, chrcolor2), length.out=nchrcolors, each=1), pal[1:ngroupcolors], "#FFFFFF", "#EBEBEB")
    }
    names(topcols) <-c(levels(factor(chrs)), levels(factor(phes)), "shade_ffffff", "shade_ebebeb")
  }
  
  p <- hudson_top(phe_dummy_data, phegroup, toptitle = NULL, colors = topcols)
  legend <- g_legend(p + theme(legend.position = legend_pos)) 
  
  plots = vector(mode = "list", length = length(datalist))
  for (i in 1:length(datalist)) {
    p = hudson_top(datalist[[i]], phegroup, tline, bline, log10, yaxis, opacity, annotate_snp, annotate_p, highlight_snp, highlight_p, highlighter, toptitle = plot_titles[i], chrcolor1, chrcolor2, groupcolors, freey, background, chrblocks, ymax, ymin, colors = topcols)
    plots[[i]] = p
  }
  p = grid_arrange_plots_and_legend(plots, position = legend_pos, legend = legend)
  
  if (!is.null(file)) {
    print(paste("Saving plot to ", file, ".png", sep=""))
    ggsave(p, filename=paste(file, ".png", sep=""), dpi=res, units="in", height=hgt, width=wi)
  }
  return(p)
}

g_legend<-function(a.gplot){
  if (!gtable::is.gtable(a.gplot))
    a.gplot <- ggplotGrob(a.gplot)
  #gtable_filter(a.gplot, 'guide-box', fixed=TRUE)
  leg <- which(sapply(a.gplot$grobs, function(x) x$name) == "guide-box")
  a.gplot$grobs[[leg]]
}


grid_arrange_plots_and_legend <- function(plots,
                                       legend,
                                       ncol = 1,
                                       nrow = length(plots),
                                       position = c("bottom", "right","top","left"),
                                       plot=TRUE)
{
  
  position <- match.arg(position)
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) {
    if (is.ggplot(x)) { x + theme(legend.position="none") } else { x }})
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "top" = arrangeGrob(
                       legend,
                       do.call(arrangeGrob, gl),
                       ncol = 1,
                       heights = unit.c(lheight, unit(1, "npc") - lheight)),
                     "bottom" = arrangeGrob(
                       do.call(arrangeGrob, gl),
                       legend,
                       ncol = 1,
                       heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "left" = arrangeGrob(
                       legend,
                       do.call(arrangeGrob, gl),
                       ncol = 2,
                       widths = unit.c(lwidth, unit(1, "npc") - lwidth)),
                     "right" = arrangeGrob(
                       do.call(arrangeGrob, gl),
                       legend,
                       ncol = 2,
                       widths = unit.c(unit(1, "npc") - lwidth, lwidth))
  )
  
  if (plot) {
    grid.newpage()
    grid.draw(combined)
  }
  
  # return gtable invisibly
  invisible(combined)
}


hudson_top <- function(top, phegroup, tline, bline, log10=TRUE, yaxis, opacity=1, annotate_snp, annotate_p, highlight_snp, highlight_p, highlighter="red", toptitle=NULL, chrcolor1="#AAAAAA", chrcolor2="#4D4D4D", groupcolors, freey=FALSE, background="variegated", chrblocks=TRUE, ymax = NULL, ymin = NULL, colors = NULL){
  
  #Sort data
  topn <- names(top)
  top$Location <- "Top"
  if(!("Shape" %in% colnames((top)))){top$Shape <- NA}
  d <- top
  d$POS <- as.numeric(as.character(d$POS))
  d$CHR <- factor(d$CHR, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"))
  if(!missing(phegroup)){
    print("Only phenotypes with grouping information will be plotted")
    d_phe <- merge(phegroup, d, by="PHE")
    names(d_phe)[names(d_phe)=="Group"] <- "Color"
  } else {
    d_phe <- d
    names(d_phe)[names(d_phe)=="PHE"] <- "Color"
  }
  d_order <- d_phe[order(d_phe$CHR, d_phe$POS), ]
  d_order$pos_index <- seq.int(nrow(d_order))
  d_order_sub <- d_order[, c("SNP", "CHR", "POS", "pvalue", "pos_index")]
  
  #Set up dataframe with color and position info
  maxRows <- by(d_order_sub, d_order_sub$CHR, function(x) x[which.max(x$pos_index),])
  minRows <- by(d_order_sub, d_order_sub$CHR, function(x) x[which.min(x$pos_index),])
  milimits <- do.call(rbind, minRows)
  malimits <- do.call(rbind, maxRows)
  lims <- merge(milimits, malimits, by="CHR")
  names(lims) <- c("Color", "snpx", "px", "posx", "posmin", "snpy", "py", "posy", "posmax")
  lims$av <- (lims$posmin + lims$posmax)/2
  lims <- lims[order(lims$Color),]
  lims$shademap <- rep(c("shade_ffffff", "shade_ebebeb"), length.out=nrow(lims), each=1)
  
  #Set up colors
  topcols <- colors
  
  #Info for y-axis
  if(log10==TRUE){
    d_order$pval <- -log10(d_order$pvalue)
    yaxislab1 <- expression(paste("-log"[10], "(p-value)", sep=""))
    yaxislab2 <- expression(paste("-log"[10], "(p-value)", sep=""))
    if(!missing(tline)) {tredline <- -log10(tline)}
    if(!missing(bline)) {bredline <- -log10(bline)}
  } else {
    d_order$pval <- d_order$pvalue
    yaxislab1 <- yaxis[1]
    yaxislab2 <- yaxis[2]
    if(!missing(tline)) {tredline <- tline}
    if(!missing(bline)) {bredline <- bline}
  }
  
  if (is.null(ymax)) {
    yaxismax1 <- ifelse(freey==FALSE, max(d_order$pval[which(d_order$pval< Inf)]), max(d_order$pval[which(d_order$pval< Inf) & d_order$Location=="Top"]))
  } else {
    yaxismax1 <- ymax
  }
  
  if (is.null(ymin)) {
    yaxismin1 <- ifelse(freey==FALSE, 0, min(d_order$pval[d_order$Location=="Top"]))
  } else {
    yaxismin1 <- ymin
  }
  
  #Theme options
  backpanel1 <- ifelse(background=="white", "NULL", "geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = yaxismin1, ymax = Inf, fill=factor(shademap)), alpha = 0.5)" )
  
  #Start plotting
  #TOP PLOT
  p1 <- ggplot() + eval(parse(text=backpanel1))
  #Add shape info if available
  if("Shape" %in% topn){
    p1 <- p1 + geom_point(data=d_order[d_order$Location=="Top",], aes(x=pos_index, y=pval, color=factor(Color), shape=factor(Shape)), alpha=opacity)
  } else {
    p1 <- p1 + geom_point(data=d_order[d_order$Location=="Top",], aes(x=pos_index, y=pval, color=factor(Color)), alpha=opacity)
  }
  p1 <- p1 + scale_x_continuous(breaks=lims$av, labels=lims$Color, expand=c(0,0))
  if(chrblocks==TRUE){
    p1 <- p1 + geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = -Inf, ymax = min(d_order$pval), fill=as.factor(Color)), alpha = 1)
  }
  #Add legend
  p1 <- p1 + scale_colour_manual(name = "Color", values = topcols) + scale_fill_manual(name = "Color", values = topcols, guides(alpha=FALSE))
  p1 <- p1 + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x=element_blank(), axis.title.x=element_blank(), legend.position="top", legend.title=element_blank())
  
  #Highlight if given
  if(!missing(highlight_snp)){
    if("Shape" %in% topn){
      p1 <- p1 + geom_point(data=d_order[d_order$SNP %in% highlight_snp & d_order$Location=="Top", ], aes(x=pos_index, y=pval, shape=Shape), colour=highlighter)
      p1 <- p1 + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p1 <- p1 + geom_point(data=d_order[d_order$SNP %in% highlight_snp & d_order$Location=="Top", ], aes(x=pos_index, y=pval), colour=highlighter)
    }
  }
  if(!missing(highlight_p)){
    if("Shape" %in% topn){
      p1 <- p1 + geom_point(data=d_order[d_order$pvalue < highlight_p[1] & d_order$Location=="Top", ], aes(x=pos_index, y=pval, shape=Shape), colour=highlighter)
      p1 <- p1 + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p1 <- p1 + geom_point(data=d_order[d_order$pvalue < highlight_p[1] & d_order$Location=="Top", ], aes(x=pos_index, y=pval), colour=highlighter)
    }
  }
  #Add pvalue threshold line
  if(!missing(tline)){
    p1 <- p1 + geom_hline(yintercept = tredline, colour="red")
  }
  #Annotate
  if(!missing(annotate_p)){
    if (!requireNamespace(c("ggrepel"), quietly = TRUE)==TRUE) {
      print("Consider installing 'ggrepel' for improved text annotation")
      p1 <- p1 + geom_text(data=d_order[d_order$pvalue < annotate_p & d_order$Location=="Top",], aes(pos_index,pval,label=SNP))
    } else {
      p1 <- p1 + ggrepel::geom_text_repel(data=d_order[d_order$pvalue < annotate_p & d_order$Location=="Top",], aes(pos_index,pval,label=SNP))
    }
  }
  if(!missing(annotate_snp)){
    if (!requireNamespace(c("ggrepel"), quietly = TRUE)==TRUE){
      print("Consider installing 'ggrepel' for improved text annotation")
      p1 <- p1 + geom_text(data=d_order[d_order$SNP %in% annotate_snp & d_order$Location=="Top",], aes(pos_index,pval,label=SNP))
    } else {
      p1 <- p1 + ggrepel::geom_text_repel(data=d_order[d_order$SNP %in% annotate_snp & d_order$Location=="Top",], aes(pos_index,pval,label=SNP))
    }
  }
  #Add title and y axis title
  p1 <- p1 + ylab(yaxislab1)
  
  if(chrblocks==TRUE){
    if(freey==TRUE){
      print("Sorry, drawing chrblocks with freey=TRUE is currently unsupported and will be ignored.")
    } else {
      p1 <- p1+theme(axis.text.x = element_text(vjust=1),axis.ticks.x = element_blank())+ylim(c(yaxismin1,yaxismax1))
    }
  } else {
    p1 <- p1+theme(axis.text.x = element_text(vjust=1),axis.ticks.x = element_blank())+ scale_y_continuous(limits=c(yaxismin1, yaxismax1),expand=expand_scale(mult=c(0,0.1)))
  }
  
  if(background=="white"){
    p1 <- p1 + theme(panel.background = element_rect(fill="white"))
  }

  if(!is.null(toptitle)) {
    p1 <- p1 + ggtitle(toptitle) + theme(plot.title = element_text(hjust = 0.5, size = 14))
  }
  return(p1)
}

