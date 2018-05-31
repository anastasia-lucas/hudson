#' gmirror
#'
#' Create mirrored Manhattan plots for GWAS
#' Dependencies: ggplot2, gridExtra
#' Suggested: RColorBrewer, ggrepel
#' @param top data frame, if not plato or plink format, must contain SNP, CHR, POS, pvalue, optional Shape
#' @param bottom data frame, if not plato or plink format, must contain SNP, CHR, POS, pvalue, optional Shape
#' @param line optional pvalue threshold to draw red line at
#' @param log10 plot -log10() of pvalue column, boolean
#' @param yaxis label for y-axis, automatically set if log10=TRUE
#' @param opacity opacity of points, from 0 to 1, useful for dense plots
#' @param annotate_snp list of RSIDs to annotate
#' @param annotate_p pvalue threshold to annotate
#' @param toptitle optional string for top plot title
#' @param bottomtitle optional string for bottom plot title
#' @param chrcolor1 first alternating color for chromosome
#' @param chrcolor2 second alternating color for chromosome
#' @param highlight_snp list of snps to highlight
#' @param highlight_p pvalue threshold to highlight
#' @param highlighter color to highlight
#' @param file file name of saved image
#' @param hgt height of plot in inches
#' @param wi width of plot in inches
#' @param res resolution of plot in pixels per inch
#' @return png image
#' @export
#' @examples
#' gmirror(top, bottom, line, log10, yaxis, opacity, annotate_snp, annotate_p, title, chrcolor1, chrcolor2, groupcolors, file, hgt, wi, res)

gmirror <- function(top, bottom, line, log10=TRUE, yaxis, opacity=1, annotate_snp, annotate_p, toptitle=NULL, bottomtitle=NULL, highlight_snp, highlight_p, highlighter="red", chrcolor1="#AAAAAA", chrcolor2="#4D4D4D", groupcolors, file="gmirror", hgt=7, wi=12, res=300 ){
  if (!requireNamespace(c("ggplot2"), quietly = TRUE)==TRUE | !requireNamespace(c("gridExtra"), quietly = TRUE)==TRUE) {
    stop("Please install ggplot2 and gridExtra to create visualization.", call. = FALSE)
  } else {
    require("ggplot2", quietly = TRUE)
    require("gridExtra", quietly = TRUE)
  }

  #Sort data
  topn <- names(top)
  bottomn <- names(bottom)
  top$Location <- "Top"
  bottom$Location <- "Bottom"
  if("Shape" %in% colnames(top) & !("Shape" %in% colnames(bottom))){
    bottom$Shape <- NA
    relhgts <- c(6/11, 5/11)
  } else if("Shape" %in% colnames(bottom) & !("Shape" %in% colnames((top)))){
      top$Shape <- NA
      relhgts <- c(5/11, 6/11)
  } else {
    relhgts <- c(1/2, 1/2)
  }
  d <- rbind(top, bottom)
  d$POS <- as.numeric(as.character(d$POS))
  d$CHR <- factor(d$CHR, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"))
  d_order <- d[order(d$CHR, d$POS), ]
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
  lims$shademap <- rep(c("shade_ffffff", "shade_ebebeb"), each=1)

  #Set up colors
  nchrcolors <- nlevels(factor(lims$Color))

  #Color by CHR
  colnames(d_order)[2] <- "Color"
  newcols <-c(rep(x=c(chrcolor1, chrcolor2), length.out=nchrcolors, each=1), "#FFFFFF", "#EBEBEB")
  names(newcols) <-c(levels(factor(lims$Color)), "shade_ffffff", "shade_ebebeb")

  #Info for y-axis
  if(log10==TRUE){
    d_order$pval <- -log10(d_order$pvalue)
    yaxislab <- expression(paste("-log"[10], "(p-value)", sep=""))
    if(!missing(line)) {redline <- -log10(line)}
  } else {
    d_order$pval <- d_order$pvalue
    yaxislab <- yaxis
    if(!missing(line)) {redline <- line}
  }
  yaxismax <- max(d_order$pval[which(d_order$pval< Inf)])

  #Start plotting
  #TOP PLOT
  p1 <- ggplot() + geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = 0, ymax = Inf, fill=factor(shademap)), alpha = 0.5)
  #Add shape info if available
  if("Shape" %in% topn){
    p1 <- p1 + geom_point(data=d_order[d_order$Location=="Top",], aes(x=pos_index, y=pval, color=factor(Color), shape=factor(Shape)), alpha=opacity)
  } else {
    p1 <- p1 + geom_point(data=d_order[d_order$Location=="Top",], aes(x=pos_index, y=pval, color=factor(Color)), alpha=opacity)
  }
  if(!missing(annotate_p)){
    if (!requireNamespace(c("ggrepel"), quietly = TRUE)==TRUE) {
      print("Consider installing 'ggrepel' for improved text annotation")
      p1 <- p1 + geom_text(data=d_order[d_order$pvalue < annotate_p & d_order$Location=="Top",], aes(pos_index,pval,label=SNP))
    } else {
      require("ggrepel", quietly = TRUE)
      p1 <- p1 + geom_text_repel(data=d_order[d_order$pvalue < annotate_p & d_order$Location=="Top",], aes(pos_index,pval,label=SNP))
    }
  }
  if(!missing(annotate_snp)){
    if (!requireNamespace(c("ggrepel"), quietly = TRUE)==TRUE){
      print("Consider installing 'ggrepel' for improved text annotation")
      p1 <- p1 + geom_text(data=d_order[d_order$SNP %in% annotate_snp & d_order$Location=="Top",], aes(pos_index,pval,label=SNP))
    } else {
      require("ggrepel", quietly = TRUE)
      p1 <- p1 + geom_text_repel(data=d_order[d_order$SNP %in% annotate_snp & d_order$Location=="Top",], aes(pos_index,pval,label=SNP))
    }
  }
  p1 <- p1 + scale_x_continuous(breaks=lims$av, labels=lims$Color, expand=c(0,0))
  p1 <- p1 + geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = -Inf, ymax = 0, fill=as.factor(Color)), alpha = 1)
  p1 <- p1 + scale_colour_manual(name = "Color", values = newcols, guides(alpha=FALSE)) + scale_fill_manual(name = "Color", values = newcols, guides(alpha=FALSE))
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
      p1 <- p1 + geom_point(data=d_order[d_order$pvalue < highlight_p & d_order$Location=="Top", ], aes(x=pos_index, y=pval, shape=Shape), colour=highlighter)
      p1 <- p1 + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p1 <- p1 + geom_point(data=d_order[d_order$pvalue < highlight_p & d_order$Location=="Top", ], aes(x=pos_index, y=pval), colour=highlighter)
    }
  }
  #Add title and y axis title
  p1 <- p1 + ylab(yaxislab)
  #Add pvalue threshold line
  if(!missing(line)){p1 <- p1 + geom_hline(yintercept = redline, colour="red")}


  #BOTTOM PLOT
  p2 <- ggplot() + geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = 0, ymax = Inf, fill=factor(shademap)), alpha = 0.5)
  #Add shape info if available
  if("Shape" %in% bottomn){
    p2 <- p2 + geom_point(data=d_order[d_order$Location=="Bottom",], aes(x=pos_index, y=pval, color=factor(Color), shape=factor(Shape)), alpha=opacity)
  } else {
    p2 <- p2 + geom_point(data=d_order[d_order$Location=="Bottom",], aes(x=pos_index, y=pval, color=factor(Color)), alpha=opacity)
  }
  if(!missing(annotate_p)){
    if (!requireNamespace(c("ggrepel"), quietly = TRUE)==TRUE) {
      print("Consider installing 'ggrepel' for improved text annotation")
      p2 <- p2 + geom_text(data=d_order[d_order$pvalue < annotate_p & d_order$Location=="Bottom",], aes(pos_index,pval,label=SNP))
    } else {
      require("ggrepel", quietly = TRUE)
      p2 <- p2 + geom_text_repel(data=d_order[d_order$pvalue < annotate_p & d_order$Location=="Bottom",], aes(pos_index,pval,label=SNP))
    }
  }
  if(!missing(annotate_snp)){
    if (!requireNamespace(c("ggrepel"), quietly = TRUE)==TRUE){
      print("Consider installing 'ggrepel' for improved text annotation")
      p2 <- p2 + geom_text(data=d_order[d_order$SNP %in% annotate_snp & d_order$Location=="Bottom",], aes(pos_index,pval,label=SNP))
    } else {
      require("ggrepel", quietly = TRUE)
      p2 <- p2 + geom_text_repel(data=d_order[d_order$SNP %in% annotate_snp & d_order$Location=="Bottom",], aes(pos_index,pval,label=SNP))
    }
  }
  p2 <- p2 + scale_x_continuous(breaks=lims$av, labels=lims$Color, expand=c(0,0))
  p2 <- p2 + geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = -Inf, ymax = 0, fill=as.factor(Color)), alpha = 1)
  p2 <- p2 + scale_colour_manual(name = "Color", values = newcols, guides(alpha=FALSE)) + scale_fill_manual(name = "Color", values = newcols, guides(alpha=FALSE))
  p2 <- p2 + theme(axis.text.x=element_text(angle=90), panel.grid.minor.x = element_blank(), panel.grid.major.x=element_blank(), axis.title.x=element_blank(), legend.position="bottom", legend.title=element_blank())
  #Highlight if given
  if(!missing(highlight_snp)){
    if("Shape" %in% bottomn){
      p2 <- p2 + geom_point(data=d_order[d_order$SNP %in% highlight_snp & d_order$Location=="Bottom", ], aes(x=pos_index, y=pval, shape=Shape), colour=highlighter)
      p2 <- p2 + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p2 <- p2 + geom_point(data=d_order[d_order$SNP %in% highlight_snp & d_order$Location=="Bottom", ], aes(x=pos_index, y=pval), colour=highlighter)
    }
  }
  if(!missing(highlight_p)){
    if("Shape" %in% bottomn){
      p2 <- p2 + geom_point(data=d_order[d_order$pvalue < highlight_p & d_order$Location=="Bottom", ], aes(x=pos_index, y=pval, shape=Shape), colour=highlighter)
      p2 <- p2 + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p2 <- p2 + geom_point(data=d_order[d_order$pvalue < highlight_p & d_order$Location=="Bottom", ], aes(x=pos_index, y=pval), colour=highlighter)
    }
  }
  #Add title and y axis title
  p2 <- p2 + ylab(yaxislab)
  #Add pvalue threshold line
  if(!missing(line)){p2 <- p2 + geom_hline(yintercept = redline, colour="red")}

  #Format
  p1 <- p1+theme(axis.text.x = element_text(vjust=1),axis.ticks.x = element_blank())+ylim(c(0,yaxismax))
  p2 <- p2+scale_y_reverse(limits=c(yaxismax,0)) + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())

  #Save
  print(paste("Saving plot to ", file, ".png", sep=""))
  p <- grid.arrange(arrangeGrob(p1, top=toptitle), arrangeGrob(p2, bottom=bottomtitle), padding=0, heights=relhgts)
  ggsave(p, filename=paste(file, ".png", sep=""), dpi=res, units="in", height=hgt, width=wi)
  return(p)
}
