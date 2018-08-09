#' phemirror
#'
#' Create mirrored Manhattan plots for PheWAS
#' Dependencies: ggplot2, gridExtra, RColorBrewer
#' Suggested: ggrepel
#' @param top data frame, if not plato or plink format, must contain PHE, SNP, CHR, POS, pvalue, columns, optional Shape
#' @param bottom data frame, if not plato or plink format, must contain PHE, SNP, CHR, POS, pvalue, columns, optional Shape
#' @param phegroup optional grouping file for phenotypes, must contain PHE and Group columns
#' @param tline optional pvalue threshold to draw red line at in top plot
#' @param bline optional pvalue threshold to draw red line at in bottom plot
#' @param log10 plot -log10() of pvalue column, boolean
#' @param yaxis label for y-axis in the format c("top", "bottom"), automatically set if log10=TRUE
#' @param opacity opacity of points, from 0-1, useful for dense plots
#' @param annotate_snp list of RSIDs to annotate
#' @param annotate_p pvalue threshold to annotate
#' @param toptitle optional string for top plot title
#' @param bottomtitle optional string for bottom plot title
#' @param chrcolor1 first alternating color for chromosome
#' @param chrcolor2 second alternating color for chromosome
#' @param highlight_snp list of snps to highlight
#' @param highlight_p pvalue threshold to highlight
#' @param highlighter color to highlight
#' @param groupcolors named list of colors where names correspond to data in 'PHE' or 'Group' column
#' @param background variegated or white
#' @param chrblocks boolean, turns on x-axis chromosome marker blocks
#' @param file file name of saved image
#' @param hgt height of plot in inches
#' @param hgtratio height ratio of plots, equal to top plot proportion
#' @param wi width of plot in inches
#' @param res resolution of plot in pixels per inch
#' @return png image
#' @export
#' @examples
#' phemirror(top. bottom, phegroup, line, log10, yaxis, opacity, annotate_snp, annotate_p, title, chrcolor1, chrcolor2, groupcolors, file, hgt, wi, res)

phemirror <- function(top, bottom, phegroup, tline, bline, log10=TRUE, yaxis, opacity=1, annotate_snp, annotate_p, highlight_snp, highlight_p, highlighter="red", toptitle=NULL, bottomtitle=NULL, chrcolor1="#AAAAAA", chrcolor2="#4D4D4D", groupcolors, background="variegated", chrblocks=FALSE, file="phemirror", hgtratio=0.5, hgt=7, wi=12, res=300 ){
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
  if("Shape" %in% colnames(top) & !("Shape" %in% colnames(bottom))){bottom$Shape <- NA}
  if("Shape" %in% colnames(bottom) & !("Shape" %in% colnames((top)))){top$Shape <- NA}
  d <- rbind(top, bottom)
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
  lims$shademap <- rep(c("shade_ffffff", "shade_ebebeb"), each=1)

  #Set up colors
  nchrcolors <- nlevels(factor(lims$Color))
  if(!missing(groupcolors)){
    dcols <- c(rep(x=c(chrcolor1, chrcolor2), length.out=nchrcolors, each=1), "#FFFFFF", "#EBEBEB")
    names(dcols) <-c(levels(factor(lims$Color)), "shade_ffffff", "shade_ebebeb")
    topcols <- c(dcols, groupcolors)
    bottomcols <- c(dcols, groupcolors)
  } else {
    if (!requireNamespace(c("RColorBrewer"), quietly = TRUE)==TRUE) {
      stop("Please install RColorBrewer to add color attribute.", call. = FALSE)
    } else {
      require("RColorBrewer", quietly=TRUE)
    }
    getPalette = colorRampPalette(brewer.pal(11, "Spectral"))
    #Top Colors
    ngroupcolors <- nlevels(factor(d_order$Color[d_order$Location=="Top"]))
    topcols <- c(rep(x=c(chrcolor1, chrcolor2), length.out=nchrcolors, each=1), getPalette(ngroupcolors), "#FFFFFF", "#EBEBEB")
    names(topcols) <-c(levels(factor(lims$Color)), levels(factor(d_order$Color[d_order$Location=="Top"])), "shade_ffffff", "shade_ebebeb")
    #Bottom Colors
    ngroupcolors <- nlevels(factor(d_order$Color[d_order$Location=="Bottom"]))
    bottomcols <- c(rep(x=c(chrcolor1, chrcolor2), length.out=nchrcolors, each=1), getPalette(ngroupcolors), "#FFFFFF", "#EBEBEB")
    names(bottomcols) <-c(levels(factor(lims$Color)), levels(factor(d_order$Color[d_order$Location=="Bottom"])), "shade_ffffff", "shade_ebebeb")
  }

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
  yaxismax <- max(d_order$pval[which(d_order$pval< Inf)])

  #Theme options
  backpanel <- ifelse(background=="white", "NULL", "geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = min(d_order$pval), ymax = Inf, fill=factor(shademap)), alpha = 0.5)" )
  
  #Start plotting
  #TOP PLOT
  p1 <- ggplot() + eval(parse(text=backpanel))
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
  
  #Start plotting
  #BOTTOM PLOT
  p2 <- ggplot() + eval(parse(text=backpanel))
  #Add shape info if available
  if("Shape" %in% bottomn){
    p2 <- p2 + geom_point(data=d_order[d_order$Location=="Bottom",], aes(x=pos_index, y=pval, color=factor(Color), shape=factor(Shape)), alpha=opacity)
  } else {
    p2 <- p2 + geom_point(data=d_order[d_order$Location=="Bottom",], aes(x=pos_index, y=pval, color=factor(Color)), alpha=opacity)
  }
  p2 <- p2 + scale_x_continuous(breaks=lims$av, labels=lims$Color, expand=c(0,0))
  if(chrblocks==TRUE){
    p2 <- p2 + geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = -Inf, ymax = min(d_order$pval), fill=as.factor(Color)), alpha = 1)
  } 
  #Add legend
  p2 <- p2 + scale_colour_manual(name = "Color", values = bottomcols) + scale_fill_manual(name = "Color", values = bottomcols, guides(alpha=FALSE))
  p2 <- p2 + theme(axis.text.x=element_text(angle=90), panel.grid.minor.x = element_blank(), panel.grid.major.x=element_blank(), axis.title.x=element_blank(), legend.position="bottom", legend.title=element_blank())
 
 #Highlight if given
  if(!missing(highlight_snp)){
    if("Shape" %in% topn){
      p1 <- p1 + geom_point(data=d_order[d_order$SNP %in% highlight_snp & d_order$Location=="Top", ], aes(x=pos_index, y=pval, shape=Shape), colour=highlighter)
      p1 <- p1 + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p1 <- p1 + geom_point(data=d_order[d_order$SNP %in% highlight_snp & d_order$Location=="Top", ], aes(x=pos_index, y=pval), colour=highlighter)
    }
    if("Shape" %in% bottomn){
      p2 <- p2 + geom_point(data=d_order[d_order$SNP %in% highlight_snp & d_order$Location=="Bottom", ], aes(x=pos_index, y=pval, shape=Shape), colour=highlighter)
      p2 <- p2 + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p2 <- p2 + geom_point(data=d_order[d_order$SNP %in% highlight_snp & d_order$Location=="Bottom", ], aes(x=pos_index, y=pval), colour=highlighter)
    }
  }
  if(!missing(highlight_p)){
    if("Shape" %in% topn){
      p1 <- p1 + geom_point(data=d_order[d_order$pvalue < highlight_p & d_order$Location=="Top", ], aes(x=pos_index, y=pval, shape=Shape), colour=highlighter)
      p1 <- p1 + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p1 <- p1 + geom_point(data=d_order[d_order$pvalue < highlight_p & d_order$Location=="Top", ], aes(x=pos_index, y=pval), colour=highlighter)
    }
    if("Shape" %in% bottomn){
      p2 <- p2 + geom_point(data=d_order[d_order$pvalue < highlight_p & d_order$Location=="Bottom", ], aes(x=pos_index, y=pval, shape=Shape), colour=highlighter)
      p2 <- p2 + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p2 <- p2 + geom_point(data=d_order[d_order$pvalue < highlight_p & d_order$Location=="Bottom", ], aes(x=pos_index, y=pval), colour=highlighter)
    }
  }
  #Add pvalue threshold line
  if(!missing(tline)){
    p1 <- p1 + geom_hline(yintercept = tredline, colour="red")
  }
  if(!missing(bline)){
    p2 <- p2 + geom_hline(yintercept = bredline, colour="red")
  }
  #Annotate
  if(!missing(annotate_p)){
    if (!requireNamespace(c("ggrepel"), quietly = TRUE)==TRUE) {
      print("Consider installing 'ggrepel' for improved text annotation")
      p1 <- p1 + geom_text(data=d_order[d_order$pvalue < annotate_p & d_order$Location=="Top",], aes(pos_index,pval,label=SNP))
      p2 <- p2 + geom_text(data=d_order[d_order$pvalue < annotate_p & d_order$Location=="Bottom",], aes(pos_index,pval,label=SNP))
    } else {
      require("ggrepel", quietly = TRUE)
      p1 <- p1 + geom_text_repel(data=d_order[d_order$pvalue < annotate_p & d_order$Location=="Top",], aes(pos_index,pval,label=SNP))
      p2 <- p2 + geom_text_repel(data=d_order[d_order$pvalue < annotate_p & d_order$Location=="Bottom",], aes(pos_index,pval,label=SNP))
    }
  }
  if(!missing(annotate_snp)){
    if (!requireNamespace(c("ggrepel"), quietly = TRUE)==TRUE){
      print("Consider installing 'ggrepel' for improved text annotation")
      p1 <- p1 + geom_text(data=d_order[d_order$SNP %in% annotate_snp & d_order$Location=="Top",], aes(pos_index,pval,label=SNP))
      p2 <- p2 + geom_text(data=d_order[d_order$SNP %in% annotate_snp & d_order$Location=="Bottom",], aes(pos_index,pval,label=SNP))
    } else {
      require("ggrepel", quietly = TRUE)
      p1 <- p1 + geom_text_repel(data=d_order[d_order$SNP %in% annotate_snp & d_order$Location=="Top",], aes(pos_index,pval,label=SNP))
      p2 <- p2 + geom_text_repel(data=d_order[d_order$SNP %in% annotate_snp & d_order$Location=="Bottom",], aes(pos_index,pval,label=SNP))
    }
  }
  #Add title and y axis title
  p1 <- p1 + ylab(yaxislab1)
  p2 <- p2 + ylab(yaxislab2)

  if(chrblocks==TRUE){
    p1 <- p1+theme(axis.text.x = element_text(vjust=1),axis.ticks.x = element_blank())+ylim(c(0,yaxismax))
    p2 <- p2+scale_y_reverse(limits=c(yaxismax, 0)) + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
    
  } else {
    p1 <- p1+theme(axis.text.x = element_text(vjust=1),axis.ticks.x = element_blank())+ scale_y_continuous(limits=c(0, yaxismax),expand=expand_scale(mult=c(0,0.1)))
    p2 <- p2+scale_y_reverse(limits=c(yaxismax,0), expand=expand_scale(mult=c(0.1,0))) + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
  }
  
  if(background=="white"){
    p1 <- p1 + theme(panel.background = element_rect(fill="white"))
    p2 <- p2 + theme(panel.background = element_rect(fill="white"))
  }

  #Save
  print(paste("Saving plot to ", file, ".png", sep=""))
  p <- grid.arrange(arrangeGrob(p1, top=toptitle), arrangeGrob(p2, bottom=bottomtitle), padding=0, heights=c(hgtratio, 1-hgtratio))
  ggsave(p, filename=paste(file, ".png", sep=""), dpi=res, units="in", height=hgt, width=wi)

  return(p)

}
