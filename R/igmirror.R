#' igmirror
#'
#' Create interactive mirrored Manhattan plots for GWAS
#' Dependencies: ggplot2, ggiraph
#' Suggested: ggrepel
#' @param top data frame, must contain SNP, CHR, POS, pvalue columns, optional Shape, Hover, Link
#'                        where 'Hover' contains tooltip information overriding the default (SNP, chr:position)
#'                        and 'Link' contains a search query, ex. "https://www.ncbi.nlm.nih.gov/snp/rs944173"
#' @param bottom data frame, must contain SNP, CHR, POS, pvalue, columns, optional Shape, Hover, Link
#'                           where 'Hover' contains tooltip information overriding the default (SNP, chr:position)
#'                           and 'Link' contains a search query, ex. "https://www.ncbi.nlm.nih.gov/snp/rs944173"
#' @param tline list of pvalues to draw red threshold lines in top plot
#' @param bline list of pvalues to draw red threshold lines in bottom plot
#' @param chroms list of chromosomes to plot in the order desired, default c(1:22, "X", "Y")
#' @param log10 plot -log10() of pvalue column, logical
#' @param yaxis title for y-axis, automatically set if log10=TRUE
#' @param opacity opacity of points, from 0-1, useful for dense plots
#' @param title figure title
#' @param chrcolor1 first alternating color for chromosome
#' @param chrcolor2 second alternating color for chromosome
#' @param chrblockmin minimum y value for chromosome block
#' @param chrblockmax maximum y value for chromosome block
#' @param highlight_snp vector of snps to highlight
#' @param highlight_p list of pvalue thresholds to annotate in the order of c(p_top, p_bottom)
#' @param highlighter color to highlight
#' @param freey allow y-axes to scale with the data
#' @param ymax override the auto set y axis limit
#' @param ymin override the auto set y axis limit
#' @param background variegated or white
#' @param chrblocks logical, turns on x-axis chromosome marker blocks
#' @param file file name of saved image
#' @param hgt height of plot in inches
#' @param wi width of plot in inches
#' @param res resolution of plot in pixels per inch
#' @return png image
#' @import ggplot2
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @export
#' @examples
#' data(gwas.t)
#' data(gwas.b)
#' gwas.t$Link <- paste0("https://www.ncbi.nlm.nih.gov/snp/", gwas.t$SNP)
#' gwas.b$Link <- paste0("https://www.ncbi.nlm.nih.gov/snp/", gwas.b$SNP)
#' igmirror(gwas.t, gwas.b)

igmirror <- function(top, bottom, tline, bline, chroms = c(1:22,"X","Y"),
                     log10=TRUE, yaxis, opacity=1, highlight_snp, 
                     highlight_p, highlighter="red", title=NULL, 
                     chrcolor1="#AAAAAA", chrcolor2="#4D4D4D", chrblockmin=-1,
                     chrblockmax=1, freey=FALSE, ymax, ymin,
                     background="variegated", chrblocks=FALSE, file="igmirror", 
                     hgt=7, wi=12, res=300 ){
  
  #Sort data
  topn <- names(top)
  bottomn <- names(bottom)
  top$Location <- "Top"
  bottom$Location <- "Bottom"
  
  #Check file formats
  if(!identical(topn, bottomn)){stop("Please ensure both inputs have the same metadata columns.")}
  
  d <- rbind(top, bottom)
  
  #Set onclick to NULL if needed
  if(!("Hover" %in% names(d))){
    d$Hover <- paste0("SNP: ", d$SNP,
                     "\nPosition: ", paste(d$CHR, d$POS, sep=":"))
  }
  if(!("Link" %in% names(d))){
    d$Link <- NULL
  } else {
    d$Link <- sprintf("window.open(\"%s\")", d$Link)
  }
  
  #Sort data
  d$CHR <- droplevels(factor(d$CHR, levels = as.character(chroms)))
  d <- d[d$CHR %in% chroms, ]
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
  lims$shademap <- rep(c("shade_ffffff", "shade_ebebeb"), length.out=nrow(lims), each=1)
  
  #Set up colors
  nchrcolors <- nlevels(factor(lims$Color))
  base_color <- c(rep(x=c(chrcolor1, chrcolor2), length.out=nchrcolors, each=1), "#FFFFFF", "#EBEBEB")
  names(base_color) <- c(levels(factor(lims$Color)), "shade_ffffff", "shade_ebebeb")
  
  #if(missing(levs)){
  #Set up colors
  nchrcolors <- nlevels(factor(lims$Color))
  
  #Color by CHR
  colnames(d_order)[2] <- "Color"
  newcols <-rep(x=c(chrcolor1, chrcolor2), length.out=nchrcolors, each=1)
  names(newcols) <-levels(factor(lims$Color))
  
  #Allow more than 6 shapes
  #3, 4 and 7 to 14 are composite symbols- incompatible with ggiraph
  if("Shape" %in% names(d)){
    allshapes <- c(16,15,17,18,0:2,5:6,19:25,33:127)
    shapevector <- allshapes[1:nlevels(as.factor(d$Shape))]
  }
  
  #Info for y-axis
  if(log10==TRUE){
    d_order$pval <- ifelse(d_order$Location=="Top", -log10(d_order$pvalue), log10(d_order$pvalue))
    yaxislab <- expression(paste("+/-log"[10], "(p-value)", sep=""))
    if(!missing(tline)) {bredline <- log10(tline)}#log10(line[2])}
    if(!missing(bline)) {tredline <- log10(bline)}#-log10(line[1])}
  } else {
    d_order$pval <- d_order$pvalue
    yaxislab <- yaxis
    if(!missing(tline)) {bredline <- tline}#log10(line[2])}
    if(!missing(bline)) {tredline <- bline}#-log10(line[1])}
  }
  
  #Theme options
  backpanel <- ifelse(background=="white", "NULL", 
                      "geom_rect(data = lims, aes(xmin = posmin, xmax = posmax, ymin = -Inf, ymax = Inf, fill=factor(shademap)), alpha = 0.3)" )
  
  #Start plotting
  p <- ggplot() + eval(parse(text=backpanel))
  #Add shape info if available
  if("Shape" %in% names(d)){
    
    p <- p + geom_point(data=d_order[is.na(d_order$Hover),],
                        aes(x=pos_index, y=pval,
                            color=factor(Color), shape=factor(Shape)),
                        alpha=opacity) +
      ggiraph::geom_point_interactive(data=d_order[!is.na(d_order$Hover),],
                                      aes(x=pos_index, y=pval, tooltip=Hover,
                                          color=factor(Color), shape=factor(Shape),
                                          onclick=Link),
                                      alpha=opacity) +
      scale_shape_manual(values=shapevector)
  } else {
    p <- p + geom_point(data=d_order[is.na(d_order$Hover),],
                        aes(x=pos_index, y=pval,
                            color=factor(Color)),
                        alpha=opacity) +
      ggiraph::geom_point_interactive(data=d_order[!is.na(d_order$Hover),],
                                      aes(x=pos_index, y=pval,
                                          color=factor(Color), tooltip=Hover,
                                          onclick=Link),
                                      alpha=opacity)
  }
  #ticklabs <- c(as.character(rep(1:15)), " ", "17", " ", "19", " ", "21", " ")
  #p <- p + scale_x_continuous(breaks=lims$av, labels=ticklabs, expand=c(0,0))
  p <- p + scale_x_continuous(breaks=lims$av, labels=lims$Color, expand=c(0,0))
  if(chrblocks==TRUE){p <- p + geom_rect(data = lims, aes(xmin = posmin, xmax = posmax, ymin = chrblockmin, ymax = chrblockmax, fill=as.factor(Color)), alpha = 1)}
  #Add legend
  p <- p + scale_colour_manual(name = "Color", values = newcols) + scale_fill_manual(values = base_color)
  p <- p + theme(panel.grid.minor.y = element_blank(),
                 panel.grid.major.y=element_blank(),
                 legend.title=element_blank(),
                 legend.position = "bottom")
  #Highlight if given
  if(!missing(highlight_snp)){
    if("Shape" %in% names(d)){
      p <- p + geom_point(data=d_order[d_order$SNP %in% highlight_snp, ], aes(x=pos_index, y=pval, shape=Shape), colour=highlighter) + scale_shape_manual(values=shapevector)
      p <- p + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p <- p + geom_point(data=d_order[d_order$SNP %in% highlight_snp, ], aes(x=pos_index, y=pval), colour=highlighter)
    }
  }
  if(!missing(highlight_p)){
    if("Shape" %in% names(d)){
      p <- p + geom_point(data=d_order[d_order$pvalue < highlight_p, ], aes(x=pos_index, y=pval, shape=Shape), colour=highlighter) + scale_shape_manual(values=shapevector)
      p <- p + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p <- p + geom_point(data=d_order[d_order$pvalue < highlight_p, ], aes(x=pos_index, y=pval), colour=highlighter)
    }
  }
  
  #Add title and y axis title
  p <- p + ggtitle(title) + ylab(yaxislab) + xlab("Chromosome")
  #Add pvalue threshold line
  if(!missing(tline)){p <- p + geom_hline(yintercept = tredline, colour="red")}
  if(!missing(bline)){p <- p + geom_hline(yintercept = bredline, colour="red")}
  
  #Theme
  if(!missing(ymax)){
    yaxismax <- ymax
  } else {
    yaxismax <- max(d_order$pval[d_order$pvalue!=0 ])
  }
  
  if(!missing(ymin)){
    yaxismin <- ymin
  } else {
    yaxismin <- min(d_order$pval[d_order$pval!=-Inf])
  }
  
  if("Shape" %in% names(d_order)){
    p <- p + scale_shape_manual(values=shapevector)
  }
  
  #Lims
  # if(equal_axis==TRUE){
  #   if(!missing(hard_limit)){
  #     true_lim <- hard_limit
  #   } else {
  #     true_lim <- max(abs(yaxismin), abs(yaxismax))
  #   }
  #   yaxismin <- -true_lim
  #   yaxismax <- true_lim
  # }
  
  if(chrblocks==TRUE){
    p <- p+ylim(c(yaxismin,yaxismax))
  } else {
    p <- p+scale_y_continuous(limits=c(yaxismin, yaxismax), expand=expand_scale(mult=c(0,0.1)))
    p <- p + geom_hline(yintercept = 0, color="black")
  }
  
  
  # annotations <- data.frame(
  #   xpos = c(Inf,Inf),
  #   ypos =  c(-Inf,Inf),
  #   annotateText = c(toptitle, bottomtitle),
  #   hjustvar = c(0,1),
  #   vjustvar = c(1,1))
  
  #if(background=="white"){p <- p + theme(panel.background = element_rect(fill="white"))}
  #p <- p + labs(caption = "some text") + 
  #  theme(plot.caption = element_text(hjust=1),
  #        plot.caption.position="panel")
  p <- p + theme(panel.background = element_rect(fill="#F8F8F8"))
  p <- p + guides(fill="none", color="none")
  p <- p + ggtitle(title)
  
  # p <- p + geom_text(data=annotations,
  #                    aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText))
  
  #Save
  print(paste0("Saving plot to ", file, ".html", sep=""))
  ip <- ggiraph::girafe(code = print(p))
  htmlwidgets::saveWidget(widget=ip, file=paste(file, ".html", sep=""))
  return("Finished")
  
}