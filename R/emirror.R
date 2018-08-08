#' emirror
#'
#' Create mirrored Manhattan plots for EWAS
#' Note: There is an issue with dev.off() if using RStudio
#' Dependencies: ggplot2, gridExtra
#' Suggested: RColorBrewer, ggrepel
#' @param top data frame, columns one and two must be Variable and pvalue; Group, Shape and Color optional
#' @param bottom data frame, columns one and two must be Variable and pvalue; Group, Shape and Color optional
#' @param tline optional pvalue threshold to draw red line at in top plot
#' @param bline optional pvalue threshold to draw red line at in bottom plot
#' @param log10 plot -log10() of pvalue column, boolean
#' @param yaxis label for y-axis, automatically set if log10=TRUE
#' @param opacity opacity of points, from 0 to 1, useful for dense plots
#' @param toptitle optional string for plot title
#' @param bottomtitle optional string for plot title
#' @param color1 first alternating color
#' @param color2 second alternating color
#' @param annotate_var list of variables to annotate
#' @param annotate_p pvalue threshold to annotate
#' @param highlight_var list of variables to highlight
#' @param highlight_p pvalue threshold to highlight
#' @param highlighter color to highlight
#' @param groupcolors named list of colors where names correspond to data in 'Color' column
#' @param rotatelabel boolean, rotate axis labels?
#' @param labelangle angle to rotate
#' @param file file name of saved image
#' @param hgtratio height ratio of plots, equal to top plot proportion
#' @param hgt height of plot in inches
#' @param wi width of plot in inches
#' @param res resolution of plot in pixels per inch
#' @return png image(s)
#' @export
#' @examples
#' emirror(top, bottom, line, file="emirror", hgt=7, wi=12, res=300 )

emirror <- function(top, bottom,  tline, bline, log10=TRUE, yaxis, opacity=1, toptitle=NULL, bottomtitle=NULL, annotate_var, annotate_p, highlight_var, highlight_p, highlighter="red", color1="#AAAAAA", color2="#4D4D4D", groupcolors, rotatelabels=FALSE, labelangle, file="emirror", hgtratio=0.5, hgt=7, wi=12, res=300){
  if (!requireNamespace(c("ggplot2"), quietly = TRUE)==TRUE | !requireNamespace(c("gridExtra"), quietly = TRUE)==TRUE) {
    stop("Please install ggplot2 and gridExtra to create visualization.", call. = FALSE)
  } else {
    require("ggplot2", quietly = TRUE)
    require("gridExtra", quietly = TRUE)
  }

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
    relhgts <- relhgts
  }
  d <- rbind(top, bottom)

  #Info for y-axis
  if(log10==TRUE){
    d$pval <- -log10(d$pvalue)
    yaxislab <- expression(paste("-log"[10], "(p-value)", sep=""))
    if(!missing(tline)) {tredline <- -log10(tline)}
    if(!missing(bline)) {bredline <- -log10(bline)}
  } else {
    d$pval <- d$pvalue
    yaxislab <- yaxis
    if(!missing(tline)) {tredline <- tline}
    if(!missing(bline)) {bredline <- bline}
  }
  yaxismax <- max(d$pval[which(d$pval< Inf)])

  #Save to merge later
  d$rowid <- seq.int(nrow(d))
  dinfo <- d[, colnames(d) %in% c("rowid", "Color", "pval", "Location", "Shape"), drop=FALSE]

  #If no group, plot raw data
  if(!"Group" %in% colnames(d)){
    print("Warning: Mirror plots contain a shared x-axis, please add grouping information if the axes are different")
    d_order <- merge(d[, colnames(d) %in% c("Variable", "pvalue", "rowid")], dinfo, by="rowid")
    if("Shape" %in% names(d)){
      if("Color" %in% names(d)){
         p1 <- ggplot() + geom_point(data=d_order[d_order$Location=="Top",], aes(x=factor(Variable), y=pval, shape=factor(Shape), color=Color), alpha=opacity)
        p2 <- ggplot() + geom_point(data=d_order[d_order$Location=="Bottom",], aes(x=factor(Variable), y=pval, shape=factor(Shape), color=Color), alpha=opacity)
      } else {
        p1 <- ggplot() + geom_point(data=d_order[d_order$Location=="Top",], aes(x=factor(Variable), y=pval, shape=factor(Shape)), alpha=opacity)
        p2 <- ggplot() + geom_point(data=d_order[d_order$Location=="Bottom",], aes(x=factor(Variable), y=pval, shape=factor(Shape)), alpha=opacity)
      }
      p1 <- p1 + theme(axis.title.x=element_blank(), legend.position="top", legend.title=element_blank())
      p2 <- p2 + theme(axis.title.x=element_blank(), legend.position="bottom", legend.title=element_blank())
    } else {
      if("Color" %in% names(d)){
        p1 <- ggplot(d_order[d_order$Location=="Top",], aes(x=factor(Variable), y=pval, color=Color)) + geom_point()
        p1 <- p1 + theme(axis.title.x=element_blank(), legend.position="top", legend.title=element_blank())
        p2 <- ggplot(d_order[d_order$Location=="Bottom",], aes(x=factor(Variable), y=pval, color=Color)) + geom_point()
        p2 <- p2 + theme(axis.title.x=element_blank(), legend.position="bottom", legend.title=element_blank())
      } else {
        p1 <- ggplot(d_order[d_order$Location=="Top",], aes(x=factor(Variable), y=pval)) + geom_point() + theme(axis.title.x=element_blank())
        p2 <- ggplot(d_order[d_order$Location=="Bottom",], aes(x=factor(Variable), y=pval)) + geom_point() + theme(axis.title.x=element_blank())
      }
    }
  } else {
    #Create position index
    subd <- d[, c("Variable", "Group", "pvalue", "rowid")]
    d_order <- subd[order(subd$Group, subd$Variable),]
    d_order$pos_index <- seq.int(nrow(d_order))

    #Set up dataframe with position info
    subd <- d_order[, colnames(d_order)!="rowid"]
    maxRows <- by(subd, subd$Group, function(x) x[which.max(x$pos_index),])
    minRows <- by(subd, subd$Group, function(x) x[which.min(x$pos_index),])
    milimits <- do.call(rbind, minRows)
    malimits <- do.call(rbind, maxRows)
    lims <- merge(milimits, malimits, by="Group")
    names(lims) <- c("Color", "Varx", "px", "posmin", "Vary", "py", "posmax")
    lims$av <- (lims$posmin + lims$posmax)/2
    lims$shademap <- rep(c("shade_ffffff","shade_ebebeb"), each=1, length.out=nrow(lims))

    #Set up colors
    nvarcolors <- nlevels(factor(lims$Color))
    if("Color" %in% names(d)){
      #Color by Color column
      if(!missing(groupcolors)){
        dcols <- c(rep(x=c(color1, color2), length.out=nvarcolors, each=1), "#FFFFFF", "#EBEBEB")
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
        ngroupcolors <- nlevels(factor(d$Color[d$Location=="Top"]))
        topcols <- c(rep(x=c(color1, color2), length.out=nvarcolors, each=1), getPalette(ngroupcolors), "#FFFFFF", "#EBEBEB")
        names(topcols) <-c(levels(factor(lims$Color)), levels(factor(d$Color[d$Location=="Top"])), "shade_ffffff", "shade_ebebeb")
        #Bottom Colors
        ngroupcolors <- nlevels(factor(d$Color[d$Location=="Bottom"]))
        bottomcols <- c(rep(x=c(color1, color2), length.out=nvarcolors, each=1), getPalette(ngroupcolors), "#FFFFFF", "#EBEBEB")
        names(bottomcols) <-c(levels(factor(lims$Color)), levels(factor(d$Color[d$Location=="Bottom"])), "shade_ffffff", "shade_ebebeb")

      }
    } else {
      #Color by Group instead
      names(d_order)[names(d_order)=="Group"] <- "Color"
      topcols <-c(rep(x=c(color1, color2), length.out=nvarcolors, each=1), "#FFFFFF", "#EBEBEB")
      names(topcols) <-c(levels(factor(lims$Color)), "shade_ffffff", "shade_ebebeb")
      bottomcols <-c(rep(x=c(color1, color2), length.out=nvarcolors, each=1), "#FFFFFF", "#EBEBEB")
      names(bottomcols) <-c(levels(factor(lims$Color)), "shade_ffffff", "shade_ebebeb")
    }

    #Start plotting
    d_order <- merge(d_order, dinfo, by="rowid")

    #TOP PLOT
    p1 <- ggplot() + geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = 0, ymax = Inf, fill=factor(shademap)), alpha = 0.5)
    #Add shape info if available
    if("Shape" %in% names(d)){
      p1 <- p1 + geom_point(data=d_order[d_order$Location=="Top",], aes(x=pos_index, y=pval, color=Color, shape=factor(Shape)), alpha=opacity)
    } else {
      p1 <- p1 + geom_point(data=d_order[d_order$Location=="Top",], aes(x=pos_index, y=pval, color=Color), alpha=opacity)
    }
    p1 <- p1 + scale_x_continuous(breaks=lims$av, labels=lims$Color, expand=c(0,0))
    p1 <- p1 + geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = -Inf, ymax = 0, fill=Color), alpha = 1)
    p1 <- p1 + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x=element_blank(), axis.title.x=element_blank(), legend.position="top", legend.title=element_blank())

    #BOTTOM PLOT
    p2 <- ggplot() + geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = 0, ymax = Inf, fill=factor(shademap)), alpha = 0.5)
    #Add shape info if available
    if("Shape" %in% bottomn){
      p2 <- p2 + geom_point(data=d_order[d_order$Location=="Bottom",], aes(x=pos_index, y=pval, color=Color, shape=factor(Shape)), alpha=opacity)
    } else {
      p2 <- p2 + geom_point(data=d_order[d_order$Location=="Bottom",], aes(x=pos_index, y=pval, color=Color), alpha=opacity)
    }
    p2 <- p2 + scale_x_continuous(breaks=lims$av, labels=lims$Color, expand=c(0,0))
    p2 <- p2 + geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = -Inf, ymax = 0, fill=Color), alpha = 1)
    p2 <- p2 + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x=element_blank(), axis.title.x=element_blank(), legend.position="bottom", legend.title=element_blank())
  }
  if("Color" %in% names(d)){
    #Add legend
    p1 <- p1 + scale_colour_manual(name = "Color", values = topcols) + scale_fill_manual(name = "Color", values = topcols, guides(alpha=FALSE))
    p2 <- p2 + scale_colour_manual(name = "Color", values = bottomcols) + scale_fill_manual(name = "Color", values = bottomcols, guides(alpha=FALSE))
  } else {
    #Don't
    p1 <- p1 + scale_colour_manual(name = "Color", values = topcols, guides(alpha=FALSE)) + scale_fill_manual(name = "Color", values = topcols, guides(alpha=FALSE))
    p2 <- p2 + scale_colour_manual(name = "Color", values = bottomcols, guides(alpha=FALSE)) + scale_fill_manual(name = "Color", values = bottomcols, guides(alpha=FALSE))
  }
  
  #Highlight if given
  if(!missing(highlight_var)){
    if("Shape" %in% topn){
      p1 <- p1 + geom_point(data=d_order[d_order$Variable %in% highlight_var & d_order$Location=="Top", ], aes(x=pos_index, y=pval, shape=Shape), colour=highlighter)
      p1 <- p1 + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p1 <- p1 + geom_point(data=d_order[d_order$Variable %in% highlight_var & d_order$Location=="Top", ], aes(x=pos_index, y=pval), colour=highlighter)
    }
    if("Shape" %in% bottomn){
      p2 <- p2 + geom_point(data=d_order[d_order$Variable %in% highlight_var & d_order$Location=="Bottom", ], aes(x=pos_index, y=pval, shape=Shape), colour=highlighter)
      p2 <- p2 + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p2 <- p2 + geom_point(data=d_order[d_order$Variable %in% highlight_var & d_order$Location=="Bottom", ], aes(x=pos_index, y=pval), colour=highlighter)
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
      p1 <- p1 + geom_text(data=d_order[d_order$pvalue < annotate_p & d_order$Location=="Top",], aes(pos_index,pval,label=Variable))
      p2 <- p2 + geom_text(data=d_order[d_order$pvalue < annotate_p & d_order$Location=="Bottom",], aes(pos_index,pval,label=Variable))
    } else {
      require("ggrepel", quietly = TRUE)
      p1 <- p1 + geom_text_repel(data=d_order[d_order$pvalue < annotate_p & d_order$Location=="Top",], aes(pos_index,pval,label=Variable))
      p2 <- p2 + geom_text_repel(data=d_order[d_order$pvalue < annotate_p & d_order$Location=="Bottom",], aes(pos_index,pval,label=Variable))
    }
  }
  if(!missing(annotate_var)){
    if (!requireNamespace(c("ggrepel"), quietly = TRUE)==TRUE){
      print("Consider installing 'ggrepel' for improved text annotation")
      p1 <- p1 + geom_text(data=d_order[d_order$Variable %in% annotate_var & d_order$Location=="Top",], aes(pos_index,pval,label=Variable))
      p2 <- p2 + geom_text(data=d_order[d_order$Variable %in% annotate_var & d_order$Location=="Bottom",], aes(pos_index,pval,label=Variable))
    } else {
      require("ggrepel", quietly = TRUE)
      p1 <- p1 + geom_text_repel(data=d_order[d_order$Variable %in% annotate_var & d_order$Location=="Top",], aes(pos_index,pval,label=Variable))
      p2 <- p2 + geom_text_repel(data=d_order[d_order$Variable %in% annotate_var & d_order$Location=="Bottom",], aes(pos_index,pval,label=Variable))
    }
  }
  #Add title and y axis title
  p1 <- p1 + ylab(yaxislab)
  p2 <- p2 + ylab(yaxislab)

  #Format
  p1 <- p1+theme(axis.text.x = element_text(vjust=1),axis.ticks.x = element_blank())+ylim(c(0,yaxismax))
  p2 <- p2+scale_y_reverse(limits=c(yaxismax,0)) + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
  if(rotatelabels==TRUE){p1 <- p1 + theme(axis.text.x = element_text(angle=labelangle))}

  #Save
  print(paste("Saving plot to ", file, ".png", sep=""))
  p <- grid.arrange(arrangeGrob(p1, top=toptitle), arrangeGrob(p2, bottom=bottomtitle), padding=0, heights=c(hgtratio, 1-hgtratio))
  ggsave(p, filename=paste(file, ".png", sep=""), dpi=res, units="in", height=hgt, width=wi)
  return(p)

}

