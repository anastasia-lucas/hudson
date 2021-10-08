#' iemirror
#'
#' Create interactive mirrored Manhattan plots for EWAS
#' Dependencies: ggplot2, ggiraph
#' Suggested: ggrepel
#' @param top data frame, columns one through three must be Variable, pvalue, and Group; Shape, Color, Hover, Link optional
#'                        where 'Hover' contains tooltip information overriding the default (Variable, Group)
#'                        and 'Link' contains a search query, ex. "https://en.wikipedia.org/wiki/"
#' @param bottom data frame, columns one and two must be Variable, pvalue, and Group; Shape, Color, Hover, Link optional
#'                           where 'Hover' contains tooltip information overriding the default (Variable, Group)
#'                           and 'Link' contains a search query, ex. "https://en.wikipedia.org/wiki/"
#' @param tline list of pvalues to draw red threshold lines in top plot
#' @param bline list of pvalues to draw red threshold lines in bottom plot
#' @param log10 plot -log10() of pvalue column, logical
#' @param yaxis title for y-axis, automatically set if log10=TRUE
#' @param xaxis title for x-axis
#' @param opacity opacity of points, from 0 to 1, useful for dense plots
#' @param toptitle optional string for plot title
#' @param bottomtitle optional string for plot title
#' @param color1 first alternating color
#' @param color2 second alternating color
#' @param highlight_var vector of variables to highlight
#' @param highlight_p list of pvalue thresholds to highlight in the order of c(p_top, p_bottom)
#' @param highlighter color to highlight
#' @param groupcolors named vector of colors where names correspond to data in 'Color' column
#' @param rotatelabels logical, rotate axis labels?
#' @param labelangle angle to rotate
#' @param freey allow y-axes to scale with data
#' @param background variegated or white
#' @param grpblocks logical, turns on x-axis group marker blocks
#' @param file file name of saved image
#' @param hgt height of plot in inches
#' @param wi width of plot in inches
#' @param res resolution of plot in pixels per inch
#' @return png image
#' @import ggplot2
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @export
#' @examples
#' data(ewas.t)
#' data(ewas.b)
#' ewas.t$Link <- paste0("https://www.google.com/search?q=", ewas.t$Variable)
#' ewas.b$Link <- paste0("https://www.google.com/search?q=", ewas.b$Variable)
#' iemirror(top=ewas.t, bottom=ewas.b)


iemirror <- function(top, bottom, tline, bline, log10=TRUE, yaxis=NULL, xaxis=NULL,
                     opacity=1, highlight_var,  highlight_p, highlighter="red", 
                     title=NULL, color1="#AAAAAA", color2="#4D4D4D", blockmin=-1,
                     blockmax=1, groupcolors, rotatelabels=FALSE, labelangle, freey=FALSE, 
                     ymax, ymin, background="variegated", grpblocks=FALSE, file="iemirror", 
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
    d$Hover <- paste0("Name: ", d$Variable,
                     "\nGroup: ", d$Group)
  }
  if(!("Link" %in% names(d))){
    d$Link <- NA
  } else {
    d$Link <- sprintf("window.open(\"%s\")", d$Link)
  }
  
  #Save to merge later
  d$rowid <- seq.int(nrow(d))
  dinfo <- d[, colnames(d) %in% c("rowid", "Color", "pval", "Location", "Shape", "Hover", "Link"), drop=FALSE]
  
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
  base_color <- c(rep(x=c(color1, color2), length.out=nvarcolors, each=1), "#FFFFFF", "#EBEBEB")
  names(base_color) <- c(levels(factor(lims$Color)), "shade_ffffff", "shade_ebebeb")
  
  if("Color" %in% names(d)){
    #if(missing(levs)){
    levs=as.character(levels(factor(d_order$Color)))
    #}
    ngroupcolors <- nlevels(factor(d_order$Color, levels=levs))
    newcols <- hudson:::Turbo(out.colors=ngroupcolors)
    names(newcols) <- levels(factor(d_order$Color, levels=levs))
  } else {
    #Color by Group instead
    names(d_order)[names(d_order)=="Group"] <- "Color"
    newcols <- rep(x=c(color1, color2), length.out=nvarcolors, each=1)
    names(newcols) <-levels(factor(lims$Color))
  }
  
  #Allow more than 6 shapes
  #3, 4 and 7 to 14 are composite symbols- incompatible with ggiraph
  if("Shape" %in% names(d)){
    allshapes <- c(16,15,17,18,0:2,5:6,19:25,33:127)
    shapevector <- allshapes[1:nlevels(as.factor(d$Shape))]
  }
  
  #Start plotting
  d_order <- merge(d_order, dinfo, by="rowid")
  
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
                                      alpha=opacity)
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
  
  if(grpblocks==TRUE){p <- p + geom_rect(data = lims, aes(xmin = posmin, xmax = posmax, ymin = blockmin, ymax = blockmax, fill=as.factor(Color)), alpha = 1)}
 
   #Add legend
  p <- p + 
        scale_colour_manual(name = "Color", values = newcols) + 
        scale_fill_manual(values = base_color) +
        theme(panel.grid.minor.y = element_blank(),
              panel.grid.major.y=element_blank(),
              legend.title=element_blank(),
              legend.position = "bottom")
  
  #Highlight if given
  if(!missing(highlight_var)){
    if("Shape" %in% topn){
      p1 <- p1 + 
              ggiraph::geom_point_interactive(data=d_order[d_order$Variable %in% highlight_var & d_order$Location=="Top", ], 
                                              aes(x=pos_index, y=pval, shape=Shape, tooltip=Hover, onclick=Link), 
                                              colour=highlighter)
      p1 <- p1 + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p1 <- p1 + 
              ggiraph::geom_point_interactive(data=d_order[d_order$Variable %in% highlight_var & d_order$Location=="Top", ], 
                                              aes(x=pos_index, y=pval, tooltip=Hover, onclick=Link), 
                                              colour=highlighter)
    }
    if("Shape" %in% bottomn){
      p2 <- p2 + 
              ggiraph::geom_point_interactive(data=d_order[d_order$Variable %in% highlight_var & d_order$Location=="Bottom", ], 
                                              aes(x=pos_index, y=pval, shape=Shape, tooltip=Hover, onclick=Link), 
                                              colour=highlighter)
      p2 <- p2 + 
              guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p2 <- p2 + 
              ggiraph::geom_point_interactive(data=d_order[d_order$Variable %in% highlight_var & d_order$Location=="Bottom", ], 
                                              aes(x=pos_index, y=pval, tooltip=Hover, onclick=Link), 
                                              colour=highlighter)
    }
  }
  if(!missing(highlight_p)){
    if("Shape" %in% topn){
      p1 <- p1 + 
              ggiraph::geom_point_interactive(data=d_order[d_order$pvalue < highlight_p[1] & d_order$Location=="Top", ], 
                                              aes(x=pos_index, y=pval, shape=Shape, tooltip=Hover, onclick=Link), 
                                              colour=highlighter)
      p1 <- p1 + 
              guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p1 <- p1 + 
              ggiraph::geom_point_interactive(data=d_order[d_order$pvalue < highlight_p[1] & d_order$Location=="Top", ], 
                                              aes(x=pos_index, y=pval, tooltip=Hover, onclick=Link), 
                                              colour=highlighter)
    }
    if("Shape" %in% bottomn){
      p2 <- p2 + 
              ggiraph::geom_point_interactive(data=d_order[d_order$pvalue < highlight_p[2] & d_order$Location=="Bottom", ], 
                                              aes(x=pos_index, y=pval, shape=Shape, tooltip=Hover, onclick=Link), 
                                              colour=highlighter)
      p2 <- p2 + 
              guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p2 <- p2 + 
              ggiraph::geom_point_interactive(data=d_order[d_order$pvalue < highlight_p[2] & d_order$Location=="Bottom", ], 
                                              aes(x=pos_index, y=pval, tooltip=Hover, onclick=Link), 
                                              colour=highlighter)
    }
  } 

  #Add title and y axis title
  p <- p + ggtitle(title) + ylab(yaxislab) + xlab(xaxis)
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
  
  if(grpblocks==TRUE){
    p <- p+ylim(c(yaxismin,yaxismax))
  } else {
    p <- p+scale_y_continuous(limits=c(yaxismin, yaxismax), expand=expansion(mult=c(0,0.1)))
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
  print(paste0("Saving plot to ", file, ".html"))

  ip <- ggiraph::girafe(code = print(p))
  htmlwidgets::saveWidget(widget=ip, file=paste0(file, ".html"))
  return("Finished")
  
}