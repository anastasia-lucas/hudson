#' iphemirror
#'
#' Create mirrored Manhattan plots for PheWAS
#' Dependencies: ggplot2, gridExtra
#' Suggested: ggrepel
#' @param top data frame, must contain: PHE, SNP, CHR, POS, pvalue, columns; optional: Shape, Hover, and Link
#' @param bottom data frame, must contain: PHE, SNP, CHR, POS, pvalue, columns; optional: Shape, Hover, and Link
#' @param phegroup optional grouping dataframe for phenotypes, must contain PHE and Group columns
#' @param tline list of pvalues to draw red threshold lines in top plot
#' @param bline list of pvalues to draw red threshold lines in bottom plot
#' @param chroms list of chromosomes to plot in the order desired, default c(1:22, "X", "Y")
#' @param log10 plot -log10() of pvalue column, logical
#' @param yaxis label for y-axis in the format c("top", "bottom"), automatically set if log10=TRUE
#' @param opacity opacity of points, from 0-1, useful for dense plots
#' @param annotate_snp vector of RSIDs to annotate
#' @param annotate_p list of pvalue thresholds to annotate in the order of c(p_top, p_bottom)
#' @param toptitle optional string for top plot title
#' @param bottomtitle optional string for bottom plot title
#' @param chrcolor1 first alternating color for chromosome
#' @param chrcolor2 second alternating color for chromosome
#' @param highlight_snp vector of snps to highlight
#' @param highlight_p list of pvalue thresholds to annotate in the order of c(p_top, p_bottom)
#' @param highlighter color to highlight
#' @param groupcolors named vector of colors where names correspond to data in 'PHE' or 'Group' column
#' @param freey allow y-axes to scale with the data
#' @param background variegated or white
#' @param chrblocks logical, turns on x-axis chromosome marker blocks
#' @param file file name of saved image
#' @param hgt height of plot in inches
#' @param hgtratio height ratio of plots, equal to top plot proportion
#' @param wi width of plot in inches
#' @return html file
#' @import ggplot2
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @export
#' @examples
#' data(phewas.t)
#' data(phewas.b)
#' phewas.t$Hover <- paste("p:", formatC(phewas.t$pvalue, format="e", digits=2))
#' phewas.b$Hover <- paste("p:", formatC(phewas.b$pvalue, format="e", digits=2))
#' phewas.t$Link <- paste0("https://www.ncbi.nlm.nih.gov/snp/", phewas.t$SNP)
#' phewas.b$Link <- paste0("https://www.ncbi.nlm.nih.gov/snp/", phewas.b$SNP)
#' iphemirror(top=phewas.t, bottom = phewas.b, 
#'            toptitle = "PheWAS Example: Data 1", bottomtitle = "PheWAS Example: Data 2")

iphemirror <- function(top, bottom, phegroup, tline, bline, chroms = c(1:22,"X","Y"),
                       log10=TRUE, yaxis, opacity=1, annotate_snp, annotate_p, highlight_snp, 
                       highlight_p, highlighter="red", toptitle=NULL, bottomtitle=NULL, 
                       chrcolor1="#AAAAAA", chrcolor2="#4D4D4D", groupcolors, freey=FALSE, 
                       background="variegated", chrblocks=TRUE, file="phemirror",
                       hgtratio=0.5, hgt=7, wi=12){
  
  #Sort data
  topn <- names(top)
  bottomn <- names(bottom)
  top$Location <- "Top"
  bottom$Location <- "Bottom"
  
  #Check file formats
  if(!identical(topn, bottomn)){stop("Please ensure both inputs have the same metadata columns.")}
  
  d <- as.data.frame(rbind(top, bottom))
  
  #Set onclick to NULL if needed
  if(!("Hover" %in% names(d))){
    d$Hover <- paste0("SNP: ", d$SNP,
                      "\nPosition: ", paste(d$CHR, d$POS, sep=":"),
                      "\nPhenotype: ", d$PHE,
                      "\np-value: ", formatC(d$pvalue, format="e", digits=2))
  }
  if(!("Link" %in% names(d))){
    d$Link <- NA
  } else {
    d$Link <- sprintf("window.open(\"%s\")", d$Link)
  }
  
  d$POS <- as.numeric(as.character(d$POS))
  d$CHR <- droplevels(factor(d$CHR, levels = as.character(chroms)))
  d <- d[d$CHR %in% chroms, ]
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
  #Chromsome & panel
  nchrcolors <- nlevels(factor(lims$Color))
  base_color <- c(rep(x=c(chrcolor1, chrcolor2), length.out=nchrcolors, each=1), "#FFFFFF", "#EBEBEB")
  names(base_color) <- c(levels(factor(lims$Color)), "shade_ffffff", "shade_ebebeb")
  
  if(!missing(groupcolors)){
    #dcols <- c(rep(x=c(chrcolor1, chrcolor2), length.out=nchrcolors, each=1), "#FFFFFF", "#EBEBEB")
    #names(dcols) <-c(levels(factor(lims$Color)), "shade_ffffff", "shade_ebebeb")
    topcols <- groupcolors
    bottomcols <- groupcolors
  } else {
    #Top Colors
    ngroupcolors <- nlevels(factor(d_order$Color[d_order$Location=="Top"]))
    if(ngroupcolors>15){
      topcols <- hudson:::Turbo(out.colors = ngroupcolors)
    } else {
      pal <- c("#009292", "#920000", "#490092", "#db6d00", "#24ff24", 
               "#ffff6d", "#000000", "#006ddb", "#004949","#924900", 
               "#ff6db6", "#6db6ff","#b66dff", "#ffb6db","#b6dbff")
      topcols <- pal[1:ngroupcolors]
    }
    names(topcols) <- levels(factor(d_order$Color[d_order$Location=="Top"]))
    #Bottom Colors
    ngroupcolors <- nlevels(factor(d_order$Color[d_order$Location=="Bottom"]))
    if(ngroupcolors>15){
      bottomcols <- hudson:::Turbo(ngroupcolors)
    } else {
      pal <- c("#009292", "#920000", "#490092", "#db6d00", "#24ff24", 
               "#ffff6d", "#000000", "#006ddb", "#004949","#924900", 
               "#ff6db6", "#6db6ff","#b66dff", "#ffb6db","#b6dbff")
      bottomcols <- pal[1:ngroupcolors]
    }
    
    names(bottomcols) <- levels(factor(d_order$Color[d_order$Location=="Bottom"]))
  }
  
  #Allow more than 6 shapes
  #3, 4 and 7 to 14 are composite symbols- incompatible with ggiraph
  if("Shape" %in% names(d)){
    allshapes <- c(16,15,17,18,0:2,5:6,19:25,33:127)
    shapevector <- allshapes[1:nlevels(as.factor(d$Shape))]
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
  
  yaxismax1 <- ifelse(freey==FALSE, max(d_order$pval[which(d_order$pval< Inf)]), max(d_order$pval[which(d_order$pval< Inf) & d_order$Location=="Top"]))
  yaxismax2 <- ifelse(freey==FALSE, max(d_order$pval[which(d_order$pval< Inf)]), max(d_order$pval[which(d_order$pval< Inf) & d_order$Location=="Bottom"]))
  yaxismin1 <- ifelse(freey==FALSE, 0, min(d_order$pval[d_order$Location=="Top"]))
  yaxismin2 <- ifelse(freey==FALSE, 0, min(d_order$pval[d_order$Location=="Bottom"]))
  
  #Theme options
  backpanel1 <- ifelse(background=="white", "NULL", "geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = yaxismin1, ymax = Inf, fill=factor(shademap)), alpha = 0.5)" )
  backpanel2 <- ifelse(background=="white", "NULL", "geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = yaxismin2, ymax = Inf, fill=factor(shademap)), alpha = 0.5)" )
  
  #Start plotting
  #TOP PLOT
  p1 <- ggplot() + eval(parse(text=backpanel1))

  #Add shape info if available
  if("Shape" %in% names(d)){
    
    p1 <- p1 + geom_point(data=d_order[is.na(d_order$Hover) & d_order$Location=="Top",],
                        aes(x=pos_index, y=pval,
                            color=factor(Color), shape=factor(Shape)),
                        alpha=opacity) +
      ggiraph::geom_point_interactive(data=d_order[!is.na(d_order$Hover) & d_order$Location=="Top",],
                                      aes(x=pos_index, y=pval, tooltip=Hover,
                                          color=factor(Color), shape=factor(Shape),
                                          onclick=Link),
                                      alpha=opacity) +
      scale_shape_manual(values=shapevector)
  } else {
    p1 <- p1 + geom_point(data=d_order[is.na(d_order$Hover) & d_order$Location=="Top",],
                        aes(x=pos_index, y=pval,
                            color=factor(Color)),
                        alpha=opacity) +
      ggiraph::geom_point_interactive(data=d_order[!is.na(d_order$Hover) & d_order$Location=="Top",],
                                      aes(x=pos_index, y=pval,
                                          color=factor(Color), tooltip=Hover,
                                          onclick=Link),
                                      alpha=opacity)
  }
  
  p1 <- p1 + scale_x_continuous(breaks=lims$av, labels=lims$Color, expand=c(0,0))
  if(chrblocks==TRUE){
    p1 <- p1 + geom_rect(data = lims, 
                         aes(xmin = posmin-.5, xmax = posmax+.5, ymin = -Inf, ymax = min(d_order$pval), fill=as.factor(Color)), 
                         alpha = 1)
  }
  #Add legend
  p1 <- p1 + 
    scale_colour_manual(name = "Color", values = topcols) + 
    scale_fill_manual(values = base_color)
  p1 <- p1 + theme(panel.grid.minor.x = element_blank(),
                   panel.grid.major.x=element_blank(), 
                   axis.title.x=element_blank(), 
                   legend.position="top", 
                   legend.title=element_blank())
  
  #Start plotting
  #BOTTOM PLOT
  p2 <- ggplot() + eval(parse(text=backpanel2))
  #Add shape info if available
  if("Shape" %in% names(d)){
    
    p2 <- p2 + geom_point(data=d_order[is.na(d_order$Hover) & d_order$Location=="Bottom",],
                          aes(x=pos_index, y=pval,
                              color=factor(Color), shape=factor(Shape)),
                          alpha=opacity) +
      ggiraph::geom_point_interactive(data=d_order[!is.na(d_order$Hover) & d_order$Location=="Bottom",],
                                      aes(x=pos_index, y=pval, tooltip=Hover,
                                          color=factor(Color), shape=factor(Shape),
                                          onclick=Link),
                                      alpha=opacity) +
      scale_shape_manual(values=shapevector)
  } else {
    p2 <- p2 + geom_point(data=d_order[is.na(d_order$Hover) & d_order$Location=="Bottom",],
                          aes(x=pos_index, y=pval,
                              color=factor(Color)),
                          alpha=opacity) +
      ggiraph::geom_point_interactive(data=d_order[!is.na(d_order$Hover) & d_order$Location=="Bottom",],
                                      aes(x=pos_index, y=pval,
                                          color=factor(Color), tooltip=Hover,
                                          onclick=Link),
                                      alpha=opacity)
  }
  
  p2 <- p2 + scale_x_continuous(breaks=lims$av, labels=lims$Color, expand=c(0,0))
  if(chrblocks==TRUE){
    p2 <- p2 + geom_rect(data = lims, 
                         aes(xmin = posmin-.5, xmax = posmax+.5, ymin = -Inf, ymax = min(d_order$pval), fill=as.factor(Color)), 
                         alpha = 1)
  } 
  #Add legend
  p2 <- p2 + scale_colour_manual(name = "Color", values = bottomcols) + scale_fill_manual(values = base_color)
  p2 <- p2 + theme(axis.text.x=element_text(angle=90), panel.grid.minor.x = element_blank(), panel.grid.major.x=element_blank(), axis.title.x=element_blank(), legend.position="bottom", legend.title=element_blank())
  
  #Highlight if given
  if(!missing(highlight_snp)){
    if("Shape" %in% topn){
      p1 <- p1 + ggiraph::geom_point_interactive(data=d_order[d_order$SNP %in% highlight_snp & d_order$Location=="Top", ], 
                                                 aes(x=pos_index, y=pval, shape=Shape, tooltip=Hover, onclick=Link), 
                                                 colour=highlighter) + 
                scale_shape_manual(values=shapevector)
      p1 <- p1 + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p1 <- p1 + ggiraph::geom_point_interactive(data=d_order[d_order$SNP %in% highlight_snp & d_order$Location=="Top", ], 
                                                 aes(x=pos_index, y=pval, tooltip=Hover, onclick=Link), 
                                                 colour=highlighter)
    }
    
    if("Shape" %in% bottomn){
      p2 <- p2 + ggiraph::geom_point_interactive(data=d_order[d_order$SNP %in% highlight_snp & d_order$Location=="Bottom", ], 
                                                 aes(x=pos_index, y=pval, shape=Shape, tooltip=Hover, onclick=Link), 
                                                 colour=highlighter) + 
                 scale_shape_manual(values=shapevector)
      p2 <- p2 + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p2 <- p2 + ggiraph::geom_point_interactive(data=d_order[d_order$SNP %in% highlight_snp & d_order$Location=="Bottom", ], 
                                                 aes(x=pos_index, y=pval, tooltip=Hover, onclick=Link), 
                                                 colour=highlighter)
    }
    
  }
  if(!missing(highlight_p)){
    if("Shape" %in% topn){
      p1 <- p1 + ggiraph::geom_point_interactive(data=d_order[d_order$pvalue < highlight_p & d_order$Location=="Top", ], 
                                                 aes(x=pos_index, y=pval, shape=Shape, tooltip=Hover, onclick=Link), 
                                                 colour=highlighter) + 
                 scale_shape_manual(values=shapevector)
      p1 <- p1 + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p1 <- p1 + ggiraph::geom_point_interactive(data=d_order[d_order$pvalue < highlight_p & d_order$Location=="Top", ], 
                                                 aes(x=pos_index, y=pval, tooltip=Hover, onclick=Link), 
                                                 colour=highlighter)
    }
    
    if("Shape" %in% bottomn){
      p2 <- p2 + ggiraph::geom_point_interactive(data=d_order[d_order$pvalue < highlight_p & d_order$Location=="Bottom", ], 
                                                 aes(x=pos_index, y=pval, shape=Shape, tooltip=Hover, onclick=Link), 
                                                 colour=highlighter) + 
                 scale_shape_manual(values=shapevector)
      p2 <- p2 + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p2 <- p2 + ggiraph::geom_point_interactive(data=d_order[d_order$pvalue < highlight_p & d_order$Location=="Bottom", ], 
                                                 aes(x=pos_index, y=pval, tooltip=Hover, onclick=Link), 
                                                 colour=highlighter)
    }
  }
  
  #Add pvalue threshold line
  if(!missing(tline)){
    for(i in 1:length(tline)){
      p1 <- p1 + geom_hline(yintercept = tredline[i], colour="red")
    }
  }
  if(!missing(bline)){
    for(i in 1:length(bline)){
      p2 <- p2 + geom_hline(yintercept = bredline[i], colour="red")
    }
  }
  #Annotate
  if(!missing(annotate_p)){
    if (!requireNamespace(c("ggrepel"), quietly = TRUE)==TRUE) {
      print("Consider installing 'ggrepel' for improved text annotation")
      p1 <- p1 + geom_text(data=d_order[d_order$pvalue < annotate_p[1] & d_order$Location=="Top",], aes(pos_index,pval,label=SNP))
      p2 <- p2 + geom_text(data=d_order[d_order$pvalue < annotate_p[2] & d_order$Location=="Bottom",], aes(pos_index,pval,label=SNP))
    } else {
      p1 <- p1 + ggrepel::geom_text_repel(data=d_order[d_order$pvalue < annotate_p[1] & d_order$Location=="Top",], aes(pos_index,pval,label=SNP))
      p2 <- p2 + ggrepel::geom_text_repel(data=d_order[d_order$pvalue < annotate_p[2] & d_order$Location=="Bottom",], aes(pos_index,pval,label=SNP))
    }
  }
  if(!missing(annotate_snp)){
    if (!requireNamespace(c("ggrepel"), quietly = TRUE)==TRUE){
      print("Consider installing 'ggrepel' for improved text annotation")
      p1 <- p1 + geom_text(data=d_order[d_order$SNP %in% annotate_snp & d_order$Location=="Top",], aes(pos_index,pval,label=SNP))
      p2 <- p2 + geom_text(data=d_order[d_order$SNP %in% annotate_snp & d_order$Location=="Bottom",], aes(pos_index,pval,label=SNP))
    } else {
      p1 <- p1 + ggrepel::geom_text_repel(data=d_order[d_order$SNP %in% annotate_snp & d_order$Location=="Top",], aes(pos_index,pval,label=SNP))
      p2 <- p2 + ggrepel::geom_text_repel(data=d_order[d_order$SNP %in% annotate_snp & d_order$Location=="Bottom",], aes(pos_index,pval,label=SNP))
    }
  }
  #Add title and y axis title
  p1 <- p1 + ylab(yaxislab1)
  p2 <- p2 + ylab(yaxislab2)
  
  if(chrblocks==TRUE){
    if(freey==TRUE){
      print("Sorry, drawing chrblocks with freey=TRUE is currently unsupported and will be ignored.")
    } else {
      p1 <- p1+theme(axis.text.x = element_text(vjust=1),axis.ticks.x = element_blank()) + ylim(c(yaxismin1,yaxismax1))
      p2 <- p2+scale_y_reverse(limits=c(yaxismax2, yaxismin2)) + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
    }
  } else {
    p1 <- p1+theme(axis.text.x = element_text(vjust=1),axis.ticks.x = element_blank())+ scale_y_continuous(limits=c(yaxismin1, yaxismax1),expand=expansion(mult=c(0,0.1)))
    p2 <- p2+scale_y_reverse(limits=c(yaxismax2,yaxismin2), expand=expansion(mult=c(0.1,0))) + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
  }
  
  if(background=="white"){
    p1 <- p1 + theme(panel.background = element_rect(fill="white"))
    p2 <- p2 + theme(panel.background = element_rect(fill="white"))
  }
  
  p1 <- p1 + guides(fill="none")
  p2 <- p2 + guides(fill="none")
  
  #Save
  print(paste0("Saving plot to ", file, ".html"))
  plot_row <- cowplot::plot_grid(p1, p2, ncol=1, rel_heights=c(hgtratio, 1-hgtratio))
  rel <- 1
  if(!is.null(toptitle)){  toptitle <- cowplot::ggdraw() + cowplot::draw_label(toptitle) ; rel <- c(0.05, rel) }
  if(!is.null(bottomtitle)){ bottomtitle <- cowplot::ggdraw() + cowplot::draw_label(bottomtitle) ; rel <- c(rel, 0.05) }
  plot_row <- list(toptitle, plot_row, bottomtitle)
  
  p <- ggiraph::girafe( ggobj = cowplot::plot_grid(plotlist = plot_row[!sapply(plot_row,is.null)],
                                                   ncol = 1, 
                                                   rel_heights = rel), 
                        width_svg = wi, height_svg = hgt)
  htmlwidgets::saveWidget(widget=p, file=paste0(file, ".html"))
  return("Finished")
  
}
