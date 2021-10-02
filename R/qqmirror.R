#' qqmirror
#'
#' Create qqplots with an assumed uniform distribution
#' @param top dataframe with at least one column, p-value; Color, Shape, and Name optional
#' @param bottom dataframe with at least one column, p-value; Color, Shape, and Name optional
#' @param CI two-sided confidence interval, default 0.95
#' @param splittop if data contains Color and/or Shape, indicate variable(s) by which the data should be subsetted for calculating CIs
#' @param splitbottom if data contains Color and/or Shape, indicate variable(s) by which the data should be subsetted for calculating CIs
#' @param opacity point opacity, default 1
#' @param toptitle title for top plot
#' @param bottomtitle title for bottom plot
#' @param groupcolors named vector of colors corresponding to data in Group column
#' @param highlight_name vector of names to highlight, dataframe must include a Name column
#' @param highlight_p list of pvalue thresholds to highlight in the order of c(p_top, p_bottom)
#' @param highlighter highlighter color
#' @param freey allow y-axes to scale with data
#' @param annotate_name vector of names to annotate, dataframe must include a Name column to use this option
#' @param annotate_p list of pvalue thresholds to annotate in the order of c(p_top, p_bottom)
#' @param tline draw a red line at (observed) pvalue threshold in top plot
#' @param bline draw a red line at (observed) pvalue threshold in bottom plot 
#' @param background can change to "white"
#' @param file filename
#' @param wi width of plot
#' @param hgt height of plot
#' @param hgtratio height ratio for plots
#' @param res resolution of plot
#' @return png image
#' @import ggplot2
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @export
#' @examples
#' data(gwas.t)
#' qqdat <- data.frame(pvalue=gwas.t$pvalue, Color=rep(c("Additive", "Dominant", "Recessive"), 
#' length.out=nrow(gwas.t), each=1))
#' qqmirror(top=qqdat, bottom=qqdat, opacity=0.6, splittop="Color", splitbottom="Color", 
#' toptitle="Example: Top", bottomtitle="Example: Bottom")

qqmirror <- function(top, bottom, CI=0.95, opacity=1, groupcolors, splittop=NULL, splitbottom=NULL, 
                     highlight_p, highlight_name, annotate_p, annotate_name, highlighter="red", 
                     freey=FALSE, tline, bline, background, toptitle=NULL, bottomtitle=NULL, 
                     file="qqmirror", wi=6, hgt=8, hgtratio=0.5, res=300){

  arglist <- list(d=top, CI=CI, opacity=opacity, splitby=splittop, highlighter=highlighter)
  if(!missing(tline)){arglist <- c(arglist, line=tline)}
  if(!missing(groupcolors)){arglist <- c(arglist, groupcolors=groupcolors)}
  if(!missing(highlight_p)){arglist <- c(arglist, highlight_p=highlight_p[2])}
  if(!missing(highlight_name)){arglist <- c(arglist, highlight_name=highlight_name)}
  if(!missing(annotate_p)){arglist <- c(arglist, annotate_p=annotate_p[1])}
  if(!missing(annotate_name)){arglist <- c(arglist, annotate_name=annotate_name)}
  if(!missing(background)){arglist <- c(arglist, background=background)}
  p1 <- do.call(qqunif, arglist)
  p1 <- p1 + theme(legend.position="top")
  p1 <- p1 + scale_x_continuous(position="top")
  
  arglist <- list(d=bottom, CI=CI, opacity=opacity, splitby=splitbottom, highlighter=highlighter, slope=-1)
  if(!missing(bline)){arglist <- c(arglist, line=bline)}
  if(!missing(groupcolors)){arglist <- c(arglist, groupcolors=groupcolors)}
  if(!missing(highlight_p)){arglist <- c(arglist, highlight_p=highlight_p[2])}
  if(!missing(highlight_name)){arglist <- c(arglist, highlight_name=highlight_name)}
  if(!missing(annotate_p)){arglist <- c(arglist, annotate_p=annotate_p[2])}
  if(!missing(annotate_name)){arglist <- c(arglist, annotate_name=annotate_name)}
  if(!missing(background)){arglist <- c(arglist, background=background)}
  p2 <- do.call(qqunif, arglist)
  
  if(freey==FALSE){
    p1 <- p1 + scale_y_continuous(limits = c(0, max(max(-log10(top$pvalue)), max(-log10(bottom$pvalue)))))
    yax <- max(max(-log10(top$pvalue)), max(-log10(bottom$pvalue)))
    p2 <- p2 + scale_y_reverse(limits=c(yax,0))
  } else {
    p2 <- p2 + scale_y_reverse()
    
  }
  
  #p1 <- p1 + theme_minimal() + theme(panel.background = element_rect(fill = NULL))
  #p2 <- p2 + theme_minimal() + theme(panel.background = element_rect(fill = NULL))
  #Save
  print(paste("Saving plot to ", file, ".png", sep=""))
  p <- grid.arrange(arrangeGrob(p1, top=toptitle), arrangeGrob(p2, bottom=bottomtitle), padding=0, heights=c(hgtratio,1-hgtratio))
  ggsave(p, filename=paste(file, ".png", sep=""), dpi=res, units="in", height=hgt, width=wi)
  return(p)

  
}
