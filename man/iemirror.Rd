% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iemirror.R
\name{iemirror}
\alias{iemirror}
\title{iemirror}
\usage{
iemirror(
  top,
  bottom,
  tline,
  bline,
  log10 = TRUE,
  yaxis,
  opacity = 1,
  toptitle = NULL,
  bottomtitle = NULL,
  annotate_var,
  annotate_p,
  highlight_var,
  highlight_p,
  highlighter = "red",
  color1 = "#AAAAAA",
  color2 = "#4D4D4D",
  groupcolors,
  rotatelabels = FALSE,
  labelangle,
  freey = FALSE,
  background = "variegated",
  grpblocks = FALSE,
  file = "emirror",
  hgtratio = 0.5,
  hgt = 7,
  wi = 12
)
}
\arguments{
\item{top}{data frame, columns one and two must be Variable, pvalue, and Group; Shape, Color, Hover, Link optional}

\item{bottom}{data frame, columns one and two must be Variable, pvalue, and Group; Shape, Color, Hover, Link optional}

\item{tline}{list of pvalues to draw red threshold lines in top plot}

\item{bline}{list of pvalues to draw red threshold lines in bottom plot}

\item{log10}{plot -log10() of pvalue column, logical}

\item{yaxis}{label for y-axis in the format c("top", "bottom"), automatically set if log10=TRUE}

\item{opacity}{opacity of points, from 0 to 1, useful for dense plots}

\item{toptitle}{optional string for plot title}

\item{bottomtitle}{optional string for plot title}

\item{annotate_var}{vector of variables to annotate}

\item{annotate_p}{list of pvalue thresholds to annotate in the order of c(p_top, p_bottom)}

\item{highlight_var}{vector of variables to highlight}

\item{highlight_p}{list of pvalue thresholds to highlight in the order of c(p_top, p_bottom)}

\item{highlighter}{color to highlight}

\item{color1}{first alternating color}

\item{color2}{second alternating color}

\item{groupcolors}{named vector of colors where names correspond to data in 'Color' column}

\item{rotatelabels}{logical, rotate axis labels?}

\item{labelangle}{angle to rotate}

\item{freey}{allow y-axes to scale with data}

\item{background}{variegated or white}

\item{grpblocks}{logical, turns on x-axis group marker blocks}

\item{file}{file name of saved image}

\item{hgtratio}{height ratio of plots, equal to top plot proportion}

\item{hgt}{height of plot in inches}

\item{wi}{width of plot in inches}
}
\value{
png image
}
\description{
Create mirrored Manhattan plots for EWAS
Dependencies: ggplot2, gridExtra
Suggested: ggrepel
}
\examples{
data(ewas.t)
data(ewas.b)
ewas.t$Link <- paste0("https://www.google.com/search?q=", ewas.t$Variable)
ewas.b$Link <- paste0("https://www.google.com/search?q=", ewas.b$Variable)
iemirror(top=ewas.t, bottom=ewas.b, annotate_p = c(0.0001, 0.0005), 
         highlight_p=c(0.0001, 0.0005), highlighter="green", 
         toptitle = "EWAS Comparison Example: Data 1", 
         bottomtitle = "EWAS Comparison Example: Data 2")
}
