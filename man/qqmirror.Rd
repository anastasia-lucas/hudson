% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qqmirror.R
\name{qqmirror}
\alias{qqmirror}
\title{qqmirror}
\usage{
qqmirror(
  top,
  bottom,
  CI = 0.95,
  opacity = 1,
  groupcolors,
  splittop = NULL,
  splitbottom = NULL,
  highlight_p,
  highlight_name,
  annotate_p,
  annotate_name,
  highlighter = "red",
  freey = FALSE,
  tline,
  bline,
  background,
  toptitle = NULL,
  bottomtitle = NULL,
  file = "qqmirror",
  type = "png",
  wi = 6,
  hgt = 8,
  hgtratio = 0.5,
  res = 300
)
}
\arguments{
\item{top}{dataframe with at least one column, p-value; Color, Shape, and Name optional}

\item{bottom}{dataframe with at least one column, p-value; Color, Shape, and Name optional}

\item{CI}{two-sided confidence interval, default 0.95}

\item{opacity}{point opacity, default 1}

\item{groupcolors}{named vector of colors corresponding to data in Group column}

\item{splittop}{if data contains Color and/or Shape, indicate variable(s) by which the data should be subsetted for calculating CIs}

\item{splitbottom}{if data contains Color and/or Shape, indicate variable(s) by which the data should be subsetted for calculating CIs}

\item{highlight_p}{list of pvalue thresholds to highlight in the order of c(p_top, p_bottom)}

\item{highlight_name}{vector of names to highlight, dataframe must include a Name column}

\item{annotate_p}{list of pvalue thresholds to annotate in the order of c(p_top, p_bottom)}

\item{annotate_name}{vector of names to annotate, dataframe must include a Name column to use this option}

\item{highlighter}{highlighter color}

\item{freey}{allow y-axes to scale with data}

\item{tline}{draw a red line at (observed) pvalue threshold in top plot}

\item{bline}{draw a red line at (observed) pvalue threshold in bottom plot}

\item{background}{can change to "white"}

\item{toptitle}{title for top plot}

\item{bottomtitle}{title for bottom plot}

\item{file}{filename}

\item{type}{plot type/extension}

\item{wi}{width of plot}

\item{hgt}{height of plot}

\item{hgtratio}{height ratio for plots}

\item{res}{resolution of plot}
}
\value{
png image
}
\description{
Create qqplots with an assumed uniform distribution
}
\examples{
data(gwas.t)
qqdat <- data.frame(pvalue=gwas.t$pvalue, Color=rep(c("Additive", "Dominant", "Recessive"), 
length.out=nrow(gwas.t), each=1))
qqmirror(top=qqdat, bottom=qqdat, opacity=0.6, splittop="Color", splitbottom="Color", 
toptitle="Example: Top", bottomtitle="Example: Bottom")
}
