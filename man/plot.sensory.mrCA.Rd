% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.sensory.mrCA.R
\name{plot.sensory.mrCA}
\alias{plot.sensory.mrCA}
\title{Plot factor plan resulting from multiple-response Correspondence Analysis (MR-CA) applied on sensory data}
\usage{
\method{plot}{sensory.mrCA}(
  x,
  axes = c(1, 2),
  alpha.total.bootstrap.test = 0.05,
  alpha.ellipse = alpha.total.bootstrap.test,
  select.desc = rownames(x$desc.coord),
  rev.x = FALSE,
  rev.y = FALSE,
  size.points = 3.5,
  size.lab = 6,
  size.head.arrow = 0.4,
  expansion = 1.25,
  title = NULL,
  ...
)
}
\arguments{
\item{x}{A list returned by \code{\link[MultiResponseR]{sensory.mrCA}}}

\item{axes}{Which dimensions of the MR-CA should be plotted?}

\item{alpha.total.bootstrap.test}{The alpha risk of the total bootstrap tests. See details}

\item{alpha.ellipse}{The alpha risk of the confidence ellipses}

\item{select.desc}{A character vector specifying the descriptors to plot. By default, all descriptors are plotted}

\item{rev.x}{Should the horizontal plotted dimension be reversed? Useful in case of map comparisons to align products}

\item{rev.y}{Should the vertical plotted dimension be reversed? Useful in case of map comparisons to align products}

\item{size.points}{The size of the points used to represent the products on the map}

\item{size.lab}{The size of the label on the map}

\item{size.head.arrow}{The size of the head of the arrows used to represent the descriptors on the map}

\item{expansion}{The factor of expansion applied to descriptors coordinates to increase readability}

\item{title}{An optional title to be added to the plot}

\item{...}{further arguments passed to or from other methods}
}
\value{
A MR-CA factor map
}
\description{
This function plots the results coming from \code{\link[MultiResponseR]{sensory.mrCA}}
}
\details{
\itemize{
\item \strong{alpha.total.bootstrap.test}: products non-significantly different at the alpha risk of \emph{alpha.total.bootstrap.test} according to the total bootstrap test are linked by a line on the plot. If these links are not required, \emph{alpha.total.bootstrap.test} can be set to 1
}
}
\examples{
data(milkchoc)

dim.sig=sensory.mr.dimensionality.test(milkchoc)$dim.sig

res=sensory.mrCA(milkchoc,nbaxes.sig=dim.sig)

plot(res)
}
