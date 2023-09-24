#' Overall analysis of multiple-response sensory data using the multiple-response chi-square framework introduced in Mahieu, Schlich, Visalli, and Cardot (2021)
#'
#' @description Successively performs \code{\link[MultiResponseR]{sensory.mr.dimensionality.test}}, \code{\link[MultiResponseR]{sensory.mrCA}} and \code{\link[MultiResponseR]{sensory.mr.sig.cell}}
#'
#' @param data A data.frame of evaluations in rows whose first two columns are factors (subject and product) and subsequent columns are binary numeric or integer, each column being a descriptor
#' @param nMC Number of Monte-Carlo simulations to consider at each step of the overall analysis
#' @param alpha The alpha risk to consider at each step of the overall analysis
#' @param cell.two.sided Logical. Should the multiple-response hypergeometric tests per cell be two-sided or not?
#'
#'
#' @return The first MR-CA factor map and the percent.derived.cont table with significant cells highlighted
#' @export
#'
#'
#' @references Mahieu, B., Schlich, P., Visalli, M., & Cardot, H. (2021). A multiple-response chi-square framework for the analysis of Free-Comment and Check-All-That-Apply data. Food Quality and Preference, 93.
#'
#'
#' @examples
#'data(milkchoc)
#'
#'sensory.overall.analysis(milkchoc)
sensory.overall.analysis=function(data,nMC=2000,alpha=0.05,cell.two.sided=TRUE){
  res.dim=sensory.mr.dimensionality.test(data,nperm=nMC,alpha=alpha)
  dim.sig=res.dim$dim.sig
  res.ca=sensory.mrCA(data,nboot=nMC,nbaxes.sig=dim.sig)
  p=plot(res.ca,alpha.total.bootstrap.test=alpha,alpha.ellipse=alpha)
  print(p)
  res.cell=sensory.mr.sig.cell(data,nsample=nMC,nbaxes.sig=dim.sig,two.sided=cell.two.sided)
  g=plot(res.cell,alpha.1=alpha,alpha.2=0)
  print(g)
}
