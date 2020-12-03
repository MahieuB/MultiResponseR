#' Overall analysis of multiple-response sensory data using the multiple-response chi-square framework introduced in Mahieu, Schlich, Visalli, and Cardot (2020)
#'
#' @description Successively performs \code{\link[MultiResponseR]{sensory.mr.dimensionality.test}}, \code{\link[MultiResponseR]{sensory.mrCA}} and \code{\link[MultiResponseR]{sensory.mr.sig.cell}}
#'
#' @param data A data.frame of evaluations in rows whose first two columns are factors (subject and product) and subsequent columns are binary numeric or integer, each column being a descriptor
#' @param nMC Number of Monte-Carlo simulations to consider at each step of the overall analysis
#' @param alpha The α risk to consider at each step of the overall analysis
#' @param cell.two.sided Logical. Should the multiple-response tests per cell be two-sided or not? By default, the tests are performed with a one-sided greater alternative hypothesis
#' @param ncores Number of cores used in the Monte-Carlo simulations. Default is 2. See details
#'
#' @details
#' \itemize{
#'   \item \strong{ncores}: The more cores are added in the process, the faster the results will be obtained. The number of available cores is accessible using \code{\link[parallel]{detectCores}}. The parallel tasks are closed once the simulations are over.
#' }
#'
#' @return The first mrCA factor map and the percent.derived.cont table with significant cells highlighted
#' @export
#'
#' @import foreach
#' @import parallel
#' @import doParallel
#' @import flextable
#' @import officer
#' @import abind
#' @importFrom candisc candisc
#' @import FactoMineR
#' @import stats
#'
#' @references Mahieu, B., Schlich, P., Visalli, M., & Cardot, H. (2020). A multiple-response chi-square framework for the analysis of Free-Comment and Check-All-That-Apply data. Manuscript submitted for publication.
#'
#'
#' @examples
#'data(milkchoc)
#'
#'parallel::detectCores()
#'
#'sensory.overall.analysis(milkchoc)
sensory.overall.analysis=function(data,nMC=2000,alpha=0.05,cell.two.sided=FALSE,ncores=2){
  res.dim=sensory.mr.dimensionality.test(data,nperm=nMC,alpha=alpha,ncores=ncores)
  dim.sig=res.dim$dim.sig
  res.ca=sensory.mrCA(data,nboot=nMC,nbaxes.sig=dim.sig,ncores=ncores)
  plt.mrCA(res.ca,alpha.total.bootstrap.test=alpha,alpha.ellipse=alpha)
  res.cell=sensory.mr.sig.cell(data,nsample=nMC,nbaxes.sig=dim.sig,two.sided=cell.two.sided,ncores=ncores)
  plt.mr.sig.cell(res.cell,alpha = alpha)
  stopImplicitCluster()
}