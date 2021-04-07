#' Multiple-response dimensionality test for sensory data
#'
#' @description Performs a multiple-response dimensionality test as defined in Mahieu, Schlich, Visalli, and Cardot (2020) using random permutations to estimate the null distribution. This test results of the transposition of the dimensionality test introduced in Mahieu, Visalli, and Schlich (2020) from the usual chi-square framework to the modified chi-square framework introduced in Loughin and Scherer (1998). The difference with \code{\link[MultiResponseR]{mr.dimensionality.test}} is that random permutations are performed within subjects rather than along all evaluations
#'
#' @param data A data.frame of evaluations in rows whose first two columns are factors (subject and product) and subsequent columns are binary numeric or integer, each column being a descriptor
#' @param nperm Number of permuted datasets to estimate the distribution of the statistic under the null hypothesis. See details
#' @param alpha The alpha risk of the test
#' @param ncores Number of cores used to estimate the null distribution. Default is 2. See details
#'
#' @details
#' \itemize{
#'   \item \strong{nperm}: The distribution of the statistic under the null hypothesis of no associations between products and descriptors is estimated using \emph{nperm} datasets generated thanks to random permutations of the response vectors along products within subjects.
#'   \item \strong{ncores}: The more cores are added in the process, the faster the results will be obtained. The number of available cores is accessible using \code{\link[parallel]{detectCores}}. The parallel tasks are closed once the \emph{nperm} datasets are generated.
#' }
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{dim.sig}{The number of significant dimensions}
#'   \item{statistics}{Observed multiple-response chi-square statistic of each dimension}
#'   \item{p.values}{P-value of the test of each dimension adjusted for closed testing procedure}
#' }
#' @references Loughin, T. M., & Scherer, P. N. (1998). Testing for Association in Contingency Tables with Multiple Column Responses. Biometrics, 54(2), 630-637.
#' @references Mahieu, B., Visalli, M., & Schlich, P. (2020). Accounting for the dimensionality of the dependence in analyses of contingency tables obtained with Check-All-That-Apply and Free-Comment. Food Quality and Preference, 83.
#' @references Mahieu, B., Schlich, P., Visalli, M., & Cardot, H. (2020). A multiple-response chi-square framework for the analysis of Free-Comment and Check-All-That-Apply data. Manuscript submitted for publication.
#' @export
#'
#' @import parallel
#' @import doParallel
#' @import foreach
#' @import iterators
#' @import stats
#'
#' @examples
#'data(milkchoc)
#'
#'parallel::detectCores()
#'
#'sensory.mr.dimensionality.test(milkchoc)
sensory.mr.dimensionality.test=function(data,nperm=2000,alpha=0.05,ncores=2){
  classe=class(data)[1]
  if (!classe%in%c("data.frame")){
    stop("data must be a data.frame")
  }
  for (j in 1:2){
    classe=class(data[,j])[1]
    if (!classe%in%c("factor")){
      stop("The first two columns of data must be factor")
    }
  }
  if(!colnames(data)[1]%in%c("sujet","subject","Sujet","Subject")){
    stop("First column name must be sujet, Sujet, subject or Subject")
  }
  if(!colnames(data)[2]%in%c("produit","product","Produit","Product")){
    stop("Second column name must be produit, product, Produit or Product")
  }
  colnames(data)[1]="sujet"
  colnames(data)[2]="produit"
  data$sujet=as.factor(as.character(data$sujet))
  data$produit=as.factor(as.character(data$produit))
  for (j in 3:ncol(data)){
    classe=class(data[,j])
    if (!classe%in%c("numeric","integer")){
      stop("Contingency data must be numeric or integer")
    }
  }
  check.bin=unique(unlist(data[,3:ncol(data)]))
  if (length(check.bin)>2){
    warning("contingency data are not composed of only ones and zeros")
  }else{
    check.un=sum(check.bin==c(0,1))
    check.deux=sum(check.bin==c(1,0))
    if (check.un!=2 & check.deux!=2){
      warning("contingency data are not composed of only ones and zeros")
    }
  }
  data=data[order(data$sujet,data$produit),]
  rownames(data)=as.character(1:nrow(data))

  org=aggregate(.~produit,data,sum)
  org$sujet=NULL
  rownames(org)=as.character(org$produit)
  org$produit=NULL

  verif.col=colSums(org)
  if (any(verif.col==0)){
    stop("Some descriptors have never been selected")
  }
  param.etendu=table(data$sujet,data$produit)
  if (length(unique(param.etendu))>1){
    warning("Data are unbalanced: products have not been evaluated a same number of times by subjects")
  }
  nplus=colSums(param.etendu)
  redui=org
  o=org
  mr=nplus
  N=sum(mr)
  mc=colSums(o)
  fij=mr%o%mc/N
  std=(((o-fij))/sqrt(fij))/(sqrt(N))
  nb.axe=min(nrow(std)-1,ncol(std))
  udv=svd(std)
  vs=udv$d[1:nb.axe]
  eig=vs^2
  chi.obs=eig
  for (j in 1:length(chi.obs)){
    chi.obs[j]=sum(eig[j:length(eig)])*N
  }

  registerDoParallel(cores = ncores)
  nperm=nperm
  sortie <- foreach(icount(nperm), .combine='rbind') %dopar% {
    virt.data=data

    for(s.perm in levels(virt.data$sujet)){
      l.s=which(virt.data$sujet==s.perm)
      loto=sample(l.s,length(l.s),replace = FALSE)
      virt.data$produit[l.s]=virt.data$produit[loto]
    }

    org=aggregate(.~produit,virt.data,sum)
    org$sujet=NULL
    rownames(org)=as.character(org$produit)
    org$produit=NULL

    param.etendu=table(virt.data$sujet,virt.data$produit)
    nplus=colSums(param.etendu)
    o=org
    mr=nplus
    N=sum(mr)
    mc=colSums(o)
    fij=mr%o%mc/N
    std=(((o-fij))/sqrt(fij))/(sqrt(N))
    nb.axe=min(nrow(std)-1,ncol(std))
    udv=svd(std)
    vs=udv$d[1:nb.axe]
    eig=vs^2

    chi.virt=eig
    for (j in 1:length(chi.virt)){
      chi.virt[j]=sum(eig[j:length(eig)])*N
    }
    return(chi.virt)
  }
  sortie=rbind(sortie,chi.obs)
  calc.pval=function(vec){
    obs=vec[length(vec)]
    virt=vec[-length(vec)]
    pval=(sum(virt>=obs)+1)/(nperm+1)
  }
  back.pval=apply(sortie, 2, calc.pval)
  for (i in 1:length(back.pval)){
    back.pval[i]=max(back.pval[1:i])
  }
  nom=paste("Dim.",1:length(chi.obs))
  names(back.pval)=names(chi.obs)=nom
  dim.sig=back.pval<=alpha
  ou=match(FALSE,dim.sig)
  if(!is.na(ou)){
    dim.sig=ou-1
  }else{
    dim.sig=length(back.pval)
  }
  axe.test=list(dim.sig=dim.sig,statistics=chi.obs,p.values=back.pval)
  back=axe.test
  stopImplicitCluster()
  return(back)
}
