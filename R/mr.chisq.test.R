#' Multiple-response chi-square test
#'
#' @description Performs a multiple-response chi-square test as defined in Loughin and Scherer (1998) using random permutations to estimate the null distribution
#'
#' @param data A data.frame of observations in rows whose first column is a factor (the categories) and subsequent columns are binary numeric or integer, each column being a response option
#' @param nperm Number of permuted datasets to estimate the distribution of the statistic under the null hypothesis. See details
#' @param ncores Number of cores used to estimate the null distribution. Default is 2. See details
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{statistic}{Observed multiple-response chi-square statistic}
#'   \item{p.value}{p-value of the test}
#' }
#' @details
#' \itemize{
#'   \item \strong{nperm}: The distribution of the statistic under the null hypothesis of no associations between categories and response options is estimated using \emph{nperm} datasets generated thanks to random permutations of the response vectors along observations. Note that this differs from the original proposition of Loughin and Scherer (1998) who used a parametric bootstrap to do so.
#'   \item \strong{ncores}: The more cores are added in the process, the faster the results will be obtained. The number of available cores is accessible using \code{\link[parallel]{detectCores}}. The parallel tasks are closed once the \emph{nperm} datasets are generated.
#' }
#' @export
#' @references Loughin, T. M., & Scherer, P. N. (1998). Testing for Association in Contingency Tables with Multiple Column Responses. Biometrics, 54(2), 630-637.
#' @references Mahieu, B., Schlich, P., Visalli, M., & Cardot, H. (2020). A multiple-response chi-square framework for the analysis of Free-Comment and Check-All-That-Apply data. Manuscript submitted for publication.
#'
#' @import parallel
#' @import doParallel
#' @import foreach
#' @import iterators
#' @import stats
#'
#' @examples
#' nb.obs=200
#' nb.response=5
#' nb.category=5
#' vec.category=paste("C",1:nb.category,sep="")
#' right=matrix(rbinom(nb.response*nb.obs,1,0.25),nb.obs,nb.response)
#' category=sample(vec.category,nb.obs,replace = TRUE)
#' dset=cbind.data.frame(category,right)
#' dset$category=as.factor(dset$category)
#'
#' parallel::detectCores()
#'
#' mr.chisq.test(dset)
mr.chisq.test=function(data,nperm=2000,ncores=2){
  classe=class(data)[1]
  if (!classe%in%c("data.frame")){
    stop("data must a data.frame")
  }
  classe=class(data[,1])
  if(!classe%in%c("factor")){
    stop("First column of data must be a factor")
  }
  for (j in 2:ncol(data)){
    classe=class(data[,j])[1]
    if (!classe%in%c("numeric","integer")){
      stop("contingency data must be integer or numeric")
    }
  }
  check.bin=unique(unlist(data[,2:ncol(data)]))
  if (length(check.bin)>2){
    warning("contingency data are not composed of only ones and zeros")
  }else{
    check.un=sum(check.bin==c(0,1))
    check.deux=sum(check.bin==c(1,0))
    if (check.un!=2 & check.deux!=2){
      warning("contingency data are not composed of only ones and zeros")
    }
  }
  colnames(data)[1]="category"
  data=data[order(data$category),]
  rownames(data)=as.character(1:nrow(data))
  original=aggregate(.~category,data,sum)
  rownames(original)=original$category
  original$category=NULL
  verif.col=colSums(original)
  if (any(verif.col==0)){
    stop("Some responses have never been selected")
  }
  nplus=table(data$category)
  nplusplus=sum(nplus)
  if (any(nplus==0)){
    stop("Some categories are not represented")
  }
  o=original
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
    loto=sample(1:nrow(virt.data),nrow(virt.data),replace = F)
    virt.data[,2:ncol(virt.data)]=virt.data[loto,2:ncol(virt.data)]

    original=aggregate(.~category,virt.data,sum)
    rownames(original)=original$category
    original$category=NULL
    nplus=table(virt.data$category)
    nplusplus=sum(nplus)
    o=original
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
  chi.test=list(statistic=chi.obs[1],p.value=back.pval[1])
  back=chi.test
  stopImplicitCluster()
  return(back)
}
