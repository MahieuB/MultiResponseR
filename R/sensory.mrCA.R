#' Multiple-response Correspondence Analysis (MR-CA) for sensory data
#'
#' @description This function performs the MR-CA of the data as well as the total bootstrap procedure (Cadoret & Husson, 2013) and the pairwise total bootstrap tests as proposed in Castura et al. (2023). The difference with \code{\link[MultiResponseR]{mrCA}} used with ellipse=TRUE is that the total bootstrap procedure is stratified with respect to subjects in \code{\link[MultiResponseR]{sensory.mrCA}}
#'
#' @param data A data.frame of evaluations in rows whose first two columns are factors (subject and product) and subsequent columns are binary numeric or integer, each column being a descriptor
#' @param nboot The number of bootstrapped panel of the total bootstrap procedure
#' @param nbaxes.sig The number of significant axes retuned by \code{\link[MultiResponseR]{sensory.mr.dimensionality.test}}. By default, all axes are considered significant. See details
#'
#' @details
#' \itemize{
#'   \item \strong{nbaxes.sig}: The number of significant axes determines the number of axes accounted for while performing the Procrustes rotations of the total bootstrap procedure (Mahieu, Schlich, Visalli, & Cardot, 2021). These same axes are accounted for the pairwise total bootstrap tests.
#' }
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{eigen}{Eigenvalues of the MR-CA and their corresponding percentages of inertia}
#'   \item{prod.coord}{Products coordinates}
#'   \item{desc.coord}{Descriptors coordinates}
#'   \item{svd}{Results of the singular value decomposition}
#'   \item{bootstrap.replicate.coord}{Coordinates of the rotated bootstrap replicates}
#'   \item{total.bootstrap.test.pvalues}{P-values of the pairwise total bootstrap tests}
#' }
#' @export
#'
#' @import stats
#' @import utils
#'
#' @references Cadoret, M., & Husson, F. (2013). Construction and evaluation of confidence ellipses applied at sensory data. Food Quality and Preference, 28(1), 106-115.
#' @references Castura, J. C., Varela, P., & Næs, T. (2023). Evaluation of complementary numerical and visual approaches for investigating pairwise comparisons after principal component analysis. Food Quality and Preference, 107.
#' @references Mahieu, B., Schlich, P., Visalli, M., & Cardot, H. (2021). A multiple-response chi-square framework for the analysis of Free-Comment and Check-All-That-Apply data. Food Quality and Preference, 93.
#'
#' @examples
#'
#'data(milkchoc)
#'
#'dim.sig=sensory.mr.dimensionality.test(milkchoc)$dim.sig
#'
#'res=sensory.mrCA(milkchoc,nbaxes.sig=dim.sig)
#'
#'plot(res)
sensory.mrCA=function(data,nboot=2000,nbaxes.sig=Inf){
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

  if (nbaxes.sig==1){
    nbaxes.proc=2
  }else if (nbaxes.sig==Inf){
    nbaxes.proc=min((nrow(org)-1),ncol(org))
    nbaxes.sig=min((nrow(org)-1),ncol(org))
  }else if (nbaxes.sig==0){
    stop("nbaxes.sig is equal to zero")
  }else{
    nbaxes.proc=nbaxes.sig
  }

  param.etendu=table(data$sujet,data$produit)
  balanced=TRUE
  if (length(unique(param.etendu))>1){
    warning("Data are unbalanced: products have not been evaluated a same number of times by subjects")
    balanced=FALSE
  }
  nplus=colSums(param.etendu)
  redui=org

  ca.actual=mrCA(data[,-1])
  actual.coord=ca.actual$row.coord[,1:nbaxes.proc]

  myProcrustes=function(tofit,target,row.weights=rep(1,nrow(target)),scaling=FALSE){
    if (!is.matrix(tofit) & !is.data.frame(tofit)){
      stop("tofit must be a matrix or a data.frame")
    }
    if (!is.matrix(target) & !is.data.frame(target)){
      stop("tofit must be a matrix or a data.frame")
    }
    X=as.matrix(tofit)
    Y=as.matrix(target)
    if (ncol(X)!=ncol(Y) | nrow(X)!=nrow(Y)){
      stop("Dimension of tofit and target are different")
    }
    if (!is.logical(scaling)){
      stop("scaling must be logical")
    }
    if (is.numeric(row.weights) | is.integer(row.weights)){
      if (length(row.weights)!=nrow(target)){
        stop("length(row.weights) must equal nrow(target)")
      }
    }else{
      stop("class(row.weights) must be numeric or integer")
    }
    vecW=row.weights/sum(row.weights)
    mX=t(as.matrix(vecW))%*%X
    X=sweep(X,2,mX,"-")
    mY=t(as.matrix(vecW))%*%Y
    loc=mY
    Y=sweep(Y,2,mY,"-")
    matW=diag(row.weights/sum(row.weights)*nrow(target))
    croise=t(X)%*%matW%*%Y
    udv=svd(croise)
    if (scaling){
      dilat=sum(udv$d)/sum(diag(crossprod(X)))
    }else{
      dilat=1
    }
    H=udv$u%*%t(udv$v)
    aligned=(dilat*X%*%H)
    aligned=sweep(aligned,2,loc,"+")
    return(aligned)
  }

  row.s=matrix(NA,nlevels(data$produit),nlevels(data$sujet))
  colnames(row.s)=levels(data$sujet)
  for (s in unique(data$sujet)){
    ou.s=which(data$sujet==s)
    row.s[1:length(ou.s),s]=ou.s
  }

  sortie=data.frame(produit=as.factor(rep(levels(data$produit),nboot)),matrix(0,nlevels(data$produit)*nboot,nbaxes.proc))
  colnames(sortie)[2:ncol(sortie)]=colnames(actual.coord)
  vec.boot=seq(1,nrow(sortie),by=nlevels(data$produit))
  pb=txtProgressBar(1,nboot,style=3)
  for (boot in 1:nboot) {
    if (balanced){
      tirage=sample(unique(data$sujet),nlevels(data$sujet),replace = T)
      vec.ligne = as.vector(row.s[,tirage])
      jdd.tirage=na.omit(data[vec.ligne,])
    }else{
      retire=TRUE
      while(retire){
        tirage=sample(unique(data$sujet),nlevels(data$sujet),replace = T)
        vec.ligne = as.vector(row.s[,tirage])
        jdd.tirage=na.omit(data[vec.ligne,])
        if (all(table(jdd.tirage$produit)>0)){
          retire=FALSE
        }
      }
    }
    jdd.tirage=jdd.tirage[order(jdd.tirage$sujet,jdd.tirage$produit),]
    rownames(jdd.tirage)=as.character(1:nrow(jdd.tirage))

    verif=colSums(jdd.tirage[,-c(1:2)])
    vire = which(verif==0)
    if(length(vire)!=0){
      vire=vire+2
      jdd.tirage=jdd.tirage[,-vire]
    }

    ca.tirage=mrCA(jdd.tirage[,-1])
    tirage.coord=ca.tirage$row.coord[,1:nbaxes.proc]
    rot = myProcrustes(tirage.coord, actual.coord,nplus)
    sortie[vec.boot[boot]:(vec.boot[boot]+nlevels(data$produit)-1),2:ncol(sortie)]=rot
    setTxtProgressBar(pb,boot)
  }

  sortie$produit=as.factor(sortie$produit)
  toellipse=sortie
  rownames(toellipse)=as.character(1:nrow(toellipse))


  coord.boot=toellipse[,1:(nbaxes.sig+1)]
  dim.sig=nbaxes.sig

  les.prod=unique(coord.boot$produit)

  diff.test=as.data.frame(matrix(1,length(les.prod),length(les.prod)))
  colnames(diff.test)=rownames(diff.test)=les.prod

  for(i in 1:(nrow(diff.test)-1)){
    for (j in (i+1):ncol(diff.test)){
      p.1=rownames(diff.test)[i]
      p.2=colnames(diff.test)[j]
      coord.p1=coord.boot[coord.boot$produit==p.1,,drop=FALSE]
      coord.p2=coord.boot[coord.boot$produit==p.2,,drop=FALSE]
      delta=as.matrix(coord.p1[,-1,drop=FALSE]-coord.p2[,-1,drop=FALSE])
      mu=colMeans(delta)
      sigma=cov(delta)
      calc.mal.sq=function(vec){
        mal.bary=t(as.matrix(vec-mu))%*%solve(sigma,tol=1e-300)%*%(as.matrix(vec-mu))
        return(as.numeric(mal.bary))
      }
      mal.sq.cloud=apply(as.matrix(delta),1,calc.mal.sq)
      mal.sq.origin=as.numeric(t(as.matrix(rep(0,length(mu))-mu))%*%solve(sigma)%*%(as.matrix(rep(0,length(mu))-mu)))
      pvalue.test=(sum(mal.sq.cloud>=mal.sq.origin)+1)/(nrow(delta)+1)
      diff.test[p.1,p.2]=pvalue.test
      diff.test[p.2,p.1]=pvalue.test
    }
  }
  back=list(eigen=ca.actual$eigen,prod.coord=ca.actual$row.coord,desc.coord=ca.actual$col.coord,
            svd=ca.actual$svd,bootstrap.replicate.coord=toellipse,total.bootstrap.test.pvalues=diff.test)
  class(back)=c("sensory.mrCA","list")
  return(back)
}
