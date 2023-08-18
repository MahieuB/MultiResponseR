#' Multiple-response tests per cell
#'
#' @description This function performs for each pair of category and response option a multiple-response hypergeometric test as defined in Mahieu, Schlich, Visalli, and Cardot (2021) using random hypergeometric samplings to estimate the null distribution
#'
#' @param data A data.frame of observations in rows whose first column is a factor (the categories) and subsequent columns are binary numeric or integer, each column being a response option
#' @param nsample Number of randomly sampled datasets to estimate the distribution of the value under the null hypothesis. See details
#' @param nbaxes.sig The number of significant axes retuned by \code{\link[MultiResponseR]{mr.dimensionality.test}}. By default, all axes are considered significant. See details
#' @param two.sided Logical. Should the tests be two-sided or not? By default, the tests are performed with a one-sided greater alternative hypothesis
#'
#' @details
#' \itemize{
#'   \item \strong{nsample}: The distribution of the value under the null hypothesis of no associations between categories and response options is estimated using \emph{nsample} datasets generated thanks to random hypergeometric samplings of the response vectors along observations.
#'   \item \strong{nbaxes.sig}: If \emph{nbaxes.sig} is lower than the total number of axes then the tests are performed on the derived contingency table corresponding to significant axes (Mahieu, Schlich, Visalli, & Cardot, 2021). This table is obtained by using the reconstitution formula of MR-CA on the first \emph{nbaxes.sig} axes.
#' }
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{original.cont}{Observed number of times each category chosen each response option}
#'   \item{percent.cont}{Within each category, percentage of observations where the response options were chosen}
#'   \item{null.cont}{Expected number of times each category chosen each response option under the null hypothesis}
#'   \item{p.values}{P-values of the tests per cell fdr adjusted by response option}
#'   \item{derived.cont}{The derived contingency table corresponding to \emph{nbaxes.sig} axes}
#'   \item{percent.derived.cont}{Within each category, percentage of observations where the response options were chosen in the derived contingency table corresponding to \emph{nbaxes.sig} axes}
#' }
#'
#'
#' @export
#'
#' @references Loughin, T. M., & Scherer, P. N. (1998). Testing for Association in Contingency Tables with Multiple Column Responses. Biometrics, 54(2), 630-637.
#' @references Mahieu, B., Schlich, P., Visalli, M., & Cardot, H. (2021). A multiple-response chi-square framework for the analysis of Free-Comment and Check-All-That-Apply data. Food Quality and Preference, 93.
#'
#' @import stats
#' @import utils
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
#' res=mr.sig.cell(dset)
#'
#' plot(res)
mr.sig.cell=function(data,nsample=2000,nbaxes.sig=Inf,two.sided=FALSE){
  classe=class(data)[1]
  if (!classe%in%c("data.frame")){
    stop("cont must be a data.frame")
  }
  classe=class(data[,1])
  if(!classe%in%c("factor")){
    stop("First column of data must be a factor")
  }
  for (j in 2:ncol(data)){
    classe=class(data[,j])[1]
    if (!classe%in%c("numeric","integer")){
      stop("data must be integer or numeric")
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

  sorted.name=sort(colnames(data[,-1]))
  d=data[,sorted.name]
  g=data[,c("category")]
  data=cbind.data.frame(g,d)
  colnames(data)[1]="category"
  data=data[order(data$category),]
  rownames(data)=as.character(1:nrow(data))


  nplus=table(data$category)
  nom=names(nplus)
  nplus=as.numeric(nplus)
  names(nplus)=nom

  nplusplus=sum(nplus)

  org=aggregate(.~category,data,sum)
  rownames(org)=as.character(org$category)
  org$category=NULL
  original=org

  if (nbaxes.sig==0){
    stop("nbaxes.sig is equal to zero")
  }else if (nbaxes.sig==Inf){
    nbaxes.sig=min(nrow(org)-1,ncol(org))
  }

  calc.cont=function(tab){
    res=mrCA(tab)
    if (nbaxes.sig==1){
      d.vs=res$svd$vs[1]
    }else{
      d.vs=diag(res$svd$vs)
      d.vs=d.vs[1:nbaxes.sig,1:nbaxes.sig]
    }
    theo=res$svd$u[,1:nbaxes.sig,drop=FALSE]%*%d.vs%*%t(res$svd$v[,1:nbaxes.sig,drop=FALSE])
    mr=nplus/nplusplus
    mc=colSums(org)/nplusplus
    e=mr%o%mc
    theo.cont=((theo*sqrt(e))+e)*nplusplus
    theo.cont=round(ifelse(theo.cont<0,0,theo.cont))
    return(theo.cont)
  }
  theo.cont=calc.cont(data)

  sortie=array(0,c(nrow(theo.cont),ncol(theo.cont),nsample))
  pb=txtProgressBar(min=0,max=nsample,style=3)
  for (ssample in 1:nsample){
    virt.data=data
    loto=sample(1:nrow(virt.data),nrow(virt.data),replace = TRUE)
    virt.data[,2:ncol(virt.data)]=virt.data[loto,2:ncol(virt.data)]

    original.virt=aggregate(.~category,virt.data,sum)
    rownames(original.virt)=original.virt$category
    original.virt$category=NULL

    verif=colSums(original.virt)
    vire = which(verif==0)
    if(length(vire)!=0){
      nom.zero=names(vire)
      original.virt=original.virt[,-vire]
      vire.extend=NULL
      for (quel.vire in names(vire)){
        ou.virer=which(colnames(virt.data)==quel.vire)
        vire.extend=c(vire.extend,ou.virer)
      }
      virt.data=virt.data[,-vire.extend]
    }

    theo.cont.virt=calc.cont(virt.data)
    sortie[,,ssample]=theo.cont.virt
    setTxtProgressBar(pb,ssample)
  }

  bb=colSums(org)/nplusplus
  aa=nplus/nplusplus

  med.mat=aa%o%bb*nplusplus

  back.pval=as.data.frame(matrix(0,nrow(original),ncol(original),dimnames=dimnames(original)))

  for (i in 1:dim(sortie)[1]){
    for (j in 1:dim(sortie)[2]){
      if(two.sided){
        obs=theo.cont[i,j]
        virt=sortie[i,j,]
        g=sum(virt<=obs)
        d=sum(virt>=obs)
        pp=((min(g,d)+1)/(nsample+1))*2
        if (pp>1){
          pp=1
        }
        back.pval[i,j]=pp
      }else{
        obs=theo.cont[i,j]
        virt=sortie[i,j,]
        pp=(sum(virt>=obs)+1)/(nsample+1)
        back.pval[i,j]=pp
      }
    }
  }

  original=as.data.frame(t(original[,sorted.name]))
  percent.cont=as.data.frame(t(as.data.frame(round(as.matrix((org)/nplus*100),2))[,sorted.name]))
  back.pval=as.data.frame(t(back.pval[,sorted.name]))
  adj.back.pval=as.data.frame(t(apply(back.pval,1,p.adjust,method="fdr")))

  back=list(original.cont=original,percent.cont=percent.cont,null.cont=as.data.frame(t(med.mat[,sorted.name])),p.value=round(adj.back.pval,4),derived.cont=as.data.frame(t(theo.cont[,sorted.name])),
            percent.derived.cont=as.data.frame(round(t(theo.cont[,sorted.name]/nplus*100),2)))


  class(back)=c("mr.sig.cell","list")
  return(back)
}


