#' Multiple-response tests per cell for sensory data
#'
#' @description This function performs for each pair of product and descriptor a multiple-response hypergeometric test as defined in Mahieu, Schlich, Visalli, and Cardot (2021) using random hypergeometric samplings to estimate the null distribution. The difference with \code{\link[MultiResponseR]{mr.sig.cell}} is that random hypergeometric samplings are performed within subjects in \code{\link[MultiResponseR]{sensory.mr.sig.cell}}
#'
#' @param data A data.frame of evaluations in rows whose first two columns are factors (subject and product) and subsequent columns are binary numeric or integer, each column being a descriptor
#' @param nsample Number of randomly sampled datasets to estimate the distribution of the value under the null hypothesis. See details
#' @param nbaxes.sig The number of significant axes retuned by \code{\link[MultiResponseR]{sensory.mr.dimensionality.test}}. By default, all axes are considered significant. See details
#' @param two.sided Logical. Should the tests be two-sided or not?
#'
#' @details
#' \itemize{
#'   \item \strong{nsample}: The distribution of the value under the null hypothesis of no associations between products and descriptors is estimated using \emph{nsample} datasets generated thanks to random hypergeometric samplings of the response vectors along products within subjects.
#'   \item \strong{nbaxes.sig}: If \emph{nbaxes.sig} is lower than the total number of axes then the tests are performed on the derived contingency table corresponding to significant axes (Mahieu, Schlich, Visalli, & Cardot, 2021) This table is obtained by using the reconstitution formula of MR-CA on the first \emph{nbaxes.sig} axes.
#' }
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{original.cont}{Observed number of times each product was described by each descriptor}
#'   \item{percent.cont}{For each product, percentage of evaluations where each descriptor was cited for this product}
#'   \item{null.cont}{Expected number of times each product was described by each descriptor under the null hypothesis}
#'   \item{p.values}{P-values of the tests per cell}
#'   \item{derived.cont}{The derived contingency table corresponding to \emph{nbaxes.sig} axes}
#'   \item{percent.derived.cont}{For each product, percentage of evaluations where each descriptor was cited for this product in the derived contingency table corresponding to \emph{nbaxes.sig} axes}
#' }
#' @export
#'
#' @import stats
#' @import utils
#'
#' @references Loughin, T. M., & Scherer, P. N. (1998). Testing for Association in Contingency Tables with Multiple Column Responses. Biometrics, 54(2), 630-637.
#' @references Mahieu, B., Schlich, P., Visalli, M., & Cardot, H. (2021). A multiple-response chi-square framework for the analysis of Free-Comment and Check-All-That-Apply data. Food Quality and Preference, 93.
#'
#'
#' @examples
#'data(milkchoc)
#'
#'dim.sig=sensory.mr.dimensionality.test(milkchoc)$dim.sig
#'
#'res=sensory.mr.sig.cell(milkchoc,nbaxes.sig=dim.sig)
#'
#'plot(res)
sensory.mr.sig.cell=function(data,nsample=2000,nbaxes.sig=Inf,two.sided=TRUE){
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

  sorted.name=sort(colnames(data[,-c(1:2)]))
  d=data[,sorted.name]
  g=data[,c("sujet","produit")]
  data=cbind.data.frame(g,d)

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
  nom=names(nplus)
  nplus=as.numeric(nplus)
  names(nplus)=nom

  nplusplus=sum(nplus)

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
    mr=as.numeric(table(tab[,1])/sum(table(tab[,1])))
    mc=as.numeric(colSums(tab[,-1])/sum(table(tab[,1])))
    e=mr%o%mc
    theo.cont=((theo*sqrt(e))+e)*sum(table(tab[,1]))
    theo.cont=round(ifelse(theo.cont<0,0,theo.cont))
    return(theo.cont)
  }
  theo.cont=calc.cont(data[,-1])

  row.s=matrix(NA,nlevels(data$produit),nlevels(data$sujet))
  colnames(row.s)=levels(data$sujet)
  for (s in unique(data$sujet)){
    ou.s=which(data$sujet==s)
    row.s[1:length(ou.s),s]=ou.s
  }
  mySample=function(vec){
    vec.retour=na.omit(vec)
    if (length(vec.retour)>1){
      vec.retour=sample(vec.retour,length(vec.retour),replace = TRUE)
    }else{
      vec.retour=vec.retour[1]
    }
    return(vec.retour)
  }

  sortie=array(0,c(nrow(theo.cont),ncol(theo.cont),nsample))
  pb=txtProgressBar(1,nsample,style = 3)
  for (ssample in 1:nsample){
    virt.data=data
    tirage=unlist(apply(row.s,2,mySample))
    virt.data[,3:ncol(virt.data)]=data[tirage,3:ncol(data)]

    verif=colSums(virt.data[,-c(1:2)])
    vire = which(verif==0)
    if(length(vire)!=0){
      nom.vire=names(vire)
      vire=vire+2
      virt.data=virt.data[,-vire]
    }

    theo.cont.virt=calc.cont(virt.data[,-1])
    if(length(vire)!=0){
      theo.zero=matrix(0,nrow = nrow(theo.cont.virt),ncol = length(nom.vire))
      rownames(theo.zero)=rownames(theo.cont.virt)
      colnames(theo.zero)=nom.vire
      theo.cont.virt=cbind(theo.cont.virt,theo.zero)
      theo.cont.virt=theo.cont.virt[,sorted.name]
    }
    sortie[,,ssample]=theo.cont.virt
    setTxtProgressBar(pb,ssample)
  }

  bb=colSums(org)/nplusplus
  aa=nplus/nplusplus

  med.mat=aa%o%bb*nplusplus

  back.pval=as.data.frame(matrix(0,nrow(org),ncol(org),dimnames=dimnames(org)))

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

  org=as.data.frame(t(org[,sorted.name]))
  percent.cont=as.data.frame(t(as.data.frame(round(t(org)/nplus*100,2))))
  back.pval=as.data.frame(t(back.pval[,sorted.name]))
  adj.back.pval=as.data.frame(t(apply(back.pval,1,p.adjust,method="fdr")))

  back=list(original.cont=org,percent.cont=percent.cont,null.cont=as.data.frame(t(med.mat[,sorted.name])),p.value=round(adj.back.pval,4),derived.cont=as.data.frame(t(theo.cont[,sorted.name])),
            percent.derived.cont=as.data.frame(round(t(theo.cont[,sorted.name]/nplus*100),2)))
  class(back)=c("sensory.mr.sig.cell","list")
  return(back)
}
