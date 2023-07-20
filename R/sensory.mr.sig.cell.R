#' Multiple-response tests per cell for sensory data
#'
#' @description This function performs for each pair of product and descriptor a multiple-response hypergeometric test as defined in Mahieu, Schlich, Visalli, and Cardot (2021) using random hypergeometric samplings to estimate the null distribution. The difference with \code{\link[MultiResponseR]{mr.sig.cell}} is that random hypergeometric samplings are performed within subjects in \code{\link[MultiResponseR]{sensory.mr.sig.cell}}
#'
#' @param data A data.frame of evaluations in rows whose first two columns are factors (subject and product) and subsequent columns are binary numeric or integer, each column being a descriptor
#' @param nsample Number of randomly sampled datasets to estimate the distribution of the value under the null hypothesis. See details
#' @param nbaxes.sig The number of significant axes retuned by \code{\link[MultiResponseR]{sensory.mr.dimensionality.test}}. By default, all axes are considered significant. See details
#' @param two.sided Logical. Should the tests be two-sided or not? By default, the tests are performed with a one-sided greater alternative hypothesis
#' @param ncores Number of cores used to estimate the null distribution. Default is 2.
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
#' @import foreach
#' @import parallel
#' @import iterators
#' @import doParallel
#' @import abind
#' @import stats
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
sensory.mr.sig.cell=function(data,nsample=2000,nbaxes.sig=Inf,two.sided=FALSE,ncores=2){
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
    mr=nplus/nplusplus
    mc=colSums(org)/nplusplus
    e=mr%o%mc
    theo.cont=((theo*sqrt(e))+e)*nplusplus
    theo.cont=round(ifelse(theo.cont<0,0,theo.cont))
    return(theo.cont)
  }

  theo.cont=calc.cont(data[,-1])

  acomb=function(...){
    abind(...,along = 3)
  }
  registerDoParallel(cores = ncores)
  nsample=nsample

  sortie <- foreach(icount(nsample), .combine='acomb',.multicombine = TRUE) %dopar% {
    virt.data=data

    for(s.perm in levels(virt.data$sujet)){
      l.s=which(virt.data$sujet==s.perm)
      if (length(l.s)==1){
        loto=l.s
      }else{
        loto=sample(l.s,length(l.s),replace = TRUE)
      }
      virt.data[l.s,3:ncol(virt.data)]=virt.data[loto,3:ncol(virt.data)]
    }

    org=aggregate(.~produit,virt.data,sum)
    org$sujet=NULL
    rownames(org)=as.character(org$produit)
    org$produit=NULL

    param.etendu=table(virt.data$sujet,virt.data$produit)

    nplus=colSums(param.etendu)
    nom=names(nplus)
    nplus=as.numeric(nplus)
    names(nplus)=nom

    nplusplus=sum(nplus)

    verif=colSums(org)
    vire = which(verif==0)
    if(length(vire)!=0){
      nom.zero=names(vire)
      org=org[,-vire]
      vire.extend=NULL
      for (quel.vire in names(vire)){
        ou.virer=which(colnames(virt.data)==quel.vire)
        vire.extend=c(vire.extend,ou.virer)
      }
      virt.data=virt.data[,-vire.extend]
    }

    ########

    mrCA=function(data,proj.row=NULL,proj.row.obs=NULL,proj.col=NULL,ellipse=FALSE,nboot=2000,nbaxes.sig=Inf,ncores=2){
      classe=class(data)[1]
      if (!classe%in%c("data.frame")){
        stop("data must be a data.frame")
      }
      classe=class(data[,1])[1]
      if (!classe%in%c("factor")){
        stop("The first column of data must be factor")
      }
      if(!colnames(data)[1]%in%c("category","Category","produit","Produit","product","Product")){
        stop("First column name must be category, Category, produit, Produit, product or Product")
      }
      colnames(data)[1]="cat"
      data$cat=as.factor(as.character(data$cat))
      for (j in 2:ncol(data)){
        classe=class(data[,j])
        if (!classe%in%c("numeric","integer")){
          stop("Contingency data must be numeric or integer")
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
      verif.col=colSums(data[,2:ncol(data)])
      if (any(verif.col==0)){
        stop("Some columns are only zeros")
      }

      cont=aggregate(.~cat,data,sum)
      rownames(cont)=cont$cat
      cont$cat=NULL

      if (nbaxes.sig==1){
        nbaxes.proc=2
      }else if (nbaxes.sig==Inf){
        nbaxes.proc=min((nrow(cont)-1),ncol(cont))
        nbaxes.sig=min((nrow(cont)-1),ncol(cont))
      }else if (nbaxes.sig==0){
        stop("nbaxes.sig is equal to zero")
      }else{
        nbaxes.proc=nbaxes.sig
      }

      row.obs=table(data$cat)
      nom.cat=names(row.obs)
      row.obs=as.numeric(row.obs)
      names(row.obs)=nom.cat

      nplus=row.obs
      nplusplus=sum(nplus)
      redui=cont
      fij=(nplus%o%colSums(redui))/nplusplus
      std=((redui-fij)/sqrt(fij))/sqrt(nplusplus)

      n=nrow(redui)
      p=ncol(redui)
      nb.axe=min(n-1,p)

      udv=svd(std)

      u=udv$u[,1:nb.axe,drop=FALSE]
      vs=udv$d[1:nb.axe]
      eig=vs^2
      v=udv$v[,1:nb.axe,drop=FALSE]
      rownames(u)=rownames(redui)
      rownames(v)=colnames(redui)

      marge.r=nplus/nplusplus
      if (nb.axe==1){
        dvs=vs
      }else{
        dvs=diag(vs)
      }

      row.coord=(u/sqrt(marge.r))%*%dvs
      col.coord=v

      if(!is.null(proj.row)){
        classe=class(proj.row)[1]
        if (!classe%in%c("matrix","data.frame")){
          stop("proj.row must be a matrix or a data.frame")
        }
        if (ncol(proj.row)!=ncol(cont)){
          stop("proj.row must have the same number of response options that data")
        }
        classe=class(proj.row.obs)[1]
        if (!classe%in%c("numeric","integer")){
          stop("proj.row.obs must be integer or numeric")
        }
        if (nrow(proj.row)!=length(proj.row.obs)){
          stop("proj.row and proj.row.obs have not the same size")
        }
        unit.supp=proj.row/proj.row.obs
        fij.sup=(colSums(redui))/nplusplus
        mat.fij.sup=as.matrix(rep(1,nrow(proj.row)))%*%fij.sup
        vec.sup=((unit.supp-mat.fij.sup)/sqrt(mat.fij.sup))
        proj.row.coord=as.matrix(vec.sup)%*%v
        rownames(proj.row.coord)=rownames(proj.row)
      }
      if(!is.null(proj.col)){
        classe=class(proj.col)[1]
        if (!classe%in%c("matrix","data.frame")){
          stop("proj.col must be a matrix or a data.frame")
        }
        if (nrow(proj.col)!=nrow(cont)){
          stop("proj.col must have the same number of categories that data")
        }
        c.supp=colSums(proj.col)/nplusplus
        o.supp=proj.col/nplusplus
        r.supp=nplus/nplusplus
        e.supp=r.supp%o%c.supp
        std.supp=(o.supp-e.supp)/sqrt(e.supp)
        proj.col.coord=t(t(u)%*%as.matrix(std.supp)/vs)
        rownames(proj.col.coord)=colnames(proj.col)
      }

      percent.eig=eig/sum(eig)*100
      cum.percent.eig=cumsum(percent.eig)
      mat.eig=cbind(eig,percent.eig,cum.percent.eig)
      colnames(mat.eig)=c("eigenvalue","percentage of inertia","cumulative percentage of inertia")
      name.dim=paste("Dim.",1:nrow(mat.eig))
      rownames(mat.eig)=colnames(col.coord)=colnames(row.coord)=name.dim
      if(!is.null(proj.row)){
        colnames(proj.row.coord)=name.dim
      }else{
        proj.row.coord=NULL
      }
      if(!is.null(proj.col)){
        colnames(proj.col.coord)=name.dim
      }else{
        proj.col.coord=NULL
      }
      colnames(u)=name.dim
      names(vs)=name.dim
      colnames(v)=name.dim

      if (ellipse){

        row.coord=row.coord[,1:nbaxes.proc,drop=FALSE]
        col.coord=col.coord[,1:nbaxes.proc,drop=FALSE]
        if(!is.null(proj.row)){
          proj.row.coord=proj.row.coord[,1:nbaxes.proc,drop=FALSE]
        }
        if(!is.null(proj.col)){
          proj.col.coord=proj.col.coord[,1:nbaxes.proc,drop=FALSE]
        }

        registerDoParallel(cores = ncores)
        nboot=nboot

        sortie <- foreach(icount(nboot), .combine='rbind') %dopar% {

          vec.ligne=NULL
          for (boot.cat in levels(data$cat)){
            les.ligne=which(data$cat==boot.cat)
            if (length(les.ligne)==1){
              loto=les.ligne
            }else{
              loto=sample(les.ligne,length(les.ligne),replace = TRUE)
            }
            vec.ligne=c(vec.ligne,loto)
          }

          jdd.tirage=data[vec.ligne,]
          rownames(jdd.tirage)=as.character(1:nrow(jdd.tirage))

          nplus.tirage=table(jdd.tirage$cat)
          nom.cat=names(nplus.tirage)
          nplus.tirage=as.numeric(nplus.tirage)
          names(nplus.tirage)=nom.cat

          cont.tirage = aggregate(.~cat,jdd.tirage,sum)
          rownames(cont.tirage)=as.character(cont.tirage$cat)
          cont.tirage$cat=NULL
          redui.tirage=cont.tirage

          verif=colSums(cont.tirage)
          vire = which(verif==0)
          if(length(vire)!=0){
            cont.tirage=cont.tirage[,-vire]
            redui.tirage=redui.tirage[,-vire]
            vire.extend=NULL
            for (quel.vire in names(vire)){
              ou.virer=which(colnames(jdd.tirage)==quel.vire)
              vire.extend=c(vire.extend,ou.virer)
            }
            jdd.tirage=jdd.tirage[,-vire.extend]
          }


          ######

          myProcrustes=function(tofit,target,scaling=FALSE){
            if (!is.matrix(tofit) & !is.data.frame(tofit)){
              stop("tofit must be a matrix or a data.frame")
            }
            X=as.matrix(tofit)
            X=sweep(X,2,colMeans(X),"-")
            if (!is.matrix(target) & !is.data.frame(target)){
              stop("tofit must be a matrix or a data.frame")
            }
            Y=as.matrix(target)
            if (!is.logical(scaling)){
              stop("scaling must be logical")
            }
            if (ncol(X)!=ncol(Y) | nrow(X)!=nrow(Y)){
              stop("Dimension of tofit and target are different")
            }
            loc=colMeans(Y)
            Y=sweep(Y,2,colMeans(Y),"-")
            croise=crossprod(X,Y)
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

          nplusplus.tirage=sum(nplus.tirage)
          fij.tirage=(nplus.tirage%o%colSums(redui.tirage))/nplusplus.tirage
          std.tirage=((redui.tirage-fij.tirage)/sqrt(fij.tirage))/sqrt(nplusplus.tirage)

          n.tirage=nrow(redui.tirage)
          p.tirage=ncol(redui.tirage)
          nb.axe.tirage=min(n.tirage-1,p.tirage)

          udv.tirage=svd(std.tirage)

          u.tirage=udv.tirage$u[,1:nb.axe.tirage]
          vs.tirage=udv.tirage$d[1:nb.axe.tirage]
          rownames(u.tirage)=rownames(redui.tirage)

          marge.r.tirage=nplus.tirage/nplusplus.tirage
          row.coord.tirage=(u.tirage/sqrt(marge.r.tirage))%*%diag(vs.tirage)
          row.coord.tirage=row.coord.tirage[,1:nbaxes.proc]
          ######

          rot = myProcrustes(row.coord.tirage, row.coord)
          rot=cbind.data.frame(rownames(rot),rot)
          colnames(rot)[1]="cat"
          return(rot)
        }

        sortie$cat=as.factor(sortie$cat)
        toellipse=sortie
        rownames(toellipse)=as.character(1:nrow(toellipse))

        coord.boot=toellipse[,1:(nbaxes.sig+1)]
        dim.sig=nbaxes.sig

        les.prod=unique(coord.boot$cat)

        diff.test=as.data.frame(matrix(1,length(les.prod),length(les.prod)))
        colnames(diff.test)=rownames(diff.test)=les.prod

        for(i in 1:(nrow(diff.test)-1)){
          for (j in (i+1):ncol(diff.test)){
            p.1=rownames(diff.test)[i]
            p.2=colnames(diff.test)[j]
            coord.p1=coord.boot[coord.boot$cat==p.1,,drop=FALSE]
            coord.p2=coord.boot[coord.boot$cat==p.2,,drop=FALSE]
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
      }else{
        toellipse=NULL
        diff.test=NULL
      }
      back=list(eigen=mat.eig,row.coord=row.coord,col.coord=col.coord,proj.row.coord=proj.row.coord,proj.col.coord=proj.col.coord,svd=list(u=u,vs=vs,v=v),bootstrap.replicate.coord=toellipse,total.bootstrap.test.pvalues=diff.test)
      class(back)=c("mrCA","list")
      return(back)
    }

    ########

    calc.cont=function(tab){
      res=mrCA(tab)
      if(length(vire)!=0){
        fun.sig=min(nbaxes.sig,nrow(res$eigen))
      }else{
        fun.sig=nbaxes.sig
      }
      if (fun.sig==1){
        d.vs=res$svd$vs[1]
      }else{
        d.vs=diag(res$svd$vs)
        d.vs=d.vs[1:fun.sig,1:fun.sig]
      }
      theo=res$svd$u[,1:fun.sig,drop=FALSE]%*%d.vs%*%t(res$svd$v[,1:fun.sig,drop=FALSE])
      mr=nplus/nplusplus
      mc=colSums(org)/nplusplus
      e=mr%o%mc
      if(length(vire)!=0){
        e=cbind(e,matrix(rep(rep(0,nrow(e)),length(nom.zero)),nrow(e),length(nom.zero)))
        colnames(e)[(ncol(e)-length(vire)+1):ncol(e)]=names(vire)
        e=e[,sorted.name]
        theo=cbind(theo,matrix(rep(rep(0,nrow(theo)),length(nom.zero)),nrow(theo),length(nom.zero)))
        colnames(theo)[(ncol(theo)-length(vire)+1):ncol(theo)]=names(vire)
        theo=theo[,sorted.name]
      }
      theo.cont=((theo*sqrt(e))+e)*nplusplus
      theo.cont=round(ifelse(theo.cont<0,0,theo.cont))
      return(theo.cont)
    }

    theo.cont=calc.cont(virt.data[,-1])
    return(theo.cont)
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

  back=list(original.cont=org,percent.cont=percent.cont,null.cont=as.data.frame(t(med.mat[,sorted.name])),p.value=round(back.pval,4),derived.cont=as.data.frame(t(theo.cont[,sorted.name])),
            percent.derived.cont=as.data.frame(round(t(theo.cont[,sorted.name]/nplus*100),2)))

  stopImplicitCluster()
  class(back)=c("sensory.mr.sig.cell","list")
  return(back)
}
