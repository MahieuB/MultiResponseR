#' Multi-response Correspondence Analysis (mrCA) for sensory data
#'
#' @description This function performs the mrCA of the data as well as the total bootstrap procedure (Cadoret & Husson, 2013) and the pairwise total bootstrap tests (Mahieu, Visalli, Thomas, & Schlich, 2020). The difference with \code{\link[MultiResponseR]{mrCA}} used with ellipse=TRUE is that the total bootstrap procedure takes into account the subject structure of sensory data in \code{\link[MultiResponseR]{sensory.mrCA}}
#'
#' @param data A data.frame of evaluations in rows whose first two columns are factors (subject and product) and subsequent columns are binary numeric or integer, each column being a descriptor
#' @param nboot The number of bootstrapped panel of the total bootstrap procedure
#' @param nbaxes.sig The number of significant axes retuned by \code{\link[MultiResponseR]{sensory.mr.dimensionality.test}}. By default, all axes are considered significant. See details
#' @param ncores Number of cores used to generate the virtual panels. Default is 2. See details
#'
#' @details
#' \itemize{
#'   \item \strong{nbaxes.sig}: The number of significant axes determines the number of axes accounted for while performing the Procrustes rotations of the total bootstrap procedure (Mahieu, Visalli, & Schlich, 2020). These same axes are accounted for the pairwise total bootstrap tests.
#'   \item \strong{ncores}: The more cores are added in the process, the faster the results will be obtained. The number of available cores is accessible using \code{\link[parallel]{detectCores}}. The parallel tasks are closed once the \emph{nboot} datasets are generated.
#' }
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{eigen}{Eigenvalues of the mrCA and their corresponding percentages of inertia}
#'   \item{row.coord}{Products coordinates}
#'   \item{col.coord}{Descriptors coordinates}
#'   \item{svd}{Results of the singular value decomposition}
#'   \item{bootstrap.replicate.coord}{Coordinates of the rotated bootstrap replicates}
#'   \item{total.bootstrap.test.pvalues}{P-values of the pairwise total bootstrap tests}
#' }
#' @export
#'
#' @references Cadoret, M., & Husson, F. (2013). Construction and evaluation of confidence ellipses applied at sensory data. Food Quality and Preference, 28(1), 106-115.
#' @references Mahieu, B., Visalli, M., Thomas, A., & Schlich, P. (2020). Free-comment outperformed check-all-that-apply in the sensory characterisation of wines with consumers at home. Food Quality and Preference, 84.
#' @references Mahieu, B., Visalli, M., & Schlich, P. (2020). Accounting for the dimensionality of the dependence in analyses of contingency tables obtained with Check-All-That-Apply and Free-Comment. Food Quality and Preference, 83.
#' @references Mahieu, B., Schlich, P., Visalli, M., & Cardot, H. (2020). A multiple-response chi-square framework for the analysis of Free-Comment and Check-All-That-Apply data. Manuscript submitted for publication.
#'
#' @import parallel
#' @import doParallel
#' @import foreach
#' @import iterators
#' @importFrom candisc candisc
#' @import stats
#'
#' @examples
#'
#'data(milkchoc)
#'
#'parallel::detectCores()
#'
#'dim.sig=sensory.mr.dimensionality.test(milkchoc)$dim.sig
#'
#'res=sensory.mrCA(milkchoc,nbaxes.sig=dim.sig)
#'
#'plt.mrCA(res)
sensory.mrCA=function(data,nboot=2000,nbaxes.sig=Inf,ncores=2){
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
  if (length(unique(param.etendu))>1){
    warning("Data are unbalanced: products have not been evaluated a same number of times by subjects")
  }
  nplus=colSums(param.etendu)
  redui=org

  ca.actual=mrCA(data[,-1])
  actual.coord=ca.actual$row.coord[,1:nbaxes.proc]

  registerDoParallel(cores = ncores)
  nboot=nboot
  nvirtual=nlevels(data$sujet)

  sortie <- foreach(icount(nboot), .combine='rbind') %dopar% {

    pioche=unique(data$sujet)
    nb.tirage=length(pioche)
    nb.panel.virtuel=nvirtual
    loto=sample(1:nb.tirage,nb.panel.virtuel,replace = TRUE)
    tirage=pioche[loto]

    vec.ligne=NULL
    for (sujet.tirage in tirage){
      nouvelle.ligne=which(data$sujet==sujet.tirage)
      vec.ligne=c(vec.ligne,nouvelle.ligne)
    }
    jdd.tirage=data[vec.ligne,]
    rownames(jdd.tirage)=as.character(1:nrow(jdd.tirage))

    param.etendu.tirage=table(jdd.tirage$sujet,jdd.tirage$produit)
    nplus.tirage=colSums(param.etendu.tirage)

    cont.tirage = aggregate(.~produit,jdd.tirage,sum)
    rownames(cont.tirage)=as.character(cont.tirage$produit)
    cont.tirage$produit=cont.tirage$sujet=NULL
    redui.tirage=cont.tirage

    ######

    procrustes <- function(amat, target, orthogonal = FALSE,
                           translate = FALSE, magnify = FALSE) {
      for (i in nrow(amat):1) {
        if (any(is.na(amat)[i, ]) | any(is.na(target)[i,
        ])) {
          amat <- amat[-i, ]
          target <- target[-i, ]
        }
      }
      dA <- dim(amat)
      dX <- dim(target)
      if (length(dA) != 2 || length(dX) != 2)
        stop("arguments amat and target must be matrices")
      if (any(dA != dX))
        stop("dimensions of amat and target must match")
      if (length(attr(amat, "tmat")))
        stop("oblique loadings matrix not allowed for amat")
      if (orthogonal) {
        if (translate) {
          p <- dX[1]
          target.m <- (rep(1/p, p) %*% target)[, ]
          amat.m <- (rep(1/p, p) %*% amat)[, ]
          target.c <- scale(target, center = target.m,
                            scale = FALSE)
          amat.c <- scale(amat, center = amat.m, scale = FALSE)
          j <- svd(crossprod(target.c, amat.c))
        }
        else {
          amat.c <- amat
          j <- svd(crossprod(target, amat))
        }
        rot <- j$v %*% t(j$u)
        if (magnify)
          beta <- sum(j$d)/sum(amat.c^2)
        else beta <- 1
        B <- beta * amat.c %*% rot
        if (translate)
          B <- B + rep(as.vector(target.m), rep.int(p,
                                                    dX[2]))
        value <- list(rmat = B, tmat = rot, magnify = beta)
        if (translate)
          value$translate <- target.m - (rot %*% amat.m)[,
          ]
      }
      else {
        b <- solve(amat, target)
        gamma <- sqrt(diag(solve(crossprod(b))))
        rot <- b * rep(gamma, rep.int(dim(b)[1], length(gamma)))
        B <- amat %*% rot
        fcor <- solve(crossprod(rot))
        value <- list(rmat = B, tmat = rot, correlation = fcor)
      }
      return(value)
    }
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
      mat.eig=cbind(eig,percent.eig)
      colnames(mat.eig)=c("eigenvalue","percentage of inertia")
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
            loto=sample(les.ligne,length(les.ligne),replace = TRUE)
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

          procrustes <- function(amat, target, orthogonal = FALSE,
                                 translate = FALSE, magnify = FALSE) {
            for (i in nrow(amat):1) {
              if (any(is.na(amat)[i, ]) | any(is.na(target)[i,
              ])) {
                amat <- amat[-i, ]
                target <- target[-i, ]
              }
            }
            dA <- dim(amat)
            dX <- dim(target)
            if (length(dA) != 2 || length(dX) != 2)
              stop("arguments amat and target must be matrices")
            if (any(dA != dX))
              stop("dimensions of amat and target must match")
            if (length(attr(amat, "tmat")))
              stop("oblique loadings matrix not allowed for amat")
            if (orthogonal) {
              if (translate) {
                p <- dX[1]
                target.m <- (rep(1/p, p) %*% target)[, ]
                amat.m <- (rep(1/p, p) %*% amat)[, ]
                target.c <- scale(target, center = target.m,
                                  scale = FALSE)
                amat.c <- scale(amat, center = amat.m, scale = FALSE)
                j <- svd(crossprod(target.c, amat.c))
              }
              else {
                amat.c <- amat
                j <- svd(crossprod(target, amat))
              }
              rot <- j$v %*% t(j$u)
              if (magnify)
                beta <- sum(j$d)/sum(amat.c^2)
              else beta <- 1
              B <- beta * amat.c %*% rot
              if (translate)
                B <- B + rep(as.vector(target.m), rep.int(p,
                                                          dX[2]))
              value <- list(rmat = B, tmat = rot, magnify = beta)
              if (translate)
                value$translate <- target.m - (rot %*% amat.m)[,
                ]
            }
            else {
              b <- solve(amat, target)
              gamma <- sqrt(diag(solve(crossprod(b))))
              rot <- b * rep(gamma, rep.int(dim(b)[1], length(gamma)))
              B <- amat %*% rot
              fcor <- solve(crossprod(rot))
              value <- list(rmat = B, tmat = rot, correlation = fcor)
            }
            return(value)
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

          rot = procrustes(row.coord.tirage, row.coord, orthogonal = T, translate = T, magnify = F)$rmat
          rot=cbind.data.frame(rownames(rot),rot)
          colnames(rot)[1]="cat"
          return(rot)
        }

        sortie$cat=as.factor(sortie$cat)
        toellipse=sortie
        rownames(toellipse)=as.character(1:nrow(toellipse))
        centre=row.coord[,1:nbaxes.proc]
        centre=cbind(as.data.frame(rownames(centre)),centre)
        colnames(centre)=colnames(toellipse)
        rownames(centre)=as.character(1:nrow(centre))
        deltapos=centre
        deltapos[,2:ncol(deltapos)]=deltapos[,2:ncol(deltapos)]-aggregate(.~cat,toellipse,mean)[,2:ncol(deltapos)]
        for (p.delta in unique(toellipse$cat)){
          choix.ligne=rownames(toellipse[toellipse$cat==p.delta,])
          choix.vec=deltapos[deltapos$cat==p.delta,][,2:ncol(deltapos)]
          choix.vec.dup=as.data.frame(as.matrix(rep(1,nboot))%*%as.matrix(choix.vec))
          toellipse[choix.ligne,2:ncol(toellipse)]=toellipse[choix.ligne,2:ncol(toellipse)]+choix.vec.dup
        }
        toellipse$cat=as.factor(toellipse$cat)
        centre$cat=as.factor(centre$cat)

        coord.boot=toellipse[,1:(nbaxes.sig+1)]
        dim.sig=nbaxes.sig

        les.prod=unique(coord.boot$cat)

        diff.test=as.data.frame(matrix(1,length(les.prod),length(les.prod)))
        colnames(diff.test)=rownames(diff.test)=les.prod

        for(i in 1:nrow(diff.test)){
          for (j in 1:ncol(diff.test)){
            if (dim.sig > 1){
              if (i != j){
                p.1=rownames(diff.test)[i]
                p.2=colnames(diff.test)[j]
                don.coord=coord.boot[coord.boot$cat==p.1 | coord.boot$cat==p.2,]
                don.coord$cat=as.factor(as.character(don.coord$cat))
                mod.man=lm(as.matrix(don.coord[,2:ncol(don.coord)])~cat,don.coord)
                res.can=candisc(mod.man,"cat",type = "3")
                proj.coord=res.can$scores
                X=proj.coord[proj.coord$cat==p.1,]
                Y=proj.coord[proj.coord$cat==p.2,]
                X$cat=Y$cat=NULL
                Z=round(X-Y,12)
                gauche=length(which(Z<=0))
                droite=length(which(Z>=0))
                pvalue.test=((min(gauche,droite)+1)/(nboot+1))*2
                if (pvalue.test>1){
                  pvalue.test=1
                }
                diff.test[p.1,p.2]=pvalue.test
              }
            }
            if (dim.sig==1){
              if (i != j){
                p.1=rownames(diff.test)[i]
                p.2=colnames(diff.test)[j]
                don.coord=coord.boot[coord.boot$cat==p.1 | coord.boot$cat==p.2,]
                don.coord$cat=as.factor(as.character(don.coord$cat))
                X=don.coord[don.coord$cat==p.1,]
                Y=don.coord[don.coord$cat==p.2,]
                X$cat=Y$cat=NULL
                Z=round(X-Y,12)
                gauche=length(which(Z<=0))
                droite=length(which(Z>=0))
                pvalue.test=((min(gauche,droite)+1)/(nboot+1))*2
                if (pvalue.test>1){
                  pvalue.test=1
                }
                diff.test[p.1,p.2]=pvalue.test
              }
            }
          }
        }


      }else{
        toellipse=NULL
        diff.test=NULL
      }
      back=list(eigen=mat.eig,row.coord=row.coord,col.coord=col.coord,proj.row.coord=proj.row.coord,proj.col.coord=proj.col.coord,svd=list(u=u,vs=vs,v=v),bootstrap.replicate.coord=toellipse,total.bootstrap.test.pvalues=diff.test)
      return(back)
    }


    ######
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

    ca.tirage=mrCA(jdd.tirage[,-1])
    tirage.coord=ca.tirage$row.coord[,1:nbaxes.proc]
    rot = procrustes(tirage.coord, actual.coord, orthogonal = T, translate = T, magnify = F)$rmat
    rot=cbind.data.frame(rownames(rot),rot)
    colnames(rot)[1]="produit"
    return(rot)
  }

  sortie$produit=as.factor(sortie$produit)
  toellipse=sortie
  rownames(toellipse)=as.character(1:nrow(toellipse))
  centre=ca.actual$row.coord[,1:nbaxes.proc]
  centre=cbind(as.data.frame(rownames(centre)),centre)
  colnames(centre)=colnames(toellipse)
  rownames(centre)=as.character(1:nrow(centre))
  deltapos=centre
  deltapos[,2:ncol(deltapos)]=deltapos[,2:ncol(deltapos)]-aggregate(.~produit,toellipse,mean)[,2:ncol(deltapos)]
  for (p.delta in unique(toellipse$produit)){
    choix.ligne=rownames(toellipse[toellipse$produit==p.delta,])
    choix.vec=deltapos[deltapos$produit==p.delta,][,2:ncol(deltapos)]
    choix.vec.dup=as.data.frame(as.matrix(rep(1,nboot))%*%as.matrix(choix.vec))
    toellipse[choix.ligne,2:ncol(toellipse)]=toellipse[choix.ligne,2:ncol(toellipse)]+choix.vec.dup
  }
  toellipse$produit=as.factor(toellipse$produit)
  centre$produit=as.factor(centre$produit)

  coord.boot=toellipse[,1:(nbaxes.sig+1)]
  dim.sig=nbaxes.sig

  les.prod=unique(coord.boot$produit)

  diff.test=as.data.frame(matrix(1,length(les.prod),length(les.prod)))
  colnames(diff.test)=rownames(diff.test)=les.prod

  for(i in 1:nrow(diff.test)){
    for (j in 1:ncol(diff.test)){
      if (dim.sig > 1){
        if (i != j){
          p.1=rownames(diff.test)[i]
          p.2=colnames(diff.test)[j]
          don.coord=coord.boot[coord.boot$produit==p.1 | coord.boot$produit==p.2,]
          don.coord$produit=as.factor(as.character(don.coord$produit))
          mod.man=lm(as.matrix(don.coord[,2:ncol(don.coord)])~produit,don.coord)
          res.can=candisc(mod.man,"produit",type = "3")
          proj.coord=res.can$scores
          X=proj.coord[proj.coord$produit==p.1,]
          Y=proj.coord[proj.coord$produit==p.2,]
          X$produit=Y$produit=NULL
          Z=round(X-Y,12)
          gauche=length(which(Z<=0))
          droite=length(which(Z>=0))
          pvalue.test=((min(gauche,droite)+1)/(nboot+1))*2
          if (pvalue.test>1){
            pvalue.test=1
          }
          diff.test[p.1,p.2]=pvalue.test
        }
      }
      if (dim.sig==1){
        if (i != j){
          p.1=rownames(diff.test)[i]
          p.2=colnames(diff.test)[j]
          don.coord=coord.boot[coord.boot$produit==p.1 | coord.boot$produit==p.2,]
          don.coord$produit=as.factor(as.character(don.coord$produit))
          X=don.coord[don.coord$produit==p.1,]
          Y=don.coord[don.coord$produit==p.2,]
          X$produit=Y$produit=NULL
          Z=round(X-Y,12)
          gauche=length(which(Z<=0))
          droite=length(which(Z>=0))
          pvalue.test=((min(gauche,droite)+1)/(nboot+1))*2
          if (pvalue.test>1){
            pvalue.test=1
          }
          diff.test[p.1,p.2]=pvalue.test
        }
      }
    }
  }

  back=list(eigen=ca.actual$eigen,row.coord=ca.actual$row.coord,col.coord=ca.actual$col.coord,
            proj.row.coord=ca.actual$proj.row.coord,proj.col.coord=ca.actual$proj.col.coord,
            svd=ca.actual$svd,bootstrap.replicate.coord=toellipse,total.bootstrap.test.pvalues=diff.test)
  stopImplicitCluster()
  return(back)
}
