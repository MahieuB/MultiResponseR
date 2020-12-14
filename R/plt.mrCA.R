#' Plot factor plan resulting from multiple-response Correspondence Analysis (mrCA)
#'
#' @description This function plots the results coming from \code{\link[MultiResponseR]{sensory.mrCA}} or \code{\link[MultiResponseR]{mrCA}}
#'
#' @param res A list returned by \code{\link[MultiResponseR]{sensory.mrCA}} or \code{\link[MultiResponseR]{mrCA}}
#' @param axes Which dimensions of the mrCA should be plotted?
#' @param alpha.total.bootstrap.test The alpha risk of the total bootstrap tests. Only useful if the mrCA was computed using \code{\link[MultiResponseR]{sensory.mrCA}} or \code{\link[MultiResponseR]{mrCA}} and ellipse=TRUE. See details
#' @param alpha.ellipse The alpha risk of the confidence ellipses. Only useful if the mrCA was computed using \code{\link[MultiResponseR]{sensory.mrCA}} or \code{\link[MultiResponseR]{mrCA}} and ellipse=TRUE
#' @param select.desc.rep A character vector specifying the descriptors/response options to plot. By default, all descriptors/response options are plotted
#' @param rev.x Should the horizontal plotted dimension be reversed? Useful in case of map comparisons to align products/categories
#' @param rev.y Should the vertical plotted dimension be reversed? Useful in case of map comparisons to align products/categories
#' @param size.points The size of the points used to represent the products/categories on the map
#' @param size.lab The size of the label on the map
#' @param expansion The factor of expansion applied to descriptors/response options coordinates to increase readability
#' @param title An optional title to be added to the plot
#'
#' @details
#' \itemize{
#'   \item \strong{alpha.total.bootstrap.test}: Products/categories non-significantly different at the alpha risk of \emph{alpha.total.bootstrap.test} according to the total bootstrap test are linked by a line on the plot. If these links are not required, \emph{alpha.total.bootstrap.test} can be set to 1
#' }
#'
#'
#' @return An mrCA factor map
#'
#' @import FactoMineR
#' @import graphics
#' @import stats
#'
#' @return
#' @export
#'
#' @examples
#' # non-sensory example
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
#' res=mrCA(dset)
#'
#' plt.mrCA(res)
#'
#' # sensory example
#'data(milkchoc)
#'
#'parallel::detectCores()
#'
#'dim.sig=sensory.mr.dimensionality.test(milkchoc)$dim.sig
#'
#'res=sensory.mrCA(milkchoc,nbaxes.sig=dim.sig)
#'
#'plt.mrCA(res)
plt.mrCA=function(res,axes=c(1,2),alpha.total.bootstrap.test=0.05,alpha.ellipse=alpha.total.bootstrap.test,
                  select.desc.rep=rownames(res$col.coord),rev.x=FALSE,rev.y=FALSE,
                  size.points=1.4,size.lab=1.4,expansion=1.25,title=NULL){
  classe=class(res)
  if (classe!="list"){
    stop("res must be a list resulting from the execution of sensory.mrCA or mrCA")
  }
  taille=length(res)
  if (taille!=8){
    stop("res must be a list resulting from the execution of sensory.mrCA or mrCA")
  }
  if (is.null(res$bootstrap.replicate.coord)){
    check.axes=nrow(res$eigen)
    if (max(axes)>check.axes){
      stop("max(axes) must be lower than or equal to the minimum between the number of categories minus one and the number of response options")
    }
  }else{
    check.axes=ncol(res$bootstrap.replicate.coord)-1
    if (max(axes)>check.axes){
      stop("max(axes) must be lower than or equal to the number of nbaxes.sig used in the mrCA")
    }
  }
  classe=class(select.desc.rep)
  if (classe!="character"){
    stop("select.desc.rep must be a character vector")
  }

  if (!is.null(res$bootstrap.replicate.coord)){
    ell=coord.ellipse(res$bootstrap.replicate.coord,axes=axes,level.conf = (1-alpha.ellipse))$res
  }

  adjusted.col.coord=res$col.coord

  if (!is.null(select.desc.rep)){
    adjusted.col.coord=adjusted.col.coord[select.desc.rep,,drop=FALSE]
  }

  if (!is.null(res$proj.row.coord)){
    max.norm.prod = max(abs(res$row.coord[,axes]),abs(res$proj.row.coord[,axes]))
  }else{
    max.norm.prod = max(abs(res$row.coord[,axes]))
  }

  if (!is.null(res$proj.col.coord)){
    max.norm.desc = max(abs(adjusted.col.coord[,axes]),abs(res$proj.col.coord[,axes]))
  }else{
    max.norm.desc = max(abs(adjusted.col.coord[,axes]))
  }

  expand = max.norm.prod/max.norm.desc*expansion
  adjusted.col.coord=adjusted.col.coord*expand
  if (!is.null(res$proj.col.coord)){
    res$proj.col.coord=res$proj.col.coord*expand
  }


  if(rev.x){
    res$row.coord[,axes[1]]=-res$row.coord[,axes[1]]
    adjusted.col.coord[,axes[1]]=-adjusted.col.coord[,axes[1]]
    ell[,2]=-ell[,2]
    if (!is.null(res$proj.col.coord)){
      res$proj.col.coord=-res$proj.col.coord[,axes[1]]
    }
    if (!is.null(res$proj.row.coord)){
      res$proj.row.coord=-res$proj.row.coord[,axes[1]]
    }
  }
  if(rev.y){
    res$row.coord[,axes[2]]=-res$row.coord[,axes[2]]
    adjusted.col.coord[,axes[2]]=-adjusted.col.coord[,axes[2]]
    ell[,3]=-ell[,3]
    if (!is.null(res$proj.col.coord)){
      res$proj.col.coord=-res$proj.col.coord[,axes[2]]
    }
    if (!is.null(res$proj.row.coord)){
      res$proj.row.coord=-res$proj.row.coord[,axes[2]]
    }
  }

  xmin=c(min(res$row.coord[,axes[1]],adjusted.col.coord[,axes[1]]))
  xmax=c(max(res$row.coord[,axes[1]],adjusted.col.coord[,axes[1]]))
  ymin=c(min(res$row.coord[,axes[2]],adjusted.col.coord[,axes[2]]))
  ymax=c(max(res$row.coord[,axes[2]],adjusted.col.coord[,axes[2]]))

  if (!is.null(res$bootstrap.replicate.coord)){
    xmin=min(xmin,min(ell[,2]))
    xmax=max(xmax,max(ell[,2]))
    ymin=min(ymin,min(ell[,3]))
    ymax=max(ymax,max(ell[,3]))
  }

  if (!is.null(res$proj.row.coord)){
    xmin=min(xmin,min(res$proj.row.coord[,axes[1]]))
    xmax=max(xmax,max(res$proj.row.coord[,axes[1]]))
    ymin=min(ymin,min(res$proj.row.coord[,axes[2]]))
    ymax=max(ymax,max(res$proj.row.coord[,axes[2]]))
  }

  if (!is.null(res$proj.col.coord)){
    xmin=min(xmin,min(res$proj.col.coord[,axes[1]]))
    xmax=max(xmax,max(res$proj.col.coord[,axes[1]]))
    ymin=min(ymin,min(res$proj.col.coord[,axes[2]]))
    ymax=max(ymax,max(res$proj.col.coord[,axes[2]]))
  }

  xmin=xmin*1.1
  xmax=xmax*1.1
  ymin=ymin*1.1
  ymax=ymax*1.1

  plot(res$row.coord[,axes],cex=size.points,col="blue",pch=16,
       xlim=c(xmin,xmax),ylim=c(ymin,ymax),
       xlab=paste("Dim ",axes[1]," (",round(res$eigen[axes[1],2],2)," %)",sep=""),
       ylab=paste("Dim ",axes[2]," (",round(res$eigen[axes[2],2],2)," %)",sep=""),
       main=title,cex.main=2,cex.lab=1.4,xaxt="n", yaxt="n")
  abline(h=0,v=0,lty="dashed")

  if (!is.null(res$bootstrap.replicate.coord)){
    for (cate in levels(ell[,1])){
      toplot=ell[ell[,1]==cate,]
      lines(toplot[,2],toplot[,3],type="l",col="blue",lwd=2)
    }
    diff.test=res$total.bootstrap.test.pvalues
    for (i in 1:nrow(diff.test)){
      for (j in i:ncol(diff.test)){
        p.1=rownames(diff.test)[i]
        p.2=colnames(diff.test)[j]
        if (diff.test[p.1,p.2]>alpha.total.bootstrap.test & i!=j){
          ac.produit.coord=as.data.frame(res$row.coord[,axes])
          p.1.coord=ac.produit.coord[p.1,]
          p.2.coord=ac.produit.coord[p.2,]
          segments(p.1.coord[,1],p.1.coord[,2],p.2.coord[,1],p.2.coord[,2],lwd=2,lty="solid",col="blue")
        }
      }
    }
  }
  if (!is.null(select.desc.rep)){
    toplot=as.matrix(adjusted.col.coord[,axes,drop=FALSE])
    for (i in 1:nrow(toplot)){
      tox=toplot[i,1]
      toy=toplot[i,2]
      arrows(0,0,tox,toy,length = 0.15,col="red",lwd=2)
    }
  }

  if(!is.null(res$proj.col.coord)){
    toplot=as.matrix(res$proj.col.coord[,axes,drop=FALSE])
    for (i in 1:nrow(toplot)){
      tox=toplot[i,1]
      toy=toplot[i,2]
      arrows(0,0,tox,toy,length = 0.15,col="red4",lwd=2)
    }
  }

  points(res$row.coord[,axes],cex=size.points,col="blue",pch=16)
  if(!is.null(res$proj.row.coord)){
    for (i in 1:nrow(res$proj.row.coord)){
      points(res$proj.row.coord[i,axes[1]],res$proj.row.coord[i,axes[2]],cex=size.points,col="darkblue",pch=16)
    }
  }

  if (!is.null(select.desc.rep)){
    lab.desc=as.matrix(adjusted.col.coord[,axes,drop=FALSE])
    rownames(lab.desc)=rownames(adjusted.col.coord)
    lab=rbind(res$row.coord[,axes],lab.desc)
    col.lab=c(rep("blue",nrow(res$row.coord[,axes])),rep("red",nrow(adjusted.col.coord)))

    if (!is.null(res$proj.col.coord)){
      lab.desc.sup=as.matrix(res$proj.col.coord[,axes,drop=FALSE])
      lab=rbind(lab,lab.desc.sup)
      col.lab=c(col.lab,rep("red4",nrow(res$proj.col.coord)))
    }

    if (!is.null(res$proj.row.coord)){
      lab=rbind(lab,res$proj.row.coord[,axes,drop=FALSE])
      col.lab=c(col.lab,rep("darkblue",nrow(res$proj.row.coord)))
    }

  }else{
    lab=res$row.coord[,axes]
    col.lab=rep("blue",nrow(res$row.coord[,axes]))

    if (!is.null(res$proj.col.coord)){
      lab.desc.sup=as.matrix(res$proj.col.coord[,axes,drop=FALSE])
      lab=rbind(lab,lab.desc.sup)
      col.lab=c(col.lab,rep("red4",nrow(res$proj.col.coord)))
    }

    if (!is.null(res$proj.row.coord)){
      lab=rbind(lab,res$proj.row.coord[,axes,drop=FALSE])
      col.lab=c(col.lab,rep("darkblue",nrow(res$proj.row.coord)))
    }

  }
  autoLab(lab,rownames(lab),col=col.lab,cex=size.lab,shadotext = TRUE)
}
