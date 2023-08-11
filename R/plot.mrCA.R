#' Plot factor plan resulting from multiple-response Correspondence Analysis (MR-CA)
#'
#' @description This function plots the results coming from \code{\link[MultiResponseR]{mrCA}}
#'
#' @param x A list returned by \code{\link[MultiResponseR]{mrCA}}
#' @param axes Which dimensions of the MR-CA should be plotted?
#' @param alpha.total.bootstrap.test The alpha risk of the total bootstrap tests. Only useful if the MR-CA was computed using \code{\link[MultiResponseR]{mrCA}} and ellipse=TRUE. See details
#' @param alpha.ellipse The alpha risk of the confidence ellipses. Only useful if the MR-CA was computed using \code{\link[MultiResponseR]{mrCA}} and ellipse=TRUE
#' @param select.rep A character vector specifying the response options to plot. By default, all response options are plotted
#' @param rev.x Should the horizontal plotted dimension be reversed? Useful in case of map comparisons to align categories
#' @param rev.y Should the vertical plotted dimension be reversed? Useful in case of map comparisons to align categories
#' @param size.points The size of the points used to represent the categories on the map
#' @param size.lab The size of the label on the map
#' @param expansion The factor of expansion applied to response options coordinates to increase readability
#' @param title An optional title to be added to the plot
#' @param ... further arguments passed to or from other methods
#'
#' @details
#' \itemize{
#'   \item \strong{alpha.total.bootstrap.test}: Categories non-significantly different at the alpha risk of \emph{alpha.total.bootstrap.test} according to the total bootstrap test are linked by a line on the plot. If these links are not required, \emph{alpha.total.bootstrap.test} can be set to 1
#' }
#'
#'
#' @return A MR-CA factor map
#'
#' @export
#'
#' @import ggplot2
#' @import ggrepel
#' @importFrom ellipse ellipse
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
#' res=mrCA(dset)
#'
#' plot(res)
plot.mrCA=function(x,
                   axes = c(1,2),
                   alpha.total.bootstrap.test = 0.05,
                   alpha.ellipse = alpha.total.bootstrap.test,
                   select.rep = rownames(x$col.coord),
                   rev.x = FALSE,
                   rev.y = FALSE,
                   size.points = 3.5,
                   size.lab = 6,
                   expansion = 1.25,
                   title = NULL,...){
  if (!inherits(x,"mrCA")){
    stop("class(x) must be mrCA")
  }
  if (is.null(x$bootstrap.replicate.coord)){
    check.axes=nrow(x$eigen)
    if (max(axes)>check.axes){
      stop("max(axes) must be lower than or equal to the minimum between the number of categories minus one and the number of response options")
    }
  }else{
    check.axes=ncol(x$bootstrap.replicate.coord)-1
    if (max(axes)>check.axes){
      stop("max(axes) must be lower than or equal to the number of nbaxes.sig used in the mrCA")
    }
  }
  classe=class(select.rep)
  if (classe!="character" & !is.null(select.rep)){
    stop("class(select.rep) must be character or NULL")
  }

  if (!is.null(x$bootstrap.replicate.coord)){
    ell=data.frame(cat=as.factor(rep(levels(x$bootstrap.replicate.coord$cat),each=100)),matrix(0,100*nlevels(x$bootstrap.replicate.coord$cat),2))
    colnames(ell)[2:3]=paste("Dim.",1:2,sep = "")
    for (c in levels(ell$cat)){
      boot.rep.c=x$bootstrap.replicate.coord[x$bootstrap.replicate.coord$cat==c,-1]
      boot.rep.c=as.matrix(boot.rep.c)
      sigma=cov(boot.rep.c)
      mu=colMeans(boot.rep.c)
      calc.mal.sq=function(vec){
        mal.bary=t(as.matrix(vec-mu))%*%solve(sigma,tol=1e-300)%*%(as.matrix(vec-mu))
        return(as.numeric(mal.bary))
      }
      mal.sq.cloud=apply(boot.rep.c,1,calc.mal.sq)
      dilat.stat=sqrt(quantile(mal.sq.cloud,1-alpha.ellipse,type=2))
      ell.c=ellipse::ellipse(sigma[c(axes[1],axes[2]),c(axes[1],axes[2])],centre=mu[c(axes[1],axes[2])],t=dilat.stat)
      ell[which(ell$cat==c),2:3]=ell.c
    }
  }

  adjusted.col.coord=x$col.coord

  if (!is.null(select.rep)){
    adjusted.col.coord=adjusted.col.coord[select.rep,,drop=FALSE]
  }

  if (!is.null(x$proj.row.coord)){
    max.norm.prod = max(abs(x$row.coord[,axes]),abs(x$proj.row.coord[,axes]))
  }else{
    max.norm.prod = max(abs(x$row.coord[,axes]))
  }

  if (!is.null(x$proj.col.coord)){
    max.norm.desc = max(abs(adjusted.col.coord[,axes]),abs(x$proj.col.coord[,axes]))
  }else{
    max.norm.desc = max(abs(adjusted.col.coord[,axes]))
  }

  expand = max.norm.prod/max.norm.desc*expansion
  adjusted.col.coord=adjusted.col.coord*expand
  if (!is.null(x$proj.col.coord)){
    x$proj.col.coord=x$proj.col.coord*expand
  }


  if(rev.x){
    x$row.coord[,axes[1]]=-x$row.coord[,axes[1]]
    adjusted.col.coord[,axes[1]]=-adjusted.col.coord[,axes[1]]
    ell[,2]=-ell[,2]
    if (!is.null(x$proj.col.coord)){
      x$proj.col.coord[,axes[1]]=-x$proj.col.coord[,axes[1]]
    }
    if (!is.null(x$proj.row.coord)){
      x$proj.row.coord[,axes[1]]=-x$proj.row.coord[,axes[1]]
    }
  }
  if(rev.y){
    x$row.coord[,axes[2]]=-x$row.coord[,axes[2]]
    adjusted.col.coord[,axes[2]]=-adjusted.col.coord[,axes[2]]
    ell[,3]=-ell[,3]
    if (!is.null(x$proj.col.coord)){
      x$proj.col.coord[,axes[2]]=-x$proj.col.coord[,axes[2]]
    }
    if (!is.null(x$proj.row.coord)){
      x$proj.row.coord[,axes[2]]=-x$proj.row.coord[,axes[2]]
    }
  }

  xmin=c(min(x$row.coord[,axes[1]],adjusted.col.coord[,axes[1]]))
  xmax=c(max(x$row.coord[,axes[1]],adjusted.col.coord[,axes[1]]))
  ymin=c(min(x$row.coord[,axes[2]],adjusted.col.coord[,axes[2]]))
  ymax=c(max(x$row.coord[,axes[2]],adjusted.col.coord[,axes[2]]))

  if (!is.null(x$bootstrap.replicate.coord)){
    xmin=min(xmin,min(ell[,2]))
    xmax=max(xmax,max(ell[,2]))
    ymin=min(ymin,min(ell[,3]))
    ymax=max(ymax,max(ell[,3]))
  }

  if (!is.null(x$proj.row.coord)){
    xmin=min(xmin,min(x$proj.row.coord[,axes[1]]))
    xmax=max(xmax,max(x$proj.row.coord[,axes[1]]))
    ymin=min(ymin,min(x$proj.row.coord[,axes[2]]))
    ymax=max(ymax,max(x$proj.row.coord[,axes[2]]))
  }

  if (!is.null(x$proj.col.coord)){
    xmin=min(xmin,min(x$proj.col.coord[,axes[1]]))
    xmax=max(xmax,max(x$proj.col.coord[,axes[1]]))
    ymin=min(ymin,min(x$proj.col.coord[,axes[2]]))
    ymax=max(ymax,max(x$proj.col.coord[,axes[2]]))
  }



  pmin=min(xmin,ymin)*1.05
  pmax=max(xmax,ymax)*1.05

  p=ggplot(as.data.frame(x$row.coord),aes(x=x$row.coord[,axes[1]],y=x$row.coord[,axes[2]]))+theme_bw()
  p=p+xlim(pmin,pmax)+ylim(pmin,pmax)+xlab(paste("Dim ",axes[1]," (",round(x$eigen[axes[1],2],2)," %)",sep=""))+ylab(paste("Dim ",axes[2]," (",round(x$eigen[axes[2],2],2)," %)",sep=""))+ggtitle(title)
  p=p+theme(axis.title.x = element_text(size = 16,face = "bold"),axis.title.y = element_text(size = 16,face = "bold"),plot.title = element_text(hjust = 0.5,face = "bold",size=20))
  p=p+geom_hline(yintercept=0,linetype="dashed",linewidth=1)+geom_vline(xintercept=0,linetype="dashed",linewidth=1)

  if (!is.null(x$bootstrap.replicate.coord)){

    p=p+geom_path(data=as.data.frame(ell),aes(x=ell[,2],y=ell[,3],group=ell[,1]),colour="blue",linewidth=1)

    diff.test=x$total.bootstrap.test.pvalues

    df.segment=NULL

    for (i in 1:nrow(diff.test)){
      for (j in i:ncol(diff.test)){
        p.1=rownames(diff.test)[i]
        p.2=colnames(diff.test)[j]
        if (diff.test[p.1,p.2]>alpha.total.bootstrap.test & i!=j){
          ac.produit.coord=as.data.frame(x$row.coord[,axes])
          p.1.coord=ac.produit.coord[p.1,]
          p.2.coord=ac.produit.coord[p.2,]
          sous.df.segment=cbind(p.1.coord,p.2.coord)
          df.segment=rbind(df.segment,sous.df.segment)
        }
      }
    }
    if (!is.null(df.segment)){
      colnames(df.segment)=as.character(1:ncol(df.segment))
      p=p+geom_segment(data=as.data.frame(df.segment),aes(x = df.segment[,1], y = df.segment[,2], xend = df.segment[,3], yend = df.segment[,4]),colour="blue",linewidth=1.3)
    }
  }

  if (!is.null(select.rep)){
    df.fleche=as.matrix(adjusted.col.coord[,axes,drop=FALSE])
    col.fleche=rep("red",length(select.rep))
    if(!is.null(x$proj.col.coord)){
      df.fleche=rbind(df.fleche,as.matrix(x$proj.col.coord[,axes,drop=FALSE]))
      col.fleche=c(col.fleche,rep("red4",nrow(x$proj.col.coord)))
    }
    p=p+geom_segment(data = as.data.frame(df.fleche), aes(x=0, y=0,xend = df.fleche[,1], yend = df.fleche[,2]), arrow=arrow(length = unit(0.4, "cm"),type = "closed"), colour=col.fleche,linewidth=1)
  }

  
  df.point=x$row.coord[,axes,drop=FALSE]
  col.point=rep("blue",nrow(x$row.coord))
  if(!is.null(x$proj.row.coord)){
    df.point=rbind(df.point,x$proj.row.coord[,axes,drop=FALSE])
    col.point=c(col.point,rep("darkblue",nrow(x$proj.row.coord)))
  }
  p=p+geom_point(data=as.data.frame(df.point),aes(x=df.point[,1],y=df.point[,2]),colour=col.point,size=size.points)

  if (!is.null(select.rep)){
    lab.desc=as.matrix(adjusted.col.coord[,axes,drop=FALSE])
    rownames(lab.desc)=rownames(adjusted.col.coord)
    lab=rbind(x$row.coord[,axes],lab.desc)
    col.lab=c(rep("blue",nrow(x$row.coord[,axes])),rep("red",nrow(adjusted.col.coord)))

    if (!is.null(x$proj.col.coord)){
      lab.desc.sup=as.matrix(x$proj.col.coord[,axes,drop=FALSE])
      lab=rbind(lab,lab.desc.sup)
      col.lab=c(col.lab,rep("red4",nrow(x$proj.col.coord)))
    }

    if (!is.null(x$proj.row.coord)){
      lab=rbind(lab,x$proj.row.coord[,axes,drop=FALSE])
      col.lab=c(col.lab,rep("darkblue",nrow(x$proj.row.coord)))
    }

  }else{
    lab=x$row.coord[,axes]
    col.lab=rep("blue",nrow(x$row.coord[,axes]))

    if (!is.null(x$proj.col.coord)){
      lab.desc.sup=as.matrix(x$proj.col.coord[,axes,drop=FALSE])
      lab=rbind(lab,lab.desc.sup)
      col.lab=c(col.lab,rep("red4",nrow(x$proj.col.coord)))
    }

    if (!is.null(x$proj.row.coord)){
      lab=rbind(lab,x$proj.row.coord[,axes,drop=FALSE])
      col.lab=c(col.lab,rep("darkblue",nrow(x$proj.row.coord)))
    }

  }
  nudge=lab*0.01
  p=p+geom_label_repel(as.data.frame(lab),mapping=aes(x=lab[,1],y=lab[,2],label=rownames(lab)),label.size=NA,colour=col.lab,size=size.lab,segment.size=1,label.padding = 0,
                       nudge_x = nudge[,1],nudge_y = nudge[,2],min.segment.length = 1)
  return(p)
}
