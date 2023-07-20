#' Plot factor plan resulting from multiple-response Correspondence Analysis (MR-CA) applied on sensory data
#'
#' @description This function plots the results coming from \code{\link[MultiResponseR]{sensory.mrCA}}
#'
#' @param x A list returned by \code{\link[MultiResponseR]{sensory.mrCA}}
#' @param axes Which dimensions of the MR-CA should be plotted?
#' @param alpha.total.bootstrap.test The alpha risk of the total bootstrap tests. See details
#' @param alpha.ellipse The alpha risk of the confidence ellipses
#' @param select.desc A character vector specifying the descriptors to plot. By default, all descriptors are plotted
#' @param rev.x Should the horizontal plotted dimension be reversed? Useful in case of map comparisons to align products
#' @param rev.y Should the vertical plotted dimension be reversed? Useful in case of map comparisons to align products
#' @param size.points The size of the points used to represent the products on the map
#' @param size.lab The size of the label on the map
#' @param expansion The factor of expansion applied to descriptors coordinates to increase readability
#' @param title An optional title to be added to the plot
#' @param ... further arguments passed to or from other methods
#'
#' @details
#' \itemize{
#'   \item \strong{alpha.total.bootstrap.test}: products non-significantly different at the alpha risk of \emph{alpha.total.bootstrap.test} according to the total bootstrap test are linked by a line on the plot. If these links are not required, \emph{alpha.total.bootstrap.test} can be set to 1
#' }
#'
#'
#' @return A MR-CA factor map
#'
#'
#' @export
#'
#' @import FactoMineR
#' @import graphics
#' @import stats
#' @import ggplot2
#' @import ggrepel
#'
#'
#' @examples
#'data(milkchoc)
#'
#'dim.sig=sensory.mr.dimensionality.test(milkchoc)$dim.sig
#'
#'res=sensory.mrCA(milkchoc,nbaxes.sig=dim.sig)
#'
#'plot(res)
plot.sensory.mrCA=function(x,
                           axes = c(1,2),
                           alpha.total.bootstrap.test = 0.05,
                           alpha.ellipse = alpha.total.bootstrap.test,
                           select.desc = rownames(x$desc.coord),
                           rev.x = FALSE,
                           rev.y = FALSE,
                           size.points = 3.5,
                           size.lab = 6,
                           expansion = 1.25,
                           title = NULL,...){
  if (!inherits(x,"sensory.mrCA")){
    stop("class(x) must be sensory.mrCA")
  }

  check.axes=ncol(x$bootstrap.replicate.coord)-1
  if (max(axes)>check.axes){
    stop("max(axes) must be lower than or equal to the number of nbaxes.sig used in the mrCA")
  }

  classe=class(select.desc)
  if (classe!="character"){
    stop("class(select.desc) must be character")
  }

  ell=coord.ellipse(x$bootstrap.replicate.coord,axes=axes,level.conf = (1-alpha.ellipse))$res

  adjusted.col.coord=x$desc.coord

  if (!is.null(select.desc)){
    adjusted.col.coord=adjusted.col.coord[select.desc,,drop=FALSE]
  }


  max.norm.prod = max(abs(x$prod.coord[,axes]))
  max.norm.desc = max(abs(adjusted.col.coord[,axes]))


  expand = max.norm.prod/max.norm.desc*expansion
  adjusted.col.coord=adjusted.col.coord*expand

  if(rev.x){
    x$prod.coord[,axes[1]]=-x$prod.coord[,axes[1]]
    adjusted.col.coord[,axes[1]]=-adjusted.col.coord[,axes[1]]
    ell[,2]=-ell[,2]
  }
  if(rev.y){
    x$prod.coord[,axes[2]]=-x$prod.coord[,axes[2]]
    adjusted.col.coord[,axes[2]]=-adjusted.col.coord[,axes[2]]
    ell[,3]=-ell[,3]
  }

  xmin=min(x$prod.coord[,axes[1]],adjusted.col.coord[,axes[1]],ell[,2])
  xmax=max(x$prod.coord[,axes[1]],adjusted.col.coord[,axes[1]],ell[,2])
  ymin=min(x$prod.coord[,axes[2]],adjusted.col.coord[,axes[2]],ell[,3])
  ymax=max(x$prod.coord[,axes[2]],adjusted.col.coord[,axes[2]],ell[,3])

  pmin=min(xmin,ymin)*1.05
  pmax=max(xmax,ymax)*1.05

  p=ggplot(as.data.frame(x$prod.coord),aes(x=x$prod.coord[,axes[1]],y=x$prod.coord[,axes[2]]))+theme_bw()
  p=p+xlim(pmin,pmax)+ylim(pmin,pmax)+xlab(paste("Dim ",axes[1]," (",round(x$eigen[axes[1],2],2)," %)",sep=""))+ylab(paste("Dim ",axes[2]," (",round(x$eigen[axes[2],2],2)," %)",sep=""))+ggtitle(title)
  p=p+theme(axis.title.x = element_text(size = 16,face = "bold"),axis.title.y = element_text(size = 16,face = "bold"),plot.title = element_text(hjust = 0.5,face = "bold",size=20))
  p=p+geom_hline(yintercept=0,linetype="dashed",linewidth=1)+geom_vline(xintercept=0,linetype="dashed",linewidth=1)
  p=p+geom_path(data=as.data.frame(ell),aes(x=ell[,2],y=ell[,3],group=ell[,1]),colour="blue",linewidth=1)

  diff.test=x$total.bootstrap.test.pvalues
  df.segment=NULL
  for (i in 1:nrow(diff.test)){
    for (j in i:ncol(diff.test)){
      p.1=rownames(diff.test)[i]
      p.2=colnames(diff.test)[j]
      if (diff.test[p.1,p.2]>alpha.total.bootstrap.test & i!=j){
        ac.produit.coord=as.data.frame(x$prod.coord[,axes])
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


  if (!is.null(select.desc)){
    df.fleche=as.matrix(adjusted.col.coord[,axes,drop=FALSE])
    col.fleche=rep("red",length(select.desc))
    if(!is.null(x$proj.col.coord)){
      df.fleche=rbind(df.fleche,as.matrix(x$proj.col.coord[,axes,drop=FALSE]))
      col.fleche=c(col.fleche,rep("red4",nrow(x$proj.col.coord)))
    }
  }

  p=p+geom_segment(data = as.data.frame(df.fleche), aes(x=0, y=0,xend = df.fleche[,1], yend = df.fleche[,2]), arrow=arrow(length = unit(0.4, "cm"),type = "closed"), colour=col.fleche,linewidth=1)

  df.point=x$prod.coord[,axes,drop=FALSE]
  col.point=rep("blue",nrow(x$prod.coord))
  p=p+geom_point(data=as.data.frame(df.point),aes(x=df.point[,1],y=df.point[,2]),colour=col.point,size=size.points)

  if (!is.null(select.desc)){
    lab.desc=as.matrix(adjusted.col.coord[,axes,drop=FALSE])
    rownames(lab.desc)=rownames(adjusted.col.coord)
    lab=rbind(x$prod.coord[,axes],lab.desc)
    col.lab=c(rep("blue",nrow(x$prod.coord[,axes])),rep("red",nrow(adjusted.col.coord)))
  }else{
    lab=x$prod.coord[,axes]
    col.lab=rep("blue",nrow(x$prod.coord[,axes]))
  }
  nudge=lab*0.01
  p=p+geom_label_repel(as.data.frame(lab),mapping=aes(x=lab[,1],y=lab[,2],label=rownames(lab)),label.size=NA,colour=col.lab,size=size.lab,segment.size=1,label.padding = 0,
                       nudge_x = nudge[,1],nudge_y = nudge[,2],min.segment.length = 1)
  return(p)
}
