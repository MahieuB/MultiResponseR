#' Plot significant cells
#'
#' @description This function plots the results coming from \code{\link[MultiResponseR]{sensory.mr.sig.cell}}
#'
#' @param x A list returned by \code{\link[MultiResponseR]{sensory.mr.sig.cell}}
#' @param alpha.1 The alpha risk to consider the tests as significant
#' @param alpha.2 The alpha risk to consider the tests as showing a trend. If trends are not to be considered, \emph{alpha.2} can be set to 0 (Default)
#' @param choice Which table from \emph{res} should be plotted? Default is percent.derived.cont
#' @param col.greater.1 The color used to highlight significant positive associations
#' @param col.lower.1 The color used to highlight significant negative associations
#' @param col.greater.2 The color used to highlight positive associations showing a trend
#' @param col.lower.2 The color used to highlight negative associations showing a trend
#' @param ... further arguments passed to or from other methods
#'
#' @return A table with cells highlighted
#'
#' @export
#'
#'
#' @import flextable
#' @import officer
#' @import stats
#' @import utils
#'
#' @examples
#'data(milkchoc)
#'
#'dim.sig=sensory.mr.dimensionality.test(milkchoc)$dim.sig
#'
#'res=sensory.mr.sig.cell(milkchoc,nbaxes.sig=dim.sig)
#'
#'plot(res)
plot.sensory.mr.sig.cell=function(x,alpha.1=0.05,alpha.2=0,choice="percent.derived.cont",col.greater.1="green3",col.lower.1="orangered",col.greater.2="lightgreen",col.lower.2="lightsalmon",...){
  if (!inherits(x,"sensory.mr.sig.cell")){
    stop("class(x) must be sensory.mr.sig.cell")
  }
  if (!choice%in%c("original.cont","percent.cont","null.cont","p.value","derived.cont","percent.derived.cont")){
    stop("choice is not valid")
  }
  ou.choice=which(names(x)==choice)
  plot.mat=x[[ou.choice]]
  mat.pval=x$p.value
  mat.exp=x$null.cont
  mat.obs=x$derived.cont

  plot.col=as.data.frame(matrix("white",nrow(plot.mat),ncol(plot.mat)))

  for (i in 1:nrow(plot.col)){
    for(j in 1:ncol(plot.col)){
      if(mat.pval[i,j]<=alpha.2 & mat.obs[i,j]>mat.exp[i,j]){
        plot.col[i,j]=col.greater.2
      }
      if(mat.pval[i,j]<=alpha.1 & mat.obs[i,j]>mat.exp[i,j]){
        plot.col[i,j]=col.greater.1
      }
      if (mat.pval[i,j]<=alpha.2 & mat.obs[i,j]<mat.exp[i,j]){
        plot.col[i,j]=col.lower.2
      }
      if (mat.pval[i,j]<=alpha.1 & mat.obs[i,j]<mat.exp[i,j]){
        plot.col[i,j]=col.lower.1
      }
    }
  }

  plot.mat=cbind.data.frame(rownames(plot.mat),plot.mat)
  colnames(plot.mat)[1]=" "
  plot.col=cbind.data.frame(rep("cyan",nrow(plot.col)),plot.col)
  plot.col=apply(plot.col, 2, as.character)
  rownames(plot.col)=rownames(plot.mat)

  sort.names=sort(rownames(plot.mat))
  plot.col=plot.col[sort.names,]
  plot.mat=plot.mat[sort.names,]

  size.line=15
  size.text=20

  p=flextable(as.data.frame(plot.mat))
  p=fontsize(p,size=size.text,part = "all")
  p=height_all(p,0.5,part="all")
  p=bg(p,bg="cyan",part="header")
  p=padding(p,padding.top = size.line,padding.bottom = size.line,part="all")
  for (i in 1:nrow(plot.col)){
    for (j in 1:ncol(plot.col)){
      p=bg(p,i,j,plot.col[i,j])
      if (j>1){
        p=align(p,j=j,align = "center",part="all")
      }
    }
  }
  p=border_inner_h(p,border = fp_border(color = "black", style = "solid", width = 1),part="all")
  p=bold(p,j=1)
  p=bold(p,part="header")
  return(p)
}
