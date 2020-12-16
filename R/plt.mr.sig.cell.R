#' Plot significant cells
#'
#' @description This function plots the results coming from \code{\link[MultiResponseR]{sensory.mr.sig.cell}} or \code{\link[MultiResponseR]{mr.sig.cell}}
#'
#' @param res A list returned by \code{\link[MultiResponseR]{sensory.mr.sig.cell}} or \code{\link[MultiResponseR]{mr.sig.cell}}
#' @param alpha The alpha risk to consider the tests as significant. Note that the tests can be slightly conservative as suggested by Appendix in Mahieu, Schlich, Visalli, and Cardot (2020)
#' @param choice Which table from \emph{res} should be plotted? Default is percent.derived.cont
#' @param col.greater The color used to highlight significant positive associations
#' @param col.lower The color used to highlight significant negative associations
#'
#'
#' @return A table with significant cells highlighted
#'
#' @export
#'
#' @references Mahieu, B., Schlich, P., Visalli, M., & Cardot, H. (2020). A multiple-response chi-square framework for the analysis of Free-Comment and Check-All-That-Apply data. Manuscript submitted for publication.
#'
#' @import flextable
#' @import officer
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
#' res=mr.sig.cell(dset)
#'
#' plt.mr.sig.cell(res)
#'
#' # sensory example
#'data(milkchoc)
#'
#'parallel::detectCores()
#'
#'dim.sig=sensory.mr.dimensionality.test(milkchoc)$dim.sig
#'
#'res=sensory.mr.sig.cell(milkchoc,nbaxes.sig=dim.sig)
#'
#'plt.mr.sig.cell(res)
plt.mr.sig.cell=function(res,alpha=0.075,choice="percent.derived.cont",col.greater="green3",col.lower="orangered"){
  classe=class(res)
  if(classe!="list"){
    stop("res must be a list resulting from the execution of sensory.mr.sig.cell or mr.sig.cell")
  }
  taille=length(res)
  if (taille!=6){
    stop("res must be a list resulting from the execution of sensory.mr.sig.cell or mr.sig.cell")
  }
  if (!choice%in%c("original.cont","percent.cont","null.cont","p.value","derived.cont","percent.derived.cont")){
    stop("choice is not valid")
  }
  ou.choice=which(names(res)==choice)
  plot.mat=res[[ou.choice]]
  mat.pval=res$p.value
  mat.exp=res$null.cont
  mat.obs=res$derived.cont

  plot.col=plot.mat

  for (i in 1:nrow(plot.col)){
    for(j in 1:ncol(plot.col)){
      if(mat.pval[i,j]<=alpha & mat.obs[i,j]>mat.exp[i,j]){
        plot.col[i,j]=col.greater
      }else if (mat.pval[i,j]<=alpha & mat.obs[i,j]<mat.exp[i,j]){
        plot.col[i,j]=col.lower
      }else{
        plot.col[i,j]="white"
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
