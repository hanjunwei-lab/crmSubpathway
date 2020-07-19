##' SubpathwayMatrix
##'
##' @title Get metabolic subpathway activity matrix
##' @description Get metabolic subpathway activity matrix.
##' @param Gexp Matrix of gene expression values (rows are genes, columns are samples).
##' @param Subpathway.list The subpathway list data. We have stored the processed human metabolic subpathway data set in variable `Metspwlist`.
##' Users can also set their own data set of subpathways.
##' @param method Method to employ in the estimation of subpathway
##' enrichment scores per sample. By default this is set to `gsva` (Hänzelmann
##' et al, 2013) and other options are `ssgsea` (Barbie et al, 2009).
##' @param kcdf Character string denoting the kernel to use during the
##' non-parametric estimation of the cumulative distribution function of
##' expression levels across samples when `spw.score.method="gsva"`. By default,
##' `kcdf="Gaussian"` which is suitable when input expression values are
##' continuous, such as microarray fluorescent units in logarithmic scale,
##' RNA-seq log-CPMs, log-RPKMs or log-TPMs. When input expression values are
##' integer counts, such as those derived from RNA-seq experiments, then this
##' argument should be set to `kcdf="Poisson"`.
##' @param min.sz Removes subpathways that contain fewer genes than `spw.min.sz` (default: 10).
##' @param max.sz Removes subpathways that contain more genes than `spw.max.sz` (default: Inf).
##' @param parallel.sz Number of processors to use when doing the calculations in
##' parallel (default value: 1). If parallel.sz=0,
##' then it will use all available core processors unless we set this argument
##' with a smaller number.
##' @details
##' Our method assesses the relative enrichment of metabolic subpathways across samples using a non-parametric approach.
##' Conceptually, this method transforms a p-gene by n-sample gene expression matrix into a g-subpathway by n-sample metabolic subpathway
##' enrichment matrix.
##' @examples
##' library(GSVA)
##' library(parallel)
##' # Get the dataset of metabolic subpathways we have processed.
##' Metspwlist<-get("Metspwlist")
##' # Get the gene expression profile of the case.
##' Geneexp<-get("Geneexp")
##' Spwmatrix<-SubpathwayMatrix(Geneexp,Metspwlist)
##' @importFrom GSVA gsva
##' @importFrom igraph V
##' @export
SubpathwayMatrix<-function(Gexp,Subpathway.list,method="gsva",kcdf="Gaussian",min.sz=1,max.sz=Inf,parallel.sz=0){
  haveGSVA <- isPackageLoaded("GSVA")
  if(haveGSVA==FALSE){
    stop("The 'GSVA' library, should be loaded first")
  }
  haveparallel <- isPackageLoaded("parallel")
  if(haveparallel==FALSE){
    stop("The 'parallel' library, should be loaded first")
  }

  if(is.vector(Subpathway.list[1])==FALSE){
  spw.l<-list()
  for(i in 1:length(Subpathway.list)){
    spw.l[[i]]<-unique(sub("hsa:","",V(Subpathway.list[[i]])$names))
  }
  names(spw.l)<-names(Subpathway.list)

  qc<-NULL
  for (i in 1:length(spw.l)) {
    if(length(spw.l[[i]])==0){
      qc<-c(qc,i)
    }
  }
  spw.l<-spw.l[-qc]

  spw_names<-names(spw.l)
  pw_names<-NULL
  for(i in 1:length(spw_names)){
    x<-strsplit(spw_names[i],split="_")
    pw_names[i]<-x[[1]][1]
  }

  qc_pw_names<-names(table(pw_names))
  qcfspw<-NULL
  for (i in 1:length(qc_pw_names)) {
    pp_index<-which(pw_names==qc_pw_names[i])
    if(length(pp_index)>1){
      x<-spw.l[pp_index]
      qc<-redup(x)
      qc<-unique(qc)
      if(length(qc)>0){
        x<-x[-qc]
        qcfspw<-c(qcfspw,x)
    }else{
        qcfspw<-c(qcfspw,x)
      }
    }else{
      qcfspw<-c(qcfspw,x)
    }
  }
  }else{
    qcfspw<-Subpathway.list
  }

  spwmatrix<-gsva(Gexp,qcfspw,method=method,kcdf=kcdf,min.sz=min.sz,max.sz=max.sz,parallel.sz=parallel.sz,mx.diff=TRUE,verbose=FALSE)
  return(spwmatrix)
}
