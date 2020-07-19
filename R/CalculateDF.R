##' CalculateDF
##'
##' @title Calculate the difference of metabolic subpathway activity
##' @description For the metabolic subpathway activity matrix, this method can calculate the difference and significance of the metabolic subpathway activity between disease and normal.
##' @param spwmatrix Metabolic subpathway activity matrix (the result of function `SubpathwayMatrix`),
##' @param smaple_class Sample phenotype vector in metabolic subpathway activity matrix.
##' @param casename Disease sample phenotype label.
##' @param controlname Control sample phenotype label.
##' @examples
##' library(limma)
##' # Get the metabolic subpathway matrix.
##' Spwmatrix<-get("Spwmatrix")
##' spwDF<-CalculateDF(Spwmatrix,colnames(Spwmatrix),"cancer","control")
##' @importFrom limma lmFit
##' @importFrom limma makeContrasts
##' @importFrom limma eBayes
##' @importFrom limma contrasts.fit
##' @importFrom limma topTable
##' @importFrom stats model.matrix
##' @export
CalculateDF<-function(spwmatrix,
                      smaple_class,
                      casename="",
                      controlname=""
                      ){
  havelimma <- isPackageLoaded("limma")
  if(havelimma==FALSE){
    stop("The 'limma' library, should be loaded first")
  }
  smaple_class1<-smaple_class
  caseindex<-which(smaple_class1==casename)
  controlindex<-which(smaple_class1==controlname)
  smaple_class1[caseindex]<-"case"
  smaple_class1[controlindex]<-"control"
  colnames(spwmatrix)<-smaple_class1
  f<-factor(smaple_class1,levels = c("case","control"))
  design<-model.matrix(~0+f)
  colnames(design) <- c("case","control")
  fit<-lmFit(spwmatrix,design)
  contrast.matrix<-makeContrasts(case-control,levels=design)
  fit2<-contrasts.fit(fit,contrast.matrix)
  fit2<-eBayes(fit2)
  DFresult<-topTable(fit2,coef = 1,number = Inf)
  colnames(spwmatrix)<-smaple_class
  result<-list(Difference_analysis=DFresult,SubpathwayMatrix=spwmatrix)
  return(result)
}
