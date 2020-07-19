##' k_clique
##'
##' @title k-clique algorithm
##' @description Mining metabolic subpathways through k-clique algorithm.
##' @param file.path Path of KGML files (XML format).
##' @param file.names A character vector of names of KGML files.
##' @param K An integer. A distance similarity parameter (default: 4).
##' @details
##' In this method, we convert the each metabolic pathways in the KEGG database into an undirected network diagram by connecting
##' two genes (enzymes) into an edge if there is a common compound in their corresponding reactions. Based on the distance
##' similarity between genes, we use the K-clique method in social network analysis to mine subpathways and the distance between
##' all genes is no exceeding K is determined as the metabolic subpathway.At the same time, our method can also be applied to other
##' non-metabolic pathways of KEGG database.
##' @examples
##' library(graph);
##' library(RBGL);
##' library(igraph);
##' library(XML);
##' file.path<-paste(system.file(package="crmSubpathway"),"/extdata/",sep="")
##' file.names<-c("hsa00010.xml","hsa00020.xml")
##' spwlist<-k_clique(file.path,file.names)
##' @export
k_clique<-function(file.path="",file.names="",K=4){
  havegraph <- isPackageLoaded("graph")
  if(havegraph==FALSE){
    stop("The 'graph' library, should be loaded first")
  }
  haveRBGL <- isPackageLoaded("RBGL")
  if(haveRBGL==FALSE){
    stop("The 'RBGL' library, should be loaded first")
  }
  haveigraph <- isPackageLoaded("igraph")
  if(haveigraph==FALSE){
    stop("The 'igraph' library, should be loaded first")
  }
  haveXML <- isPackageLoaded("XML")
  if(haveXML==FALSE){
    stop("The 'XML' library, should be loaded first")
  }
  Metabolicxml<-getPathway(file.path,file.names)
  MetabolicGraph<-getMetabolicGraph(Metabolicxml)
  MetabolicUGraph<-getUGraph(MetabolicGraph)
  metagl<-filterNode(MetabolicUGraph,nodeType=c("map"))
  simMeta<-simplifyGraph(metagl,nodeType="geneProduct")
  expMeta<-expandNode(simMeta)
  Subpathway<-getKcSubiGraph(k=K,expMeta)
  return(Subpathway)
}
