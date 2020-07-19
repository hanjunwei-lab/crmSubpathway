## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(crmSubpathway)

## ----eval=FALSE, include=FALSE------------------------------------------------
#  library(graph);
#  library(RBGL);
#  library(igraph);
#  library(XML);
#  file.path<-paste(system.file(package="crmSubpathway"),"/inst/extdata/",sep="")
#  file.names<-c("hsa00010.xml","hsa00020.xml")
#  spwlist<-k_clique(file.path,file.names)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(GSVA)
library(parallel)
# Get the dataset of metabolic subpathways we have processed.
Metspwlist<-get("Metspwlist")
# Get the gene expression profile of the case.
Geneexp<-get("Geneexp")
Spwmatrix<-SubpathwayMatrix(Geneexp,Metspwlist)
head(Spwmatrix)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(limma)
# Get the metabolic subpathway matrix.
Spwmatrix<-get("Spwmatrix")
spwDF<-CalculateDF(Spwmatrix,colnames(Spwmatrix),"cancer","control")

## ----fig.height=7, fig.width=7, message=FALSE, warning=FALSE------------------
library(igraph)
plotNetGraph(spwid="00010_1")

## ----fig.height=7, fig.width=7, message=FALSE, warning=FALSE------------------
library(pheatmap)
DFspw<-get("DFspw")
plotspwheatmap(DFspw,cluster_rows = TRUE,show.colnames=FALSE)

