# Internal function
getKcSubiGraph<-function(k=4,graphList){
      if(k<1)
            stop("k can't <1 ")
	  require(graph,quietly = TRUE)
	  require(RBGL,quietly = TRUE)
      subGraphIndex<-0
      subGraph<-list()
      subNames<-character()
	  #graphList<-getSimpleGraph(graphList)
	  if(length(graphList)>0){
      for(i in 1:length(graphList)){
	     if(vcount(graphList[[i]])>0){
	     graph1<-igraph.to.graphNEL(graphList[[i]])
         kc<-kCliques(graph1)
         if(length(kc)>0){
            if(k<=length(kc)){
                  kkc<-kc[[k]]
            }
            else{
                  kkc<-kc[[length(kc)]]
            }
            node_name<-V(graphList[[i]])$name
			V_seq<-seq(1,vcount(graphList[[i]]))  #new!important
            for(j in 1:length(kkc)){
                  subGraphIndex<-subGraphIndex+1
				  #index<-V(graphList[[i]])[V(graphList[[i]])$name %in% kkc[[j]]]
				  #print(kkc[[j]])
				  index<-V_seq[node_name %in% kkc[[j]]]
				  #print(index)
                  subGraph[[subGraphIndex]]<-induced.subgraph(graphList[[i]],index)
                  subGraph[[subGraphIndex]]$number<-paste(subGraph[[subGraphIndex]]$number,j,sep="_")
				  names(subGraph)[subGraphIndex]<-subGraph[[subGraphIndex]]$number
            }
          }
		  }
      }
	  }
      return(subGraph)
}
