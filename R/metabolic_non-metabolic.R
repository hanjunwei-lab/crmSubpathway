# Internal function
####################################################
#expandNode filterNode merge getSimpleGraph getUGraph mapNode simplifyGraph
######################################################

expandNode<-function(graphList,nodeType=c("ortholog","enzyme","gene","compound","map")){
    expandgraphList<-list()
    graphListLength<-length(graphList)
    for(i in 1:graphListLength){
    ###############################################################
       name<-c();id<-c();names<-c();type<-c();reaction<-c();link<-c()
       graphics_name<-c();graphics_fgcolor<-c();graphics_bgcolor<-c()
       graphics_type<-c();graphics_x<-c();graphics_y<-c()
       graphics_width<-c();graphics_height<-c();graphics_coords<-c()
       Vcount<-vcount(graphList[[i]])
       if(Vcount>0){
         for(j in 1:Vcount){
		 #all new!
		  if(get.vertex.attribute(graphList[[i]],"type",j) %in% nodeType){
           vlist<-strsplit(get.vertex.attribute(graphList[[i]],"names",j)," ")
           for(k in 1:length(vlist[[1]])){
                name<-c(name,paste(get.vertex.attribute(graphList[[i]],"name",j),"_",k,sep=""))
                id<-c(id,paste(get.vertex.attribute(graphList[[i]],"id",j),"_",k,sep=""))
                names<-c(names,vlist[[1]][k])
                type<-c(type,get.vertex.attribute(graphList[[i]],"type",j))
                reaction<-c(reaction,get.vertex.attribute(graphList[[i]],"reaction",j))
				newlink<-paste("http://www.kegg.jp/dbget-bin/www_bget?",vlist[[1]][k],sep="")
				link<-c(link,newlink)
                #link<-c(link,get.vertex.attribute(graphList[[i]],"link",j))
				graphics_name<-c(graphics_name,vlist[[1]][k])
                #graphics_name<-c(graphics_name,get.vertex.attribute(graphList[[i]],"graphics_name",j))
                graphics_fgcolor<-c(graphics_fgcolor,get.vertex.attribute(graphList[[i]],"graphics_fgcolor",j))
                graphics_bgcolor<-c(graphics_bgcolor,get.vertex.attribute(graphList[[i]],"graphics_bgcolor",j))
                graphics_type<-c(graphics_type,get.vertex.attribute(graphList[[i]],"graphics_type",j))
                graphics_x<-c(graphics_x,get.vertex.attribute(graphList[[i]],"graphics_x",j))
                graphics_y<-c(graphics_y,get.vertex.attribute(graphList[[i]],"graphics_y",j))
                graphics_width<-c(graphics_width,get.vertex.attribute(graphList[[i]],"graphics_width",j))
                graphics_height<-c(graphics_height,get.vertex.attribute(graphList[[i]],"graphics_height",j))
                graphics_coords<-c(graphics_coords,get.vertex.attribute(graphList[[i]],"graphics_coords",j))
           }
          }
		  else{
                name<-c(name,get.vertex.attribute(graphList[[i]],"name",j))
                id<-c(id,get.vertex.attribute(graphList[[i]],"id",j))
                names<-c(names,get.vertex.attribute(graphList[[i]],"names",j))
                type<-c(type,get.vertex.attribute(graphList[[i]],"type",j))
                reaction<-c(reaction,get.vertex.attribute(graphList[[i]],"reaction",j))
                link<-c(link,get.vertex.attribute(graphList[[i]],"link",j))
                graphics_name<-c(graphics_name,get.vertex.attribute(graphList[[i]],"graphics_name",j))
                graphics_fgcolor<-c(graphics_fgcolor,get.vertex.attribute(graphList[[i]],"graphics_fgcolor",j))
                graphics_bgcolor<-c(graphics_bgcolor,get.vertex.attribute(graphList[[i]],"graphics_bgcolor",j))
                graphics_type<-c(graphics_type,get.vertex.attribute(graphList[[i]],"graphics_type",j))
                graphics_x<-c(graphics_x,get.vertex.attribute(graphList[[i]],"graphics_x",j))
                graphics_y<-c(graphics_y,get.vertex.attribute(graphList[[i]],"graphics_y",j))
                graphics_width<-c(graphics_width,get.vertex.attribute(graphList[[i]],"graphics_width",j))
                graphics_height<-c(graphics_height,get.vertex.attribute(graphList[[i]],"graphics_height",j))
                graphics_coords<-c(graphics_coords,get.vertex.attribute(graphList[[i]],"graphics_coords",j))
		  }
         }
       }
      entry1<-c();entry2<-c();str<-list()
      Ecount<-ecount(graphList[[i]])
      if(Ecount>0){
         edgesAttr<-list.edge.attributes(graphList[[i]])
         edgesAttrLenngth<-length(edgesAttr)
         if(edgesAttrLenngth>0){
            for(j in 1:edgesAttrLenngth){
               str[[j]]<-c("first")
            }
         }
         for(j in 1:Ecount){

           e<-ends(graphList[[i]],j)

           l_1<-strsplit(get.vertex.attribute(graphList[[i]],"names",e[1])," ")
           l1<-length(l_1[[1]])
           l_2<-strsplit(get.vertex.attribute(graphList[[i]],"names",e[2])," ")
           l2<-length(l_2[[1]])

           e1<-c()
		  if(get.vertex.attribute(graphList[[i]],"type",e[1]) %in% nodeType){
           for(m in 1:l1){
               e1<-c(e1,paste(get.vertex.attribute(graphList[[i]],"id",e[1]),"_",m,sep=""))
           }
		  }
		  else{
		    e1<-get.vertex.attribute(graphList[[i]],"id",e[1])
			l1<-1
		  }
		  if(get.vertex.attribute(graphList[[i]],"type",e[2]) %in% nodeType){
           for(n in 1:l2){
               entry1<-c(entry1,e1)
               entry2<-c(entry2,rep(paste(get.vertex.attribute(graphList[[i]],"id",e[2]),"_",n,sep=""),l1))
           }
		  }
		  else{
		    entry1<-c(entry1,e1)
		    entry2<-c(entry2,rep(get.vertex.attribute(graphList[[i]],"id",e[2]),l1))
			l2<-1
		  }
           ############

           if(edgesAttrLenngth>0){
               for(k in 1:edgesAttrLenngth){
                  str[[k]]<-c(str[[k]],rep(get.edge.attribute(graphList[[i]],edgesAttr[k],j),l1*l2))
               }
           }
        }
     }

 	 #print(entry1);print(entry2)
	 #print(name)
   ###############################################################################

     if(Ecount>0){
         V<-data.frame(name=name,id=id,names=names,type=type,reaction=reaction,link=link,graphics_name=graphics_name,
                    graphics_fgcolor=graphics_fgcolor,graphics_bgcolor=graphics_bgcolor,graphics_type=graphics_type,
                    graphics_x=graphics_x,graphics_y=graphics_y,graphics_width=graphics_width,graphics_height=graphics_height,
                    graphics_coords=graphics_coords)
        edgesAttrLenngth<-length(edgesAttr)
        if(edgesAttrLenngth>0){
           E<-data.frame(entry1,entry2)
           for(k in 1:edgesAttrLenngth){
              E<-data.frame(E,str[[k]][2:length(str[[k]])])
           }
           names(E)<-c("entry1","entry2",edgesAttr)
        }
        else{
           E<-data.frame(entry1=entry1,entry2=entry2)
        }
        expandgraphList[[i]]<-graph.data.frame(E,directed=is.directed(graphList[[i]]),V)
     }else{
        expandgraphList[[i]]<-graph.empty(n=0,directed=is.directed(graphList[[i]]))
        expandgraphList[[i]]<-add.vertices(expandgraphList[[i]],length(name),name=name,id=id,names=names,type=type,reaction=reaction,link=link,graphics_name=graphics_name,
                graphics_fgcolor=graphics_fgcolor,graphics_bgcolor=graphics_bgcolor,graphics_type=graphics_type,
                graphics_x=graphics_x,graphics_y=graphics_y,graphics_width=graphics_width,graphics_height=graphics_height,
                graphics_coords=graphics_coords)
     }
        expandgraphList[[i]]<-set.graph.attribute(expandgraphList[[i]],"name",get.graph.attribute(graphList[[i]],"name"))
        expandgraphList[[i]]<-set.graph.attribute(expandgraphList[[i]],"number",get.graph.attribute(graphList[[i]],"number"))
        expandgraphList[[i]]<-set.graph.attribute(expandgraphList[[i]],"org",get.graph.attribute(graphList[[i]],"org"))
        expandgraphList[[i]]<-set.graph.attribute(expandgraphList[[i]],"title",get.graph.attribute(graphList[[i]],"title"))
        expandgraphList[[i]]<-set.graph.attribute(expandgraphList[[i]],"image",get.graph.attribute(graphList[[i]],"image"))
        expandgraphList[[i]]<-set.graph.attribute(expandgraphList[[i]],"link",get.graph.attribute(graphList[[i]],"link"))
   }
   names(expandgraphList)<-names(graphList)
   return(expandgraphList)
}

#################################################################################################3
###c("ortholog","enzyme","gene","compound","map")
filterNode<-function(graphList,nodeType=c("map")){
       graphListLength<-length(graphList)
	   if(graphListLength>0){
       for(i in 1:graphListLength){
           vCount<-vcount(graphList[[i]])
           deleteId<-c()
		   if(vCount>0){
           for(j in 1:vCount){
               #if(get.vertex.attribute(graphList[[i]],"type",j-1) %in% nodeType){
			   if(get.vertex.attribute(graphList[[i]],"type",j) %in% nodeType){
                   deleteId<-c(deleteId,j)
               }
           }
           graphList[[i]]<-delete.vertices(graphList[[i]],deleteId)
		   }
       }
	   }
       return(graphList)
}


#################################################################################################3

mergeNode<-function(graphList,simpleGraph=TRUE){
    mergegraphList<-list()
    graphListLength<-length(graphList)
	if(graphListLength>0){
    for(i in 1:graphListLength){


       Vcount<-vcount(graphList[[i]])
	   middle<-list();mid<-c()
       if(Vcount>0){
         for(j in 1:Vcount){
           #temp_name<-get.vertex.attribute(graphList[[i]],"names",j-1)
		   #new!
           temp_name<-get.vertex.attribute(graphList[[i]],"names",j)
		   matched_name_index<-match(temp_name,mid)
           if(is.na(matched_name_index)){
                mid<-c(mid,temp_name)
                #middle[[length(mid)]]<-j-1
				middle[[length(mid)]]<-j
           }
		   else{
 		        #middle[[matched_name_index]]<-c(middle[[matched_name_index]],j-1)

				middle[[matched_name_index]]<-c(middle[[matched_name_index]],j)
		   }
         }
       }
	 middleLength<-length(middle)
       name<-c();id<-c();names<-c();type<-c();reaction<-c();link<-c()
       graphics_name<-c();graphics_fgcolor<-c();graphics_bgcolor<-c()
       graphics_type<-c();graphics_x<-c();graphics_y<-c()
       graphics_width<-c();graphics_height<-c();graphics_coords<-c()
     if(middleLength>0){
	   for(j in 1:middleLength){

                name[j]<-get.vertex.attribute(graphList[[i]],"names",middle[[j]][1])
                id[j]<-get.vertex.attribute(graphList[[i]],"id",middle[[j]][1])
                names[j]<-get.vertex.attribute(graphList[[i]],"names",middle[[j]][1])
                type[j]<-get.vertex.attribute(graphList[[i]],"type",middle[[j]][1])
                reaction[j]<-get.vertex.attribute(graphList[[i]],"reaction",middle[[j]][1])
                link[j]<-get.vertex.attribute(graphList[[i]],"link",middle[[j]][1])
                graphics_name[j]<-get.vertex.attribute(graphList[[i]],"graphics_name",middle[[j]][1])
                graphics_fgcolor[j]<-get.vertex.attribute(graphList[[i]],"graphics_fgcolor",middle[[j]][1])
                graphics_bgcolor[j]<-get.vertex.attribute(graphList[[i]],"graphics_bgcolor",middle[[j]][1])
                graphics_type[j]<-get.vertex.attribute(graphList[[i]],"graphics_type",middle[[j]][1])
                graphics_x[j]<-get.vertex.attribute(graphList[[i]],"graphics_x",middle[[j]][1])
                graphics_y[j]<-get.vertex.attribute(graphList[[i]],"graphics_y",middle[[j]][1])
                graphics_width[j]<-get.vertex.attribute(graphList[[i]],"graphics_width",middle[[j]][1])
                graphics_height[j]<-get.vertex.attribute(graphList[[i]],"graphics_height",middle[[j]][1])
                graphics_coords[j]<-get.vertex.attribute(graphList[[i]],"graphics_coords",middle[[j]][1])
            middleJLength<-length(middle[[j]])
            if(middleJLength>1){
             for(k in 2:middleJLength){
                id[j]<-paste(id[j],get.vertex.attribute(graphList[[i]],"id",middle[[j]][k]),sep=";")
                reaction[j]<-paste(reaction[j],get.vertex.attribute(graphList[[i]],"reaction",middle[[j]][k]),sep=";")
             }
            }
	   }
     }

      entry1<-c();entry2<-c();str<-list()
      Ecount<-ecount(graphList[[i]])
      if(Ecount>0){
         edgesAttr<-list.edge.attributes(graphList[[i]])
         edgesAttrLenngth<-length(edgesAttr)
         if(edgesAttrLenngth>0){
            for(j in 1:edgesAttrLenngth){
               str[[j]]<-c("first")
            }
         }
         for(j in 1:Ecount){
           ####  e?д??ŵ??ǵ?j???ߵĶ???id
           #e<-ends(graphList[[i]],j-1)
           e<-ends(graphList[[i]],j)
           entry1<-c(entry1,get.vertex.attribute(graphList[[i]],"names",e[1]))
           entry2<-c(entry2,get.vertex.attribute(graphList[[i]],"names",e[2]))

           if(edgesAttrLenngth>0){
               for(k in 1:edgesAttrLenngth){
                  #str[[k]]<-c(str[[k]],get.edge.attribute(graphList[[i]],edgesAttr[k],j-1))
                  str[[k]]<-c(str[[k]],get.edge.attribute(graphList[[i]],edgesAttr[k],j))
               }
           }
        }
     }


     if(Ecount>0){
         V<-data.frame(name=name,id=id,names=names,type=type,reaction=reaction,link=link,graphics_name=graphics_name,
                    graphics_fgcolor=graphics_fgcolor,graphics_bgcolor=graphics_bgcolor,graphics_type=graphics_type,
                    graphics_x=graphics_x,graphics_y=graphics_y,graphics_width=graphics_width,graphics_height=graphics_height,
                    graphics_coords=graphics_coords)
        edgesAttrLenngth<-length(edgesAttr)
        if(edgesAttrLenngth>0){
           E<-data.frame(entry1,entry2)
           for(k in 1:edgesAttrLenngth){
              E<-data.frame(E,str[[k]][2:length(str[[k]])])
           }
           names(E)<-c("entry1","entry2",edgesAttr)
        }
        else{
           E<-data.frame(entry1=entry1,entry2=entry2)
        }
       mergegraphList[[i]]<-graph.data.frame(E,directed=is.directed(graphList[[i]]),V)
     }
     else{
        mergegraphList[[i]]<-graph.empty(n=0,directed=is.directed(graphList[[i]]))
        mergegraphList[[i]]<-add.vertices(mergegraphList[[i]],length(name),name=name,id=id,names=names,type=type,reaction=reaction,link=link,graphics_name=graphics_name,
                graphics_fgcolor=graphics_fgcolor,graphics_bgcolor=graphics_bgcolor,graphics_type=graphics_type,
                graphics_x=graphics_x,graphics_y=graphics_y,graphics_width=graphics_width,graphics_height=graphics_height,
                graphics_coords=graphics_coords)
     }
	    mergegraphList[[i]]<-set.graph.attribute(mergegraphList[[i]],"name",get.graph.attribute(graphList[[i]],"name"))
        mergegraphList[[i]]<-set.graph.attribute(mergegraphList[[i]],"number",get.graph.attribute(graphList[[i]],"number"))
        mergegraphList[[i]]<-set.graph.attribute(mergegraphList[[i]],"org",get.graph.attribute(graphList[[i]],"org"))
        mergegraphList[[i]]<-set.graph.attribute(mergegraphList[[i]],"title",get.graph.attribute(graphList[[i]],"title"))
        mergegraphList[[i]]<-set.graph.attribute(mergegraphList[[i]],"image",get.graph.attribute(graphList[[i]],"image"))
        mergegraphList[[i]]<-set.graph.attribute(mergegraphList[[i]],"link",get.graph.attribute(graphList[[i]],"link"))
   }
   }
   names(mergegraphList)<-names(graphList)
   if(simpleGraph==TRUE){
      mergegraphList<-getSimpleGraph(mergegraphList)
   }
   return(mergegraphList)
}
###################################################################################################
getSimpleGraph<-function(graphList){
     pathwayListLength<-length(graphList)
     pathwayList<-graphList
	 if(pathwayListLength>0){
     for(i in 1:pathwayListLength){
         deleteEdges<-c()
		 EMultiple<-c()
         eCount<-ecount(pathwayList[[i]])
         if(eCount>0){
		 if(any(is.loop(pathwayList[[i]]))||any(is.multiple(pathwayList[[i]]))){
		   #index<-seq(0,eCount-1)
		   #new!
            index<-seq(1,eCount)
		   deleteEdges<-c(deleteEdges,index[is.loop(pathwayList[[i]])])

		   EMultiple<-index[count.multiple(pathwayList[[i]])>1]
           EMultipleLength<-length(EMultiple)

           if(EMultipleLength>0){
               mulEdges<-get.edges(pathwayList[[i]],EMultiple)
               a<-unique(mulEdges)
               mid<-list()
               xLength<-length(a[,1])
			   a_string<-apply(a,1,function(x) paste(x,collapse=";"))
			   mulEdges_string<-apply(mulEdges,1,function(x) paste(x,collapse=";"))
               for(j in 1:xLength){
				  matched_mulEdges_string_index<-match(mulEdges_string,a_string[j])
                  mid[[j]]<-EMultiple[!is.na(matched_mulEdges_string_index)]
                       for(n in 2:length(mid[[j]])){
                          deleteEdges<-c(deleteEdges,mid[[j]][n])
                       }
               }####for(j in 1:xLength)
			   #print(Sys.time())
               edgeAttr<-list.edge.attributes(pathwayList[[i]])

			 edgeAttrLength<-length(edgeAttr)
             if(edgeAttrLength>0){
             for(k in 1:edgeAttrLength){
                 edgeAtt.k<-get.edge.attribute(pathwayList[[i]],edgeAttr[k])
                 str_k<-c()
                 str_k_index<-c()
                 for(j in 1:xLength){
                       str<-edgeAtt.k[mid[[j]][1]+1]
                       for(m in 2:length(mid[[j]])){
                          str<-paste(str,edgeAtt.k[mid[[j]][m]+1],sep=";")
                       }
					   str_k<-c(str_k,str)
					   str_k_index<-c(str_k_index,mid[[j]][1])
                       #pathwayList[[i]]<-set.edge.attribute(pathwayList[[i]],edgeAttr[k],mid[[j]][1],str)
                }####for(j in 1:xLength)

				pathwayList[[i]]<-set.edge.attribute(pathwayList[[i]],edgeAttr[k],str_k_index,str_k)
            }###for(k in 1:edgeAttrLength)
            }###if(edgeAttrLength>0)
           }###if(EMultipleLength>0)
		}####if(any(is.loop
        }####if(eCount+1>0)
        pathwayList[[i]]<-delete.edges(pathwayList[[i]],deleteEdges)
     }####for(i in 1:pathwayListLength)
	 }
     return(pathwayList)
}


#####################getUGraph##########################
getUGraph<-function(graphList,simpleGraph=TRUE){
    graphListLength<-length(graphList)
	if(graphListLength>0){
         for(i in 1:graphListLength){
              graphList[[i]]<-as.undirected(graphList[[i]],mode = "each")
        }
	}
	if(simpleGraph==TRUE)
	  graphList<-getSimpleGraph(graphList)
    return(graphList)
}

#####################mapNode##########################
mapNode<-function(graphList){
    newgraphList<-list()
    graphListLength<-length(graphList)
	if(graphListLength>0){
    for(i in 1:graphListLength){
	     org<-graphList[[i]]$org
		 nodeType<-""
		 if(org=="ec"){
		     nodeType<-"enzyme"
		 }else if(org=="ko"){
		     nodeType<-"ortholog"
		 }else{
		     nodeType<-"gene"
		 }
		 name<-c();id<-c();names<-c();type<-c();reaction<-c();link<-c()
         graphics_name<-c();graphics_fgcolor<-c();graphics_bgcolor<-c()
         graphics_type<-c();
         graphics_width<-c();graphics_height<-c();graphics_coords<-c()
		 Vcount<-vcount(graphList[[i]])
		 if(Vcount>0){
		 for(j in 1:Vcount){
		     #if(get.vertex.attribute(graphList[[i]],"type",j-1)==nodeType){
		     if(get.vertex.attribute(graphList[[i]],"type",j)==nodeType){
		         #node_name<-get.vertex.attribute(graphList[[i]],"names",j-1)
		         node_name<-get.vertex.attribute(graphList[[i]],"names",j)
		         expand_node_names<-unlist(strsplit(node_name,"[ ;]"))
				 if(org=="ko"){
	                 genes<-getGeneFromKO(expand_node_names)
				 }else if(org=="ec"){
	                 genes<-getGeneFromEnzyme(expand_node_names)
				 }else{
	                 genes<-getGeneFromKGene(expand_node_names)
				 }
				 if(length(genes)>0){
		             new_node_names<-paste(genes,collapse=" ")
				     #print(new_node_names)
                     names<-c(names,new_node_names)
                     type<-c(type,"gene")
                     link<-c(link,"unknow")
					 if(length(genes)>1){
                         graphics_name<-c(graphics_name,paste(genes[1],"...",sep=""))
					 }else{
					     graphics_name<-c(graphics_name,genes[1])
					 }
                     graphics_fgcolor<-c(graphics_fgcolor,"#000000")
                     graphics_bgcolor<-c(graphics_bgcolor,"#BFFFBF")
                }else{
                     #names<-c(names,get.vertex.attribute(graphList[[i]],"names",j-1))
                     #type<-c(type,get.vertex.attribute(graphList[[i]],"type",j-1))
                     #link<-c(link,get.vertex.attribute(graphList[[i]],"link",j-1))
                     #graphics_name<-c(graphics_name,get.vertex.attribute(graphList[[i]],"graphics_name",j-1))
                     names<-c(names,get.vertex.attribute(graphList[[i]],"names",j))
                     type<-c(type,get.vertex.attribute(graphList[[i]],"type",j))
                     link<-c(link,get.vertex.attribute(graphList[[i]],"link",j))
                     graphics_name<-c(graphics_name,get.vertex.attribute(graphList[[i]],"graphics_name",j))
                     graphics_fgcolor<-c(graphics_fgcolor,"#000000")
                     graphics_bgcolor<-c(graphics_bgcolor,"#FFFFFF")
                }

			}else{
                 #names<-c(names,get.vertex.attribute(graphList[[i]],"names",j-1))
                 #type<-c(type,get.vertex.attribute(graphList[[i]],"type",j-1))
                 #link<-c(link,get.vertex.attribute(graphList[[i]],"link",j-1))
                 #graphics_name<-c(graphics_name,get.vertex.attribute(graphList[[i]],"graphics_name",j-1))
                names<-c(names,get.vertex.attribute(graphList[[i]],"names",j))
                 type<-c(type,get.vertex.attribute(graphList[[i]],"type",j))
                 link<-c(link,get.vertex.attribute(graphList[[i]],"link",j))
                 graphics_name<-c(graphics_name,get.vertex.attribute(graphList[[i]],"graphics_name",j))
                 graphics_fgcolor<-c(graphics_fgcolor,"#000000")
                 graphics_bgcolor<-c(graphics_bgcolor,"#FFFFFF")
			}
		 }#end for(j in 1:Vcount)
		 }#end if(Vcount>0)
		 newgraphList[[i]]<-set.vertex.attribute(graphList[[i]],"names",value=names)
		 newgraphList[[i]]<-set.vertex.attribute(newgraphList[[i]],"type",value=type)
		 newgraphList[[i]]<-set.vertex.attribute(newgraphList[[i]],"link",value=link)
		 newgraphList[[i]]<-set.vertex.attribute(newgraphList[[i]],"graphics_name",value=graphics_name)
 		 newgraphList[[i]]<-set.vertex.attribute(newgraphList[[i]],"graphics_fgcolor",value=graphics_fgcolor)
		 newgraphList[[i]]<-set.vertex.attribute(newgraphList[[i]],"graphics_bgcolor",value=graphics_bgcolor)

		number<-get.graph.attribute(graphList[[i]],"number")
		new_org<-getOrgAndIdType()[1]
		new_idType<-getOrgAndIdType()[2]
	    newgraphList[[i]]<-set.graph.attribute(newgraphList[[i]],"name",paste("path:",new_org,number,sep=""))
        newgraphList[[i]]<-set.graph.attribute(newgraphList[[i]],"org",paste(new_org,new_idType,sep=";"))
        newgraphList[[i]]<-set.graph.attribute(newgraphList[[i]],"image",paste("http://www.genome.jp/kegg/pathway/",
		  new_org,"/",new_org,number,".png",sep=""))
        newgraphList[[i]]<-set.graph.attribute(newgraphList[[i]],"link",paste("http://www.genome.jp/kegg-bin/show_pathway?",
		  new_org,number,sep=""))
        names(newgraphList)[i]<-names(graphList)[i]
	}
	}
	return (newgraphList)
}


########################################################################################
simplifyGraph<-function(graphList,nodeType="geneProduct",directEdge=TRUE,verbose=FALSE){
     if(nodeType!="geneProduct"&&nodeType!="compound") stop("nodeType must be compound or geneProduct.")
     #if(!is.igraph(graphList)) stop("graphList should is a igraph class.")
     pathwayList<-graphList
     pathwayListLength<-length(pathwayList)
     graphList<-list()
     if(pathwayListLength>0){
    for(i in 1:pathwayListLength){
         if(class(pathwayList[[i]])!="igraph") stop(paste("graph ",i," must belong to a igraph class in graphList",sep=""))
         #if(is.directed(pathwayList[[i]])==FALSE) stop(paste("graph ",i," must be directed in graphList",sep=""))
   	     if(verbose==TRUE)
           print(paste("deal with the graph ",i," in ",pathwayListLength, " graphs",sep=""))
		 directed<-is.directed(pathwayList[[i]])
	     #Vcount<-vcount(pathwayList[[i]]) #new!
	     #if(Vcount>0){
		 enzymeType<-""
	     if(get.graph.attribute(pathwayList[[i]],"org")=="ec"){
	         enzymeType<-"enzyme"
	     }else if(get.graph.attribute(pathwayList[[i]],"org")=="ko"){
	         enzymeType<-"ortholog"
	     }else{
	         enzymeType<-"gene"
	     }
		 temp_nodeType<-""
		 if(nodeType=="geneProduct"){
		     temp_nodeType<-enzymeType
		 }else{
		     temp_nodeType<-"compound"
		 }
         vCount<-vcount(pathwayList[[i]])-1
		 #vCount<-vcount(pathwayList[[i]])
         enzyme_index<-c()#

         ID<-c();id<-c();names<-c();type<-c();reaction<-c();link<-c()
         graphics_name<-c();graphics_fgcolor<-c();graphics_bgcolor<-c()
         graphics_type<-c();graphics_x<-c();graphics_y<-c()
         graphics_width<-c();graphics_height<-c();graphics_coords<-c()
        if(vCount>0){
        for(j in 0:vCount){
        #for(j in 1:vCount){
            if(get.vertex.attribute(pathwayList[[i]],"type",j+1)==temp_nodeType){
                 enzyme_index<-c(enzyme_index,j+1)
                 id<-c(id,get.vertex.attribute(pathwayList[[i]],"id",j+1))
                 names<-c(names,get.vertex.attribute(pathwayList[[i]],"names",j+1))
                 type<-c(type,get.vertex.attribute(pathwayList[[i]],"type",j+1+1))
                 reaction<-c(reaction,get.vertex.attribute(pathwayList[[i]],"reaction",j+1))
                 link<-c(link,get.vertex.attribute(pathwayList[[i]],"link",j+1))
                 graphics_name<-c(graphics_name,get.vertex.attribute(pathwayList[[i]],"graphics_name",j+1))
                 graphics_fgcolor<-c(graphics_fgcolor,get.vertex.attribute(pathwayList[[i]],"graphics_fgcolor",j+1))
                 graphics_bgcolor<-c(graphics_bgcolor,get.vertex.attribute(pathwayList[[i]],"graphics_bgcolor",j+1))
                 graphics_type<-c(graphics_type,get.vertex.attribute(pathwayList[[i]],"graphics_type",j+1))
                 graphics_x<-c(graphics_x,get.vertex.attribute(pathwayList[[i]],"graphics_x",j+1))
                 graphics_y<-c(graphics_y,get.vertex.attribute(pathwayList[[i]],"graphics_y",j+1))
                 graphics_width<-c(graphics_width,get.vertex.attribute(pathwayList[[i]],"graphics_width",j+1))
                 graphics_height<-c(graphics_height,get.vertex.attribute(pathwayList[[i]],"graphics_height",j+1))
                 graphics_coords<-c(graphics_coords,get.vertex.attribute(pathwayList[[i]],"graphics_coords",j+1))
            }
          ID<-id
        }##for(j in 1:vCount)
        }
        ###
        start<-c();terminate<-c();middle<-list()
        eCount<-ecount(pathwayList[[i]])
        if(eCount>0){
		     flag<-enzymeType
             if(nodeType=="geneProduct"){
                 flag<-"compound"
             }
             VindexLength<-length(enzyme_index)
            if(VindexLength>0){
            for(j in 1:VindexLength){
                 outIndex<-c()#
                 outIndex<-as.integer(neighbors(pathwayList[[i]],enzyme_index[j],"out"))
                 outIndexLength<-length(outIndex)
			     Out<-c();In<-c();Mid<-list()
			     #
                 if(outIndexLength>0){
                     for(k in 1:outIndexLength){
				        outType<-get.vertex.attribute(pathwayList[[i]],"type",outIndex[k])
                        if(outType==flag){
                             Index<-as.integer(neighbors(pathwayList[[i]],outIndex[k],"out"))
                             IndexLength<-length(Index)
                             if(IndexLength>0){
                                 for(m in 1:IndexLength){

						             if((directed==TRUE&&Index[m]!=enzyme_index[j])||
								       (directed==FALSE&&Index[m]>enzyme_index[j])){
                                        if(unlist(strsplit(get.vertex.attribute(pathwayList[[i]],"id",Index[m]),"_"))[1]!=unlist(strsplit(get.vertex.attribute(pathwayList[[i]],"id",enzyme_index[j]),"_"))[1]&&get.vertex.attribute(pathwayList[[i]],"type",Index[m])==temp_nodeType){
										matched_enzyme_index<-match(Index[m],In)
										     if(is.na(matched_enzyme_index)){
	                                             Out<-c(Out,enzyme_index[j])
                                                 In<-c(In,Index[m])
										         Mid[[length(In)]]<-outIndex[k]
										     }else{
											     Mid[[matched_enzyme_index]]<-c(Mid[[matched_enzyme_index]],outIndex[k])
										     }
                                        }
							        }
                                }#end for(m in 1:IndexLength)
                            }###end if(IndexLength>0)
                        }else if(directEdge==TRUE&&outType==temp_nodeType){

						        if((directed==TRUE&&outIndex[k]!=enzyme_index[j])||
								 (directed==FALSE&&outIndex[k]>enzyme_index[j])){
                                    if(unlist(strsplit(get.vertex.attribute(pathwayList[[i]],"id",outIndex[k]),"_"))[1]!=unlist(strsplit(get.vertex.attribute(pathwayList[[i]],"id",enzyme_index[j]),"_"))[1]){
										   matched_enzyme_index<-match(outIndex[k],In)
										   if(is.na(matched_enzyme_index)){
	                                             Out<-c(Out,enzyme_index[j])
                                                 In<-c(In,outIndex[k])
												 Mid[[length(In)]]<--1
										   }
										   else{
											Mid[[matched_enzyme_index]]<-c(Mid[[matched_enzyme_index]],-1)
										   }
                                    }
							    }
                        }#end if(outType==flag)
                    }###for(k in 1:outIndexLength)
                }###if(outIndexLength>0)
			    start<-c(start,Out);terminate<-c(terminate,In);middle<-c(middle,Mid)
            }##for(j in 1:VindexLength)
            }###if(VindexLength>0)
        }##if(eCount>0)


        startLength<-length(start)
        if(startLength>0){
            entry1<-c();entry2<-c();Eid<-c();Enames<-c();Etype<-c();Ereaction<-c();Egraphics_name<-c()
            for(j in 1:startLength){
                 entry1[j]<-get.vertex.attribute(pathwayList[[i]],"id",start[j])
                 entry2[j]<-get.vertex.attribute(pathwayList[[i]],"id",terminate[j])

			    if(middle[[j]][1]!=-1){
			         Eid[j]<-get.vertex.attribute(pathwayList[[i]],"id",middle[[j]][1])
                     Enames[j]<-get.vertex.attribute(pathwayList[[i]],"names",middle[[j]][1])
                     Etype[j]<-get.vertex.attribute(pathwayList[[i]],"type",middle[[j]][1])
                     Ereaction[j]<-get.vertex.attribute(pathwayList[[i]],"reaction",middle[[j]][1])
			         Egraphics_name[j]<-get.vertex.attribute(pathwayList[[i]],"graphics_name",middle[[j]][1])
			    }else{
			         Eid[j]<-"unknow"
                     Enames[j]<-"unknow"
                     Etype[j]<-"unknow"
                     Ereaction[j]<-"unknow"
			         Ereaction[j]<-"unknow"
			         Egraphics_name[j]<-"unknow"
			    }
                middleJLength<-length(middle[[j]])
                if(middleJLength>1){
                     for(k in 2:middleJLength){
			             if(middle[[j]][k]!=-1){#20120709revised,
               #print(paste(middle[[j]][k],"ddd"))
               Eid[j]<-paste(Eid[j],";",get.vertex.attribute(pathwayList[[i]],"id",middle[[j]][k]),sep="")
               Enames[j]<-paste(Enames[j],";",get.vertex.attribute(pathwayList[[i]],"names",middle[[j]][k]),sep="")
               Etype[j]<-paste(Etype[j],";",get.vertex.attribute(pathwayList[[i]],"type",middle[[j]][k]),sep="")
               Ereaction[j]<-paste(Ereaction[j],";",get.vertex.attribute(pathwayList[[i]],"reaction",middle[[j]][k]),sep="")
			   Egraphics_name[j]<-paste(Egraphics_name[j],";",get.vertex.attribute(pathwayList[[i]],"graphics_name",middle[[j]][k]),sep="")
			            }else{
               Eid[j]<-paste(Eid[j],";","unknow",sep="")
               Enames[j]<-paste(Enames[j],";","unknow",sep="")
               Etype[j]<-paste(Etype[j],";","unknow",sep="")
               Ereaction[j]<-paste(Ereaction[j],";","unknow",sep="")
			   Egraphics_name[j]<-paste(Egraphics_name[j],";","unknow",sep="")
			            }

                    }
                }
            }##for(j in 1:startLength)
        }##end if(startLength>0)

#######################################
		if(startLength>0){
             vertex<-data.frame(ID=id,id=id,names=names,type=type,reaction=reaction,link=link,
		     graphics_name=graphics_name,graphics_fgcolor=graphics_fgcolor,graphics_bgcolor=graphics_bgcolor,
		     graphics_type=graphics_type,graphics_x=graphics_x,graphics_y=graphics_y,graphics_width=graphics_width,
		     graphics_height=graphics_height,graphics_coords=graphics_coords)
             edges<-data.frame(entry1=entry1,entry2=entry2,id=Eid,names=Enames,type=Etype,reaction=Ereaction,graphics_name=Egraphics_name)
		#edges<-unique(edges)
             graphList[[i]]<-graph.data.frame(edges,directed=directed,vertex)
        }else{
             graphList[[i]]<-graph.empty(n=0,directed=directed)
             graphList[[i]]<-add.vertices(graphList[[i]],length(id),name=id,id=id,names=names,type=type,
		         reaction=reaction,link=link,graphics_name=graphics_name,graphics_fgcolor=graphics_fgcolor,
		         graphics_bgcolor=graphics_bgcolor,graphics_type=graphics_type,graphics_x=graphics_x,
		         graphics_y=graphics_y,graphics_width=graphics_width,graphics_height=graphics_height,
                 graphics_coords=graphics_coords)
        }##else
	    graphList[[i]]<-set.graph.attribute(graphList[[i]],"name",get.graph.attribute(pathwayList[[i]],"name"))
        graphList[[i]]<-set.graph.attribute(graphList[[i]],"number",get.graph.attribute(pathwayList[[i]],"number"))
        graphList[[i]]<-set.graph.attribute(graphList[[i]],"org",get.graph.attribute(pathwayList[[i]],"org"))
        graphList[[i]]<-set.graph.attribute(graphList[[i]],"title",get.graph.attribute(pathwayList[[i]],"title"))
        graphList[[i]]<-set.graph.attribute(graphList[[i]],"image",get.graph.attribute(pathwayList[[i]],"image"))
        graphList[[i]]<-set.graph.attribute(graphList[[i]],"link",get.graph.attribute(pathwayList[[i]],"link"))

    }##for(i in 1:pathwayListLength)
    }
    names(graphList)<-names(pathwayList)
    return(graphList)
}
