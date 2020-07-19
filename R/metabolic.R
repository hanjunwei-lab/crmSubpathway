# Internal function
getMetabolicGraph<-function(pathwayList,verbose=FALSE){
   if(!is.list(pathwayList)) stop("users must input a list type of data which include pathways")
   pathwayListLength<-length(pathwayList)
   graphList<-list()
   Na<-c()##Na
   if(pathwayListLength>0){
   for(t in 1:pathwayListLength){
   	      if(verbose==TRUE)
              print(paste("deal with the pathwayList ",t," in ",pathwayListLength, " pathwayLists",sep=""))
     ######
           #if(pathwayList[[t]]$pathwayAttrs$number
           Na[t]<-pathwayList[[t]]$pathwayAttrs$number
           entrylength<-0
           entrylength<-length(pathwayList[[t]]$entry)

           ID<-c();id<-c();names<-c();type<-c();reaction<-c();link<-c()
           graphics_name<-c();graphics_fgcolor<-c();graphics_bgcolor<-c()
           graphics_type<-c();graphics_x<-c();graphics_y<-c()
           graphics_width<-c();graphics_height<-c();graphics_coords<-c()

           entry1<-c();entry2<-c()
           idsize<-0;entry1size<-0;entry2size<-0
           unenzymeSize<-0;unenzyme_reaction<-c();enzyme_Id<-c()
           unsubstrateSize<-0;unsubstrate_name<-c();substrate_Id<-c()
           unproductSize<-0;unproduct_name<-c();product_Id<-c()
          #############################################################################################
		  j<-0
		  k<-0
		  map_id<-c()
          for(i in 1:entrylength){
		  	if(pathwayList[[t]]$entry[[i]]$type!="group"){
			   j<-j+1
               id[j]<-pathwayList[[t]]$entry[[i]]$id
               names[j]<-pathwayList[[t]]$entry[[i]]$name
               type[j]<-pathwayList[[t]]$entry[[i]]$type
               reaction[j]<-pathwayList[[t]]$entry[[i]]$reaction
               link[j]<-pathwayList[[t]]$entry[[i]]$link
               graphics_name[j]<-pathwayList[[t]]$entry[[i]]$graphics$name
               if(is.null(pathwayList[[t]]$entry[[i]]$graphics$fgcolor)){
                         graphics_fgcolor[j]<-"unknow"
               }
               else{
                         graphics_fgcolor[j]<-pathwayList[[t]]$entry[[i]]$graphics$fgcolor
               }
               if(is.null(pathwayList[[t]]$entry[[i]]$graphics$bgcolor)){
                         graphics_bgcolor[j]<-"unknow"
               }
               else{
                         graphics_bgcolor[j]<-pathwayList[[t]]$entry[[i]]$graphics$bgcolor
               }
               if(is.null(pathwayList[[t]]$entry[[i]]$graphics$type)){
                         graphics_type[j]<-"unknow"
               }
               else{
                         graphics_type[j]<-pathwayList[[t]]$entry[[i]]$graphics$type
               }
               graphics_x[j]<-pathwayList[[t]]$entry[[i]]$graphics$x
               graphics_y[j]<-pathwayList[[t]]$entry[[i]]$graphics$y
               if(is.null(pathwayList[[t]]$entry[[i]]$graphics$width)){
                         graphics_width[j]<-"unknow"
               }
               else{
                         graphics_width[j]<-pathwayList[[t]]$entry[[i]]$graphics$width
               }
               if(is.null(pathwayList[[t]]$entry[[i]]$graphics$height)){
                         graphics_height[j]<-"unknow"
               }
               else{
                         graphics_height[j]<-pathwayList[[t]]$entry[[i]]$graphics$height
               }
               graphics_coords[j]<-pathwayList[[t]]$entry[[i]]$graphics$coords

			   if(pathwayList[[t]]$entry[[i]]$type=="map"){
			   		       k<-k+1
			             map_id[k]<-pathwayList[[t]]$entry[[i]]$id
			   }
			 }#end "group"
          }#end entry
		  #entrylength_revised<-j
###############################################################
#######relationrelationrelationrelationrelation################

       relationlength<-0
       relationlength<-length(pathwayList[[t]]$relation)
       if(relationlength>0){
	         for(i in 1:relationlength){
               if(pathwayList[[t]]$relation[[i]]$type=="maplink"){

                     if(pathwayList[[t]]$relation[[i]]$entry1 %in% map_id){
                        entry1<-c(entry1,pathwayList[[t]]$relation[[i]]$entry1)
                        entry2<-c(entry2,pathwayList[[t]]$relation[[i]]$subtype[[1]]$value)
                     }
                     if(pathwayList[[t]]$relation[[i]]$entry2 %in% map_id){
                        entry1<-c(entry1,pathwayList[[t]]$relation[[i]]$subtype[[1]]$value)
                        entry2<-c(entry2,pathwayList[[t]]$relation[[i]]$entry2)
                    }
               }
			 }
        }
##################################################################################################################
####
    reactionlength<-0
    reactionlength<-length(pathwayList[[t]]$reaction)
    if((reactionlength>1)||(pathwayList[[t]]$reaction[[1]]$id!="unknow")){#revised
        for(i in 1:reactionlength){

            enzymeId<-c();enzymeSize<-0;
            for(j in 1:entrylength){
                 if(pathwayList[[t]]$entry[[j]]$id==pathwayList[[t]]$reaction[[i]]$id){
                      enzymeSize<-enzymeSize+1
                      enzymeId[enzymeSize]<-pathwayList[[t]]$entry[[j]]$id
                 }
            }
            #
            if(enzymeSize==0){
   			    enzymeSize<-enzymeSize+1
                unenzymeSize<-unenzymeSize+1
                enzymeId[enzymeSize]<-paste("ue",unenzymeSize,sep="")
                enzyme_Id[unenzymeSize]<-paste("ue",unenzymeSize,sep="")
                unenzyme_reaction[unenzymeSize]<-pathwayList[[t]]$reaction[[i]]$name
            }
        ##########################################################################################

            tempSubId<-c();subId<-c()
            substrateSize<-length(pathwayList[[t]]$reaction[[i]]$substrate)
            Snumber<-0
            S_number<-0
            if(pathwayList[[t]]$reaction[[i]]$substrate[[1]]$id=="unknow"){	#revised

                          unsubstrateSize<-unsubstrateSize+1
                          S_number<-S_number+1
                          tempSubId[S_number]<-paste("us",unsubstrateSize,sep="")
                          substrate_Id[unsubstrateSize]<-paste("us",unsubstrateSize,sep="")
                          unsubstrate_name[unsubstrateSize]<-pathwayList[[t]]$reaction[[i]]$substrate[[1]]$name

			}
            else{
                  for(j in 1:substrateSize){
                      Snumber<-Snumber+1
                      subId[Snumber]<-pathwayList[[t]]$reaction[[i]]$substrate[[j]]$id
                  }
            }

           substrateId<-c()
           substrateId<-c(subId,tempSubId)

        #################################################################################################

           tempProId<-c();proId<-c()
           productSize<-length(pathwayList[[t]]$reaction[[i]]$product)
           Pnumber<-0
           P_number<-0
           if(pathwayList[[t]]$reaction[[i]]$product[[1]]$id=="unknow"){#revised
                          unproductSize<-unproductSize+1
                          P_number<-P_number+1
                          tempProId[P_number]<-paste("up",unproductSize,sep="")
                          product_Id[unproductSize]<-paste("up",unproductSize,sep="")
                          unproduct_name[unproductSize]<-pathwayList[[t]]$reaction[[i]]$product[[1]]$name
		   }
           else{
                for(j in 1:productSize){
                    Pnumber<-Pnumber+1
                    proId[Pnumber]<-pathwayList[[t]]$reaction[[i]]$product[[j]]$id
                }
           }

          productId<-c()#
          productId<-c(proId,tempProId)
        #####################################################################################################
          #
          length_substrateId<-length(substrateId)
          length_productId<-length(productId)
          if(pathwayList[[t]]$reaction[[i]]$type=="reversible"){
                if(length_substrateId!=0 & length_productId!=0){
                     for(j in 1:enzymeSize){
                          entry1<-c(entry1,substrateId)
                          entry2<-c(entry2,rep(enzymeId[j],length_substrateId))
                          entry1<-c(entry1,rep(enzymeId[j],length_productId))
                          entry2<-c(entry2,productId)

                          entry1<-c(entry1,productId)
                          entry2<-c(entry2,rep(enzymeId[j],length_productId))
                          entry1<-c(entry1,rep(enzymeId[j],length_substrateId))
                          entry2<-c(entry2,substrateId)
                     }
                }
                if(length_substrateId!=0 & length_productId==0){
                      for(j in 1:enzymeSize){
                         entry1<-c(entry1,substrateId)
                         entry2<-c(entry2,rep(enzymeId[j],length_substrateId))

                         entry1<-c(entry1,rep(enzymeId[j],length_substrateId))
                         entry2<-c(entry2,substrateId)
                      }
                }
                if(length_substrateId==0 & length_productId!=0){
                     for(j in 1:enzymeSize){
                         entry1<-c(entry1,rep(enzymeId[j],length_productId))
                         entry2<-c(entry2,productId)

                         entry1<-c(entry1,productId)
                         entry2<-c(entry2,rep(enzymeId[j],length_productId))
                     }
                }
         }###if(pathwayList$reaction[[i]]$type=="reversible")
         else{
               if(length_substrateId!=0 & length_productId!=0){
                     for(j in 1:enzymeSize){
                        entry1<-c(entry1,substrateId)
                        entry2<-c(entry2,rep(enzymeId[j],length_substrateId))
                        entry1<-c(entry1,rep(enzymeId[j],length_productId))
                        entry2<-c(entry2,productId)
                     }
                }
                if(length_substrateId!=0 & length_productId==0){
                     for(j in 1:enzymeSize){
                        entry1<-c(entry1,substrateId)
                        entry2<-c(entry2,rep(enzymeId[j],length_substrateId))
                     }
                }
                if(length_substrateId==0 & length_productId!=0){
                     for(j in 1:enzymeSize){
                        entry1<-c(entry1,rep(enzymeId[j],length_productId))
                        entry2<-c(entry2,productId)
                     }
               }
        }###??if(pathwayList$reaction[[i]]$type=="reversible")
      }###for(i in 1:reactionlength)
##############################################################################################################
          id<-c(id,enzyme_Id,substrate_Id,product_Id)
          names<-c(names,rep("unknow",unenzymeSize),unsubstrate_name,unproduct_name)
          if(pathwayList[[t]]$pathwayAttrs$org=="ec"){
               type<-c(type,rep("enzyme",unenzymeSize),rep("compound",unsubstrateSize),rep("compound",unproductSize))
          }
          else if(pathwayList[[t]]$pathwayAttrs$org=="ko"){
               type<-c(type,rep("ortholog",unenzymeSize),rep("compound",unsubstrateSize),rep("compound",unproductSize))
          }
		  else{
               type<-c(type,rep("gene",unenzymeSize),rep("compound",unsubstrateSize),rep("compound",unproductSize))
          }
          reaction<-c(reaction,unenzyme_reaction,rep("unknow",unsubstrateSize),rep("unknow",unproductSize))
          link<-c(link,rep("unknow",unenzymeSize),rep("unknow",unsubstrateSize),rep("unknow",unproductSize))
          graphics_name<-c(graphics_name,rep("unknow",unenzymeSize),rep("unknow",unsubstrateSize),rep("unknow",unproductSize))
          graphics_fgcolor<-c(graphics_fgcolor,rep("#000000",unenzymeSize),rep("#000000",unsubstrateSize+unproductSize))
          graphics_bgcolor<-c(graphics_bgcolor,rep("#BFBFFF",unenzymeSize),rep("#FFFFFF",unsubstrateSize+unproductSize))
          graphics_type<-c(graphics_type,rep("rectangle",unenzymeSize),rep("circle",unsubstrateSize+unproductSize))
          graphics_x<-c(graphics_x,rep("unknow",unenzymeSize+unsubstrateSize+unproductSize))
          graphics_y<-c(graphics_y,rep("unknow",unenzymeSize+unsubstrateSize+unproductSize))
          graphics_width<-c(graphics_width,rep("46",unenzymeSize),rep("8",unsubstrateSize+unproductSize))
          graphics_height<-c(graphics_height,rep("17",unenzymeSize),rep("8",unsubstrateSize+unproductSize))
          graphics_coords<-c(graphics_coords,rep("unknow",unenzymeSize+unsubstrateSize+unproductSize))
          #component_id<-c(component_id,rep("unknow",unenzymeSize+unsubstrateSize+unproductSize))
    }##if(reactionlength>1)
##############################################################################################################################
######
     if(length(entry1)>0){
             vertices<-data.frame()

             vertices<-data.frame(ID=id,id=id,names=names,type=type,reaction=reaction,link=link,graphics_name=graphics_name,
                graphics_fgcolor=graphics_fgcolor,graphics_bgcolor=graphics_bgcolor,graphics_type=graphics_type,
                graphics_x=graphics_x,graphics_y=graphics_y,graphics_width=graphics_width,graphics_height=graphics_height,
                graphics_coords=graphics_coords)

             edges<-data.frame()
             edges<-data.frame(entry1,entry2)
             edges<-unique(edges)
             graphList[[t]]<-graph.data.frame(edges,directed=TRUE,vertices)
    }
    else{
             graphList[[t]]<-graph.empty(n=0,directed=TRUE)
             graphList[[t]]<-add.vertices(graphList[[t]],length(id),name=id,id=id,names=names,type=type,reaction=reaction,link=link,graphics_name=graphics_name,
             graphics_fgcolor=graphics_fgcolor,graphics_bgcolor=graphics_bgcolor,graphics_type=graphics_type,
             graphics_x=graphics_x,graphics_y=graphics_y,graphics_width=graphics_width,graphics_height=graphics_height,
             graphics_coords=graphics_coords)
    }
	         graphList[[t]]<-set.graph.attribute(graphList[[t]],"name",pathwayList[[t]]$pathwayAttrs$name)
             graphList[[t]]<-set.graph.attribute(graphList[[t]],"number",pathwayList[[t]]$pathwayAttrs$number)
             graphList[[t]]<-set.graph.attribute(graphList[[t]],"org",pathwayList[[t]]$pathwayAttrs$org)
             graphList[[t]]<-set.graph.attribute(graphList[[t]],"title",pathwayList[[t]]$pathwayAttrs$title)
             graphList[[t]]<-set.graph.attribute(graphList[[t]],"image",pathwayList[[t]]$pathwayAttrs$image)
             graphList[[t]]<-set.graph.attribute(graphList[[t]],"link",pathwayList[[t]]$pathwayAttrs$link)
  }#end pathwayList
  }
     names(graphList)<-Na
	 return(graphList)
}





