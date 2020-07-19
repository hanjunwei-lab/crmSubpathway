##' Internal function
##' @title internal function
##' @export
####################################################################
getPathway<-function(path,filelist=list.files(path),verbose=FALSE){
  #begintime<-Sys.time()
  mapList<-list()
  filelist_length<-length(filelist)
  jj<-0
  if(filelist_length>=1){
  for(j in 1:filelist_length){
	  if(verbose==TRUE)
         print(paste("deal with the pathway ",j," in ",filelist_length, " pathways",sep=""))
      #print(filelist[[j]])


	  top_temp<-tryCatch(xmlTreeParse(paste(path,filelist[[j]],sep=""),error=NULL),error=function(e) "error")
    if(class(top_temp)[1]=="character"){
	   if(verbose==TRUE){
         cat(paste("warning:the pathway ",filelist[j]," don't exist or has errors.\n Therefore, it is deleted from lists.\n",sep=""))
       }
	}
	else{
      top<-xmlRoot(top_temp)
      pathwaylist<-list()
      pathwayAttrs<-list()
      pathwayAttrs[[1]]<-xmlGetAttr(top,"name","unknow")
      pathwayAttrs[[2]]<-xmlGetAttr(top,"number","unknow")
      pathwayAttrs[[3]]<-xmlGetAttr(top,"org","unknow")
      pathwayAttrs[[4]]<-xmlGetAttr(top,"title","unknow")
      pathwayAttrs[[5]]<-xmlGetAttr(top,"image","unknow")
      pathwayAttrs[[6]]<-xmlGetAttr(top,"link","unknow")
      names(pathwayAttrs)<-c("name","number","org","title","image","link")
      entry<-list()
      relation<-list()
      reaction<-list()
      entryK<-0
      relationK<-0
      reactionK<-0
	  #tag1<-TRUE;tag2<-TRUE;tag3<-TRUE
      for(i in 1:length(top)){
	     #print(paste("entry",Sys.time(),sep=""))
		 top_i<-top[[i]]
		 xml_name<-xmlName(top_i)
         if(xml_name=="entry"){

             entryK<-entryK+1
             entry[[entryK]]<-getEntry(top_i)
			 #if(tag1==TRUE) {print(paste("entry",Sys.time(),sep=""));tag1=FALSE}
          }
         else if(xml_name=="relation"){
             relationK<-relationK+1
             relation[[relationK]]<-getRelation(top_i)
             #if(tag2==TRUE) {print(paste("lat",Sys.time(),sep=""));tag2=FALSE}

          }
          else if(xml_name=="reaction"){
		     #if(tag3==TRUE) {print(paste("reaction",Sys.time(),sep=""));tag3=FALSE}
             reactionK<-reactionK+1
             reaction[[reactionK]]<-getReaction(top_i)
          }
      }
      if(entryK==0){entry[[1]]<-"this file don't have entrys"}
      if(relationK==0){relation[[1]]<-getUnknowRelation()}
      if(reactionK==0){reaction[[1]]<-getUnknowReaction()}
      pathwaylist[[1]]<-pathwayAttrs
      pathwaylist[[2]]<-entry
      pathwaylist[[3]]<-relation
      pathwaylist[[4]]<-reaction
      names(pathwaylist)<-c("pathwayAttrs","entry","relation","reaction")
	  jj<-jj+1
	  mapList[[jj]]<-pathwaylist
	  names(mapList)[jj]<-filelist[j]
	 }
  }#end for filelist
  }
  else{
     stop("should at least input one pathway file.")
  }
  #print(Sys.time()-begintime)
  return(mapList)
}


######################################################
getGraphics<-function(graphics){
         graphicslist<-list()
         graphicslist[[1]]<-xmlGetAttr(graphics,"name","unknow")
         graphicslist[[2]]<-xmlGetAttr(graphics,"fgcolor","unknow")
         graphicslist[[3]]<-xmlGetAttr(graphics,"bgcolor","unknow")
         graphicslist[[4]]<-xmlGetAttr(graphics,"type","unknow")
         graphicslist[[5]]<-xmlGetAttr(graphics,"x","unknow")
         graphicslist[[6]]<-xmlGetAttr(graphics,"y","unknow")
         graphicslist[[7]]<-xmlGetAttr(graphics,"width","unknow")
         graphicslist[[8]]<-xmlGetAttr(graphics,"height","unknow")
         graphicslist[[9]]<-xmlGetAttr(graphics,"coords","unknow")
         names(graphicslist)<-c("name","fgcolor","bgcolor","type","x","y","width","height","coords")
         return(graphicslist)
}

###########################################################
getEntry<-function(entry){
     Enlist<-list()
     Childrenlength<-0
     componentId<-c()
     Enlist[[1]]<-xmlGetAttr(entry,"id","unknow")
     Enlist[[2]]<-xmlGetAttr(entry,"name","unknow")
     Enlist[[3]]<-xmlGetAttr(entry,"type","unknow")
     Enlist[[4]]<-xmlGetAttr(entry,"reaction","unknow")
     Enlist[[5]]<-xmlGetAttr(entry,"link","unknow")
     xml_children<-xmlChildren(entry)
     Childrenlength<-length(xml_children)
     Enlist[[6]]<-getGraphics(xml_children[[1]])
     if(Childrenlength>=2){
          for(i in 2:Childrenlength){
               componentId[i-1]<-xmlGetAttr(xml_children[[i]],"id","unknow")
          }
     }
     if(length(componentId)==0){componentId[1]<-"unknow"}
     Enlist[[7]]<-componentId
     names(Enlist)<-c("id","name","type","reaction","link","graphics","component")
     return(Enlist)
}

##############################################################
getSubtype<-function(subtype){
     subtypelist<-list()
     subtypelist[[1]]<-xmlGetAttr(subtype,"name","unknow")
     subtypelist[[2]]<-xmlGetAttr(subtype,"value","unknow")
     names(subtypelist)<-c("name","value")
     return(subtypelist)
}
##############################################################
getUnknowSubtype<-function(){
       subtypelist<-list("unknow","unknow")
       names(subtypelist)<-c("name","value")
       return(subtypelist)
}
##############################################################
getRelation<-function(relation){
     relationlist<-list()
     subtype<-list()
     subtypelength<-0
     relationlist[[1]]<-xmlGetAttr(relation,"entry1","unknow")
     relationlist[[2]]<-xmlGetAttr(relation,"entry2","unknow")
     relationlist[[3]]<-xmlGetAttr(relation,"type","unknow")
     xml_children<-xmlChildren(relation)
	 subtypelength<-length(xml_children)
	 if(subtypelength>=1){
         for(i in 1:subtypelength){
             subtype[[i]]<-getSubtype(xml_children[[i]])
         }
	 }
	 else{
	     subtype[[1]]<-getUnknowSubtype()
	 }
     relationlist[[4]]<-subtype
     names(relationlist)<-c("entry1","entry2","type","subtype")
     return(relationlist)
}

##############################################################
getUnknowRelation<-function(){
     relationlist<-list("unknow","unknow","unknow","unknow")
     names(relationlist)<-c("entry1","entry2","type","subtype")
     return(relationlist)
}

##############################################################
getReaction<-function(reaction){
     reactionlist<-list()
     substrate<-list()
     product<-list()
     childrenlength<-0
     sk<-0
     pk<-0
	 xml_children<-xmlChildren(reaction)
     childrenlength<-length(xml_children)
     if(childrenlength>=1){
         for(i in 1:childrenlength){
		      xml_children_i<-xml_children[[i]]
		      xml_name<-xmlName(xml_children_i)
              if(xml_name=="substrate"){
                  sk<-sk+1
				  substrate[[sk]]<-getSubstrate(xml_children_i)
              }
              if(xml_name=="product"){
                  pk<-pk+1
				  product[[pk]]<-getProduct(xml_children_i)
              }
          }
		  if(sk<1){substrate[[1]]<-getUnknowSubstrate()}
		  if(pk<1){product[[1]]<-getUnknowProduct()}
      }
	 else{
	     substrate[[1]]<-getUnknowSubstrate()
		 product[[1]]<-getUnknowProduct()
	 }
	 reactionlist[[1]]<-xmlGetAttr(reaction,"id","unknow")
     reactionlist[[2]]<-xmlGetAttr(reaction,"name","unknow")
     reactionlist[[3]]<-xmlGetAttr(reaction,"type","unknow")
     reactionlist[[4]]<-substrate
     reactionlist[[5]]<-product
     names(reactionlist)<-c("id","name","type","substrate","product")
     return(reactionlist)
}
##############################################################
getUnknowReaction<-function(){
     reactionlist<-list("unknow","unknow","unknow","unknow","unknow")
     names(reactionlist)<-c("id","name","type","substrate","product")
     return(reactionlist)
}
##############################################################
getSubstrate<-function(substrate){
     substratelist<-list()
     substratelist[[1]]<-xmlGetAttr(substrate,"id","unknow")
     substratelist[[2]]<-xmlGetAttr(substrate,"name","unknow")
     names(substratelist)<-c("id","name")
     return(substratelist)
}
##############################################################
getProduct<-function(product){
     productlist<-list()
     productlist[[1]]<-xmlGetAttr(product,"id","unknow")
     productlist[[2]]<-xmlGetAttr(product,"name","unknow")
     names(productlist)<-c("id","name")
     return(productlist)
}
##############################################################
getUnknowSubstrate<-function(){
       substratelist<-list("unknow","unknow")
       names(substratelist)<-c("id","name")
       return(substratelist)
}
##############################################################
getUnknowProduct<-function(){
       productlist<-list("unknow","unknow")
       names(productlist)<-c("id","name")
       return(productlist)
}


