library(R.matlab)
library(entropy)
library(data.table)
library(RColorBrewer)
library(scatterplot3d)
library(igraph)
library(lmtest)
#library(networkDynamic)
#library(ndtv)

# housekeeping
options(scipen = 999)

# set working directoy
setwd("/Users/mlr/Documents/git-projects/mamem-cascades/src/")

# helpers
degree.distribution <- function (graph, cumulative = FALSE, ...) 
{
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  cs <- degree(graph, ...)
  hi <- hist(cs, -1:max(cs), plot = FALSE)$count
  if (!cumulative) {
    res <- hi
  }
  else {
    res <- rev(cumsum(rev(hi)))
  }
  res
}

# parameters

patients <- c('001','002','003','004','005','006','007','008','009','010','011')
suffices <- c('ai','aii','bi','bii','ci','cii','di','dii','ei','eii') 

levels <- c(1,10,50,100) # configure the threshold levels to inspect here

nChannels <- 14 # configure the number of channels of the EEG measurement here

# preprocessing

stats <- data.frame(matrix(vector(), 6, 0))
for(lev in levels){
  dir.create(file.path("../output/", lev), showWarnings = FALSE)
  
  for(pat in patients){
    for(suf in suffices){
  
      #replace with path to EEG data
      mat_dat <- readMat(paste("../data/U",pat,suf,".mat",sep=''))
      
      tmp <- t(mat_dat$eeg)
      
      events<-list()
      specs<-list()
      for(i in 1:99){
        for(j in 1:nChannels){
          ds1 <- tmp[((i*178)-177):(i*178),j]
          specs[[j]]<-stats::spectrum(ds1, plot=FALSE)
          #specs[[j]]<-ds1
        }
        
        events[[paste(i,sep='')]]<-c()
        
        for(nChan in 1:nChannels){
          events[[paste(i,sep='')]] <- cbind(events[[paste(i,sep='')]],round(specs[[nChan]]$spec,1))
          #events[[paste(i,sep='')]] <- cbind(events[[paste(i,sep='')]],specs[[nChan]])
        }
      }
      
      nodes <- c()
      links <- c()
      rootsl <- list()
      
      matched<-list()
      
      for(j in 1:length(events)){
        interact<-list()
        lastInteract <- list()
        if(j>1){
          for(k in j-1:1){
            for(l in 1:nChannels){
              if(is.null(interact[[paste0(l)]])) interact[[paste0(l)]]<-l
              if(k<j && is.null(matched[[paste(j,l,sep='_')]])){
                #euclidian distance of the spectral densities
                #diff<-dist(rbind(unlist(events[[j]][,l]),unlist(events[[k]][,l])), method = "euclidean")
                #diff <- round(grangertest(unlist(events[[j]][,l]),unlist(events[[k]][,l]))$F[2])
                #print(diff)
                #mean of the raw delta of the spectral densities
                diff<-round(sum(data.frame(events[j])[,l]-data.frame(events[k])[,l]),0)
                
                #threshold for euclidian
                #if(abs(diff)<lev){
                #threshold for raw delta
                if(diff<=lev){
                  if(length(tail(which(links[,3]==l),1)) > 0){
                    if(links[tail(which(links[,3]==l),1),2]!=as.character(k)){
                      links<-rbind(links,c(links[tail(which(links[,3]==l),1),2],k,paste0(l),0.2))
                    }
                  }
                  links<-rbind(links,c(k,j,paste0(l),1))
                  rootsl[[paste0(l)]] <- c(k,paste0(l))
                  matched[[paste(j,l,sep='_')]]<-1
                  lastInteract[[paste0(l)]]<-j
                }
              }
            }
          }
        }
        nodes <- rbind(nodes, c(as.numeric(j),paste(interact,collapse=', '),as.numeric(j)))
      }
      
      roots<-do.call(rbind.data.frame, rootsl)
      colnames(nodes) <- c('node_id','tags','dpub')
      
      if(length(links)>0){
      #   colnames(links) <- c('source','target','tag')
      #   colnames(roots) <- c('root_node_id','tag')
      #   
      #   nd <- as.networkDynamic(network.initialize(0))
      #   set.network.attribute(nd,"vertex.pid","vertex.names")
      #   set.network.attribute(nd,"edge.pid","edge.names")
      #   for(i in 1:nrow(links)){
      #     fromN <- get.vertex.id(nd,unlist(links[i,1]))
      #     if(is.na(fromN)){
      #       add.vertices(nd,nv=1,vertex.pid=c(unlist(links[i,1])))
      #       fromN <- get.vertex.id(nd,unlist(links[i,1]))
      #       set.vertex.attribute(nd,'content',unlist(nodes[which(nodes[,1]==unlist(links[i,1])),2]),v=c(fromN))
      #       set.vertex.attribute(nd,'content2',paste(unlist(words300B[unlist(nodes[which(nodes[,1]==unlist(links[i,1])),1])]),collapse=' '),v=c(fromN))
      #       set.vertex.attribute(nd,'step',unlist(nodes[which(nodes[,1]==unlist(links[i,1])),1]),v=c(fromN))
      #       activate.vertices(nd,onset=as.numeric(unlist(links[i,2])),terminus=max(as.numeric(unlist(links[,2])))+1,v=c(fromN))
      #     }
      #     
      #     toN <- get.vertex.id(nd,unlist(links[i,2]))
      #     if(is.na(toN)){
      #       add.vertices(nd,nv=1,vertex.pid=c(unlist(links[i,2])))
      #       toN <- get.vertex.id(nd,unlist(links[i,2]))
      #       set.vertex.attribute(nd,'content',unlist(nodes[which(nodes[,1]==unlist(links[i,2])),2]),v=c(toN))
      #       set.vertex.attribute(nd,'content2',paste(unlist(words300B[unlist(nodes[which(nodes[,1]==unlist(links[i,2])),1])]),collapse=' '),v=c(toN))
      #       set.vertex.attribute(nd,'step',unlist(nodes[which(nodes[,1]==unlist(links[i,2])),1]),v=c(toN))
      #       activate.vertices(nd,onset=as.numeric(unlist(links[i,2])),terminus=max(as.numeric(unlist(links[,2])))+1,v=c(toN))
      #     }
      #     edgeID <- which(get.edge.attribute(nd,'ident')==paste(unlist(links[i,1]),unlist(links[i,2]),sep='-'))
      #     if(length(edgeID)==0){
      #       add.edges.active(nd,onset=as.numeric(unlist(links[i,2])), terminus=max(as.numeric(unlist(links[,2])))+1,head=toN,tail=fromN,names.eval=list(list('set','ident')),vals.eval=list(list(links[i,3][[1]],paste(unlist(links[i,1]),unlist(links[i,2]),sep='-'))))
      #     } else{
      #       linkLabel <- paste(get.edge.attribute(nd,'set',unlist=FALSE)[[edgeID]],unlist(links[i,3]),sep=", ")
      #       set.edge.attribute(nd, attrname='set', value=linkLabel, e=c(edgeID))
      #     }
      #   }
      #   
      #   compute.animation(nd, animation.mode = "kamadakawai", chain.direction=c('forward'),weight.dist=T,default.dist=3)
      #   
      #   #interactive
      #   render.d3movie(nd, filename=paste("../dynamic-network.html",sep=''),launchBrowser=T, 
      #                  displaylabels = T, label=nd %v% "vertex.names",
      #                  vertex.col="white",edge.col="darkgray",label.cex=.6,
      #                  vertex.cex = function(slice){ degree(slice)/10 }, vertex.border="#000000",
      #                  vertex.tooltip = paste("<span style='font-size: 10px;'><b>Slice:</b>", (nd %v% "step") , "<br>","<b>Matched characters:</b>", (nd %v% "content"), "<br>"),
      #                  edge.tooltip = paste("<b>Link:</b>", (nd %e% "set"),"</span>" ))
        
        colnames(links) <- c("id1","id2","label","weight")
        g <- graph.data.frame(links,directed=TRUE)
        V(g)$frame.color <- "white"
        V(g)$color <- "orange"
        
        E(g)$weight <- links[,4]
        
        qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
        col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
        
        colrs<-sample(col_vector, 14)
        E(g)$color <- colrs[as.numeric(E(g)$label)]
        E(g)$width <- E(g)$weight
        lay <- layout_on_grid(g)
        lay <- norm_coords(lay, ymin=-1, ymax=1, xmin=-1, xmax=1)
        pdf(paste("../output/",lev,"/",pat,suf,"_grid_net.pdf",sep=''),10,10)
        plot(g,rescale=F,vertex.size=5,vertex.label.cex=0.6,edge.label.cex=0.2,layout=lay,edge.arrow.size=0.1)
        dev.off()
        
        lay <- layout_with_dh(g)
        lay <- norm_coords(lay, ymin=-1, ymax=1, xmin=-1, xmax=1)
        
        pdf(paste("../output/",lev,"/",pat,suf,"_dh_net.pdf",sep=''),10,10)
        plot(g,rescale=F,vertex.size=5,vertex.label.cex=0.6,edge.label.cex=0.2,layout=lay,edge.arrow.size=0.1)
        dev.off()
        
        jpeg(paste("../output/",lev,"/",pat,suf,"_degree_distri.jpg",sep=''))
        plot(degree.distribution(g))
        dev.off()
        wtc <- cluster_walktrap(g)
        gstat <- data.frame(c(diameter(g),min(degree.distribution(g)),max(degree.distribution(g)),mean(degree.distribution(g)),edge_density(g),modularity(wtc)))
        
        tri <- mean(as.numeric(substr(matrix(triangles(g), nrow=3), nchar(matrix(triangles(g), nrow=3))-1+1, nchar(matrix(triangles(g), nrow=3)))))
        nods <- mean(as.numeric(substr(V(g)$name, nchar(V(g)$name)-1+1, nchar(V(g)$name))))
        
        
        jpeg(paste("../output/",lev,"/",pat,suf,"_triangles_nodes.jpg",sep=''))
        triangles <- cbind(tri,nods)
        colnames(triangles)<-c('triangles','all nodes')
        midpoints <- barplot(triangles, xlab="Mean unit value of node id")
        text(midpoints, 2, labels=triangles)
        dev.off()
        
        jpeg(paste("../output/",lev,"/",pat,suf,"_links_source_nrow.jpg",sep=''))
        plot(links[,2],c(1:nrow(links)),pch=".")
        dev.off()
        
        jpeg(paste("../output/",lev,"/",pat,suf,"_links_source_target.jpg",sep=''))
        plot(links[,2],links[,1],pch=".")
        dev.off()
        
        jpeg(paste("../output/",lev,"/",pat,suf,"_links_targets.jpg",sep=''))
        plot(links[,2],as.numeric(links[,3]),pch=".")
        for(ab in 1:20){
          if(ab %% 2 != 0) abline(v=(5*ab), col = "lightgray", lty = "dotted")
          else abline(v=(5*ab), col = "lightgray", lty = "dashed")
        }
        dev.off()
        write.csv(cbind(links[,2],as.numeric(links[,3])),file=paste('../output/',lev,'/',pat,suf,'_targets.txt',sep=''))
        
        jpeg(paste("../output/",lev,"/",pat,suf,"_links_sources.jpg",sep=''))
        plot(links[,1],as.numeric(links[,3]),pch=".")
        for(ab in 1:20){
          if(ab %% 2 != 0) abline(v=(5*ab), col = "lightgray", lty = "dotted")
          else abline(v=(5*ab), col = "lightgray", lty = "dashed")
        }
        dev.off()
        write.csv(cbind(links[,1],as.numeric(links[,3])),file=paste('../output/',lev,'/',pat,suf,'_sources.txt',sep=''))
        
        casc <- c()
        inter <- c()
        ent <- c()
        wien <- c()
        colnames(links) <- c('source','target','tag','weight')
        
        coordinates <- c()
        spec <- list()
        div=1
        
        for(z in 1:nrow(nodes)){
            #entropy
          if(z==1){
            ent <- rbind(ent,c(0,1,0,1))
          }
          if(length(links[which(links[,2]==nodes[z,1]),3])>0){
            
            inter <- rbind(inter,paste(sort(as.numeric(links[which(links[,2]==nodes[z,1]),3])), collapse=', '))
            nextI <- paste(sort(as.numeric(links[which(links[,2]==nodes[z,1]),3])), collapse=', ')
            
            if(length(spec)==0){
              coordinates <- rbind(coordinates,c(nodes[z,1],0,0))
              spec[[nextI]] <- c(0,0)
            }
            else{
              if(is.null(spec[[nextI]])){
                spec[[nextI]] <- c(0,div)
                coordinates <- rbind(coordinates,c(nodes[z,1],spec[[nextI]][1],div))
                div <- div+1
              }else{
                spec[[nextI]] <- c(spec[[nextI]][1]+1,spec[[nextI]][2])
                coordinates <- rbind(coordinates,c(nodes[z,1],spec[[nextI]][1],spec[[nextI]][2]))
              }
            }
            
            interact <- list()
            for(v in 1:nrow(inter)){
              interactions <- strsplit(unlist(inter[v,1]),', ')
              for( m in 1:length(interactions)){
                if(is.null(interact[[interactions[[1]][m]]])) interact[[interactions[[1]][m]]] <- 1
                else interact[[interactions[[1]][m]]] <- interact[[interactions[[1]][m]]] + 1
              }
            }
            df <- data.frame(unlist(interact))
            tmp<-df[,1]/colSums(df)
            df$loga<-log(tmp)
            df$piloga<-tmp*log(tmp)
            if(is.nan((-1*(colSums(df)[3]))/log(nrow(df)))){
              ent <- rbind(ent,c(entropy.empirical(df[,1], unit="log2"),1,-1*(colSums(df)[3]),1))
            } else{
              ent <- rbind(ent,c(entropy.empirical(df[,1], unit="log2"),(-1*(colSums(df)[3]))/log2(nrow(df)),-1*(colSums(df)[3]),(-1*(colSums(df)[3]))/log(nrow(df))))
            }
          }
            colnames(ent)<-c('empEntropy','evenness_log2','entropy','evenness')
        }
        colnames(coordinates) <- c("t","specificity","diversity")
        jpeg(paste("../output/",lev,"/",pat,suf,"_coordinates.jpg",sep=''))
        scatterplot3d(coordinates[,2],coordinates[,1],coordinates[,3],pch=16, highlight.3d=TRUE,type="h",xlab="Specificity",ylab="Node index",zlab="Diversity")
        dev.off()
        
        write.csv(ent,file=paste('../output/',lev,'/eeg_temp_stats_',pat,suf,'_entropy.txt',sep=''))
        jpeg(paste("../output/",lev,"/",pat,suf,"_entropy.jpg",sep=''))
        plot(ent[,1],type="l")
        for(ab in 1:20){
          if(ab %% 2 != 0) abline(v=(5*ab), col = "lightgray", lty = "dotted")
          else abline(v=(5*ab), col = "lightgray", lty = "dashed")
        }
        dev.off()
        jpeg(paste("../output/",lev,"/",pat,suf,"_evenness.jpg",sep=''))
        plot(ent[,2],type="l")
        for(ab in 1:20){
          if(ab %% 2 != 0) abline(v=(5*ab), col = "lightgray", lty = "dotted")
          else abline(v=(5*ab), col = "lightgray", lty = "dashed")
        }
        dev.off()
      } else{
        gstat <- data.frame(c('-','-','-','-','-','-'))
      }
      colnames(gstat) <- c(paste(pat,suf,sep=""))
      
      stats <- cbind(stats,gstat)
      rownames(stats) <- c('diameter','min degree','max degree','avg degree','density','modularity')
    }
  }
  
  write.csv2(stats,paste("../output/",lev,"/gstats.csv",sep=''),sep=';',row.names = F,col.names = T)
}