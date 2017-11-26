library(R.matlab)
library(entropy)
library(data.table)
library(RColorBrewer)
library(scatterplot3d)
library(igraph)
library(lmtest)
library(tuneR)
library(poweRlaw)
library(plyr)

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

#patients <- c('S001','S002','S003','S004','S005','S006','S007','S008','S009','S010','S011')
patients <- c('S001')
suffices <- c('am')

levels <- c(50) # configure the threshold levels to inspect here
cohlev <- 0.99 #this is for spectral coherence 

nChannels <- 256 # configure the number of channels of the EEG measurement here
sRate <- 250
secs <- 2 # configure how long the slices shall be
# preprocessing

stats <- data.frame(matrix(vector(), 6, 0))
for(lev in levels){
  dir.create(file.path("../output/", lev), showWarnings = FALSE)
  
  for(pat in patients){
    for(suf in suffices){
      
      #replace with path to EEG data
      mat_dat <- readMat(paste("../data/",pat,suf,".mat",sep=''))
      
      tmp <- t(mat_dat$val)
      
      dPoints <- length(tmp[,1])
      
      events<-list()
      specs<-list()
      for(i in 1:floor(dPoints/(sRate*secs))){ #todo -> fix the slicing to make it generic based on sampling rate and signal length
        
        for(j in 1:nChannels){
          ds1 <- tmp[((i*sRate*secs)-((sRate*secs)-1)):(i*sRate*secs),j] # todo -> see above for generic slicing todo
          
          #for power spectra comparison use the power spectra of the slices as initialised below
          #wv<-Wave(ds1,samp.rate=sRate)
          #s1 <- wv@left
          #s1 <- s1 / 2^(wv@bit -1)
          #n <- length(s1)
          #p <- fft(s1)
          #nUniquePts <- ceiling((n+1)/2)
          #p <- p[1:nUniquePts] #select just the first half since the second half 
          #p <- abs(p)  #take the absolute value, or the magnitude 
          #p <- p / n #scale by the number of points so that
          #p <- p^2  # square it to get the power 
          #if (n %% 2 > 0){
          #  p[2:length(p)] <- p[2:length(p)]*2 # we've got odd number of points fft
          #} else {
          #  p[2: (length(p) -1)] <- p[2: (length(p) -1)]*2 # we've got even number of points fft
          #}
          #freqArray <- (0:(nUniquePts-1)) * (wv@samp.rate / n) #  create the frequency array 
          #plot(freqArray/1000, 10*log10(p), type='l', col='black', xlab='Frequency (kHz)', ylab='Power (dB)')
          #specs[[j]]<-10*log10(p)
          #old
          #specs[[j]]<-stats::spectrum(ds1, plot=FALSE)
          
          #for raw signal comparison use the raw data as initialized below
          specs[[j]]<-ds1
        }
        
        events[[paste(i,sep='')]]<-c()
        
        for(nChan in 1:nChannels){
          #for power spectra comparison use the power spectra of the slices as initialised below, old
          #events[[paste(i,sep='')]] <- cbind(events[[paste(i,sep='')]],round(specs[[nChan]]$spec,1))
          
          #for raw signal comparison use the raw data as initialized below
          events[[paste(i,sep='')]] <- cbind(events[[paste(i,sep='')]],specs[[nChan]])
        }
      }
      
      nodes <- c()
      links <- c()
      rootsl <- list()
      
      matched<-list()
      matchedSource<-list()
      
      for(j in 1:length(events)){
        interact<-list()
        if(j>1){
          for(k in 1:(j-1)){
            
            #setup parallel backend to use many processors
            #cores=detectCores()
            #cl <- makeCluster(cores[1]-1)
            #registerDoParallel(cl)
            #foreach(i=1:nChannels) %dopar% {
            for(l in 1:nChannels){
              if(is.null(matched[[paste(j,l,sep='_')]]) && is.null(matchedSource[[paste(k,l,sep='_')]])){
                #euclidian distance of the spectral densities or the raw signals
                #diff<-dist(rbind(unlist(events[[j]][,l]),unlist(events[[k]][,l])), method = "euclidean")
                
                #spectral coherence
                #diff<-stats::spectrum(cbind(unlist(events[[j]][,l]),unlist(events[[k]][,l])), plot=FALSE,spans=c(3,5))$coh
                diff<-max(ccf(unlist(events[[j]][,l]),unlist(events[[k]][,l]),plot = F)$acf)
                
                if(is.na(diff[1])){
                  diff <- 1
                }
                #compare the power spectra with the granger test
                #diff <- round(grangertest(unlist(events[[j]][,l]),unlist(events[[k]][,l]))$F[2])
                
                #sum of the raw delta of the spectral densities
                #diff<-round(sum(data.frame(events[j])[,l]-data.frame(events[k])[,l]),0)
                
                #threshold for euclidian
                #if(abs(diff)<lev){
                
                #threshold for granger test -> we are interested in F score around 1 for similarity
                #if(round(diff)==1){
                # for sectral coherence
                if(mean(diff)>=cohlev){
                  
                  #print(mean(diff))
                  interact[[paste0(l)]]<-l
                  
                  #if(length(tail(which(links[,3]==l),1)) > 0){
                  #  if(links[tail(which(links[,3]==l),1),2]!=as.character(k)){
                  #    links<-rbind(links,c(links[tail(which(links[,3]==l),1),2],k,paste0(l),0.2))
                  #  }
                  #}
                  links<-rbind(links,c(k,j,paste0(l),0.4))
                  rootsl[[paste0(l)]] <- c(k,paste0(l))
                  matched[[paste(j,l,sep='_')]]<-1
                  matchedSource[[paste(k,l,sep='_')]]<-1
                }
              }
            }
            #stop cluster
            #stopCluster(cl)
          }
        }
        nodes <- rbind(nodes, c(as.numeric(j),paste(interact,collapse=', '),as.numeric(j)))
      }
      
      roots<-do.call(rbind.data.frame, rootsl)
      colnames(nodes) <- c('node_id','tags','dpub')
      
      if(length(links)>0){
        links <- data.frame(source=unlist(links[,1]),target=unlist(links[,2]),tag=unlist(links[,3]),weight=unlist(links[,4]),stringsAsFactors = F)
        colnames(links) <- c('source','target','tag','weight')
        linksDelta <- as.integer(links$target)-as.integer(links$source)
        jpeg(paste0("../output/",lev,"/",pat,suf,"_links_delta.jpg"))
        plot(linksDelta,type='l')
        abline(h=mean(linksDelta),col="red")
        abline(h=median(linksDelta),col="blue")
        dev.off()
        
        linksDelta.count <- count(linksDelta)
        jpeg(paste0("../output/",lev,"/",pat,suf,"_links_delta_distri.jpg"))
        plot(linksDelta.count$x,linksDelta.count$freq,pch=20)
        dev.off()
        
        jpeg(paste0("../output/",lev,"/",pat,suf,"_links_delta_distri_loglog.jpg"))
        plot(linksDelta.count$x,linksDelta.count$freq,pch=20,log="xy")
        dev.off()
        
        #power law?
        
        m_bl = displ$new(linksDelta)
        est = estimate_xmin(m_bl)
        m_bl$setXmin(est)
        m_ln = dislnorm$new(linksDelta)
        est = estimate_xmin(m_ln)
        m_ln$setXmin(est)
        m_pois = dispois$new(linksDelta)
        est = estimate_xmin(m_pois)
        m_pois$setXmin(est)
        
        jpeg(paste0("../output/",lev,"/",pat,suf,"_links_delta_distri_plaw.jpg"))
        plot(m_bl, ylab="CDF")
        text(100,0.15,bquote(x[min] ~ .(paste0("=")) ~ .(m_bl$xmin) ~ .(paste0(", ")) ~ alpha ~ .(paste0("=")) ~ .(m_bl$pars)))
        lines(m_bl, col=2)
        lines(m_ln, col=3)
        lines(m_pois, col=4)
        dev.off()
        
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
        
        colrs<-sample(col_vector, 256, replace = T)
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
        dev.off()
        write.csv(cbind(links[,2],as.numeric(links[,3])),file=paste('../output/',lev,'/',pat,suf,'_targets.txt',sep=''))
        
        jpeg(paste("../output/",lev,"/",pat,suf,"_links_sources.jpg",sep=''))
        plot(links[,1],as.numeric(links[,3]),pch=".")
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
        dev.off()
        jpeg(paste("../output/",lev,"/",pat,suf,"_evenness.jpg",sep=''))
        plot(ent[,2],type="l")
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