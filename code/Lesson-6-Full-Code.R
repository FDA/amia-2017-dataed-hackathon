library(dendextend)
library(dendextendRcpp)
library(NMF)


getwd()
orig_data <- read.delim("Suppl_Figure6.cdt", 
                        header = FALSE, 
                        as.is = TRUE)

expressions <- orig_data[4:nrow(orig_data), 5:ncol(orig_data)]

mat <- matrix(as.numeric(unlist(expressions)), nrow = nrow(expressions),
              ncol = ncol(expressions))

dimnames(mat) <- list(orig_data[4:nrow(orig_data), 3],
                      orig_data[1, 5:ncol(orig_data)])

gweight <- as.numeric(orig_data[4:nrow(orig_data), 4])

eweight <- as.numeric(orig_data[3, 5:ncol(orig_data)])


sim<-function(x,y,w){
  
  if(sum(is.na(x+y))>0){
    complete <-cbind(x,y,w)[-which(is.na(x+y)),]
  }
  else{
    complete <-cbind(x,y,w)
  }
  
  return(cov.wt(complete[,1:2], wt= complete[,3],center=FALSE,cor=TRUE)$cor[1,2])
}


a.sim.mat<-matrix(ncol=ncol(mat), nrow=ncol(mat))
dimnames(a.sim.mat)<-list(colnames(mat),colnames(mat))

for(i in 1:ncol(mat)){
  for(j in 1:ncol(mat)){
    if(is.na(a.sim.mat[i,j])){
      a.sim.mat[i,j]<-sim(mat[,i], mat[,j], gweight)
      a.sim.mat[j,i]<-a.sim.mat[i,j]
    }
  }
} 

a.dist<-1-as.dist(a.sim.mat)


g.sim.mat<-matrix(ncol=nrow(mat), nrow=nrow(mat))
dimnames(g.sim.mat)<-list(rownames(mat),rownames(mat))

for(i in 1:nrow(mat)){
  for(j in 1:nrow(mat)){
    if(is.na(g.sim.mat[i,j])){
      g.sim.mat[i,j]<-sim(mat[i,], mat[j,], eweight)
      g.sim.mat[j,i]<-g.sim.mat[i,j]
    }
  }
} 

g.dist<-1-as.dist(g.sim.mat)


a.sim.mat[which(rownames(a.sim.mat)=="NormBreast2"), which(rownames(a.sim.mat)=="NormBreast3")]


g.sim.mat[which(rownames(g.sim.mat)=="116219 GATA3 GATA binding protein 3 Hs.169946 H72474 "), which(rownames(g.sim.mat)=="101778 GATA3 GATA binding protein 3 Hs.169946 R31441 ")]


a.hclust<-hclust(a.dist, method="average")

a.dend.plot<-as.dendrogram(rotate(a.hclust, order=match(1:122,a.hclust$order)))

g.hclust<-hclust(g.dist, method="average")

g.dend.plot <- as.dendrogram(rotate(g.hclust, order=match(1:552,g.hclust$order)))


col.arr <-rep("black", 122)
col.arr[which(labels(a.dend.plot)=="Norway 83-BE")]<-"red"
col.arr[which(labels(a.dend.plot)=="Norway FU02-BE")]<-"red"

png("Lesson6_arrDend.png", width=1920, height=960)
a.dend.plot %>% color_branches(col = col.arr) %>% color_labels(col = col.arr) %>% plot(main="Array Cluster Dendrogram")
dev.off()


col.genes <-rep("black", 552)
col.genes[which(labels(g.dend.plot)=="108924  **Homo sapiens cDNA FLJ34425 fis, clone HHDPC2008297 Hs.120638 N81017 ")]<-"red"
col.genes[which(labels(g.dend.plot)=="120444 N33 Putative prostate cancer tumor suppressor Hs.71119 H13424 ")]<-"red"

g.dend.plot <- set(g.dend.plot,"labels", substr(labels(g.dend.plot), start=1, stop=10))
png("Lesson6_geneDend.png", width=3600, height=1000)
g.dend.plot %>% assign_values_to_leaves_nodePar(c(rep(0.7, 552)), "lab.cex") %>% color_branches(col = col.genes)%>% color_labels(col = col.genes) %>% hang.dendrogram() %>% plot(main="Gene Cluster Dendrogram")
dev.off()


png("Lesson6_heatmap.png", width=960, height=2400)
aheatmap(mat, Rowv=g.dend.plot, Colv=a.dend.plot, col=c("green", "black", "red"))
dev.off()
