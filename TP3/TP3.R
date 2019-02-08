library(miic)
library(cluster)
library(fpc)
library(Matrix)
library("tidyverse")
library("igraph")
library("ggraph")
library(bnlearn)
library(pcalg)
library(graph)
library(BiocManager)
library(qgraph)
library(dplyr)






?cosmicCancer
data("cosmicCancer");typeof(cosmicCancer)
colnames(cosmicCancer)
dim(cosmicCancer)
df=data.matrix(cosmicCancer) #on transforme les donn??es str en int pour la suite
head(df)
summary(df)






####Exploration
library(DataExplorer)
create_report(cosmicCancer,"/Users/jzk/Documents/M2/m??thodo/cosmicCancer.html")
plot_str(cosmicCancer, type="r")
introduce(cosmicCancer)
dev.off()
plot_intro(cosmicCancer)
plot_missing(cosmicCancer)
plot_bar(cosmicCancer)
plot_histogram(cosmicCancer)
plot_boxplot(cosmicCancer, by = "arr_delay")

config <- list(
  "introduce" = list(),
  "plot_str" = list(
    "type" = "diagonal",
    "fontSize" = 35,
    "width" = 1000,
    "margin" = list("left" = 350, "right" = 250)
  ),
  "plot_missing" = list(),
  "plot_histogram" = list(),
  "plot_qq" = list(sampled_rows = 1000L),
  "plot_bar" = list(),
  "plot_correlation" = list("cor_args" = list("use" = "pairwise.complete.obs")),
  #  "plot_prcomp" = list(),
  "plot_boxplot" = list(),
  "plot_scatterplot" = list(sampled_rows = 1000L)
)

create_report(dt,"cosmicCancerDt.htmlxxx",y="Ploidy")

library(corrplot)
mcor <- cor(dt, method = "pearson")
corrplot(mcor, type="upper", order="hclust", tl.col="black", tl.srt=45)


col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = mcor, col = col, symm = TRUE)








#####Premi??re Approche: HC

df2=hc(cosmicCancer) ###Marche pas, cherche pk les NA g??nent
df1=cosmicCancer[complete.cases(cosmicCancer),]
head(df1);is.na(df1)
df1$Ploidy <- sapply(df1$Ploidy, as.factor) ### On change le type de Ploidy qui pouvait faire bugger
dim(df1)
setdiff(rownames(cosmicCancer), rownames(df1)) ###on regarde quelles cellules ??taient vides
df1=df1[,apply(data.matrix(df1), 2, var, na.rm=TRUE) != 0] ###On vire les colonnes constantes
df2=hc(df1) #on applique enfin notre hillclimbing
adj2=amat(df2)
?qgraph.layout.fruchtermanreingold
net2 <- graph_from_adjacency_matrix(adj2)
Isolated = which(igraph::degree(net2)==0); net2= igraph::delete.vertices(net2, Isolated) # on vire les noeuds sans edges
e <- get.edgelist(net2,names=FALSE)
l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(net2),area=15*(vcount(net2)^1.9),repulse.rad=(vcount(net2)^3.1))
cluster=cbind(colnames(df1),as.integer(colnames(df1) == toupper(colnames(df1))))
cluster[162,2]=2;rownames(cluster)=cluster[,1]; cluster=cluster[,-c(1)];cluster

V(net2)$color[V(net2)$name %in% names(cluster[cluster == 0])] <- "#ff2afb" # Bleu
V(net2)$color[V(net2)$name %in% names(cluster[cluster == 1])] <- "#ff3141" # Rouge
V(net2)$color[V(net2)$name %in% names(cluster[cluster == 2])] <- "#4f42ff" # pink
V(net2)$size <- 10
E(net2)$size <- 15
E(net2)$arrow.size <- 0.2
E(net2)$width <- 1
E(net2)$weight <- 1
plot(net2,layout=l,vertex.size=4)




hub=hub_score(net2)
sort(hub$vector)
E(net2) [ to("Ploidy") ]
E(net2) [ from("Ploidy") ]
E(net2) [ from("tp53") ]
eb=edge_betweenness(net2, e = E(net2), directed = TRUE, weights = NULL)
index <- which(eb >= sort(eb, decreasing=T)[10], arr.ind=TRUE)
get.edgelist(net2)[index,]

vb=betweenness(net2, v = V(net2))
index <- which(vb >= sort(vb, decreasing=T)[10], arr.ind=TRUE)
V(net2)[index]







####Deuxieme Approche: PC

dt <- data.matrix(df1)-1
View(dt)
V <- colnames(dt)
summary(dt)

vect = c()
for(i in 1:162){
  vect = cbind(vect, length(unique(dt[,i])))
}
suffStat <- list(dm = dt, nlev = vect, adaptDF = FALSE)




pc.D <- pcalg::pc(suffStat,indepTest = disCItest, alpha = 0.01, labels = V, verbose = TRUE)
pc.D <- pcalg::pc(suffStat,indepTest = disCItest, alpha = 0.1, labels = V, verbose = TRUE)
pc.D <- pcalg::pc(suffStat,indepTest = disCItest, alpha = 0.0001, labels = V, verbose = TRUE)
pc.D <- pcalg::pc(suffStat,indepTest = disCItest, alpha = 0.9, labels = V, verbose = TRUE)
pc.D <- pcalg::pc(suffStat,indepTest = disCItest, alpha = 0.3, labels = V, verbose = TRUE)

amat=as(pc.D, "amat")
net2=graph_from_adjacency_matrix(amat)

Isolated = which(igraph::degree(net2)==0); net2= igraph::delete.vertices(net2, Isolated) # on vire les noeuds sans edges
e <- get.edgelist(net2,names=FALSE)
l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(net2),area=15*(vcount(net2)^1.9),repulse.rad=(vcount(net2)^3.1))
cluster=cbind(colnames(dt),as.integer(colnames(dt) == toupper(colnames(dt))))
cluster[162,2]=2;rownames(cluster)=cluster[,1]; cluster=cluster[,-c(1)];cluster

V(net2)$color[V(net2)$name %in% names(cluster[cluster == 0])] <- "#ff2afb" # Bleu
V(net2)$color[V(net2)$name %in% names(cluster[cluster == 1])] <- "#ff3141" # Rouge
V(net2)$color[V(net2)$name %in% names(cluster[cluster == 2])] <- "#4f42ff" # pink
V(net2)$size <- 10
E(net2)$size <- 15
E(net2)$arrow.size <- 0.2
E(net2)$width <- 1
E(net2)$weight <- 1
plot(net2,layout=l,vertex.size=4)




hub=hub_score(net2)
sort(hub$vector)
E(net2) [ to("Ploidy") ]
E(net2) [ from("Ploidy") ]
E(net2) [ from("tp53") ]
E(net2) [ to("tp53") ]

eb=edge_betweenness(net2, e = E(net2), directed = TRUE, weights = NULL)
index <- which(eb >= sort(eb, decreasing=T)[10], arr.ind=TRUE)
get.edgelist(net2)[index,]

vb=betweenness(net2, v = V(net2))
index <- which(vb >= sort(vb, decreasing=T)[10], arr.ind=TRUE)
V(net2)[index]






####Troisi??me question: MIIC
net2 = miic(inputData = cosmicCancer, categoryOrder = cosmicCancer_stateOrder, latent = TRUE,
                confidenceShuffle = 100, confidenceThreshold = 0.001)

miic.res = miic(inputData = cosmicCancer, categoryOrder = cosmicCancer_stateOrder, latent = TRUE,
                confidenceShuffle = 100, confidenceThreshold = 0.001)
miic.plot(net2)

miic.write.network.cytoscape(net2,file="/Users/jzk/Documents/M2/m??thodo/dumpgraph")

miic.write.style.cytoscape(file="/Users/jzk/Documents/M2/m??thodo/dumpstyle2")

amat=net2$adjMatrix


net2=graph_from_adjacency_matrix(amat)
net2 <- simplify(net2, remove.multiple = T, remove.loops = T)

Isolated = which(igraph::degree(net2)==0); net2= igraph::delete.vertices(net2, Isolated) # on vire les noeuds sans edges
e <- get.edgelist(net2,names=FALSE)
l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(net2),area=15*(vcount(net2)^1.9),repulse.rad=(vcount(net2)^3.1))
cluster=cbind(colnames(dt),as.integer(colnames(dt) == toupper(colnames(dt))))
cluster[162,2]=2;rownames(cluster)=cluster[,1]; cluster=cluster[,-c(1)];cluster

V(net2)$color[V(net2)$name %in% names(cluster[cluster == 0])] <- "#ff2afb" # Bleu
V(net2)$color[V(net2)$name %in% names(cluster[cluster == 1])] <- "#ff3141" # Rouge
V(net2)$color[V(net2)$name %in% names(cluster[cluster == 2])] <- "#4f42ff" # pink
V(net2)$size <- 10
E(net2)$size <- 15
E(net2)$arrow.size <- 0.2
E(net2)$width <- 1
E(net2)$weight <- 1
plot(net2,layout=l,vertex.size=4)
hub=hub_score(net2)
sort(hub$vector)
E(net2) [ to("Ploidy") ]
E(net2) [ from("Ploidy") ]
E(net2) [ from("tp53") ]
E(net2) [ to("tp53") ]

eb=edge_betweenness(net2, e = E(net2), directed = TRUE, weights = NULL)
index <- which(eb >= sort(eb, decreasing=T)[10], arr.ind=TRUE)
get.edgelist(net2)[index,]

vb=betweenness(net2, v = V(net2))
index <- which(vb >= sort(vb, decreasing=T)[10], arr.ind=TRUE)
V(net2)[index]

###################################
##########AUTRES METHODES##########
###################################


###FCI et RFCI
library(minet)
library(parmigene)

pc.D=fci(suffStat, indepTest = disCItest, alpha=0.8, labels= V,verbose = TRUE)
pc.D=fci(suffStat, indepTest = disCItest, alpha=0.1, labels= V,verbose = TRUE)
amat=as(pc.D, "amat")
net2=graph_from_adjacency_matrix(amat)

Isolated = which(igraph::degree(net2)==0); net2= igraph::delete.vertices(net2, Isolated) # on vire les noeuds sans edges
e <- get.edgelist(net2,names=FALSE)
l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(net2),area=15*(vcount(net2)^1.9),repulse.rad=(vcount(net2)^3.1))
cluster=cbind(colnames(dt),as.integer(colnames(dt) == toupper(colnames(dt))))
cluster[162,2]=2;rownames(cluster)=cluster[,1]; cluster=cluster[,-c(1)];cluster

V(net2)$color[V(net2)$name %in% names(cluster[cluster == 0])] <- "#ff2afb" # Bleu
V(net2)$color[V(net2)$name %in% names(cluster[cluster == 1])] <- "#ff3141" # Rouge
V(net2)$color[V(net2)$name %in% names(cluster[cluster == 2])] <- "#4f42ff" # pink
V(net2)$size <- 10
E(net2)$size <- 15
E(net2)$arrow.size <- 0.2
E(net2)$width <- 1
E(net2)$weight <- 1
plot(net2,layout=l,vertex.size=4)
hub=hub_score(net2)
sort(hub$vector)
E(net2) [ to("Ploidy") ]
E(net2) [ from("Ploidy") ]
E(net2) [ from("tp53") ]
E(net2) [ to("tp53") ]

eb=edge_betweenness(net2, e = E(net2), directed = TRUE, weights = NULL)
index <- which(eb >= sort(eb, decreasing=T)[10], arr.ind=TRUE)
get.edgelist(net2)[index,]

vb=betweenness(net2, v = V(net2))
index <- which(vb >= sort(vb, decreasing=T)[10], arr.ind=TRUE)
V(net2)[index]


####Aracne

data(syn.data)
mim <- build.mim(syn.data, estimator="spearman")
net <- mrnet(mim)
net=graph_from_adjacency_matrix(net)
plot(net)

mi =minet(dt, estimator="mi.empirical" )
mi  <- parmigene::knnmi.all(t(dt))
amat=mrnet(mi)
amat=aracne.a(mi, eps=10^-5)
net2=graph_from_adjacency_matrix(amat)



#####D??tection de communaut??s
fc <- fastgreedy.community(as.undirected(g))

#make colors for different communities
V(g)$color <- ifelse(membership(fc)==1,"red","blue")
plot(g)


print("LABEL PROPAGATION")
fc<-cluster_label_prop(g)

print("Leading Eigen")
fc<-cluster_leading_eigen(g)

print("SpinGlass")
fc<-cluster_spinglass(g, stop.temp = 0.05)

print("walktrap")
fc<-cluster_walktrap(g, steps=4)

library(MCL)
print("MCL")
adj<-get.adjacency(g)
fc<-mcl(adj,addLoops=TRUE)

