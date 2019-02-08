library(miic)
library(cluster)
library(fpc)
library(NMF)
library(skmeans)
library(factoextra)
library(foreach)
library(doParallel)
library(fastICA)
library(wordcloud)
library(topicmodels)
library(Matrix)
library("tidyverse")
library("corrr")
library("igraph")
library("ggraph")


data("hematoData")
df=hematoData
View(df)
dim(df)
colnames(df)
rownames(df)

sparsity <-function(df){
  return (sum(df == 0)/(dim(df)[1]*dim(df)[2]))
}

sparsity(df)

###Definition tfidf ####

term.frequency=function(row){ ###Calcul TF
  row/sum(row)
}

inverse.doc.frequency=function(col){ ####Calcul IDF
  corpus.size=length(col)
  doc.count=length(which(col > 0))
  
  log10(corpus.size/(1+doc.count))
}

tf.idf=function(tf,idf){
  tf * idf
}

mat.df <- apply(df, 1, term.frequency)
dim(mat.df)
mat.idf <- apply(df, 2, inverse.doc.frequency)
dim(mat.idf)
mat.tfidf <-  apply(mat.df, 2, tf.idf, idf = mat.idf)
dim(mat.tfidf)
mat.tfidf=data.matrix(df)
set.seed(1337)
k.max <- 15
wss <- sapply(2:k.max, 
              function(k){skmeans(mat.tfidf, k,control = list(verbose = TRUE))$value})
wss
plot(2:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Value of the criterion")
findElbow(wss)

weight=Matrix(rep(1,dim(mat.tfidf)[1]*dim(mat.tfidf)[2]),nrow=dim(mat.tfidf)[1]);dim(weight)
sk=skmeans(mat.tfidf,3,control = list(verbose = TRUE));dim(sk$prototypes) # Skmeans sur le meilleur k et affichage de la dimension des prototypes

summary(sk$cluster)
quantile(sk$cluster);length(sk$cluster);length(sk$cluster[sk$cluster==1]);length(sk$cluster[sk$cluster==2]);length(sk$cluster[sk$cluster==3])

init <- nmfModel(3,mat.tfidf, W=1, H=sk$prototypes)###creation de la methode de seeding via les r??sultats du spherical kmeans
nmfSpherical=nmf(mat.tfidf,3,method="ls-nmf",.options="vt",seed=init,weight=as.matrix(weight)) ####NmfSpherical
res=nmfSpherical
snmf=extractFeatures(nmfSpherical,method="max")
plot(nmfSpherical)

res.coef <- coef(nmfSpherical)####on r??cup??re H
res.bas <- basis(nmfSpherical)####on r??cup??re W

##recuperation des basismap + consensusmap des matrices ?? r??aliser pour chaque nmf
png(filename = "/Users/jzk/Documents/M2/m??thodo/dumpR/test2.png", width = 450, height = 450,
    pointsize = 12, bg = "white",  res = NA)
par( mfrow = c( 1, 2 ) )
basismap(res)
coefmap(res)
dev.off()

#### T-SNE
library("Rtsne")
colors = rainbow(length(unique(sk$cluster)));colors
names(colors) = unique(sk$cluster);names
tsne <- Rtsne(mat.tfidf, dims = 2, perplexity=50, verbose=TRUE, max_iter = 500, check_duplicates = FALSE)
exeTimeTsne<- system.time(Rtsne(mat.tfidf, dims = 2, perplexity=50, verbose=TRUE, max_iter = 500, check_duplicates = FALSE))
plot(tsne$Y, t='n', main="tsne")
text(tsne$Y,labels=sk$cluster,col=colors[sk$cluster])
km = kmeans(t(mat.tfidf),3)


weight=Matrix(rep(1,dim(mat.tfidf)[1]*dim(mat.tfidf)[2]),nrow=dim(mat.tfidf)[1]);dim(weight)
sk=skmeans(t(mat.tfidf),3,control = list(verbose = TRUE));dim(sk$prototypes)


df=hematoData
df_cors <- df %>% correlate(method = "pearson") %>% stretch() %>% filter(abs(r) > .6)
df_cors$color <- "#00a9ff" #NEGATIF - bleu
df_cors$color[df_cors$r > 0 ] <- "#9c9c9c"

net <- graph_from_data_frame(df_cors,directed=TRUE)

net <- igraph::simplify(net, remove.loops = TRUE, remove.multiple = TRUE )
E(net)$size <- 20
E(net)$arrow.size <- 0
E(net)$color <- df_cors$color
E(net)$width <- (1 + df_cors$r)^2
E(net)$weight <- 1

V(net)$color[V(net)$name %in% names(sk$cluster[sk$cluster == 1])] <- "#ff2afb" # Bleu
V(net)$color[V(net)$name %in% names(sk$cluster[sk$cluster == 2])] <- "#ff3141" # Rouge
V(net)$color[V(net)$name %in% names(sk$cluster[sk$cluster == 3])] <- "#4f42ff" # pink
V(net)$size <- 10


plot(net)



library(ppcor)

df2=hematoData
#corr_df <- pcor(df2, method = "pearson")
#df_cors2 <- matrix(corr_df$estimate) %>% stretch() %>% filter(abs(r) > .4)

tesst <- cov(df2, method = "pearson")
diag(tesst) <- diag(tesst) + 0.1
m_o <- solve(tesst)
#m <- m_o / max(m_o)
m <- m_o * -1
#partiel_cor_temp <- data.frame(row=rownames(m)[row(m)], col=colnames(m)[col(m)], corr=c(m))
partiel_cor_temp <- data.frame(row=rownames(m)[row(m)[upper.tri(m)]],
                               col=colnames(m)[col(m)[upper.tri(m)]],
                               corr=m[upper.tri(m)])

# Normalisation
#partiel_cor_temp <- partiel_cor_temp %>% filter(corr < 2.5)
partiel_cor_temp$corr <- partiel_cor_temp$corr / max(partiel_cor_temp$corr)
#a <- partiel_cor_temp$corr
#partiel_cor_temp$corr <- 2*((a - min(a))/(max(a) - min(a))) - 1


boxplot(partiel_cor_temp$corr)
partiel_cor <- partiel_cor_temp %>% filter(abs(corr) > 0.4)
partiel_cor$color <- "#00a9ff" #NEGATIF - bleu
partiel_cor$color[partiel_cor$corr > 0 ] <- "#9c9c9c"

net <- graph_from_data_frame(partiel_cor,directed=TRUE)
#net <- igraph::simplify(net, remove.loops = TRUE, remove.multiple = TRUE )

E(net)$size <- 20
E(net)$arrow.size <- 0
E(net)$color <- partiel_cor$color
#E(net)$width <- (1 + partiel_cor$corr)^2
E(net)$width <- 1
E(net)$weight <- 1

V(net)$color[V(net)$name %in% names(sk$cluster[sk$cluster == 1])] <- "#ff2afb" # Bleu
V(net)$color[V(net)$name %in% names(sk$cluster[sk$cluster == 2])] <- "#ff3141" # Rouge
V(net)$color[V(net)$name %in% names(sk$cluster[sk$cluster == 3])] <- "#4f42ff" # pink
V(net)$size <- 10

plot(net)

#boxplot(m / 9.730427)





a <- partiel_cor_temp$corr

qsdsss <- 2*((a - min(a))/(max(a) - min(a))) - 1

boxplot(qsdsss)







library(corrplot)
mcor <- cor(df, method = "pearson")
corrplot(mcor, type="upper", order="hclust", tl.col="black", tl.srt=45)


col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = mcor, col = col, symm = TRUE)



oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(3)

sparsity <-function(df){
  return (sum(df == 0)/(dim(df)[1]*dim(df)[2]))
}





#####Miic
miic.res = miic(inputData = hematoData, latent = TRUE,
                confidenceShuffle = 10, confidenceThreshold = 0.001)
miic.plot(miic.res)
