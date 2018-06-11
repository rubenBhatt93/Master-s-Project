library(TDAmapper)
library(networkD3)
library(locfit)

liver <- read.table("C:/Users/Ujjal Bhattacharya/Desktop/TDA dataset/Liver.csv", 
                    header=TRUE, sep=",")
liver$Albumin_and_Globulin_Ratio = ifelse(is.na(liver$Albumin_and_Globulin_Ratio),
                                          ave(liver$Albumin_and_Globulin_Ratio, 
                                              FUN = function(x)mean(x, na.rm = TRUE)),
                                          liver$Albumin_and_Globulin_Ratio)
liver <- subset(liver,Gender=="Male")
rownames(liver)= 1:nrow(liver)
liver <- liver[-2]

#Building backward elimination
regressor = lm(formula= Dataset~Age+Total_Bilirubin+Direct_Bilirubin+Alkaline_Phosphotase
               +Alamine_Aminotransferase+Aspartate_Aminotransferase+Total_Protiens+
                 Albumin+Albumin_and_Globulin_Ratio,data=liver)
summary(regressor)

regressor = lm(formula= Dataset~Age+Direct_Bilirubin+Alkaline_Phosphotase
               +Alamine_Aminotransferase+Aspartate_Aminotransferase+Total_Protiens+
                 Albumin+Albumin_and_Globulin_Ratio,data=liver)
summary(regressor)

regressor = lm(formula= Dataset~Age+Direct_Bilirubin+Alkaline_Phosphotase
               +Alamine_Aminotransferase+Aspartate_Aminotransferase+Total_Protiens+
                 Albumin,data=liver)
summary(regressor)

regressor = lm(formula= Dataset~Age+Direct_Bilirubin+Alkaline_Phosphotase
               +Alamine_Aminotransferase+Total_Protiens+
                 Albumin,data=liver)
summary(regressor)

regressor = lm(formula= Dataset~Age+Direct_Bilirubin+Alkaline_Phosphotase
               +Alamine_Aminotransferase+
                 Albumin,data=liver)
summary(regressor)
#Direct_Bilirubin the most important followed by Age,Alamine,Alkaline and Albumin
#New data after removing features based in BE
liver = liver[-2]
liver = liver[-5]
liver = liver[-5]
liver = liver[-6]


liver$Dataset= ifelse(liver$Dataset==1,"Diseased","Non-Diseased")
liverBoth=liver
liver <- subset(liver,Dataset=="Diseased")

rownames(liver)= 1:nrow(liver)
liver<-liver[]
liverBeforeScaling = liver
#liver[1:5] = scale(liver[1:5])
mmnorm <-
  function (data,minval=0,maxval=1) 
  {
    
    d=dim(data)
    c=class(data)
    cnames=colnames(data)
    
    #remove classes from dataset
    classes=data[,d[2]]
    data=data[,-d[2]]
    
    minvect=apply(data,2,min)
    maxvect=apply(data,2,max)
    rangevect=maxvect-minvect
    zdata=scale(data,center=minvect,scale=rangevect)
    
    newminvect=rep(minval,d[2]-1)
    newmaxvect=rep(maxval,d[2]-1)
    newrangevect=newmaxvect-newminvect
    zdata2=scale(zdata,center=FALSE,scale=(1/newrangevect))
    zdata3=zdata2+newminvect
    
    zdata3=cbind(zdata3,classes)
    
    if (c=="data.frame") zdata3=as.data.frame(zdata3)
    colnames(zdata3)=cnames
    return(zdata3)
  }

liver=mmnorm(liver)
liver <- liver[-6]
indx <- sapply(liver, is.factor)
liver[indx] <- lapply(liver[indx], function(x) as.numeric(as.character(x)))

distMatrix = dist(liver,method="euclidean")
clusters <- hclust(distMatrix)
plot(clusters)


library(dynamicTreeCut)
cutreeDynamicTree(clusters)

clusterCut=cutree(clusters,3)
#liver$cluster_label = clusters$cluster
liverBeforeScaling=liverBeforeScaling[-6]
liverBeforeScaling$cluster_label = clusterCut
liver$cluster_label = clusterCut
liver$cluster_label= ifelse(liver$cluster_label==1,"Cluster 1",
                            ifelse(liver$cluster_label==2,"Cluster 2","Cluster 3"))
pcaData=princomp(liver[1:5])
plot(pcaData, type = "l")

#liverDf <- data.frame(pcaData$scores, "cluster" = factor(carsClusters))
#liverDf <- transform(liverDf, cluster_name = paste("Cluster",carsClusters))

library(devtools)
install_github("ggbiplot", "vqv")

library(ggbiplot)
g <- ggbiplot(pcaData, obs.scale = 1, var.scale = 1, 
              groups = liver$cluster_label, ellipse = TRUE, 
              circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)









# library(plotly)
# 
# p <- plot_ly(as.data.frame(pcaData), x = Comp.1 , y = Comp.2, text = rownames(pcaData),
#              mode = "markers", color = cluster_label, marker = list(size = 11)) 
# 
# p <- layout(p, title = "PCA Clusters from Hierarchical Clustering of Liver Patient Data", 
#             xaxis = list(title = "PC 1"),
#             yaxis = list(title = "PC 2"))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# pal <- c("red", "blue", "green")
# pal <- setNames(pal, c("cluster1", "cluster2", "cluster3"))
# 
# plot_ly(data = liverBeforeScaling, x = ~Direct_Bilirubin, y = ~Alamine_Aminotransferase, color = ~cluster_label, colors = pal)
# 
# plot_ly(data = liverBeforeScaling, x = ~Alamine_Aminotransferase, y = ~Albumin, color = ~cluster_label, colors = pal)
# 
# plot_ly(data = liverBeforeScaling, x = ~Alamine_Aminotransferase, y = ~Alkaline_Phosphotase, color = ~cluster_label, colors = pal)
# 
# plot_ly(data = liverBeforeScaling,x = ~Age, y = ~Alkaline_Phosphotase, color = ~cluster_label, colors = pal)
# 
# plot_ly(data = liverBeforeScaling,x = ~Age, y = ~Albumin, color = ~cluster_label, colors = pal)
# 
# plot_ly(data = liverBeforeScaling,x = ~Age, y = ~Alamine_Aminotransferase, color = ~cluster_label, colors = pal)
# 
# plot_ly(data = liverBeforeScaling,x = ~Age, y = ~Direct_Bilirubin, color = ~cluster_label, colors = pal)
# 
# plot_ly(data = liverBeforeScaling,x = ~Alkaline_Phosphotase, y = ~Albumin, color = ~cluster_label, colors = pal)
# 
# plot_ly(data = liverBeforeScaling,x = ~Alkaline_Phosphotase, y = ~Direct_Bilirubin, color = ~cluster_label, colors = pal)
# 
# plot_ly(data = liverBeforeScaling,x = ~Albumin, y = ~Direct_Bilirubin, color = ~cluster_label, colors = pal)
# 

#Finding mean of features of each cluster
#Cluster1
cluster1=subset(liverBeforeScaling,cluster_label==1)
cluster1=cluster1[-6]
cluster2=subset(liverBeforeScaling,cluster_label==2)
cluster2=cluster2[-6]
cluster3=subset(liverBeforeScaling,cluster_label==3)
cluster3=cluster3[-6]