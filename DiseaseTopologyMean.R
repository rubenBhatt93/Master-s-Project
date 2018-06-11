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

liver.dist = dist(liver[,1:5],method = "manhattan")

library(ks)

filter.kde <- kde(liver[,1:5],H=diag(1,nrow = 5),
                  eval.points = liver[,1:5])$estimate

library(igraph)
#Good practice to choose number of intervals as square root of number of data
liver.mapper = mapper(
  dist_object = liver.dist,
  filter_values = filter.kde,
  num_intervals = 7,
  percent_overlap = 75,
  num_bins_when_clustering = 10)
liver.graph <- graph.adjacency(liver.mapper$adjacency, mode="undirected")
plot(liver.graph)

l = length(V(liver.graph))

Age.maj.vertex <-c()
Direct_Bilirubin.maj.vertex <- c() #c() returns a vector, a 1D array
Alkaline_Phosphotase.maj.vertex<- c()
Alamine_Aminotransferase.maj.vertex <- c()
Albumin.maj.vertex<-c()

filter.kde.vertex <- c()

# for (i in 1:l){
#   points.in.vertex <- liver.mapper$points_in_vertex[[i]]
#   Mean.in.vertex <- mean(liver$Direct_Bilirubin[points.in.vertex])
#   Direct_Bilirubin.maj.vertex <- c(Direct_Bilirubin.maj.vertex,
#                                    as.character(Mean.in.vertex))
#   Mean.in.vertex <- mean(liver$Age[points.in.vertex])
#   Age.maj.vertex <- c(Age.maj.vertex,
#                       as.character(Mean.in.vertex))
#   Mean.in.vertex<- mean(liver$Alamine_Aminotransferase[points.in.vertex])
#   Alamine_Aminotransferase.maj.vertex <- c(Alamine_Aminotransferase.maj.vertex,
#                                            as.character(Mean.in.vertex))
#   Mean.in.vertex<- mean(liver$Alkaline_Phosphotase[points.in.vertex])
#   Alkaline_Phosphotase.maj.vertex <- c(Alkaline_Phosphotase.maj.vertex,
#                                        as.character(Mean.in.vertex))
#   Mean.in.vertex<- mean(liver$Albumin[points.in.vertex])
#   Albumin.maj.vertex <- c(Albumin.maj.vertex,
#                           as.character(Mean.in.vertex))
#   
#   filter.kde.vertex <- c(filter.kde.vertex,mean(filter.kde[points.in.vertex]))
# }

for (i in 1:l){
  points.in.vertex <- liver.mapper$points_in_vertex[[i]]
  
  Mean.in.vertex <- mean(liverBeforeScaling$Direct_Bilirubin[points.in.vertex])
  Direct_Bilirubin.maj.vertex <- c(Direct_Bilirubin.maj.vertex,
                                   as.character(Mean.in.vertex))
  Mean.in.vertex <- mean(liverBeforeScaling$Age[points.in.vertex])
  Age.maj.vertex <- c(Age.maj.vertex,
                      as.character(Mean.in.vertex))
  Mean.in.vertex<- mean(liverBeforeScaling$Alamine_Aminotransferase[points.in.vertex])
  Alamine_Aminotransferase.maj.vertex <- c(Alamine_Aminotransferase.maj.vertex,
                                           as.character(Mean.in.vertex))
  Mean.in.vertex<- mean(liverBeforeScaling$Alkaline_Phosphotase[points.in.vertex])
  Alkaline_Phosphotase.maj.vertex <- c(Alkaline_Phosphotase.maj.vertex,
                                       as.character(Mean.in.vertex))
  Mean.in.vertex<- mean(liverBeforeScaling$Albumin[points.in.vertex])
  Albumin.maj.vertex <- c(Albumin.maj.vertex,
                          as.character(Mean.in.vertex))
  
  filter.kde.vertex <- c(filter.kde.vertex,mean(filter.kde[points.in.vertex]))
}

vertex.size <- rep(0,l)

for (i in 1:l){
  points.in.vertex <- liver.mapper$points_in_vertex[[i]]
  vertex.size[i] <- length((liver.mapper$points_in_vertex[[i]]))
}

MapperNodes <- mapperVertices(liver.mapper, 1:nrow(liver))
MapperNodes$Age.maj.vertex <- as.factor(Age.maj.vertex)
MapperNodes$Direct_Bilirubin.maj.vertex <- as.factor(Direct_Bilirubin.maj.vertex)
MapperNodes$Alkaline_Phosphotase.maj.vertex <- as.factor(Alkaline_Phosphotase.maj.vertex)
MapperNodes$Alamine_Aminotransferase.maj.vertex <- as.factor(Alamine_Aminotransferase.maj.vertex)
MapperNodes$Albumin.maj.vertex <- as.factor(Albumin.maj.vertex)
MapperNodes$filter.kde <- filter.kde.vertex
MapperNodes$Nodesize <- vertex.size

MapperLinks <- mapperEdges(liver.mapper)

#ColourScale <- 'd3.scaleOrdinal()
#            .domain(["Diseased", "Non-Diseased"])
#.range(["red", "black"]);'

forceNetwork(Nodes = MapperNodes, Links = MapperLinks, 
             Source = "Linksource", Target = "Linktarget",
             Value = "Linkvalue", NodeID = "Nodename",
             Group = "filter.kde", opacity = 1, 
             linkDistance = 55, charge = -8,legend = TRUE,
             Nodesize = "Nodesize")

forceNetwork(Nodes = MapperNodes, Links = MapperLinks, 
             Source = "Linksource", Target = "Linktarget",
             Value = "Linkvalue", NodeID = "Nodename",
             Group = "Direct_Bilirubin.maj.vertex", opacity = 1, 
             linkDistance = 55, charge = -8,legend = TRUE,
             Nodesize = "Nodesize")

forceNetwork(Nodes = MapperNodes, Links = MapperLinks, 
             Source = "Linksource", Target = "Linktarget",
             Value = "Linkvalue", NodeID = "Nodename",
             Group = "Alkaline_Phosphotase.maj.vertex", opacity = 1, 
             linkDistance = 55, charge = -8,legend = TRUE,
             Nodesize = "Nodesize")

forceNetwork(Nodes = MapperNodes, Links = MapperLinks, 
             Source = "Linksource", Target = "Linktarget",
             Value = "Linkvalue", NodeID = "Nodename",
             Group = "Alamine_Aminotransferase.maj.vertex", opacity = 1, 
             linkDistance = 55, charge = -8,legend = TRUE,
             Nodesize = "Nodesize")

forceNetwork(Nodes = MapperNodes, Links = MapperLinks, 
             Source = "Linksource", Target = "Linktarget",
             Value = "Linkvalue", NodeID = "Nodename",
             Group = "Albumin.maj.vertex", opacity = 1, 
             linkDistance = 55, charge = -8,legend = TRUE,
             Nodesize = "Nodesize")

forceNetwork(Nodes = MapperNodes, Links = MapperLinks, 
             Source = "Linksource", Target = "Linktarget",
             Value = "Linkvalue", NodeID = "Nodename",
             Group = "Age.maj.vertex", opacity = 1, 
             linkDistance = 55, charge = -8,legend = TRUE,
             Nodesize = "Nodesize")

  vertexNo <- liver.mapper$points_in_vertex[[10]]
  new$Alamine_Aminotransferase=liverBeforeScaling$Alamine_Aminotransferase[vertexNo]
  new$Age=liverBeforeScaling$Age[vertexNo]
  new$Albumin=liverBeforeScaling$Albumin[vertexNo]
  new$Direct_Bilirubin=liverBeforeScaling$Direct_Bilirubin[vertexNo]
  new$Alkaline_Phosphotase=liverBeforeScaling$Alkaline_Phosphotase[vertexNo]
  new=as.data.frame(new)
  new=new[12:16]
  new=as.numeric(new)
  
  
  
  newLog=log(new)
  # print(AvgAge <- mean(liverBeforeScaling$Age[vertexNo]))
  # print(AvgDb <- mean(liverBeforeScaling$Direct_Bilirubin[vertexNo]))
  # print(AVgAlkaline<-mean(liverBeforeScaling$Alkaline_Phosphotase[vertexNo]))
  # print(AVgAlanine<-mean(liverBeforeScaling$Alamine_Aminotransferase[vertexNo]))
  # print(AVgAlbumin<-mean(liverBeforeScaling$Albumin[vertexNo]))

logFunct <- function(x) log10(x)
dim(data.frame(logFunct(liverBeforeScaling$Direct_Bilirubin[vertexNo])))
