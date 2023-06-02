rm(list=ls())
###############################
# Intersection over Union (IoU) metric
###############################
IOU<-function(tablecv){
  vf<-c(rep(NA,nclass))
  for (j in 1:nclass){
    TP<-diag(tablecv)[j]
    FP<-sum(tablecv[j,])-TP
    FN<-sum(tablecv[,j])-TP
    iou<-TP/(TP+FP+FN)
    vf[j]<-iou
  }
  return(vf)
}
###############################
IOU2<-function(tablecv){
  vf2<-vf<-c(rep(NA,nclass))
  TP<-diag(tablecv)
  for (j in 1:nclass){
    FP<-sum(tablecv[j,])-TP[j]
    FN<-sum(tablecv[,j])-TP[j]
    #iou<-TP[j]/(TP[j]+FP+FN)
    #vf[j]<-TP[j]
    vf2[j]<-(TP[j]+FP+FN)
  }
  return(sum(TP)/sum(vf2))
}
###############################

features<-c("Eigenvalue sum","Omnivariance","Eigenentropy","Anisotropy","Planarity","Linearity","PCA1","PCA2","Surface Variation","Sphericity","Verticality",
            "Nx","Ny","Nz")
#[2]<-  "Height"
#######################

# ldata list with:
# "df" data.frame with 834809 registers observed (rows) in 846 (columns)
# 14 features (fdata objects) 834809 registers observed (rows) in 60 scales (columns)
load(file="ldata.pinar.rda")
class(ldf.train)
dim(ldf.train$df) # 834809 846
names(ldf.train$df) # Multivariate data 
names(ldf.train) # functional curves
argsvals<-seq(1,60,by=2)
for (i in 1:14){
  ldf.train[[i+1]]<-ldf.train[[i+1]][,argsvals]
}
ldf.train$f1$argvals
dim(ldf.train$f1$data)
ldat <- ldf.train
# Muestra pequeña load(file="./data/ldf.train.RData")
######################
names(table(ldat$df$class))
# df <- ldat$df[,1:22]
# write.csv(df,file="pinar834809.csv")
# save(df,file="pinar834809.rda")
# #list.train <- list.test <- list()
# save(list.train,file="list.train.rda")
# save(list.test,file="list.test.rda")
# set.seed(1:5)

# CPU TIME
# RF feature selection fincor
# RF+BSP
# RF+BSP+VS
# Pointnet
# DGCNN
# library(fda.usc)
library(fda.usc.devel) # leer del 
##############################################
source("R/classif.ML.R")
source("R/classif.ML.vs.R")
classifKgroups <- fda.usc.devel:::classifKgroups
fdata2model<-fda.usc.devel:::fdata2model


library(randomForest)
classif <- "classif.randomForest"
tab <- table(ldat$df$parcela)
iparcela <- c(0,2,6)
sum(tab[-iparcela]) # 100k,250k,500k
prop.table(tab)
ii <-1:400000
plot(ldat$df[ii,1:2],col=ldat$df$parcela[ii]+1,asp=T)
ii <- ldat$df$parcela %in% c(0,2,6) 
points(ldat$df[ii,1:2],col=1)

#ii <- ldat$df$parcela %in% c(0,2,3,6) 
ii <- ldat$df$parcela %in% c(0,2,6) 
ldf.train <- subset(ldat,ii==T)
table(ldf.train$df$class)
lev <- levels(ldf.train$df$class)
lev

# sapply(ldf.train,dim)
ii <- ldat$df$parcela %in% c(1,4,5,3) 
ldf.test <- subset(ldat,ii==T)
# sapply(ldf.test,dim)

##############################################
# multiclass random forest
######################################
mdl <- c("RF","RF.BSP","RF.VS.BSP","RF.f4f9f11BSP")
dc <- 0.05
nclass <- 4   
labels <- 1:4 # etiquetas clases
nclass <- 4
grid <- expand.grid(c(1e4,2e4),c(1e5,2e5,4e5))
grid <- expand.grid(c(5e3,1e4,2e4),c(1e5,2e5,4e5))
grid <- expand.grid(c(5e3,1e4),c(1e5,2e5))

grid
ngrid <- nrow(grid) 
nbsp <- 7
nmodel <- length(mdl)
balaccu<-array(NA,dim=c(nclass,ngrid,nmodel)) #balanced accuracy
dimnames(balaccu)<-list(labels,1:ngrid,mdl)

dfiou<-array(NA,dim=c(nclass,ngrid,nmodel)) #balanced accuracy
dimnames(dfiou)<-list(labels,1:ngrid,mdl)
overall    <-  array(NA,dim=c(3,ngrid,nmodel)) #balanced accuracy
dimnames(overall)<-list(c(" Accuracy","Kappa","IOU"),1:ngrid,mdl)

tt3 <- tt4 <- matrix(NA,ngrid,nmodel)
tt3
vsBSP <- matrix(NA,ngrid,nmodel)
rownames(vsBSP)<-1:ngrid

nam.f <- names(ldf.train)[-1]
nam.vs <- nam.f[c(3,4,5,7,9,11)]
j<-4
set.seed(j)
ipc <- 4
#####################################################
classif <- "classif.randomForest"
library(caret)
lord <- list("VS.BSP"=NULL)
nlev <- length(lev)
list.train <- list.test <- list()
library(tictoc)
for (i in 1:ngrid){
  # i <- 1
  print("iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii")
  print(i)
  table(ldf.train$df$parcela)
  #ldf.train$df$class <- factor(ldf.train$df$class)
  
  ntrain <- grid[i,1]
  ntest <- grid[i,2]

   set.seed(1) # Faltaría un for para repetirlo para otras semillas
  #ltrain <- ldf.train[sample(nrow(ldf.train$df),ntrain),row=T]
  # ltest <- ldf.test[sample(nrow(ldf.test$df),ntest),row=T]
  ij <- NULL
  for (j in 1:nlev)
    ij <- c(ij,sample(which(ldf.train$df$class==lev[j]),ntrain/nlev))
  list.train[[i]] <- ij
  ltrain <- ldf.train[ij,row=T]
  #ij <- NULL
  # for (j in 1:nlev)
  #   ij <- c(ij,sample(which(ldf.test$df$class==lev[j]),ntest/nlev))
  set.seed(1)
  ij <- sample(nrow(ldf.test$df),ntest)
  ltest <- ldf.test[ij,row=T]
  list.test[[i]] <- ij

  sum(table(ltest$df$class))
  ######################### SIN DERIVADAS
  nfunc <- length(nam.f)
  nam.f <- names(ltrain)[-1]
  
  
  #####  Tiempos multivariantes
  tt1<-tictoc::tic()
  predictors <- names(ltrain$df[-c(1,2,3,5,6,831)])
  # dfcor <- cor(ltrain$df[,c(predictors)])
  # diag(dfcor)<-0
  # cpu_time_mult <- TRUE
  # # searches through a correlation matrix and returns a vector of integers 
  # # corresponding to columns to remove to reduce pair-wise correlations 
  # # cutoff 0.85 (value for the pair-wise absolute correlation cutoff)
  # tt1<-tic()
  # hc <- findCorrelation(dfcor, cutoff=0.85, names=TRUE)
  # predictors_nc <- setdiff(predictors, hc) 
  # rm(dfcor)
  train <- ltrain$df
  test <- ltest$df
  table(train$class)
  #set.seed(j)
  mod_rf<- randomForest(x = ltrain$df[,predictors_nc], 
                        y = ltrain$df$class)
  tt2<- tictoc::toc()
  tt3[i,1]<-tt2$toc-tt2$tic
  
  # varImpPlot(mod_rf)
  tt1<-tic()
  pred <- predict(mod_rf,test)
  tt2<- tictoc::toc()
  tt4[i,1]<-tt2$toc-tt2$tic
  
  tab_rf<-table(pred,ltest$df$class)
  result<-confusionMatrix(tab_rf)
  balaccu[,i,1]<-data.frame(result$byClass)$Balanced.Accuracy
  dfiou  [,i,1]<-IOU(tab_rf)
  overall[,i,1] <- c(result$overall[1:2],"IOU"=IOU2(tab_rf))
  
  lpc3 <- lbsp <- list()
  rtt <- ltrain$f1$rangeval
  tt <- ltrain$f1$argvals
  
  salto <- length(nam.f)
  for (ii in 1:salto){
    #lbsp[[nam.f[ii]]] <-       create.pc.basis(ltrain[[nam.f[ii]]],ipc)
    lbsp[[nam.f[ii]]] <-       create.bspline.basis(rtt,nbsp)
  }
  covar <- nam.f
  form <- as.formula(paste0("class~altura+",paste0(nam.f[1:14],collapse="+")))
  tt1<-tic()
  mod <- classif.randomForest(form,data=ltrain,basis.x=lbsp)
  tt2<- tictoc::toc()
  tt3[i,2]<-tt2$toc-tt2$tic
  
  tt1<-tic()
  pred <- predict(mod, ltest)
  tt2<- tictoc::toc()
  tt4[i,2]<-tt2$toc-tt2$tic
  tab<-table(pred,ltest$df$class)
  result<-confusionMatrix(tab)
  balaccu[,i,2]<-data.frame(result$byClass)$Balanced.Accuracy
  dfiou[,i,2]<-IOU(tab)
  overall[,i,2] <- c(result$overall[1:2],"IOU"=IOU2(tab))
  rm(mod)
  #############
  # cpu time
  ldist <- as.list(ltrain$df[,c("class","altura")])
  #ldist <- c(ldist,ltrain[-1])
  ldist <- c(ldist,ltrain[nam.vs])
  ldist <- fda.usc.devel:::dist.list(ldist)
  
  # set.seed(j)
  tt1<-tic()
  # form <- as.formula(paste0("class~altura+",paste0(nam.f[c(3,5,7,9,11)],collapse="+")))
  #res1 <- classif.randomForest(form,data=ltrain,basis.x=lbsp)
  res1 <- classif.ML.vs(ltrain,"class",c("altura",nam.vs),
                         classif=classif, xydist=ldist
                         ,basis.x = lbsp,dcor.min = dc)
  tt2<- tictoc::toc()
  tt3[i,3]<-tt2$toc-tt2$tic
  tt1<-tic()
  pred <- predict(res1, ltest)
  tt2<- tictoc::toc()
  tt4[i,3]<-tt2$toc-tt2$tic
  tab<-table(pred,ltest$df$class)
  result<-confusionMatrix(tab)
  balaccu[,i,3]<-data.frame(result$byClass)$Balanced.Accuracy
  dfiou[,i,3]<-IOU(tab)
  overall[,i,3] <- c(result$overall[1:2],"IOU"=IOU2(tab))
  #vsBSP[i, ] <- res1$i.predictor
  #lord[["VS.BSP"]]<-c(lord[["VS.BSP"]],res1$ipredictor)
  
  #set.seed(j) 
  form <- as.formula("class~altura+f9+f11")
  tt1<-tic()
  mod <- classif.randomForest(form,data=ltrain,basis.x=lbsp)
  tt2<- tictoc::toc()
  tt3[i,4]<-tt2$toc-tt2$tic
  tt1<-tic()
  pred <- predict(mod, ltest)
  tt2<- tictoc::toc()
  tt4[i,4]<-tt2$toc-tt2$tic
  tab<-table(pred,ltest$df$class)
  result<-confusionMatrix(tab)
  # result$overall[1:2]
  balaccu[,i,4]<-data.frame(result$byClass)$Balanced.Accuracy
  dfiou[,i,4]<-IOU(tab)
  overall[,i,4] <- c(result$overall[1:2],"IOU"=IOU2(tab))
  print(overall[,i,])
  print(tt3[,4])
  print(tt4[,4])
}
#colnames(vs)<-nam.f
#####################################################
# Functional PCA Basis not included
# Derivatives of Functional Covariates not included
#####################################################
colMeans(tt3[c(1,4,7),])
colMeans(tt3[c(1,4,7)+1,])
colMeans(tt3[c(1,4,7)+2,])
colMeans(tt4[c(1:3),])
colMeans(tt4[c(1:3)+1,])
colMeans(tt4[c(1:3)+2,])
 
colMeans(tt3[c(1,3),])
colMeans(tt3[c(1,3)+1,])
colMeans(tt4[c(1:2),])
colMeans(tt4[c(1:2)+2,])

round(apply(balaccu,c(1,3),mean,na.rm=T),2)
round(apply(dfiou,c(1,3),mean,na.rm=T),2)
round(apply(dfiou,c(2,3),mean,na.rm=T),2)
colMeans(tt3)
colMeans(tt4)
tt3[2,]/tt3[1,]
tt4[2,]/tt4[1,]
apply(overall[1,,],2,mean);apply(overall[2,,],2,mean);apply(overall[3,,],2,mean)
round(overall[1,,],2);round(overall[2,,],2);round(overall[3,,],2)
out <- rbind(apply(overall,c(1,3),mean,na.rm=T),apply(balaccu,c(1,3),mean,na.rm=T),apply(dfiou,c(1,3),mean,na.rm=T))
round(out*100,1)
tab.nam<- paste0("output/",classif,"-ntrain",ntrain,"-dc",dc,"-4PC7BSP",j,"seed.txt")
tab.nam
# write.table(round(out,3),tab.nam,sep="\t",row.names=TRUE)
# 
# vs <- rbind(vsPC,colMeans(vsPC))
# rownames(vs)[5]<-"Means"
# vs <- rbind(vs,vsBSP,colMeans(vsBSP))
# rownames(vs)[10]<-"Means"
# tab.nam<- paste0("output/","VS",classif,"-ntrain",ntrain,"-dc",dc,"-4PC7BSP",j,"seed.txt")
# tab.nam
# colnames(vs)<-c("Altura",nam.f[1:14])
# vs
# write.table(vs,tab.nam,sep="\t",row.names=TRUE)
# 
# vs <- rbind(vsPCd,colMeans(vsPCd))
# rownames(vs)[5]<-"Means"
# vs <- rbind(vs,vsBSPd,colMeans(vsBSPd))
# rownames(vs)[10]<-"Means"
# tab.nam<- paste0("output/","VS",classif,"-ntrain",ntrain,"-dc",dc,"-4PC7BSPd",j,"seed.txt")
# tab.nam
# colnames(vs)<-c("Altura",nam.f)
# write.table(vs,tab.nam,sep="\t",row.names=TRUE)
# 
# barplot(table(unlist(lord)))



