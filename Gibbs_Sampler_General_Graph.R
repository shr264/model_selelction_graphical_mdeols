
####----This code works for replicate 1 under scenario I, lines that need to be 
####----modified to run another replicate/scenario are commented accordingly

#setwd('') 

#####################################################################################################
#####################################################################################################
########------------------------Program to read data from QMSim-----------------------------#########
#####################################################################################################
#####################################################################################################

#--To read genotypes

#replace "Pop1_qtl_001.txt" by "Pop1_qtl_00j.txt", j=2,3,4 to run another replicate 

QTL1=read.table("Pop1_qtl_001.txt",header=FALSE) 
                     
QTL1=as.matrix(QTL1[,-1])
QTL.Matrix=function(QTL){
	template=seq(1,ncol(QTL)-1,2)
	A=matrix(0,nrow=nrow(QTL),ncol=ncol(QTL))
	for(i in 1:nrow(QTL)){
		for(j in template){
			A[i,j]=QTL[i,j]+QTL[i,j+1]-3
			        }
			     }
	A=A[,template]
	return(A)
				}



###---Reading pedigree, phenotypes and total additive effects     

#replace "Pop1_data_001.txt" by "Pop1_data_00j.txt", j=2,3,4 to run another replicate                                      
data1=read.table("Pop1_data_001.txt",header=T)    

W=QTL.Matrix(QTL=QTL1)
nloci=ncol(W)

#########################################################
###--------Create training and testing sets--------######
#########################################################

data1.Training=data1[-which(data1$G==3), ]
W.Train=W[-which(data1$G==3), ]
n.train=nrow(W.Train)

data1.Test=data1[which(data1$G==3),]

###############################################################
########---Functions to build the adjacency matrices----#######
###############################################################

Adj=function(Nodes,nloci){
  library(Matrix)
  Adj=as(matrix(0,nrow=nloci,ncol=nloci),"dgTMatrix")
  for(i in 1:nrow(Nodes)){
    Adj[Nodes[i,1],Nodes[i,2]]=1
  }
  return(Adj)
}


###########################################
######----Read adjacency matrix-----#######
###########################################

Adjo=read.csv("selgraph1.csv",header=FALSE)


############################################################################
###---Function to sample from full cond. in Concentration graph models---### 
#############----------------General graphs----------------------###########
############################################################################

library(BDgraph)
SampG=function(Delta,U,Adjo,g){
as.matrix(round(rgwish(n=1,adj.g=Adjo,b=Delta+1,D=tcrossprod(g)+U)[,,1],6))
}

################################################################################
################################################################################
###########------Gibbs sampler for decomposable graphs---------#################
################################################################################              
################################################################################


BibiGraph.W=function(y,nloci,Nsim,a,U,Tau,V,W,Adjo,Burnin){
          
   n=nrow(W)
   
   #####---Create arrays to save samples

   g=matrix(0,nrow=Nsim,ncol=nloci)
   norms=matrix(0,nrow=Nsim)
   gnorms=matrix(0,nrow=Nsim)
   res.var=matrix(0,nrow=Nsim)

  #####---Initial values of parameters

   library(mvtnorm) 
   g[1, ]=rmvnorm(1,mean= matrix(0,nrow=nloci),sigma=diag(nloci),method="chol")

   library(pscl)
   res.var[1]=rigamma(1,alpha=V/2,beta=Tau/2)
   
  
   ##########################################
   #####-----Start algorithm-----############
   ##########################################

   WTW=crossprod(W)
   WTy=crossprod(W,y)
   Ginv.hat=matrix(0,nrow=nloci,ncol=nloci)
   for(i in 2:Nsim){

     Ginv=SampG(Delta=a,U=U,Adjo=Adjo,g=g[i-1, ])
     norms[i]=norm(Ginv,type="F")
     cov=chol2inv(chol(Ginv + WTW/res.var[i-1]))
     g[i, ]=rmvnorm(1, mean=crossprod(cov,WTy)/res.var[i-1],sigma=as.matrix(cov),method="chol")
     gnorms[i]=norm(matrix(g[i, ]))
     res.var[i]=rigamma(1,alpha=(V+n)/2,beta=(Tau+crossprod(y-W%*%g[i, ]))/2)
     if(i > Burnin){
                 Ginv.hat=(Ginv.hat*(i-Burnin-1)+Ginv)/(i-Burnin)
                                     }
                          }
   
   PostmeanMar=colMeans(g[(Burnin+1):Nsim, ])
   
   ###--In the following lines, replace "1" by "j", j=2,3,4 when running a different replicate

   write.table(matrix(res.var[(Burnin+1):Nsim]),col.names=FALSE,row.names=FALSE,
   file="ResVar1.csv", sep=",")

   write.table(norms[(Burnin+1):Nsim],col.names=FALSE,row.names=FALSE,
   file="norms1.csv", sep=",")

   write.table(gnorms[(Burnin+1):Nsim],col.names=FALSE,row.names=FALSE,
   file="gnorms1.csv", sep=",")

   write.table(PostmeanMar,col.names=FALSE,row.names=FALSE,
   file="PostmeanMark1.csv", sep=",")

   write.table(as.matrix(Ginv.hat),col.names=FALSE,row.names=FALSE,
   file="OmegaHat1.csv", sep=",")
  
   return(list(PostmeanMar,res.var[Burnin:Nsim],norms[(Burnin+1):Nsim],Ginv.hat))

                               
}

U=diag(nloci)

Test=BibiGraph.W(y=data1.Training$Phen,nloci=nloci,Nsim=50000,a=10,U=U,Tau=7000,V=72,W=W.Train,Adjo=Adjo,Burnin=20000)

#--Predicted Breeding values
W.Test=W[which(data1$G==4), ]
PredBVTrainGEMB=W.Train%*%as.matrix(Test[[1]])
PredBVTestGEMB=W.Test%*%as.matrix(Test[[1]])

##################################################
###--------Predictive abilities---------------####
##################################################

Predabil.GB=cor(data1.Test$Phen,PredBVTestGEMB,method="pearson")

#########################################
###--correlation of breeding values---###
#########################################

corrBVGB=cor(data1.Test$QTL,PredBVTestGEMB,method="pearson")

#########################################################
###---Correlations of breeding values in training----####
#########################################################

corrBVTrainGB=cor(data1.Training$QTL,PredBVTrainGEMB,method="pearson")

Performance=matrix(c(Predabil.GB,corrBVGB,corrBVTrainGB))

#Replace Bayesperform1.csv by Bayesperformj.csv, j=2,3,4 when running another replicate.
#Moreover, ad "II" at the end when running scenario 2.

write.table(Performance,col.names=FALSE,row.names=FALSE,file="Bayesperform1.csv", sep=",")


################################################################
################################################################
######---------------Evaluating performance----------------#####
################################################################
################################################################


####-----Randomly setting  zeros------####
###--This is a slightly different version of the function above--###

Rand.Sparse=function(nloci,thrsld,delta,U,seed1,seed2){
  Nodes=Def.Graph.Band(band.size=nloci,nloci=nloci)[[1]]
  set.seed(seed2)
  matrix=rwish(n=1,p=nloci,b=nloci+1,D=diag(nloci))
  matrix=matrix(matrix,nrow=nloci,byrow=TRUE)
  tmp1=cbind(Nodes,as.matrix(matrix[lower.tri(matrix)]))
  nonzeros=tmp1[which(abs(tmp1[,3])>thrsld),1:2] ##--Nodes1 
  Adja=Adj(Nodes=nonzeros,nloci=nloci)
  set.seed(seed1)
  Omega=round(rgwish(n=1,adj.g=Adja,b=delta,D=U),4)
  Omega=matrix(Omega,nrow=nloci,byrow=TRUE)
  tmp2=cbind(Nodes,as.matrix(Omega[lower.tri(Omega)]))
  nonzeros2=tmp2[which(abs(tmp2[,3])>0),1:2] ##--Nodes2 
  sparsity.pct=100*(nrow(nonzeros2)/(nloci*(nloci-1)/2))
  return(list(Omega,nonzeros2,sparsity.pct))            
}


nloci=300
thrsld=35
U=diag(nloci)

##--In the following line

##--For scenario 1 use: 
##--seed1=40, seed2=250 for replicate 1
##--seed1=20, seed2=250 for replicate 2
##--seed1=10, seed2=250 for replicate 3
##--seed1=50, seed2=250 for replicate 4

##--For scenario 2 use: 
##--seed1=500, seed2=2500 for replicate 1
##--seed1=501, seed2=2500 for replicate 2
##--seed1=502, seed2=2500 for replicate 3
##--seed1=600, seed2=2500 for replicate 4

Om1=Rand.Sparse(nloci=nloci,thrsld=thrsld,delta=10,U=U,seed1=40,seed2=250)

#############################################################
####---Comparisson with true concentration matrix------######
#############################################################

NORM.DIFF=norm(Test[[4]]-Om1[[1]],type="F")/norm(Om1[[1]],type="F")

##################################################
############----Save estimated Omega-----#########
##################################################

write.table(Test[[4]],row.names=FALSE,col.names=FALSE,file="OmegaHatBayes1.csv",sep=",")

############################################
####---Get true adjacency matrices----######
############################################

Adj1=Adj(Nodes=Om1[[2]],nloci)
Adj1=as.matrix(Adj1[upper.tri(Adj1)])


############################################
##---Get estimated adjacency matrices---####
############################################

dj.hat=Adj(selgraph1.csv,nloci)
Adj.hat=as.matrix(Adj.hat[upper.tri(Adj.hat)])

###----Differences---#####

Dif=as.matrix(Adj1-Adj.hat)

####----Computing measures of performance----#######
P.template=which(Adj1==1)
P=length(P.template)
N.template=which(Adj1==0)
N=length(N.template)
TN=length(which(Adj.hat[N.template]==0))
N=nrow(Adj1)-P
TP=length(which(Adj.hat[P.template]==1))
FP=N-TN
FN=P-TP
Sensitivity=TP/P
Specificity=TN/N
FPR=FP/N
RCCCP=length(which(Dif==0))/length(Dif)

