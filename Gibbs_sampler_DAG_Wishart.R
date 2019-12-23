
#setwd("")

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
###-----Create training and testing population-----######
#########################################################

data1.Training=data1[-which(data1$G==3), ]
W.Train=W[-which(data1$G==3), ]
n.train=nrow(W.Train)

data1.Test=data1[which(data1$G==3), ]

###############################################################
########---Function to build the adjacency matrices-----#######
###############################################################

Adj=function(Nodes,nloci){
  library(Matrix)
  Adj=as(matrix(0,nrow=nloci,ncol=nloci),"dgTMatrix")
  for(i in 1:nrow(Nodes)){
    Adj[Nodes[i,1],Nodes[i,2]]=1
  }
  return(Adj)
}

####---Function to define sets of parent nodes and create file "Nodes"----######

##--Parents have smaller indices

##sel.G = selected graph, replace selDAG1.csv by selDAGj.csv, j=2,3,4 when running another replicate.

Sel.G=read.csv("selDAG1.csv",header=FALSE)

pai=function(i){matrix(which(sel.G[i,1:(i-1)]!=0))}

parents=apply(matrix(seq(2,nrow(sel.G),1)),1,pai)

getnodes=function(parents){
 seed=matrix(0,ncol=2)
  for(i in 1:length(parents)){
   npa=length(parents[[i]])
   if(npa==0)
   next
   b=cbind(rep(i+1,npa),parents[[i]])
   seed=rbind(seed,b)
  }
 return(seed[-1,])
}  


Nodes=getnodes(parents)

pa=matrix(0,nrow=nloci)
pa[1]=0
for(i in 2:nloci){pa[i]=length(parents[[i-1]])}
alpha=as.matrix((rep((max(pa)+3),nloci)))


############################################################################
###---Function to sample from full cond. under a DAG-Wishart prior-------### 
############################################################################

SampLD.DAG=function(alpha,U,nloci,Nodes,g){
 alpha=alpha+1
 S=tcrossprod(g)
 postm=S+U
 D=matrix(0,nrow=nloci)
 X=list(c(rep(NA,nloci-1)))
 L=diag(nloci)
 library(pscl)
 library(mvtnorm) 
 library(Matrix)
 for(j in 2:(nloci)){
  pointer=Nodes[which(Nodes[ ,1]==j),2]
  if(length(pointer)>0){
    Post1=postm[pointer,pointer]
    cov=solve(Post1) ##slower but less problematic
    postmj=postm[j,pointer]
    D[j]=rgamma(1,shape=(alpha[j]+length(pointer)+2)/2,rate=(postm[j,j]-postmj%*%(cov%*%postmj))/2)
    X[[j]]=rmvnorm(1,mean=-cov%*%postmj,sigma=cov/D[j],method="chol")
    L[j,pointer]=X[[j]]
    }
 } 
      D[1]=rgamma(1,shape=(alpha+2)/2,rate=U[1,1]/2)
      OMEGA=tcrossprod(L%*%Diagonal(length(D),sqrt(D)))
      
      return(list(OMEGA,matrix(do.call(c,X)),D))
	                                        }



################################################################################
################################################################################
###########----------------------Gibbs sampler-----------------#################
################################################################################              
################################################################################


BibiDAG.W=function(y,nloci,Nsim,alpha,U,Tau,V,W,Nodes,Burnin){
          
   n=nrow(W)
   X.size=nrow(Nodes)

   #####---Create arrays to save samples

   g=matrix(0,nrow=Nsim,ncol=nloci)
   norms=matrix(0,nrow=Nsim)
   gnorms=matrix(0,nrow=Nsim)
   res.var=matrix(0,nrow=Nsim)

      #####---Initial values of parameters

   library(mvtnorm) 
   g[1, ]=rmvnorm(1,mean= matrix(0,nrow=nloci),sigma=U,method="chol")

   library(pscl)
   res.var[1]=rigamma(1,alpha=V/2,beta=Tau/2)
   
  
   ##########################################
   #####-----Start algorithm-----############
   ##########################################

   WTW=crossprod(W)
   WTy=crossprod(W,y)
   Ginv.hat=matrix(0,nrow=nloci,ncol=nloci)
   for(i in 2:Nsim){

     Ginv=SampLD.DAG(alpha=alpha,U=U,nloci=nloci,Nodes=Nodes,g=g[i-1, ])
     norms[i]=norm(Ginv[[1]],type="F")
     cov=chol2inv(chol(Ginv[[1]] + WTW/res.var[i-1]))
     g[i, ]=rmvnorm(1, mean=crossprod(cov,WTy)/res.var[i-1],sigma=as.matrix(cov),method="chol")
     gnorms[i]=norm(matrix(g[i, ]))
     res.var[i]=rigamma(1,alpha=(V+n)/2,beta=(Tau+crossprod(y-W%*%g[i, ]))/2)
     if(i > Burnin){
                 Ginv.hat=(Ginv.hat*(i-Burnin-1)+Ginv[[1]])/(i-Burnin)
                                     }
                          }
   
   PostmeanMar=colMeans(g[(Burnin+1):Nsim, ])
   
   ###--In the following lines, replace "1" by "j", j=2,3,4 when running a different replicate
   
   write.table(matrix(res.var[(Burnin+1):Nsim]),col.names=FALSE,row.names=FALSE,
   file="Phen1ResVar1.csv", sep=",")

   write.table(gnorms[(Burnin+1):Nsim],col.names=FALSE,row.names=FALSE,
   file="Phen1gnorms1.csv", sep=",")
   
   write.table(norms[(Burnin+1):Nsim],col.names=FALSE,row.names=FALSE,
   file="Phen1norms1.csv", sep=",")
   
   write.table(PostmeanMar,col.names=FALSE,row.names=FALSE,
   file="Phen1PostmeanMar1.csv", sep=",")

   write.table(as.matrix(Ginv.hat),col.names=FALSE,row.names=FALSE,
   file="OmegaHat1.csv", sep=",")

   return(list(PostmeanMar,res.var[Burnin:Nsim],norms[(Burnin+1):Nsim],Ginv.hat))

                   
                                                          }

U=diag(nloci)
Test=BibiDAG.W(y=data1.Training$Phen,nloci=nloci,Nsim=50000,alpha=alpha,U=U,Tau=7000,V=72,W=W.Train,Nodes=Nodes,Burnin=20000)

###---Build the estimated adjacency matrix---######

Adj.hat=Test[[4]]

for(i in 1:nrow(Adj.hat)){
 for(j in 1:ncol(Adj.hat)){
  Adj.hat[i,j]=as.numeric(Adj.hat[i,j]!=0)
 }
}


Nod=Def.Graph.Band(band.size=nloci,nloci=nloci)[[1]]
tmp1=cbind(Nod,as.matrix(Adj.hat[lower.tri(Adj.hat)]))
nonzeros=tmp1[which(abs(tmp1[,3])>0),1:2] 
sparsity.pct=100*(nrow(nonzeros)/(nloci*(nloci-1)/2))

Adj.hat=Adj(nonzeros,nloci)
Adj.hat=as.matrix(Adj.hat[upper.tri(Adj.hat)])


################################################################
################################################################
######---Evaluating the performance of graph selection-----#####
################################################################
################################################################


####-----Randomly setting  zeros------####

Rand.Sparse=function(nloci,thrsld,delta,U,seed1,seed2){
      Nodes=Def.Graph.Band(band.size=nloci,nloci=nloci)[[1]]
       set.seed(seed2)
	  matrix=rwish(n=1,p=nloci,b=nloci+1,D=diag(nloci))
	   matrix=matrix(matrix,nrow=nloci,byrow=TRUE)
	    tmp1=cbind(Nodes,as.matrix(matrix[lower.tri(matrix)]))
	     nonzeros=tmp1[which(abs(tmp1[,3])>thrsld),1:2]
	      Adja=Adj(Nodes=nonzeros,nloci=nloci)
	      set.seed(seed1)
	     Omega=round(rgwish(n=1,adj.g=Adja,b=delta,D=U),4)
          Omega=matrix(Omega,nrow=nloci,byrow=TRUE)
         tmp2=cbind(Nodes,as.matrix(Omega[lower.tri(Omega)]))
        nonzeros2=tmp2[which(abs(tmp2[,3])>0),1:2] 
       sparsity.pct=100*(nrow(nonzeros2)/(nloci*(nloci-1)/2))
      return(list(Omega,nonzeros2,sparsity.pct))            
}

nloci=300
thrsld=35
#thrsld=50 #for scenario 2
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

############################################
####---Get true adjacency matrices----######
############################################

Adj1=Adj(Nodes=Om1[[2]],nloci)
Adj1=as.matrix(Adj1[upper.tri(Adj1)])

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


#####################################################################
############----Compute accuracies and predictive ability-----#######
#####################################################################


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

#Replace BayesDAGWperform1.csv by BayesDAGWperformj.csv, j=2,3,4 when running another replicate
write.table(Performance,col.names=FALSE,row.names=FALSE,file="BayesDAGWperform1.csv", sep=",")

#############################################################
####---Comparisson with true concentration matrix------######
#############################################################

NORM.DIFF=norm(Test[[4]]-Om1[[1]],type="F")/norm(Om1[[1]],type="F")




