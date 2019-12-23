
####----This code works for replicate 1 under scenario II, lines that need to be 
####----modified to run another replicate/scenario are commented accordingly
####----Recall that Bayes DAG-Sel only works when n > m, which means that this code 
####----is only used in scenario 1. 

#setwd("")

#####################################################################################################
#####################################################################################################
########------------------------Program to read data from QMSim-----------------------------#########
#####################################################################################################
#####################################################################################################

#--To read genotypes

#replace "Pop1_qtl_001.txt" by "Pop1_qtl_00j.txt", j=2,3,4 to run another replicate 

QTL1=read.table("Pop1_qtl_001.txt",header=F) 

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

###---Reading phenotypes and total additive effects     

#replace "Pop1_data_001.txt" by "Pop1_data_00j.txt", j=2,3,4 to run another replicate

data1=read.table("Pop1_data_001.txt",header=T)                                     

W=QTL.Matrix(QTL=QTL1)
nloci=ncol(W)

#########################################################
###-----Create training and testing population-----######
#########################################################

data1.Training=data1[-which(data1$G==4), ]
W.Train=W[-which(data1$G==4), ]
n.train=nrow(W.Train)
data1.Test=data1[which(data1$G==4),]

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

###############----Read multiple files with the DAG's from CSCS_EM-------###############

temp=list.files(pattern="*.txt")
Graphs=lapply(temp, read.table)
dim(Graphs[[1]])

ad=function(k){
 for(i in 1:nloci){
  for(j in 1:nloci){
    Graphs[[k]][i,j]=as.numeric(Graphs[[k]][i,j]!=0)
   }
  }
 return(Graphs[[k]])
}

Graphs=apply(matrix(c(1:length(Graphs))),1,ad)

################################################################################
################################################################################
##########-----Function to get posterior score based on DAG-Wishart n>m-----####
################################################################################
################################################################################


WTWInv=solve(crossprod(W.Train)+diag(0.01,nloci))
Ay=tcrossprod(WTWInv,W.Train)%*%data1.Training$Phen
U2=tcrossprod(Ay)+U

####---Function to define sets of parent nodes----######

##--Parents have smaller indices

pai=function(i){
 matrix(which(L[i,1:(i-1)]!=0))
}

getnodes=function(parents){
 seed=matrix(0,ncol=2)
  for(i in 1:length(parents)){
   if(length(parents[[i]])>0){
      b=cbind(rep(i+1,nrow(parents[[i]])),parents[[i]])
      e=rbind(seed,b)
     }
   seed=e
  } 
 return(e[-1,])
} 
 

###---Get the score----#####

zs=function(i,U2,pa,parents,alpha){
alf=(alpha[i]-pa[i]-2)*0.5
  num=gamma(alf)*(2^(alf+(pa[i]*0.5)))*(sqrt(pi)^pa[i])*(det(as.matrix(U2[parents[[i-1]],parents[[i-1]]]))^(alf-0.5))
  den=(det(as.matrix(U2[c(parents[[i-1]],i),c(parents[[i-1]],i)])))^alf
 z=num/den
return(z)
}


get.score2=function(alpha,parents,pa){
 alpha=alpha+1
  z1=(gamma((alpha[1]-2)*0.5)*(2^(alpha[1]*0.5-1))*det(U2))/(U2[1,1]^((alpha[1]-2)*0.5))
    Scores=0
   for(i in 2:nloci){Scores=Scores+log(zs(i,U2,pa,parents,alpha))}
  score=log(z1)+Scores
 return(score)
}


DAG.score2=function(L){
 parents=list(NA)
  for(i in 2:nrow(L)){parents[[i-1]]=matrix(which(L[i,1:(i-1)]!=0))}
   pa=matrix(0,nrow=nloci)
    pa[1]=0
    for(i in 2:nloci){pa[i]=length(parents[[i-1]])}
   alpha=as.matrix((rep((max(pa)+3),nloci)))
  score=get.score2(alpha,parents,pa)
 score
}



###############################################################################################################
###############################################################################################################
######---------Functions to perform the SSS algorithm of Ben-David et al. (2015)----###########################
###############################################################################################################
###############################################################################################################


getNgraphs <- function(D,N){
    N1 <- list()
    dimD <- dim(D)
    for(i in 1:N){
        N1[[i]] <- D
        x <- floor(runif(1,1,dimD[1]))
        y <- floor(runif(1,1,x))
        N1[[i]][x,y] = ifelse(N1[[i]][x,y]==0,1,0)
    }
    N1
}

 
 
getscores <- function(graphs){
    scores <- rep(0,length(graphs))
    for(i in 1:length(graphs)){
        scores[i]<- DAG.score2(L=graphs[[i]])
        }
    scores
}

 
getDnew <- function(gamma, scores2, graphs2){
    pvec <- exp(scores2^gamma)
    pvec <- pvec/sum(pvec)
    cpvec <- cumsum(pvec)
    u <- runif(1,0,1)
    Dind <- min(sum(u > cpvec) + 1,length(scores2))
    graphs2[[Dind]]
}

 
getLk <- function(gamma,M,N1,D0){
 
    graphs <- getNgraphs(D0,N1)
    scores <- getscores(graphs)
    D0 <- getDnew(gamma,scores,graphs)
    for(i in 2:M){
        graphs2 <- getNgraphs(D0,N1)
        scores2 <- getscores(graphs2)
        D0 <- getDnew(gamma,scores2,graphs2)
        graphs <- append(graphs,graphs2)
        scores <- append(scores,scores2)
    }
    return(list(graphs=graphs, scores=scores))
}

 
 
Largestscoregraph <- function(gamma,D,M,N1){
    L <- list()
    for(i in 1:length(D)){
        newL <- getLk(gamma,M,N1,D[[i]])
        L$graphs <- append(newL$graphs,L$graphs)
        L$scores <- append(newL$scores,L$scores)
    print(i)
    }
    L$graphs[[which.max(L$scores)]]
}


Sel.G=Largestscoregraph(gamma=0.5,D=Graphs,M=3,N1=10)

#Replace selDAG1.csv by selDAGj.csv when running a different replicate, j=2,3,4
write.table(sel.G,col.names=FALSE,row.names=FALSE,file="selDAG1.csv", sep=",")

##--This file (selDAG1.csv) is an input in program Gibbs_Sampler_General_Graph
