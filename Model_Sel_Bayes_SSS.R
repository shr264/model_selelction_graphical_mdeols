
####----This code works for replicate 1 under scenario I, lines that need to be 
####----modified to run another replicate/scenario are commented accordingly

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

####---Reading graph files----#####

#--As explained in the paper, for each replicate, the SSS algortihm is run for 
#--for each one of the 15 graphs.

#X=1,2,3,4    
#Y=1,2,...,15

GraphX_Y=read.table("GraphX Y .txt",header=F)  

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

#########----Function to do full MC integration (i.e.,NO Laplace approximation)-------#########

MCint2=function(y,W,Tau=210,V=4.1,nsamples=2000,Adjac){
 nloci=ncol(W)
  library(mvtnorm)
   library(pscl)
    n=length(y)
     I=diag(n)
      WTW=crossprod(W)
	 yt=t(y)
	  Wt=t(W)
	   funct=matrix(0,nrow=nsamples,ncol=1)
	    for(i in 1:nsamples){
         resvar=rigamma(1,alpha=V/2,beta=Tau/2)
        Omega=matrix(round(rgwish(n=1,adj.g=Adjac,b=10,D=diag(nloci)),8),ncol=nloci,byrow=TRUE)
       Sigma=as.matrix(forceSymmetric(chol2inv(chol(Omega))))
      g=matrix(rmvnorm(1,mean=matrix(0,nrow=nloci,ncol=1),sigma=Sigma,method="svd"))
      funct[i]=dmvnorm(yt,mean=W%*%g,sigma=diag(resvar,nrow=n),log=TRUE)/100 
    }
   estimate=mean(funct)
  var.estimate=var(funct)*(nsamples-1)/(nsamples^2)
 return(list(estimate,var.estimate,funct))        
}


#########-------Function to compute the score when n > m-------#########

norm.const <- function(y,W,Adjac,delta){
  library(BDgraph)
  nloci=ncol(W)
  WTW=crossprod(W)
  ghat=chol2inv(chol(WTW))%*%(t(W)%*%y)
  S=tcrossprod(ghat)
  score=gnorm(adj.g=Adjac,b=delta+1,D=diag(nloci)+S,iter=niter)
}

#####################################################################################################################
#####################################################################################################################
######--------Functions to perform the SSS algorithm of Ben-David et al. (2015)-----------###########################
#####################################################################################################################
#####################################################################################################################

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
        scores[i] <-MCint2(y=data1.Training$Phen,W=W.Train,Tau=210,V=4.1,nsamples=2000,Adjac=graphs[[i]])[[1]]
        #scores[i] <- norm.const(y=data1.Training$Phen,W=W.Train,Adjac=graphs[[i]],delta=10) #-For m<n case only, scenario 2
        }
    scores
}

getDnew <- function(gamma, scores2, graphs2){
    pvec <- exp(scores2*gamma) 
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

###---Get the score for a given Y, Y=1,2,...,15 
 
Adjac=Adj(Nodes=GraphX_Y,nloci=nloci)

#X=1,2,3,4 
#Y=1,2,...,15

ScoreX_Y=getLk(gamma=0.5,M=3,N1=10,D0=Adjac)

####--So,this gives a list of 30 graphs and scores for a given Y, Y=1,2,...,15. 

###---Pick the graph with highest score for this particular Y----####

Scores=matrix(ScoreX_Y$scores)
Graph.number=which.max(Scores)
Max.score_Y=Scores[Graph.number]
Max.score_Y 

#--Thus, for each replicate, after running this program for each of the 15 graphs from
#--CONCORD_EMyou will have Max.score_1,...,Max.score_15, pick the largest one and run the 
#--the follwing command in the corresponding file to get the selected graph.
#--This graph will be used as input file in program Gibbs_Sampler_General_Graph

sel.graph1=ScoreX_Y$graphs[[Graph.number]]

#Replace selgraph1.csv by selgraphj.csv, j=2,3,4 when running a different replicate

write.table(sel.graph,col.names=FALSE,row.names=FALSE,file="selgraph1.csv", sep=",")
