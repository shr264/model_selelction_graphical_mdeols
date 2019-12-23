
####----This code works for replicate 1 under scenario I, lines that need to be 
####----modified to run another replicate/scenario are commented accordingly

#setwd("")

#####################################################################################################
#####################################################################################################
########-------------------------------To read data ----------------------------------------#########
#####################################################################################################
#####################################################################################################

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


######################################################################
####-----------Auxiliary Function to be used later---------------#####
######################################################################

Def.Graph.Band=function(band.size,nloci){
      win.size=band.size-1
	size=band.size*(nloci-band.size+1)+sum(seq(2,max(band.size-1,2),1))-nloci+1  
	Connect.nodes=matrix(nrow=size,ncol=2)
      template=sort(rep(seq(1,nloci-win.size,1),win.size))
      for(i in 0:(nloci-win.size-1)){
      	Connect.nodes[which(template==i+1),2]=seq(2+i,2+i+win.size-1,1)
      	seq(2+i,2+i+win.size-1,1)
                                     }
      Connect.nodes[1:length(template),1]=template
      counter=seq(2,max(2,win.size-1),1)
      links=list(rep(NA,win.size-1))
      for(i in 1:(win.size-1)){
		links[[i]]=rep(nloci-i,i)
                             }
	Connect.nodes[(length(template)+1):size,1]=sort(do.call(c,links),decreasing=FALSE)
	links2=list(rep(NA,win.size-1))
      for(i in 1:(win.size-1)){
		links2[[i]]=seq(nloci-win.size+i+1,nloci,1)
                             }
	Connect.nodes[(length(template)+1):size,2]=do.call(c,links2)	
      return(list(Connect.nodes,size))	
                                        }

#######################################################################################################################
#######################################################################################################################
##############---------------Function to run Concord using S as an argument---------###################################
#######################################################################################################################
#######################################################################################################################

custConcord1 = function(S,n,rmax=100,eps=10^(-5),lambda){
(p=dim(S)[1])  
(omegahatcurrent = diag(p))
(r = 1)
(converged  = FALSE)
(maxdiff = eps/10)
while((converged == FALSE)&&(r<rmax)){
    #cat('Iter:',r,maxdiff,'\n')
    (maxdiff = eps/10)
    (omegahatold = omegahatcurrent)
    for (i in 1:(p-1)){
    for (j in (i+1):p){
    if(i!=j){
        #cat('r=',r,'i=',i,'j=',j,'i!=j\n')
        (x = -t(omegahatcurrent[i,!(1:p == j)])%*%(S[j,!(1:p == j)])
         - t(omegahatcurrent[j,!(1:p == i)])%*%(S[i,!(1:p == i)]))
        (omegahatcurrent[i,j] = sign(x)*max(abs(x)-(lambda/n),0))
        (omegahatcurrent[i,j] = omegahatcurrent[i,j]/(S[i,i]+S[j,j]))
        (omegahatcurrent[j,i] = omegahatcurrent[i,j])
        (maxdiff = max(maxdiff, abs(omegahatcurrent[i,j]-omegahatold[i,j])))
    }}}
for(i in 1:p)
{
#cat('r=',r,'i=',i,'i==j\n')
    (omegahatcurrent[i,i] = -t(omegahatcurrent[i,!(1:p == i)])%*%(S[i,!(1:p == i)]) + sqrt((t(omegahatcurrent[i,!(1:p == i)])%*%(S[i,!(1:p == i)]))^2+4*S[i,i]))
    (omegahatcurrent[i,i] = omegahatcurrent[i,i]/(2*S[i,i]))
    (maxdiff = max(maxdiff, abs(omegahatcurrent[i,i]-omegahatold[i,i])))
}
if(maxdiff<eps){
converged = TRUE
}
else{
    r = r+1
}
}
return(list(Omega = omegahatcurrent,Iter = r))
}


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

#---This second version builds a symmetric matrix instead of an upper triangular one

Adj2=function(Nodes,nloci){
	library(Matrix)
	Adj=as(diag(nloci),"dgTMatrix")
      for(i in 1:nrow(Nodes)){
		Adj[Nodes[i,1],Nodes[i,2]]=1
            Adj[Nodes[i,2],Nodes[i,1]]=1
                              }
	return(Adj)
			       }


#######################################
###---Grouping data by families----####
#######################################

Sires.Order=matrix(sort(data1.Training$Sire))
Number.Sires=1 ##Stars in 1 to count 0, i.e., unknown sire
Temp=matrix(-1,nrow=nrow(Sires.Order))
for(i in 1:(nrow(Sires.Order)-1)){
                         if(Sires.Order[i]!=Sires.Order[i+1]){
                                            Number.Sires=Number.Sires+1
                                            Temp[i]=Sires.Order[i+1]
                                                              }
                                  }

f=Number.Sires

Sires.List=matrix(c(0,Temp[which(Temp!=-1)])) ##Does include 0 which means unknown sire, these are grouped as founders
length(Sires.List)

templ=list(rep(NA,f))
for(i in 1:f){
            templ[[i]]=matrix(which(data1.Training$Sire==Sires.List[i]))
              }

sizes=matrix(0,nrow=f)
for(i in 1:f){sizes[i]=length(templ[[i]])}

column=list(rep(NA,f))
for(i in 1:f){column[[i]]=matrix(rep(i,sizes[i]))}
column=do.call(rbind,column)

#-Col1=Family member position, Col2=Original Family code, Col3=Consecutive family code
template=cbind(do.call(rbind,templ),Sires.Order,column)

#################################
#####----Define Omega0-----######
#################################


####-----Randomly setting  zeros------####

library(BDgraph)

Rand.Sparse=function(nloci,thrsld,delta,U,seed){
      Nodes=Def.Graph.Band(band.size=nloci,nloci=nloci)[[1]]
       set.seed(250)
	  matrix=rwish(n=1,p=nloci,b=nloci+1,D=diag(nloci))
	   matrix=matrix(matrix,nrow=nloci,byrow=TRUE)
	    tmp1=cbind(Nodes,as.matrix(matrix[lower.tri(matrix)]))
	     nonzeros=tmp1[which(abs(tmp1[,3])>thrsld),1:2] 
	     sparsity.pct=100*(nrow(nonzeros)/(nloci*(nloci-1)/2))
	     Adja=Adj(Nodes=nonzeros,nloci=nloci)
	   set.seed(seed)
	  Omega=round(rgwish(n=1,adj.g=Adja,b=delta,D=U),8)
	 Omega=matrix(Omega,nrow=nloci,byrow=TRUE)
      return(list(Omega,nonzeros,sparsity.pct))            
                                              }

thrsld=35
U=diag(nloci)

SparseSample=Rand.Sparse(nloci=nloci,thrsld=thrsld,delta=10,U=U,seed=40)
Omega0=SparseSample[[1]]
Init.sparse=SparseSample[[3]]


Bibi.GEMCONCORD=function(y,W,Omega0,res.var0,epsilon=1e-04,Maxiter=100,template,sizes,rho=0.001,thr=0.001){

	library(glasso)
      library(Matrix)
	N=length(y)
      f=max(template[ ,3])
	nloci=ncol(W)
      npos=nloci*(nloci+1)*0.5
	iter=0
      SIGMA0=chol2inv(chol(as.matrix(Omega0)))
      temp0=(W%*%SIGMA0)%*%t(W)+(res.var0*diag(N))
	Omega=as.matrix(Omega0)
	res.var=res.var0
      delta=epsilon*2
      Im=diag(nloci)
      
    while(epsilon<delta && iter< Maxiter){
      Oldres.var = res.var
      OldOmega=Omega
	Sg.sum=matrix(0,nrow=nloci,ncol=nloci)
	trace=0
      inner=0
	SIGMA=chol2inv(chol(Omega))

###---Expectation step, Using Formula 1

	for(i in 1:f){
		indices = template[which(template[,3]==i),1]
            WindSig = W[indices, ]%*%SIGMA
		Vinv=chol2inv(chol((WindSig)%*%t(W[indices, ])+(diag(sizes[i])*res.var)))
		summand=SIGMA%*%(Im - (t(W[indices, ])%*%(Vinv-crossprod(t(y[indices])%*%Vinv)))%*%WindSig)
		Sg.sum=Sg.sum+summand
		       }

	Sg.expect=Sg.sum/f 
	Vin=chol2inv(chol(W%*%(SIGMA%*%t(W))+diag(N)*res.var))
	inner.res.expect=as.numeric(res.var*(N-res.var*(sum(diag(Vin))-crossprod(Vin%*%y))))

##---Maximization step

	res.var=inner.res.expect/N
	Omega=custConcord1(S=Sg.expect,n=f,rmax=50,eps=thr,lambda=rho)[[1]]+diag(0.1,nloci)
       
##---Evaluate convergence

	#delta=max(abs(max(Omega-OldOmega)),abs(Oldres.var-res.var))
	delta=(npos*mean(abs(Omega-OldOmega))+abs(Oldres.var-res.var))/(npos+1)
     
	iter=iter+1
                                   
                                     } #end of the while loop

      print(iter)
	if(delta>epsilon){print("Convergence criteria was not met")}
      Omega=as.matrix(Omega)
      Nodes=Def.Graph.Band(band.size=nloci,nloci=nloci)[[1]]
      tmp1=cbind(Nodes,as.matrix(Omega[lower.tri(Omega)]))
	nonzeros=tmp1[which(abs(tmp1[,3])>0),1:2] ##--i.e.Nodes
	sparsity.pct=100*(nrow(nonzeros)/(nloci*(nloci-1)/2))
     
	
	return(list(Omega,res.var,iter,delta,sparsity.pct,nonzeros))
	
	}


###################################################################################################
###################################################################################################
####----------------------------SOLVING MME TO GET g_hat---------------------------------##########
###################################################################################################
###################################################################################################

Solve.MME=function(y,W,OMEGA,Resvar){
	RHS=crossprod(W,y)
	LHS=crossprod(W)+Resvar*OMEGA
	ghat=solve(LHS)%*%RHS
	return(ghat)
                                    }

############################################################################################################
############################################################################################################
#####------Program to write GLasso-EM or CONCORD-EM over a path of the tuning parameter----#################
############################################################################################################
############################################################################################################

#lambda.grid=matrix(seq(3.3,6.2,length.out=15)) ##--Use this one for scenario 2
lambda.grid=matrix(seq(1.1,2.4,length.out=15))
library(mvtnorm)
n.Train=nrow(W.Train)
Wt=t(W.Train)
W.Test=W[which(data1$G==4), ]

F=function(lambda){
 Run=Bibi.GEMCONCORD(y=data1.Training$Phen,W=W.Train,Omega0=Omega0,res.var0=105.24,epsilon=0.001,Maxiter=50,template=template,
  sizes=sizes,rho=lambda,thr=0.001)
   OMEGA=Run[[1]]
    SIGMA=solve(OMEGA)  
     ghat=Solve.MME(y=data1.Training$Phen,W=W.Train,OMEGA=OMEGA,Resvar=Run[[2]])
      PredBVTrainGEMB=W.Train%*%ghat
       PredBVTestGEMB=W.Test%*%ghat
        Predabil.GB=cor(data1.Test$Phen,as.numeric(PredBVTestGEMB),method="pearson")
         corrBVGB=cor(data1.Test$QTL,as.numeric(PredBVTestGEMB),method="pearson")
         corrBVTrainGB=cor(data1.Training$QTL,as.numeric(PredBVTrainGEMB),method="pearson")
        V=as.matrix(forceSymmetric(W.Train%*%(SIGMA%*%Wt)+diag(Run[[2]],n.Train)))
       tuning.criteriaLIK=dmvnorm(data1.Training$Phen, mean = rep(0, n.Train), sigma = V, log = TRUE)
      tuning.criteriaBIC=-2*tuning.criteriaLIK+log(n.Train)*nrow(Run[[6]])
     tuning.criteriaRES=sum((data1.Test$Phen-as.numeric(PredBVTestGEMB))^2)
   PAR1=matrix(c(tuning.criteriaBIC,tuning.criteriaLIK,tuning.criteriaRES,Predabil.GB,
  corrBVGB,corrBVTrainGB,as.numeric(Run[[5]])))
 return(list(PAR1,Run[[6]]))
}

Out=apply(lambda.grid,1,F)

summary.tune=matrix(nrow=length(lambda.grid),ncol=8)

summary.tune[,1]=lambda.grid
for(i in 1:length(lambda.grid)){summary.tune[i,2:8]=Out[[i]][[1]]}

selected.graphLIK=Out[[which.max(summary.tune[,3])]][[2]]
selected.graphSPARSITY=Out[[6]][[2]]
selected.graphBIC=Out[[which.min(summary.tune[,2])]][[2]]
selected.graphRES=Out[[which.min(summary.tune[,4])]][[2]]
selected.graphPA=Out[[which.max(summary.tune[,5])]][[2]]

#################################################################################################
#############-------------------------Saving outputs---------------------########################
#################################################################################################

## In the following lines, replace the last digit of the output file name by j,
#j=2,3,4, in order to avoid confusion. 

write.table(summary.tune,row.names=FALSE,col.names=c("Lambda","BIC","LIK","RES","PA","ACCV","ACCT","SPARSITY"),
file="Tuning1CONCORD.csv",sep=",")

write.table(selected.graphBIC,row.names=FALSE,col.names=FALSE,file="EdgesBICCONCORD1.csv",sep=",")

write.table(selected.graphLIK,row.names=FALSE,col.names=FALSE,file="EdgesLIKCONCORD1.csv",sep=",")

write.table(selected.graphRES,row.names=FALSE,col.names=FALSE,file="EdgesRESCONCORD1.csv",sep=",")

write.table(selected.graphSPARSITY,row.names=FALSE,col.names=FALSE,file="EdgesSPACONCORD1.csv",sep=",")

write.table(selected.graphPA,row.names=FALSE,col.names=FALSE,file="EdgesPACONCORD1.csv",sep=",")

##############################################################################
############-----Save the 15 Graphs---------##################################
##############################################################################

#In the following line, replace "GRAPH1" by "GRAPHj", j=2,3,4 when analyzing another replicate.
#Moreover, for scenario 2, add "II" at the end of filename such that it has the form "GraphjII"
#j=1,2,3,4.

for(i in 1:length(Out)){
  write.table(as.matrix(Out[[i]][[3]]),row.names=FALSE,col.names=FALSE,file=paste("Graph1",toString(i),".txt",sep=" "))
}


################################################################
################################################################
######---Evaluating the performance of graph selection-----#####
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


############################################
##---Get estimated adjacency matrices---####
############################################

Adj.hat=Adj(selected.graphLIK,nloci)
Adj.hat=as.matrix(Adj.hat[upper.tri(Adj.hat)])

Adj.hat=Adj(selected.graphSPARSITY,nloci)
Adj.hat=as.matrix(Adj.hat[upper.tri(Adj.hat)])

Adj.hat=Adj(selected.graphPA,nloci)
Adj.hat=as.matrix(Adj.hat[upper.tri(Adj.hat)])

Adj.hat=Adj(selected.graphBIC,nloci)
Adj.hat=as.matrix(Adj.hat[upper.tri(Adj.hat)])

Adj.hat=Adj(selected.graphRES,nloci)
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

####---The following part uses the selected graph under each criterion
####---which correspond to objects selected.graphLIK, selected.graphSPARSITY, 
####---selected.graphBIC, selected.graphRES and selected.graphPA defined above.
####---It uses each graph as an input to estimate Omega. This step is required
####---because CONCORD does not guarantee positive definite estimators
####---for further details see the manuscript

#####################################################################################################
##########---------Function to implement the IPF algorithm of Speed and Kiiveri----------############
#####################################################################################################


ipf <- function(S,amat,tol= 10^(-4)){
    cat('...IPF')
    p = dim(S)[2]
    k = ncol(S)
    W0 = W = theta = S
    it = 0
    converge = FALSE
    while( !converge ) {
        it = it+1
        for (j in 1:k){
            W11 = W[-j,-j, drop=FALSE]     
            w12 = W[-j,j]     
            s12 = S[-j,j, drop=FALSE]
            s22 = S[j,j]
            paj = amat[j,] == 1; # neighbors
            paj = paj[-j]
            beta = rep(0, k-1)
            if (all(!paj)){
                w = rep(0, k-1)  
            }
            else{
                beta[paj] <- solve(W11[paj, paj], s12[paj, ])
                w = W11 %*% beta
            }
            W[-j, j] = w
            W[j, -j] = w
            theta[j,j] = 1/(s22-sum(s12*beta))
            theta[-j,j] = -beta*theta[j,j]
        }
        di <- norm(W0-W)      
        if (di < tol){
            converge = TRUE
        }
        else {
            W0 <- W 
          }
    }
    return(list(Omegahat = theta, it=it))
}


####---Function implementing an algortihm combining IPF and EM to estimate Omega

Bibi.IPF=function(y,W,Omega0,res.var0,epsilon=1e-04,Maxiter=100,template,sizes,thr=0.001,Nodes){

	library(Matrix)
	N=length(y)
      f=max(template[ ,3])
	nloci=ncol(W)
      npos=nloci*(nloci+1)*0.5
	iter=0
      SIGMA0=chol2inv(chol(as.matrix(Omega0)))
      temp0=(W%*%SIGMA0)%*%t(W)+(res.var0*diag(N))
	Omega=as.matrix(Omega0)
	res.var=res.var0
      delta=epsilon*2
      Im=diag(nloci)
      Adja=Adj2(Nodes=Nodes,nloci=nloci)+diag(nloci)
	      
    while(epsilon<delta && iter< Maxiter){
      Oldres.var = res.var
      OldOmega=Omega
	Sg.sum=matrix(0,nrow=nloci,ncol=nloci)
	trace=0
      inner=0
	SIGMA=chol2inv(chol(Omega))

###---Expectation step, Using Formula 1

	for(i in 1:f){
		indices = template[which(template[,3]==i),1]
            WindSig = W[indices, ]%*%SIGMA
		Vinv=chol2inv(chol((WindSig)%*%t(W[indices, ])+(diag(sizes[i])*res.var)))
		summand=SIGMA%*%(Im - (t(W[indices, ])%*%(Vinv-crossprod(t(y[indices])%*%Vinv)))%*%WindSig)
		Sg.sum=Sg.sum+summand
		       }

	Sg.expect=Sg.sum/f 
	Vin=chol2inv(chol(W%*%(SIGMA%*%t(W))+diag(N)*res.var))
	inner.res.expect=as.numeric(res.var*(N-res.var*(sum(diag(Vin))-crossprod(Vin%*%y))))

##---Maximization step

	res.var=inner.res.expect/N
	
      Omega=ipf(S=Sg.expect,amat=Adja,tol=thr)[[1]] 
##---Evaluate convergence

	delta=(npos*mean(abs(Omega-OldOmega))+abs(Oldres.var-res.var))/(npos+1)

     
	iter=iter+1
                                   
                                     } #end of the while loop

      print(iter)
	if(delta>epsilon){print("Convergence criteria was not met")}
      Omega=as.matrix(Omega)
           
	
   	return(list(Omega,res.var,iter,delta))

  }



SK.LIK=Bibi.IPF(y=data1.Training$Phen,W=W.Train,Omega0=Omega0,res.var0=83.37,epsilon=0.001,Maxiter=150,template=template,
  sizes=sizes,thr=0.001,Nodes=selected.graphLIK)

SK.SPA=Bibi.IPF(y=data1.Training$Phen,W=W.Train,Omega0=Omega0,res.var0=83.37,epsilon=0.001,Maxiter=150,template=template,
  sizes=sizes,thr=0.001,Nodes=selected.graphSPARSITY)

SK.PA=Bibi.IPF(y=data1.Training$Phen,W=W.Train,Omega0=Omega0,res.var0=83.37,epsilon=0.001,Maxiter=150,template=template,
  sizes=sizes,thr=0.001,Nodes=selected.graphPA)

SK.RES=Bibi.IPF(y=data1.Training$Phen,W=W.Train,Omega0=Omega0,res.var0=83.37,epsilon=0.001,Maxiter=150,template=template,
sizes=sizes,thr=0.001,Nodes=selected.graphRES)


##################################################
############----Save estimated Omega-----#########
##################################################

write.table(SK.LIK[[1]],row.names=FALSE,col.names=FALSE,file="OmegaHatLIK1.csv",sep=",")
write.table(SK.SPA[[1]],row.names=FALSE,col.names=FALSE,file="OmegaHatSPA1.csv",sep=",")
write.table(SK.PA[[1]],row.names=FALSE,col.names=FALSE,file="OmegaHatPA1.csv",sep=",")
write.table(SK.RES[[1]],row.names=FALSE,col.names=FALSE,file="OmegaHatRES1.csv",sep=",")

#############################################################
####---Comparisson with true concentration matrix------######
#############################################################

NORM.DIFFLIK=norm(SK.LIK[[1]]-Om1[[1]],type="F")/norm(Om1[[1]],type="F")
NORM.DIFFSPA=norm(SK.SPA[[1]]-Om1[[1]],type="F")/norm(Om1[[1]],type="F")
NORM.DIFFPA=norm(SK.PA[[1]]-Om1[[1]],type="F")/norm(Om1[[1]],type="F")
NORM.DIFFRES=norm(SK.RES[[1]]-Om1[[1]],type="F")/norm(Om1[[1]],type="F")


#######-----Evaluate predictive performance-------########

ghatLIK=Solve.MME(y=data1.Training$Phen,W=W.Train,OMEGA=SK.LIK[[1]],Resvar=SK.LIK[[2]])
ghatSPA=Solve.MME(y=data1.Training$Phen,W=W.Train,OMEGA=SK.SPA[[1]],Resvar=SK.SPA[[2]])
ghatPA=Solve.MME(y=data1.Training$Phen,W=W.Train,OMEGA=SK.PA[[1]],Resvar=SK.PA[[2]])
ghatRES=Solve.MME(y=data1.Training$Phen,W=W.Train,OMEGA=SK.RES[[1]],Resvar=SK.RES[[2]])

W.Test=W[which(data1$G==3), ]

PredBVTrainLIK=W.Train%*%ghatLIK
PredBVTestLIK=W.Test%*%ghatLIK
Predabil.LIK=cor(data1.Test$Phen,as.numeric(PredBVTestLIK),method="pearson")
corrBVLIK=cor(data1.Test$QTL,as.numeric(PredBVTestLIK),method="pearson")
corrBVTrainLIK=cor(data1.Training$QTL,as.numeric(PredBVTrainLIK),method="pearson")
Performance=matrix(c(Predabil.LIK,corrBVLIK,corrBVTrainLIK))
write.table(Performance,row.names=FALSE,col.names=FALSE,file="CONCORDPERFORMLIK.csv",sep=",")


PredBVTrainSPA=W.Train%*%ghatSPA
PredBVTestSPA=W.Test%*%ghatSPA
Predabil.SPA=cor(data1.Test$Phen,as.numeric(PredBVTestSPA),method="pearson")
corrBVSPA=cor(data1.Test$QTL,as.numeric(PredBVTestSPA),method="pearson")
corrBVTrainSPA=cor(data1.Training$QTL,as.numeric(PredBVTrainSPA),method="pearson")
Performance=matrix(c(Predabil.SPA,corrBVSPA,corrBVTrainSPA))
write.table(Performance,row.names=FALSE,col.names=FALSE,file="CONCORDPERFORMSPA.csv",sep=",")

PredBVTrainPA=W.Train%*%ghatPA
PredBVTestPA=W.Test%*%ghatPA
Predabil.PA=cor(data1.Test$Phen,as.numeric(PredBVTestPA),method="pearson")
corrBVPA=cor(data1.Test$QTL,as.numeric(PredBVTestPA),method="pearson")
corrBVTrainPA=cor(data1.Training$QTL,as.numeric(PredBVTrainPA),method="pearson")
Performance=matrix(c(Predabil.PA,corrBVPA,corrBVTrainPA))
write.table(Performance,row.names=FALSE,col.names=FALSE,file="CONCORDPERFORMPA.csv",sep=",")

PredBVTrainRES=W.Train%*%ghatRES
PredBVTestRES=W.Test%*%ghatRES
Predabil.RES=cor(data1.Test$Phen,as.numeric(PredBVTestRES),method="pearson")
corrBVRES=cor(data1.Test$QTL,as.numeric(PredBVTestRES),method="pearson")
corrBVTrainRES=cor(data1.Training$QTL,as.numeric(PredBVTrainRES),method="pearson")
Performance=matrix(c(Predabil.RES,corrBVRES,corrBVTrainRES))
write.table(Performance,row.names=FALSE,col.names=FALSE,file="CONCORDPERFORMRES.csv",sep=",")


