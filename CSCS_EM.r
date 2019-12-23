
####----This code works for replicate 1 under scenario II, lines that need to be 
####----modified to run another replicate/scenario are commented accordingly

#setwd("")

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


######################################################################
####----------Auxiliary function to be used later----------------#####
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
	     nonzeros=tmp1[which(abs(tmp1[,3])>thrsld),1:2] ##--Nodes
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


#####################################
#####---Function to run CSCS---######
#####################################

CSCS2 = function(S, n, lambda, L=NULL, maxitr=100, tol=1e-4, warmstart=FALSE) {

  #### inputs
  ## Y: n by p matrix of data
  ## lambda: l1-penalty parameter
  ## maxitr: maximum number of iterations allowed (diagonal/off-diagonal optimization)
  ## tol: if maximum difference of iterates is less than tol, consider converged
  ## warmstart: warmstarting actually made the runs slower (overhead may be too expensive)

  #### outputs
  ## L: lower triangular matrix of cholesky factor computed with CSCS algorithm
  ## itr_log: (p-1)-vector of number of iterations 
  ## eps_log: (p-1)-vector of number maximum difference for considering convergence
  
  p = ncol(S)

  if (is.null(L)) L = diag(p)

  itr_log = eps_log = NULL

  L[1, 1] = 1/sqrt(S[1,1])
  for (k in 2:p){ ## Loop in Algorithm 2

    ## nu_k vector (equation 2.5)
    nuk_old = nuk_new = c(rep(0, k-1), 1) 
    r = 0

    repeat {      ## Loop in Algorithm 1

      r = r + 1
      
      km1_ = 1:(k-1)    ## 1, ..., k-1 is off-diagonal elements indices

      ## Update off-diagonal terms
      hk = lassoshooting(XtX    =  S[km1_, km1_, drop=FALSE], 
			 Xty    = -nuk_old[k] * S[km1_, k],
			 lambda =  0.5*lambda)
      nuk_new[km1_] = hk$coefficients

      ## Update diagonal term
      sumterm = sum(nuk_new[km1_] * S[k, km1_])
      nuk_new[k] = (-sumterm + sqrt(sumterm^2 + 4*S[k, k]))/(2*S[k, k])

      ## Check convergence
      maxdiff = max(abs(nuk_new - nuk_old))
      if (maxdiff < tol || r >= maxitr){
	L[k, 1:k] = nuk_new
	eps_log = c(eps_log, maxdiff)
	itr_log = c(itr_log, r)
	break
      } else {
	nuk_old = nuk_new
      }
    }
  }

  list(L=L, itr=itr_log, eps=eps_log)

}


####################################################################################
########---Function implementing the combination of EM algorithm and CSCS----#######
####################################################################################

Bibi.GEMDAG=function(y,W,Omega0,res.var0,epsilon=1e-04,Maxiter=100,template,sizes,rho=0.01,thr=0.001){
	library(glasso)
      library(Matrix)
      library(lassoshooting)
	N=length(y)
      f=max(template[ ,3])
	nloci=ncol(W)
      npos=nloci*(nloci+1)*0.5
	iter=0
      
	Omega=as.matrix(Omega0)
      SIGMA=chol2inv(chol(Omega))
	res.var=res.var0
      delta=epsilon*2
      Im=diag(nloci)
      
    while(epsilon<delta && iter< Maxiter){
      Oldres.var = res.var
      OldOmega=Omega
	Sg.sum=matrix(0,nrow=nloci,ncol=nloci)
	trace=0
      inner=0
	    Block=list(rep(NA,f))

###---Expectation step, Using Formula 1

	for(i in 1:f){
		indices = template[which(template[,3]==i),1]
            WindSig = W[indices, ]%*%SIGMA
		Vinv=chol2inv(chol((WindSig)%*%t(W[indices, ])+(diag(sizes[i])*res.var)))
            Block[[i]]=Vinv
		summand=SIGMA%*%(Im - (t(W[indices, ])%*%(Vinv-crossprod(t(y[indices])%*%Vinv)))%*%WindSig)
		Sg.sum=Sg.sum+summand
		       }

	Sg.expect=Sg.sum/f 
      Vin=as.matrix(bdiag(Block))
	inner.res.expect=as.numeric(res.var*(N-res.var*(sum(diag(Vin))-crossprod(Vin%*%y))))

##---Maximization step

	res.var=inner.res.expect/N
	run=CSCS2(S=Sg.expect, n=N, lambda=rho, L=NULL, maxitr=50, tol=thr) 
      Omega=crossprod(run$L)
      SIGMA=tcrossprod(solve(run$L))
       
##---Evaluate convergence

	delta=(npos*mean(abs(Omega-OldOmega))+abs(Oldres.var-res.var))/(npos+1)
     
	iter=iter+1
                                   
                                     } #end of the while loop

      print(iter)
	if(delta>epsilon){print("Convergence criteria was not met")}
      Omega=as.matrix(Omega)
      SIGMA=as.matrix(SIGMA)
      L=as.matrix(run$L)
      Nodes=Def.Graph.Band(band.size=nloci,nloci=nloci)[[1]]
      tmp1=cbind(Nodes,as.matrix(Omega[lower.tri(Omega)]))
	nonzeros=tmp1[which(abs(tmp1[,3])>0),1:2] ##--i.e.Nodes
	sparsity.pct=100*(nrow(nonzeros)/(nloci*(nloci-1)/2))
     
	 	return(list(Omega,res.var,iter,delta,sparsity.pct,nonzeros,SIGMA,L))

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

lambda.grid=matrix(seq(0.002,0.004,length.out=15)) 
#lambda.grid=matrix(seq(0.0028,0.0045,length.out=15)) #use this for scenario 2

library(mvtnorm)
n.Train=nrow(W.Train)
Wt=t(W.Train)
W.Test=W[which(data1$G==4), ]

F=function(lambda){
 Run=Bibi.GEMDAG(y=data1.Training$Phen,W=W.Train,Omega0=Omega0,res.var0=93.23,epsilon=0.001,Maxiter=50,template=template,
  sizes=sizes,rho=lambda,thr=0.001)
   OMEGA=Run[[1]]
    SIGMA=Run[[7]]  
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
 return(list(PAR1,Run[[6]],Run[[8]]))
}

Out=apply(lambda.grid,1,F)

summary.tune=matrix(nrow=length(lambda.grid),ncol=8)
summary.tune[,1]=lambda.grid
for(i in 1:length(lambda.grid)){summary.tune[i,2:8]=Out[[i]][[1]]}
selected.graphLIK=Out[[which.max(summary.tune[,3])]][[2]]
selected.graphSPARSITY=Out[[8]][[2]]
selected.graphBIC=Out[[which.min(summary.tune[,2])]][[2]]
selected.graphRES=Out[[which.min(summary.tune[,4])]][[2]]
selected.graphPA=Out[[which.max(summary.tune[,5])]][[2]]

#################################################################################################
#############-------------------------Saving outputs---------------------########################
#################################################################################################

#-In the following lines replace "1" by "j", j=2,3,4 when runninganother replicate
#-In addition, add "II" at the end of each file's name when running scenario 2.

write.table(summary.tune,row.names=FALSE,col.names=c("Lambda","BIC","LIK","RES","PA","ACCV","ACCT","SPARSITY"),
file="Tuning1IIDAG.csv",sep=",")

write.table(selected.graphBIC,row.names=FALSE,col.names=FALSE,file="EdgesBICDAG1.csv",sep=",")

write.table(selected.graphLIK,row.names=FALSE,col.names=FALSE,file="EdgesLIKDAG1.csv",sep=",")

write.table(selected.graphRES,row.names=FALSE,col.names=FALSE,file="EdgesRESDAG1.csv",sep=",")

write.table(selected.graphSPARSITY,row.names=FALSE,col.names=FALSE,file="EdgesSPADAG1.csv",sep=",")

write.table(selected.graphPA,row.names=FALSE,col.names=FALSE,file="EdgesPADAG1.csv",sep=",")

###################################################################################################
############-----Save the 15 DAGS (only in scenario 2)---------####################################
###################################################################################################

#-Saving these 15 files is only necessary for scenario 2 because Bayes DAG-Sel only works when n>m
#-In the following line, replace "GRAPH1" by "GRAPHj", j=2,3,4 when analyzing another replicate. 

# for(i in 1:length(Out)){
#   write.table(as.matrix(Out[[i]][[3]]),row.names=FALSE,col.names=FALSE,file=paste("DAG1II",toString(i),".txt",sep=" "))
# }


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

############################################
####---Get true adjacency matrices----######
############################################

Adj1=Adj(Nodes=Om1[[2]],nloci)
Adj1=as.matrix(Adj1[upper.tri(Adj1)])

####--The following lines have to be run each time tha Adj.hat is defined

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


####---The following part requires setting the tuned values of the penalty parameter "rho"
####---which correspond to objects selected.graphLIK, selected.graphSPARSITY, 
####---selected.graphBIC, selected.graphRES and selected.graphPA defined above.


#############################################################
####---Comparisson with true concentration matrix------######
#############################################################

SPA=Bibi.GEMDAG(y=data1.Training$Phen,W=W.Train,Omega0=Omega0,res.var0=83.36,epsilon=0.001,Maxiter=50,template=template,
sizes=sizes,rho=0.0027,thr=0.001)

PA=Bibi.GEMDAG(y=data1.Training$Phen,W=W.Train,Omega0=Omega0,res.var0=83.36,epsilon=0.001,Maxiter=50,template=template,
sizes=sizes,rho=0.004,thr=0.001)

RES=LIK=PA #--This is so because in this particular case the three approaches selected the same graph.

NORM.DIFFSPA=norm(SPA[[1]]-Om1[[1]],type="F")/norm(Om1[[1]],type="F")
NORM.DIFFPA=norm(PA[[1]]-Om1[[1]],type="F")/norm(Om1[[1]],type="F")

##################################################
############----Save estimated Omega-----#########
##################################################

write.table(SPA[[1]],row.names=FALSE,col.names=FALSE,file="OmegaHatSPA1.csv",sep=",")
write.table(PA[[1]],row.names=FALSE,col.names=FALSE,file="OmegaHatPA1.csv",sep=",")

