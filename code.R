
rm(list=ls())
## Load Packages ################################################################
library(R.matlab)
library(pscl)
library(LaplacesDemon)
library(MASS)
library(gdata)
library(GIGrvg)
library(mvtnorm)
library(MCMCpack)
library(parallel)
library(doParallel)
library(foreach)

################################################################################
## Scenario - Simulated Data ########################################################
################################################################################

N <- 70                 ## Sample size
V <- 30                 ## No. of nodes
q <- V*(V-1)/2       ## No. of elements in the upper triangle of the symmetric tensor
R.generate <- 2    ## true dim of latent variables
pi <- 0.4                ## true node sparsity
umean <- 0.8        ## mean of true PARAFAC modes
usd <- 0.5             ## s.d. of the true PARAFAC modes
n.rep <- 50            ## no. of replications

response <- list()
predictor <- list()
indic.coef1 <- list()
coef1 <- list()
coef0 <- numeric()
coef2 <- numeric()
coef3 <- numeric()

for(rep in 1:n.rep){
   X1 <- rnorm(N) ## Simulated predictor of interest
   X2 <- rnorm(N) ## Simulated auxiliary predictor
   X3 <- rnorm(N) ## Simulated auxiliary predictor

   ## Create true beta0, beta2, beta3, sigma #############################################
   true.b0  <- 0.2   ## true intercept
   true.b2  <- 0.4   ## true auxiliary predictor coefficient
   true.b3  <- -0.1  ## true auxiliary predictor coefficient
   sigma.tr <- 1      ## true error variance

   seq.indx <- c(0,0)
   for(i in 1:(V-1)){
      for(j in (i+1):V){
       seq.indx <- rbind(seq.indx,c(i,j))
      }
   }
   seq.indx <- seq.indx[-1,]

   ## Create true beta1 using tensor structure #############################################

   U.gen <- matrix(NA,nrow=V,ncol=R.generate)
   W.gen <- numeric()

   ## Sample Locations for 1

   num.1          <- floor(pi*V)
   loc.1          <- sample(1:V,num.1,replace = F)
   bin.gen        <- c(rep(0,V))
   bin.gen[loc.1] <- 1   ## which tensor nodes are active

   for(i in 1:V){
         if (bin.gen[i] == 1) {U.gen[i,] <- rnorm(R.generate,umean,usd)
    	} else {U.gen[i,] <- c(rep(0,R.generate))}
   }

   for(i in 1:nrow(seq.indx)){  
     W.gen[i] <- sum(U.gen[seq.indx[i,1],]*U.gen[seq.indx[i,2],])
   }

   true.b1 <- W.gen  ## true coefficient corresponding to the predictor of interest

   #################################################################################
   ## Generate binary response which are elements in the network matrix ########################
   #################################################################################
   
   pred <- matrix(NA,nrow=N,ncol=q)   ## predictor, each row for an individual
   y  <- matrix(NA,nrow=N,ncol=q)       ## response, each for an individual
   
   for (count in 1:N){
      pred[count,] <- true.b0 + true.b1*X1[count] + true.b2*X2[count] + true.b3*X3[count]
      y[count,] <- rnorm(q,pred[count,],sqrt(sigma.tr))
   }
   response[[rep]]  <- y                         ## response 
   predictor[[rep]] <- cbind(X1,X2,X3)  ## predictors
   coef1[[rep]] <- true.b1
   indic.coef1[[rep]] <- bin.gen
   coef0[rep] <- true.b0
   coef2[rep] <- true.b2
   coef3[rep] <- true.b3
}

####################################################################################
## Running replications of the model in parallel ##############################################
####################################################################################

replicated_2d_cont <- function(rep){
    library(R.matlab)
    library(pscl)
    library(LaplacesDemon)
    library(MASS)
    library(gdata)
    library(GIGrvg)
    library(mvtnorm)
    library(MCMCpack)

    y     <- response[[rep]]
    X1  <- predictor[[rep]][,1]
    X2  <- predictor[[rep]][,2]
    X3  <- predictor[[rep]][,3] 
    N    <- length(X1)           ## Sample size
    q    <- dim(y)[2]
    V    <- (1+sqrt(1+8*q))/2
    R    <- 4                          ## Maximum dimension of the latent vector
    a    <- 1
    b    <- 1
    a.wish  <- 10
    S.scale <- diag(rep(1,R))
    niter   <- 15000               ## No. of Iterations
    eta     <- 1.01
    a.ghap  <- numeric()
    b.ghap  <- numeric()

    for(r in 1:R){
       a.ghap[r] <- 1
       b.ghap[r] <- r^(eta)
    }

    ## Initialize ##

    beta0 <- c(rep(NA,niter))
    W1 <- matrix(NA,nrow=niter,ncol=q)
    beta2 <- c(rep(NA,niter))
    beta3 <- c(rep(NA,niter))
    lambda <- c(rep(NA,niter))
    ulist <- list()
    for (i in 1:niter){
	  ulist[[i]] <- matrix(NA,nrow=V,ncol=R)  ## u1,...,uV for each iteration.
    }
    rn.gen <- matrix(NA,niter,V)   ## Denotes whether a node is active or not (0 or 1 in every iteration) 
    kappa <- matrix(NA,niter,R)   ## 0 or 1 corresponding to whether a dimension is active or not.
    pi.ghap <- matrix(NA,niter,R) ## Prob. corresponding to each kappa
    pi <- numeric()
    M <- list()   ## Denotes M matrix
    sigma <- numeric()

    ####################################################################################
    ## Gibbs Sampler  #####################################################################
    ## Initial Values ## The list counts actually count the iterate (not the observation) ###################
    ####################################################################################

    for(r in 1:R){
       pi.ghap[1,r] <- rbeta(1,a.ghap[r],b.ghap[r])
       kappa[1,r] <- rbinom(1,1,pi.ghap[1,r])
    }
    beta0[1] <- 0.20
    beta2[1] <- true.b2
    beta3[1] <- true.b3
    pi[1] <- 0.5
    M[[1]] <- riwish(a.wish,S.scale)
    
    ## Each of the 'niter' no. of matrices in the list represents u1,u2,u3,...,uV. (Each of the u's is of dimension R.)
    ulist[[1]] <- matrix(rnorm(V*R,0,1),nrow=V,ncol=R)
    Byatha <- diag(kappa[1,])
    U.first <- ulist[[1]]
    for(i in 1:nrow(seq.indx)){  
       W1[1,i] <- sum(U.first[seq.indx[i,1],]*U.first[seq.indx[i,2],])
    }

   rn.gen[1,] <- rbinom(V,1,0.5)
   sigma[1] <- 1

   ########################################################################################
   ## To find out which edges are associated with a certain node ######################################
   ## Each row in "whichedges" gives the edge numbers associated with a certain node "k" #################
   ########################################################################################

   whichedges <- matrix(NA,V,(V-1))
   indx.rest <- list()
   for (k in 1:V){
      type1.indx <- which(seq.indx[,1] == k)
      type2.indx <- which(seq.indx[,2] == k)
      whichedges[k,] <- c(type1.indx,type2.indx)
      indx1 <- seq.indx[which(seq.indx[,1]==k),2]
      indx2 <- seq.indx[which(seq.indx[,2]==k),1]
      indx.rest[[k]] <- c(indx1,indx2)
   }

   #########################################################################################
   ## Run Gibbs iterations ######################################################################

   for (i in 2:niter){
        
       ## (1) Update beta0

       mean.beta0.denom <- (N*q)/sigma[i-1] + 1
       mean.beta0.num   <- sum(sapply(1:N,function(ll){
                           sum(y[ll,] - W1[i-1,]*X1[ll] - 
                              beta2[i-1]*X2[ll] - beta3[i-1]*X3[ll])}))/sigma[i-1]
    
       mean.beta0 <- mean.beta0.num/mean.beta0.denom
       var.beta0  <- 1/mean.beta0.denom   
       beta0[i]   <- rnorm(1,mean.beta0,sqrt(var.beta0)) ## update intercept
        
       ## (2) Update sigmasq
 
       sum_sigma <- 0      
       for(count in 1:N){
    	    par.omega <- beta0[i] + W1[i-1,]*X1[count] + beta2[i-1]*X2[count] + beta3[i-1]*X3[count]
    	    sum_sigma <- sum_sigma+sum((y[count,]-par.omega)^2)
       }
       a.sig    <- 1
       b.sig    <- 1
       sigma[i] <- rigamma(1,a.sig+N*q/2,b.sig+sum_sigma/2) ## update error variance
    
       U.first <- ulist[[i-1]]
       for(ii in 1:nrow(seq.indx)){  
           W1[i-1,ii] <- sum(U.first[seq.indx[ii,1],]*U.first[seq.indx[ii,2],]*kappa[i-1,])
       }


      ## (3) Update beta2 and beta3 together

      X23       <- cbind(X2,X3)
      var.beta  <- matrix(0,2,2)
      mean.beta <- rep(0,2)
      for (j in 1:q){
    	    Z_tilda <- y[,j] - beta0[i] - matrix(X1)*W1[i-1,j]
          var.beta  <- var.beta+t(X23)%*%X23/sigma[i]
          mean.beta <- mean.beta+c(t(X23)%*%c(Z_tilda))/sigma[i]
      }
      var.b23  <- solve(var.beta+diag(2))
      mean.b23 <- var.b23%*%mean.beta
      beta23   <- rmvnorm(1,mean.b23,var.b23)
      beta2[i] <- beta23[1]                   ## Update auxiliary coefficient 1
      beta3[i] <- beta23[2]                   ## Update auxiliary coefficient 2
    
      
      ## (4) Update u_k (Each u_k for k=1,...,V is an R-dim vector)
      ## Update ulist[[i]].Has V rows for u1,u2,u3,...,uV, and R columns for R-dim uk's.

    
      X10 <- matrix(X1)
      Z_tilda2_all <- matrix(NA,N,q)
      for (j in 1:q){
     	  Z_tilda2_all[,j] <- y[,j] - beta0[i] - matrix(X2)*beta2[i] - matrix(X3)*beta3[i]
      }    
   
      for(k in 1:V){
         
         ## Find Zu

         Kappa_big <- matrix(rep(kappa[i-1,],(V-1)),(V-1),R,byrow=T)
         ubar      <- ulist[[i-1]][indx.rest[[k]],]*Kappa_big                 
         Zu        <- list()
         Z_tilda2  <- list()    
        
        Zu[[k]] <- matrix(NA,1,R)
        for (w in 1:(V-1)){
    	        Zu[[k]] <- rbind(Zu[[k]],X10%*%ubar[w,])
        }
        Zu[[k]] <- Zu[[k]][-1,]
            
        ## Find Z_tilda2
        
        Z_tilda2[[k]] <- c(NA)
        for (j in 1:(V-1)){
            Z_tilda2[[k]] <- c(c(Z_tilda2[[k]]),Z_tilda2_all[,whichedges[k,j]])
        }
        Z_tilda2[[k]] <- Z_tilda2[[k]][-1]   
      
        ## Find Omega_k

        cov.uk      <- solve(t(Zu[[k]])%*%Zu[[k]]/sigma[i]+solve(M[[i-1]]))
        mean.uk     <- cov.uk%*%t(Zu[[k]])%*%c(Z_tilda2[[k]])/sigma[i]
        first.term  <- cov.uk 
        sec.term    <- c(t(Zu[[k]])%*%c(Z_tilda2[[k]]))/sigma[i]
        third.term  <- sum(Z_tilda2[[k]]^2)/sigma[i]
        exponent    <- (-1/2)*(third.term-crossprod(sec.term,c(first.term%*%sec.term)))
        det.log     <- (-(N*(V-1))/2)*log(2*3.14159)+(1/2)*log(det(first.term))-((N*(V-1))/2)*log(sigma[i])-(1/2)*log(det(M[[i-1]]))
        pi1         <- c(exponent)+det.log
        pi2         <- sum(dnorm(Z_tilda2[[k]],0,sqrt(sigma[i]),log=T))
        piup        <- exp(pi2-pi1)
    	rn.gen[i,k] <- rbinom(1,1,pi[i-1]/(pi[i-1]+(1-pi[i-1])*piup))        ## Whether the k-th node is active
        if(rn.gen[i,k]==1){uk <- mvrnorm(1,mean.uk, cov.uk)            ## if active, draw from the slab component
        }else{uk <- c(rep(0,R))}                                                          ## else draw from the spike component
        ulist[[i]][k,] <- uk      
        }

    U.i <- ulist[[i]]
    for(r in 1:R){
        kappa.temp1    <- kappa[i-1,]
        kappa.temp2    <- kappa[i-1,]
        kappa.temp1[r] <- 1
        kappa.temp2[r] <- 0
        W1.temp1 <- sapply(1:nrow(seq.indx),function(ii){sum(U.i[seq.indx[ii,1],]*U.i[seq.indx[ii,2],]*kappa.temp1)})
        W1.temp2 <- sapply(1:nrow(seq.indx),function(ii){sum(U.i[seq.indx[ii,1],]*U.i[seq.indx[ii,2],]*kappa.temp2)})     
        prob.up <- sum(sapply(1:N,function(count){sum(dnorm(y[count,],beta0[i] + W1.temp1*X1[count] + beta2[i]*X2[count] + beta3[i]*X3[count] ,sqrt(sigma[i]),log=T))}))
        prob.lo  <- sum(sapply(1:N,function(count){sum(dnorm(y[count,],beta0[i] + W1.temp2*X1[count] + beta2[i]*X2[count] + beta3[i]*X3[count],sqrt(sigma[i]),log=T))}))
        puplo  <- exp(prob.lo-prob.up)
        kappa[i,r] <- rbinom(1,1,pi.ghap[i-1,r]/(pi.ghap[i-1,r]+(1-pi.ghap[i-1,r])*puplo))  ## Decide if r-th latent dimension is informative
    }
    
    for(ii in 1:nrow(seq.indx)){  
       W1[i,ii] <- sum(U.i[seq.indx[ii,1],]*U.i[seq.indx[ii,2],]*kappa[i,])   ## Update coefficient corresponding to predictor of interest
    }

    for(r in 1:R){
        pi.ghap[i,r] <- rbeta(1,a.ghap[r]+kappa[i,r],b.ghap[r]+1-kappa[i,r])   ## Update probability corresponding to the r-th latent dimension
    }

    #### (5) Update pi

    pi[i] <- rbeta(1,(a+sum(rn.gen[i,])),(b+V-sum(rn.gen[i,])))  ## Update probability of the nonzero mixture component in the node selection

    #### (4) Update matrix M

    len.nonzero.ui <- length(which(rn.gen[i,]==1))
    U.U            <- matrix(0,R,R)
    for(j in 1:V){ 
    	U.U <- U.U+tcrossprod(ulist[[i]][j,])
    }
    M[[i]] <- riwish(a.wish+len.nonzero.ui,S.scale+U.U)  
    }  
    hh <- list(W1,kappa,beta0,beta2,beta3,sigma,rn.gen)
    names(hh) <- c("W1","kappa","beta0","beta2","beta3","sigma","indic.latvar")
    return(hh)  
}

cl <- makeCluster(n.rep)  
registerDoParallel(cl)

## Start time

strt<-Sys.time()

## Parallelized computation of replications in different cores

obj <- foreach(rep =1:n.rep) %dopar% replicated_2d_cont(rep)  

## Total time for parallelized inference

final.time <- Sys.time()-strt  
stopCluster(cl)

### MSE - b1

R <- 4
MSE1 <- numeric()
MSE2 <- numeric()
MSE3 <- numeric()
avg.cov <- numeric()
avg.length <- numeric()
node_detect <- matrix(NA,V,2*n.rep)
eff_dim <- matrix(NA,(R+1),n.rep)
niter <- 15000
for(rep in 1:n.rep){
   ee                                    <- obj[[rep]]$W1
   ee2                                  <- obj[[rep]]$beta2
   ee3                                  <- obj[[rep]]$beta3
   ee4                                  <- obj[[rep]]$kappa
   MSE                                <- mean((coef1[[rep]] - colMeans(ee[seq(5001,niter,by=2),]))^2)
   MSE                                <- MSE/mean(coef1[[rep]]^2)     
   cov.method                     <- cbind(apply(ee[seq(5001,niter,by=2),],2,quantile,0.975),coef1[[rep]],apply(ee[seq(5001,niter,by=2),],2,quantile,0.025))
   node_detect[,(2*(rep-1)+1):(2*rep)] <- cbind(indic.coef1[[rep]],colMeans((obj[[rep]]$indic.latvar)[seq(5001,niter,by=2),]))
   MSE1[rep]                       <- MSE                                                                            ## Scaled MSE for predictor of interest
   MSE2[rep]                       <- (coef2[rep] - mean(ee2[seq(5001,niter,by=2)]))^2/coef2[rep]^2  ## scaled MSE for auxiliary predictor 1
   MSE3[rep]                       <- (coef3[rep] - mean(ee3[seq(5001,niter,by=2)]))^2/coef3[rep]^2  ## scaled MSE for auxiliary predictor 2
   avg.cov[rep]                    <- length(intersect(which(cov.method[,3]<=cov.method[,2]),
                                                        which(cov.method[,1]>=cov.method[,2])))/q        ## Avg. coverage of 95% CI for the pred. coeff. of interest
   avg.length[rep]                <- mean(cov.method[,1]-cov.method[,3])                          ## Avg. lengfth of 95% CI for the pred. coeff. of interest
   for(r in 1:(R+1)){
      eff_dim[r,rep] <- length(which(rowSums(ee4[seq(5001,niter,by=2),])==(r-1)))/((niter-5000)/2)       ## Probability dist. of effective dimensionality
   }
}

####################################################################################
## True and False Positive Rates for node identification with different choices of cut-off ###############
####################################################################################

vec.q <- c(0,0.01,0.25,0.50,0.75,0.99,1)
TPR <- matrix(NA,n.rep,length(vec.q))
FPR <- matrix(NA,n.rep,length(vec.q))
for(rep in 1:n.rep){
   for(jj in 1:length(vec.q)){
     TPR[rep,jj] <- length(intersect(which(node_detect[,(2*(rep-1)+1)]==1),which(node_detect[,2*rep]>=vec.q[jj])))/
                    length(which(node_detect[,(2*(rep-1)+1)]==1))
     FPR[rep,jj] <- length(intersect(which(node_detect[,(2*(rep-1)+1)]==0),which(node_detect[,2*rep]>=vec.q[jj])))/
                    length(which(node_detect[,(2*(rep-1)+1)]==0)) 
  }
}



