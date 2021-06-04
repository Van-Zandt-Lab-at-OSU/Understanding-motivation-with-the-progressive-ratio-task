
##########################################
# Figure 5: Simulation

# function to plot whiskers
posterior_bars <- function(param,tru,name,sim,I,cond_num,prior.var){
   p <- rnorm(10000,0,sqrt(prior.var))
   plot.new()
   plot.window(xlim=c(-1,I+10),ylim=c(min(param),max(param)))
   axis(1)
   axis(2)
   quan <- quantile(p,probs=c(0.025,0.5,0.975))
   lines(x=c(0,0),y=c(quan[1],quan[3]),col='orange')  
   points(0,quan[2],col='purple')
   for(i in 1:I){
      quan <- quantile(param[,,i],probs=c(0.025,0.5,0.975))
      lines(x=c(i,i),y=c(quan[1],quan[3]),col='blue')
      points(i,quan[2],col='blue')
      if(sim==1){
         points(i,tru[i],col='red')
      }
      abline(v=(cond_num[2]+cond_num[3])/2,col="red")
      abline(v=(cond_num[4]+cond_num[5])/2,col="red")
   }
   abline(h=0,col="light grey")
   title(name,xlab="Participants",ylab="Value",line=2.5)
}

# simulated data
load("simulation_data_90.Rdata")

# result from simulation
load("simulation_result_90.Rdata")

I <- dat$I
num_id <- dat$num_id

par(mfrow=c(4,2),mai=c(0.5,0.5,0.5,0.5))

# plots. dat$t_param are true parameter values used to generate the simulated data

posterior_bars(fit$beta_ic,dat$t_param$beta,expression(~beta[ic]),1,I,num_id,3*25)
posterior_bars(fit$mu_10_ic,dat$t_param$mu_10,expression({~mu^{"10"}} [ic]),1,I,num_id,3*25)
posterior_bars(fit$mu_01_ic,dat$t_param$mu_01,expression({~mu^{"01"}} [ic]),1,I,num_id,3*25)
posterior_bars(fit$gamma_p1_ic,dat$t_param$gamma_p1,expression({~gamma^{"p"}} [ic]),1,I,num_id,3*9)
posterior_bars(fit$gamma_p2_ic,dat$t_param$gamma_q1,expression({~gamma^{"q"}} [ic]),1,I,num_id,3*9)
posterior_bars(fit$eta_p1_ic,dat$t_param$eta_p1,expression({~eta^{"p"}} [ic]),1,I,num_id,3*9)
posterior_bars(fit$eta_p2_ic,dat$t_param$eta_q1,expression({~eta^{"q"}} [ic]),1,I,num_id,3*9)




###################################################
######## Figure 8


load("simulation_data_90.Rdata")

I <- dat$I

post_predictive_start <- function(mu,mu_10,mu_01,psi,eta_p1,eta_q1,eta_10,eta_01,
                                  gamma_p1,gamma_q1,gamma_10,gamma_01,
                                  V_dum1,V_dum2,M,tim){
   L <- length(V_dum1)
   X <- rep(NA,L)
   for(l in 1:L){
      V1 <- (mu + mu_10*V_dum1[l] + mu_01*V_dum2[l] + eta_p1*M[l]
             + gamma_p1*tim[l]  + psi*tim[l]*M[l]
             + (eta_10*V_dum1[l]+eta_01*V_dum2[l])*M[l]
             + (gamma_10*V_dum1[l]+gamma_01*V_dum2[l])*tim[l]) 
      p1 <- 1/(1+exp(-V1))
      X[l] <- rbinom(1,1,p1)
   }
   return(X)
}

post_predictive_middle <- function(mu,mu_10,mu_01,psi,eta_p1,eta_q1,eta_10,eta_01,
                                   gamma_p1,gamma_q1,gamma_10,gamma_01,
                                   V_dum1,V_dum2,M,tim){
   L <- length(V_dum1)
   X <- rep(NA,L)
   for(l in 1:L){
      V1 <- (mu + mu_10*V_dum1[l] + mu_01*V_dum2[l] + eta_p1*M[l]
             + gamma_p1*tim[l]  + psi*tim[l]*M[l]
             + (eta_10*V_dum1[l]+eta_01*V_dum2[l])*M[l]
             + (gamma_10*V_dum1[l]+gamma_01*V_dum2[l])*tim[l]) 
      p1 <- 1/(1+exp(-V1))
      V2 <- (mu + eta_q1*M[l]
             + gamma_q1*tim[l]  + psi*tim[l]*M[l])
      p2 <- 1/(1+exp(-V2))
      X[l] <- rbinom(1,1,p2*(1-p1))
   }
   return(X)
}

load("simulation_result_90.Rdata")

I <- dat$I

dat[["cond_num"]] <- NULL
dat[["I"]] <- NULL
dat[["N_RT"]] <- NULL
dat[["num_id"]] <- NULL
dat[["N"]] <- NULL
dat[["rtlabel"]] <- NULL
dat[["RT"]] <- NULL
dat[["t_param"]] <- NULL # only for simulated


df <- data.frame(matrix(unlist(dat), nrow = max(lengths(dat)), byrow = F))
names(df) <- c("X","Y","V0","V_dum1","V_dum2","M","tim","ilabel")


b_break <- array(-5,dim=c(3000,I,21))
m_break <- array(-5,dim=c(3000,I,21))
for(i in 1:I){
   sub <- subset(df,ilabel==i)
   L <- dim(sub)[1]
   for(j in 1:3000){
      b_break[j,i,(1:L)] <- post_predictive_start(fit$beta_ic[j,1,i],fit$mu_10_ic[j,1,i],fit$mu_01_ic[j,1,i],
                                                  fit$psi_ic[j,1,i],fit$eta_p1_ic[j,1,i],fit$eta_q1_ic[j,1,i],
                                                  fit$eta_10_ic[j,1,i],fit$eta_01_ic[j,1,i],fit$gamma_p1_ic[j,1,i],
                                                  fit$gamma_q1_ic[j,1,i],fit$gamma_10_ic[j,1,i],fit$gamma_01_ic[j,1,i],
                                                  sub$V_dum1,sub$V_dum2,log(sub$M),log(sub$tim))
      m_break[j,i,(1:L)] <- post_predictive_middle(fit$beta_ic[j,1,i],fit$mu_10_ic[j,1,i],fit$mu_01_ic[j,1,i],
                                                   fit$psi_ic[j,1,i],fit$eta_p1_ic[j,1,i],fit$eta_q1_ic[j,1,i],
                                                   fit$eta_10_ic[j,1,i],fit$eta_01_ic[j,1,i],fit$gamma_p1_ic[j,1,i],
                                                   fit$gamma_q1_ic[j,1,i],fit$gamma_10_ic[j,1,i],fit$gamma_01_ic[j,1,i],
                                                   sub$V_dum1,sub$V_dum2,log(sub$M),log(sub$tim))
   }
#   for(j in 1:1500){
#      b_break[j+1500,i,(1:L)] <- post_predictive_start(fit$beta_ic[j,2,i],fit$mu_10_ic[j,2,i],fit$mu_01_ic[j,2,i],
#                                                       fit$psi_ic[j,2,i],fit$eta_p1_ic[j,2,i],fit$eta_q1_ic[j,2,i],
#                                                       fit$eta_10_ic[j,2,i],fit$eta_01_ic[j,2,i],fit$gamma_p1_ic[j,2,i],
#                                                       fit$gamma_q1_ic[j,2,i],fit$gamma_10_ic[j,2,i],fit$gamma_01_ic[j,2,i],
#                                                       sub$V_dum1,sub$V_dum2,log(sub$M),log(sub$tim))
#      m_break[j+1500,i,(1:L)] <- post_predictive_middle(fit$beta_ic[j,2,i],fit$mu_10_ic[j,2,i],fit$mu_01_ic[j,2,i],
#                                                        fit$psi_ic[j,2,i],fit$eta_p1_ic[j,2,i],fit$eta_q1_ic[j,2,i],
#                                                        fit$eta_10_ic[j,2,i],fit$eta_01_ic[j,2,i],fit$gamma_p1_ic[j,2,i],
#                                                        fit$gamma_q1_ic[j,2,i],fit$gamma_10_ic[j,2,i],fit$gamma_01_ic[j,2,i],
#                                                        sub$V_dum1,sub$V_dum2,log(sub$M),log(sub$tim))
#   }
   print(i)
}


# plot it. make the never-reached sets empty

par(mfrow=c(4,4),mai=c(0.4,0.4,0.4,0.4))

for(i in 1:I){
   sub <- subset(df,ilabel==i)
   # start of session plot
   L <- rep(0,21)
   value <- c(50,25,10)
   for(c in 1:3){
      L[(7*c-6):(7*c-7+sum(sub$V0==value[c]))] <- 1
   }
   b.mean <- rep(0,21)
   m.mean <- rep(0,21)
   X.0 <- rep(NA,21)
   Y.0 <- rep(NA,21)
   k <- 1
   for(c in 1:21){
      if(L[c]==1){
         b.mean[c] <- mean(b_break[,i,k])
         m.mean[c] <- mean(m_break[,i,k])
         X.0[c] <- sub$X[k]
         Y.0[c] <- sub$Y[k]
         k <- k+1
      }
   }
   a <- barplot(b.mean,ylim=c(0,1.1),main=paste("participant",i,"quit at initialization"))
   text(a,b.mean+0.1,X.0)
   title(xlab="sets",mgp=c(1,1,0), cex.lab=1.2)
   title(ylab="quit probability",mgp=c(2,1,0), cex.lab=1.2)
   #  mid-session plot
   a <- barplot(m.mean,ylim=c(0,1.1),main=paste("participant",i,"quit within"),col="blue")
   text(a,m.mean+0.1,Y.0)
   title(xlab="sets",mgp=c(1,1,0), cex.lab=1.2)
   title(ylab="quit probability",mgp=c(2,1,0), cex.lab=1.2)
}




#####################################
############ Figure 11


load("simulation_data_90.Rdata")

load("simulation_result_90.Rdata")

I <- dat$I
cond_num <- dat$cond_num

M <- dat$M
tim <- dat$tim

beta_c <- rbind(fit$beta_c[,1,])
mu_10_c <- rbind(fit$mu_10_c[,1,])
mu_01_c <- rbind(fit$mu_01_c[,1,])
eta_p1_c <- rbind(fit$eta_p1_c[,1,])
eta_q1_c <- rbind(fit$eta_q1_c[,1,])
eta_10_c <- rbind(fit$eta_10_c[,1,])
eta_01_c <- rbind(fit$eta_01_c[,1,])
gamma_p1_c <- rbind(fit$gamma_p1_c[,1,])
gamma_q1_c <- rbind(fit$gamma_q1_c[,1,])
gamma_10_c <- rbind(fit$gamma_10_c[,1,])
gamma_01_c <- rbind(fit$gamma_01_c[,1,])
psi_c <- rbind(fit$psi_c[,1,])

########

# k: iteration

posteriors_group <- function(K,I=cond_num,V_dum1,V_dum2,M,tim){
   # the probabilities
   p1 <- array(dim=c(K,I))
   q1 <- array(dim=c(K,I))
   for(i in 1:I){
      v <- (beta_c[,i]  + mu_10_c[,i]*V_dum1 + mu_01_c[,i]*V_dum2
            + eta_p1_c[,i]*M + (eta_10_c[,i]*V_dum1 + eta_01_c[,i]*V_dum2)*M
            + gamma_p1_c[,i]*tim + (gamma_10_c[,i]*V_dum1 + gamma_01_c[,i]*V_dum2)*tim
            + psi_c[,i]*M*tim)
      p1[,i] <- 1/(1+exp(-v))
      v <- (beta_c[,i]
            + eta_q1_c[,i]*M
            + gamma_q1_c[,i]*tim
            + psi_c[,i]*M*tim) 
      q1[,i] <- 1/(1+exp(-v))
   }
   result <- list(p1=p1,q1=q1)
   return(result)
}

# The plots
# only time varies
K <- 60000
C <- 6
L <- 6

p_group <- array(dim=c(K,C,cond_num,L))
q_group <- array(dim=c(K,C,cond_num,L))

tim_seq <- quantile(log(tim),probs=seq(0.1,0.9,length.out=C))
M_seq <- log(c(12,26,45,100,167,500))
M_seq_2 <- c(12,26,45,100,167,500)

for(l in 1:length(M_seq)){
   for(c in 1:C){
      res2 <- posteriors_group(K,cond_num,V_dum1=0,V_dum2=0,M=M_seq[l],tim=tim_seq[c])
      p_group[,c,,l] <- res2$p1
      q_group[,c,,l] <- res2$q1 
      print(c)
   }
}

p <- p_group + (1-p_group)*q_group

# plot the 3 group levels
par(mfrow=c(8,1),mai=c(0.02,0.8,0.02,0.8))


plot.new()
plot.window(xlim=c(0,20),ylim=c(0,10))
text(0.5,1,"time")
for(c in 1:C){
   text(3*c-1,1,paste(trunc(exp(tim_seq[c])),"s"))
}
for(l in 1:L){
   plot.new()
   plot.window(xlim=c(0,20),ylim=c(0,1))
   axis(2)
   for(c in 1:C){
      boxplot(p[,c,1,l],at=3*c-2,col="grey",outcol="light grey",outpch=20,add=T)
      boxplot(p[,c,2,l],at=3*c-1,col="blue",outcol="light grey",outpch=20,add=T)  
      boxplot(p[,c,3,l],at=3*c,col="red",outcol="light grey",outpch=20,add=T) 
      if(c<C) abline(v=3*c+0.5,col="green")
   }
   text(19.5,0.5,paste((M_seq_2[l]),"trials"))
   title(ylab="quit probability",mgp=c(2,1,0))
}

plot.new()
plot.window(xlim=c(0,20),ylim=c(0,10))
for(c in 1:C){
   text(3*c-2,9,"con")
   text(3*c-1,9,"rel")
   text(3*c,9,"sch")
}






