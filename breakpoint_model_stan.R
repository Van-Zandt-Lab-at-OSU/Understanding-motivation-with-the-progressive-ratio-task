
Model_breakpoint <- '
functions {
   real binom_dist_lpmf(int X, int Y,real p, real q){
      real dens;
// finished all trials
      if(X+Y==0){
         dens =  log(1-p) + log(1-q);
      }
// quit at initialization screen
      else if(X==1){
         dens = log(p);
      }
// quit after initialization, within the set
      else{
         dens = log(1-p) + log(q);
      }
      return(dens);
   }
}
data {
   int<lower=0> I;                   //number of subjects
   int<lower=0> N;                   //number of responses from all subjects
   int<lower=0> cond_num;            //number of groups (3)
   int<lower=0> num_id[2*cond_num];  // identify the start and end ID of each group, 
                                     // without schizophrenia: 1-37, relatives: 38-77, schizophrenia: 78-118
   int X[N];                         //initialization reponses. 1=quit, 0=proceed.
   int Y[N];                         //within set reponses. 1=quit somewhere, 0=finish.
   int V_dum1[N];                    //reward value dummy variable, =1 when reward is 25  
   int V_dum2[N];                    //reward value dummy variable, =1 when reward is 10  
   real M[N];                        //effort at current set    
   real tim[N];                      //elapsed time at current set
   int ilabel[N];                    //subject label of the response
}

//log transformation to add discounting effects
transformed data {
   real M_main[N];              
   real tim_main[N];

   for(n in 1:N){
      M_main[n] = log(M[n]);
      tim_main[n] = log(tim[n]);
   }
}

parameters {

// parameters for quit probability
// use p1 and q1 instead of p and q in parameters, because they are somehow not allowed!!!
   real eta_p1;
   real eta_q1;
   real eta_10;
   real eta_01;
   real gamma_p1;
   real gamma_q1;
   real gamma_10;
   real gamma_01;
   real beta;
   real mu_10;
   real mu_01;
   real psi;

// group level

   real eta_p1_c[cond_num];
   real eta_q1_c[cond_num];
   real eta_10_c[cond_num];
   real eta_01_c[cond_num];
   real gamma_p1_c[cond_num];
   real gamma_q1_c[cond_num];
   real gamma_10_c[cond_num];
   real gamma_01_c[cond_num];
   real beta_c[cond_num];
   real mu_10_c[cond_num];
   real mu_01_c[cond_num];
   real psi_c[cond_num];
   
// individual level

   real gamma_p1_ic[I];
   real gamma_q1_ic[I];
   real beta_ic[I];
   real eta_p1_ic[I];
   real eta_q1_ic[I];
   real mu_10_ic[I];
   real mu_01_ic[I];
}

transformed parameters {
   real<lower=0,upper=1> p[N];
   real<lower=0,upper=1> q[N];
   real eta_10_ic[I];
   real eta_01_ic[I];
   real gamma_10_ic[I];
   real gamma_01_ic[I];
   real psi_ic[I]; 
   real V;
   for(c in 1:cond_num){
      for(i in num_id[2*c-1]:num_id[2*c]){
        psi_ic[i] = psi_c[c];
        eta_10_ic[i] = eta_10_c[c];
        eta_01_ic[i] = eta_01_c[c];
        gamma_10_ic[i] = gamma_10_c[c];
        gamma_01_ic[i] = gamma_01_c[c];
      }
   }
   for(n in 1:N){
      V = beta_ic[ilabel[n]]+ mu_10_ic[ilabel[n]]*V_dum1[n] + mu_01_ic[ilabel[n]]*V_dum2[n]
          + (eta_p1_ic[ilabel[n]] + eta_10_ic[ilabel[n]]*V_dum1[n] + eta_01_ic[ilabel[n]]*V_dum2[n])*M_main[n] 
          + (gamma_p1_ic[ilabel[n]] + gamma_10_ic[ilabel[n]]*V_dum1[n] + gamma_01_ic[ilabel[n]]*V_dum2[n])*tim_main[n]
          + psi_ic[ilabel[n]]*M_main[n]*tim_main[n];
      p[n] = inv_logit(V);
      V = beta_ic[ilabel[n]]
          + (eta_q1_ic[ilabel[n]])*M_main[n] 
          + (gamma_q1_ic[ilabel[n]])*tim_main[n]
          + psi_ic[ilabel[n]]*M_main[n]*tim_main[n];
      q[n] = inv_logit(V);
   }      
}
model {  

   eta_p1 ~ normal(0.,3.);
   eta_q1 ~ normal(0.,3.);
   eta_10 ~ normal(0.,3.);
   eta_01 ~ normal(0.,3.);
   gamma_p1 ~ normal(0.,3.);
   gamma_q1 ~ normal(0.,3.);
   gamma_10 ~ normal(0.,3.);
   gamma_01 ~ normal(0.,3.);
   psi ~ normal(0.,3.);
   beta ~ normal(0.,5.);
   mu_10 ~ normal(0.,5.);
   mu_01 ~ normal(0.,5.);
   
   for(c in 1:cond_num){

      eta_p1_c[c] ~ normal(eta_p1,3.);
      eta_q1_c[c] ~ normal(eta_q1,3.);
      eta_10_c[c] ~ normal(eta_10,3.); 
      eta_01_c[c] ~ normal(eta_01,3.); 
      gamma_p1_c[c] ~ normal(gamma_p1,3.);
      gamma_q1_c[c] ~ normal(gamma_q1,3.);
      gamma_10_c[c] ~ normal(gamma_10,3.);
      gamma_01_c[c] ~ normal(gamma_01,3.);
      beta_c[c] ~ normal(beta,5.); 
      mu_10_c[c] ~ normal(mu_10,5.);
      mu_01_c[c] ~ normal(mu_01,5.);
      psi_c[c] ~ normal(psi,3.); 

      for(i in num_id[2*c-1]:num_id[2*c]){
          eta_p1_ic[i] ~ normal(eta_p1_c[c],3.);
          eta_q1_ic[i] ~ normal(eta_q1_c[c],3.);
          beta_ic[i] ~ normal(beta_c[c],5.);
          gamma_p1_ic[i] ~ normal(gamma_p1_c[c],3.);
          gamma_q1_ic[i] ~ normal(gamma_q1_c[c],3.);
          mu_10_ic[i] ~ normal(mu_10_c[c],5.);
          mu_01_ic[i] ~ normal(mu_01_c[c],5.);
      }
    }
  
    for(n in 1:N){
         target += binom_dist_lpmf(X[n] | Y[n], p[n],q[n]);
    }
}
'



require(rstan)

# load data file
load("empirical_dummy_updated.Rdata")

# initial values for two chains
I <- dat$I
C <- dat$cond_num

initial_string <- list(beta_ic=rep(-5,I),
                       mu_10_ic=rep(0,I),
                       mu_01_ic=rep(0,I),
                       eta_p1_ic=rep(0,I),
                       eta_q1_ic=rep(0,I),
                       gamma_p1_ic=rep(0,I),
                       gamma_q1_ic=rep(0,I),
                       eta_10_c=rep(0,C),
                       eta_01_c=rep(0,C),
                       gamma_10_c=rep(0,C),
                       gamma_01_c=rep(0,C),
                       psi_c=rep(0,C))

init <- list(initial_string,initial_string)

# seed
set.seed(11117)

# fit with two chains
fit <- stan(model_code=Model_breakpoint,data=dat,chains=2,iter=35000,warmup=5000,init=init)

# extract the parameters to save independently, to save some space

lp <- extract(fit,par=c("lp__"),permuted=F)

beta <- extract(fit,par=c("beta"),permuted=F)
mu_10 <- extract(fit,par=c("mu_10"),permuted=F)
mu_01 <- extract(fit,par=c("mu_01"),permuted=F)
eta_p1 <- extract(fit,par=c("eta_p1"),permuted=F)
eta_q1 <- extract(fit,par=c("eta_q1"),permuted=F)
gamma_p1 <- extract(fit,par=c("gamma_p1"),permuted=F)
gamma_q1 <- extract(fit,par=c("gamma_q1"),permuted=F)
eta_10 <- extract(fit,par=c("eta_10"),permuted=F)
eta_01 <- extract(fit,par=c("eta_01"),permuted=F)
psi <- extract(fit,par=c("psi"),permuted=F)
gamma_10 <- extract(fit,par=c("gamma_10"),permuted=F)
gamma_01 <- extract(fit,par=c("gamma_01"),permuted=F)

mu_c <- extract(fit,par=c("mu_c"),permuted=F)
mu_10_c <- extract(fit,par=c("mu_10_c"),permuted=F)
mu_01_c <- extract(fit,par=c("mu_01_c"),permuted=F)
eta_p1_c <- extract(fit,par=c("eta_p1_c"),permuted=F)
eta_q1_c <- extract(fit,par=c("eta_q1_c"),permuted=F)
gamma_p1_c <- extract(fit,par=c("gamma_p1_c"),permuted=F)
gamma_q1_c <- extract(fit,par=c("gamma_q1_c"),permuted=F)
eta_10_c <- extract(fit,par=c("eta_10_c"),permuted=F)
eta_01_c <- extract(fit,par=c("eta_01_c"),permuted=F)
psi_c <- extract(fit,par=c("psi_c"),permuted=F)
gamma_10_c <- extract(fit,par=c("gamma_10_c"),permuted=F)
gamma_01_c <- extract(fit,par=c("gamma_01_c"),permuted=F)

eta_p1_ic <- extract(fit,par=c("eta_p1_ic"),permuted=F)
eta_q1_ic <- extract(fit,par=c("eta_q1_ic"),permuted=F)
gamma_p1_ic <- extract(fit,par=c("gamma_p1_ic"),permuted=F)
gamma_q1_ic <- extract(fit,par=c("gamma_q1_ic"),permuted=F)
eta_10_ic <- extract(fit,par=c("eta_10_ic"),permuted=F)
eta_01_ic <- extract(fit,par=c("eta_01_ic"),permuted=F)
psi_ic <- extract(fit,par=c("psi_ic"),permuted=F)
gamma_10_ic <- extract(fit,par=c("gamma_10_ic"),permuted=F)
gamma_01_ic <- extract(fit,par=c("gamma_01_ic"),permuted=F)
mu_ic <- extract(fit,par=c("mu_ic"),permuted=F)
mu_10_ic <- extract(fit,par=c("mu_10_ic"),permuted=F)
mu_01_ic <- extract(fit,par=c("mu_01_ic"),permuted=F)


fit <- list(lp=lp,mu=mu,mu_10=mu_10,mu_01=mu_01,mu_c=mu_c,
            mu_10_c=mu_10_c,mu_ic=mu_ic,
            mu_01_c=mu_01_c,mu_10_ic=mu_10_ic,mu_01_ic=mu_01_ic,
            gamma_p1=gamma_p1,gamma_q1=gamma_q1,gamma_10=gamma_10,gamma_01=gamma_01,
            gamma_p1_c=gamma_p1_c,gamma_q1_c=gamma_q1_c,gamma_10_c=gamma_10_c,
            gamma_01_c=gamma_01_c,gamma_p1_ic=gamma_p1_ic,gamma_q1_ic=gamma_q1_ic,
            gamma_10_ic=gamma_10_ic,gamma_01_ic=gamma_01_ic,
            eta_p1=eta_p1,eta_q1=eta_q1,eta_10=eta_10,eta_01=eta_01,
            eta_p1_c=eta_p1_c,eta_q1_c=eta_q1_c,eta_10_c=eta_10_c,
            eta_01_c=eta_01_c,
            eta_p1_ic=eta_p1_ic,eta_q1_ic=eta_q1_ic,eta_10_ic=eta_10_ic,eta_01_ic=eta_01_ic,
            psi=psi,psi_c=psi_c,psi_ic=psi_ic)

# save results to Rdata file

save(fit,file="empirical_breakpoint_result_updated.Rdata")





