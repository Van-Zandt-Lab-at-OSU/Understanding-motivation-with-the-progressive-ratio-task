# add a wald distribution to model RT


Model_RT <- '
functions {
// self-defined functions for the inverse-Gaussian distribution
  real IG_lpdf(real t, real v,real b) {
      real density1;
      real density2;
      real density;
      density1 = log(b) - 0.5*(log(2)+log(pi())+3*log(t));
      density2 = -0.5*(v*t-b)^2/t;
      density = density1 + density2;
      return(density);
  }
  real IG_cdf(real t,real v,real b) {
      real cdf1;
      real cdf2;
      real cdf;
      cdf1 = normal_cdf((v*t-b)/sqrt(t),0,1);
      cdf2 = exp(2*v*b)*normal_cdf((-v*t-b)/sqrt(t),0,1);
      cdf = cdf1 + cdf2;
      return(cdf);
  }
  real minimum_lpdf(real t, real b,real a{ 
    real density;
    real d;
    d = 1-IG_cdf(t,b,a);
    density = log(2) + IG_lpdf(t | b,a) + log(d);
    return(density);
  }
}
data {
   int<lower=0> I;                   // number of subjects
   int<lower=0> cond_num;            // number of conditions
   int<lower=0> num_id[2*cond_num];  // identify the start and end of each group: 1-37, 38-77, 78-122
   int<lower=0> N_RT;                // number of RTs
   real<lower=0> RT[N_RT];           // response time in seconds
   int rtlabel[N_RT];                // identify the subject for each response time
}
parameters {
// parameters for response times
   real<lower=0> a;
   real<lower=0> b;
   real<lower=0> a_c[cond_num];
   real<lower=0> b_c[cond_num];
   real<lower=0> a_ic[I];
   real<lower=0> b_ic[I];
}
model {  

   a ~ gamma(1.,1.);
   b ~ gamma(1.,1.);
   
   for(c in 1:cond_num){

      a_c[c] ~ gamma(1.,1/a);
      b_c[c] ~ gamma(1.,1/b);

      for(i in num_id[2*c-1]:num_id[2*c]){

          a_ic[i] ~ gamma(1.,1/a_c[c]);
          b_ic[i] ~ gamma(1.,1/b_c[c]);
      }
    }

    for(n in 1:N_RT){
          target += minimum_lpdf(RT[n] | b_ic[rtlabel[n]],
                                        a_ic[rtlabel[n]]);
    }
}
'



require(rstan)

# load data file
load("empirical_dummy_updated.Rdata")

# initial value
I <- dat$I
C <- dat$cond_num

init <- list(list(a_ic=rep(1,I),
                  b_ic=rep(1,I)),
             list(a_ic=rep(1,I),
                  b_ic=rep(1,I)))

# seed
set.seed(11117)

# fit model
fit <- stan(model_code=Model_RT,data=dat,chains=2,iter=6000,warmup=1000,init=init)



# extract relevant values and save independently to save some space

lp <- extract(fit,par=c("lp__"),permuted=F)

a <- extract(fit,par=c("a"),permuted=F)
b <- extract(fit,par=c("b"),permuted=F)

a_c <- extract(fit,par=c("a_c"),permuted=F)
b_c <- extract(fit,par=c("b_c"),permuted=F)

a_ic <- extract(fit,par=c("a_ic"),permuted=F)
b_ic <- extract(fit,par=c("b_ic"),permuted=F)

fit <- list(lp=lp,a=a,b=b,
            a_c=a_c,b_c=b_c,
            a_ic=a_ic,b_ic=b_ic)

save(fit,file="empirical_rt_result.Rdata")

