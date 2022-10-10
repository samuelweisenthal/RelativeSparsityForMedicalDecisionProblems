# simulates data. used in MC experiments

library(mvtnorm)
source('utils.R')

gen.data.3arg = function(b0,n,init.state.mean,init.state.var,trans.var,tau,R,seed){

  set.seed(seed)
  K = dim(init.state.mean)[1]
  eps = list()
  for (i in 1:n){
    
    sp = t(mvtnorm::rmvnorm(n=1,mean=init.state.mean,
                            sigma=init.state.var))
    As = list()
    Ss = list()
    Rs = list()
    

      s =  sp
      Ss = s
      
      a = rbinom(1,1,prob=expit(t(b0)%*%s))
     
      
      sp = t(mvtnorm::rmvnorm(n=1,mean=s+matrix(tau*s*rep(a,dim(s)[1])),
                              sigma=trans.var))
      As = a
      Rs = R(s,a,sp)

    eps[[i]] = list(Ss=Ss,As=As,Rs=Rs)
  }
  eps
  
}
