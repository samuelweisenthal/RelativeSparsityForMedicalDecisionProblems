
# a stable expit
expit = function(x){
    # k goes away, but stabilizes
    k  = max(0,x)#(0>x)*0+(x>=0)*x#max(0,x)
    as.vector(exp(x-k)/(exp(0-k)+exp(x-k)))
}

# penalty
pen = function(beta,b,betahat,bhat,sigmas){
  sum(abs(beta-b))
}

# return (here, just reward)
G=function(Rs,u){
  
  g=Rs
  
  g
}

# pi(A=1|S=s). eta is old argument, not used
ps1 = function(theta,eta,s){
  # issue with categorical cov. maybe after scaling?
  expit(t(theta)%*%(s))
}

# pi(A=a|S=s)
psa = function(theta,eta,s,a){
  if(is.na(matrix(s)[1,])){
    1 # for censoring. not used now.
    }else{
  ps1(theta,eta,s)^a*(1-ps1(theta,eta,s))^(1-a)
  }
}

# wrapper for pi over time if needed
kappa = function(theta,eta,Ss,As){

  k=slog(psa(theta,eta,Ss,As))
  exp(k)
}

logkappa = function(theta,eta,Ss,As){
  slog(kappa(theta,eta,Ss,As))
}

# importance anonymouspling for one observation
vi = function(beta,b,Ss,As,Rs,u){

  log.kappab = slog(psa(b,b,Ss,As))
  log.kappabeta = +slog(psa(beta,b,Ss,As))
  G=G(Rs,u)
  
  exp(log.kappabeta)/exp(log.kappab)*G
}

# sum of vis
vn = function(beta,b,eps,u){
  n=length(eps)
  v=0
  for (i in 1:n){
    v = v + vi(beta,b,eps[[i]]$Ss,eps[[i]]$As,eps[[i]]$Rs,u)
  }
  v/n
}

# variance of vn
sigma.vn = function(beta,b,eps,u){
  # try w mn
  vvn = vn(beta,b,eps,u)
  n=length(eps)
  sigma2 = 0
  for (i in 1:n){
    sigma2 = sigma2 + (vi(beta,b,eps[[i]]$Ss,eps[[i]]$As,eps[[i]]$Rs,u) - vvn)^2
  }
  # it's only n, because we want variance of sqrt(n)V_n, not of V_n
  sqrt(sigma2/(n))
}

# old function
mn = function(beta,b,eps,u,sigmas){

   vn(beta,b,eps,u) 
}

# jn
wn = function(beta,b,eps,u,betahat,bhat,lambda,sigmas){
  n = length(eps) #can rem
  mn(beta,b,eps,u,sigmas) - lambda*pen(beta,b,betahat,bhat,sigmas)
}



