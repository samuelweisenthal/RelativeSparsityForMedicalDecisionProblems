source('center.scale.R')
# general helpers.
# some not used for experiments in main text

# log with epsilon
slog = function(x){log(x+1e-50)}

# make a scalar a vector
tile = function(r,dim){
  matrix(r,nrow=dim[1],ncol=dim[2],byrow=TRUE)
}

# another log with epsilon
# eps.log = function(x,EPS=1e-15){
#   log(x+EPS)
# }

# reshape elements of list
unlist.reshape=function(theta.stars){
  theta.stars.m = matrix(unlist(theta.stars),
                         nrow=length(theta.stars),
                         ncol=length(theta.stars[[1]]),
                         byrow=TRUE)
  theta.stars.m
}

# reshape elements of list
unlist.reshape.eps.ss.sumTsxK = function(ss){

  n=length(ss)
  if (is.null(ss[[1]])){
    print("Null ss")
    browser()
  }
  K=dim(ss[[1]])[1]
  matrix(unlist(ss),nrow=n,ncol=K,byrow=TRUE)
}

get.sparse.combn = function(cont){
  # each variable plotted only 1 to 2 times
  # for A,B,C,D,E
  # gives just (A,B),(B,C),(C,D),(D,E),(A,E)
  combs = t(combn(c(cont),2))
  combs.red=rbind(combs[!duplicated(combs[,1]),],c(cont[1],cont[length(cont)]))
  split(combs.red, seq(nrow(combs.red)))
}

# convert vectors/dataframes of episodes in real data to list, processed by modeling functions
get.eps.ms = function(s,a,r){
  
  n=dim(s)[1]
  a = as.matrix(a)
  r = as.matrix(r)
  eps=list()
  
  for (i in 1:n){
      Ss=matrix(as.vector(s[i,])) 
      As=a[i,1]
      Rs=r[i,1] 
    eps[[i]] = list(Ss=Ss,As=As,Rs=Rs)
  }
  eps
}




