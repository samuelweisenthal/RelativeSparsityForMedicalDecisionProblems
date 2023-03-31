# various functions to run and analyze Monte-Carlo results
#source('liops.R')
library(latex2exp)
library(glmnet)
library(predtools)
library(DynTxRegime)

# processes episodes - ie computes betanlambda, jn, var vn, etc
process.eps = function(eps,b0,
                       u,
                       lambda,plt.obj,scale.s,
                       split,other.args=NULL){



  use.diff = other.args$use.diff
  
  n = length(eps)

  K = dim(eps[[1]]$Ss)[1]
  
  # maybe make as part of mimic.R
  check.p = 0
  if (check.p){
    check.pos(eps,cov.of.int=1:K,pen.b=1,train.test=0)
    browser()
  }
  
  hf = as.integer(n/2)
  n.test = NA#length(eps.test)
  n.train = NA#length(eps)
  
  # split data into test and training
  if(split){
    # need to deal with this
    eps.test = eps[(hf+1):n]
    n.test = length(eps.test)
    eps = eps[1:hf]
    n.train = length(eps)
  }else{
    eps.test=eps
  }

  if(scale.s){
    so = center.scale.s(eps)
    sigmas = rep(1,K) # we don t do scale factor in penalty bc we scale state var
    so.test = center.scale.s(eps.test,means=so$means,sds=so$sds)
    
    eps.test = so.test$e
    
  }else{
    # don t center and scale
    so = center.scale.s(eps,means=rep(0,K),sds = rep(1,K)) #not centering and scaling
    # but we still need the sds for the scaling factors.
    soc = center.scale.s(eps)
    sigmas = soc$sds
    # generally, we would do the above with scale factors. instead
    # override scale factors by setting sigmas to 1. we are not scaling for simulations. 
    # we just generate data with variance 1
    #we just scale for real data
    sigmas = rep(1,K)
    
    ### eps.test no longer defined. I never don't scale.
    # if do scale, need to fix
    print("Need to define eps.test!")
    browser()
    
  }
  # a not affected by scaling
  aa = unlist(lapply(eps,function(x){x$As}))
  sm = so$sm.cs
  eps=so$e
  rs = unlist(lapply(eps,function(x){x$Rs}))
 
  
  
  ######


  # epsq = so$e.raw
  # browser()
  # quantmax = 5
  # smq = so$sm
  # #smq[,1]=as.factor(round(sm[,1],0))#cut(sm[,1],c(-Inf,seq(-quantmax,quantmax,1),Inf))
  # #smq[,2]=as.factor(round(sm[,2],0))#cut(sm[,2],c(-Inf,seq(-quantmax,quantmax,1),Inf))
  # mq = matrix(c(1:100),nrow=10,ncol=10)
  # sq = matrix(NA,nrow=dim(smq)[1],ncol=1)
  # spq =  matrix(NA,nrow=dim(smq)[1],ncol=1)
  # for (i in 1:dim(smq)[1]){
  #   print(i)
  #  sq[i,]= mq[round(epsq[[i]]$Ss[1,],0),round(epsq[[i]]$Ss[2,],0)]   
  #  spq[i,] = mq[round(epsq[[i]]$Sp[1,],0),round(epsq[[i]]$Sp[2,],0)]   
  # }
  # browser()

  
  # fit nuisance
  # if split, eps is already of length hf, so don't really need this

  #if (split){
  #  bm=glm(aa[1:hf]~sm[1:hf,]-1,family='binomial')
    
  #}else{
  
  # this will be overwritten, but you can set value to be this and penalty to be
  # cv.glmnet below. However decoupling not easy. best for absolute sparsity, 
  # if decide to use it

  #}
  

  #bn = b.optim$par
  pen.b=other.args$pen.b
    if (pen.b){

      bmp = cv.glmnet(sm,aa,intercept = FALSE,alpha=0,family='binomial',
                      type.measure=other.args$type.measure,
                      lambda=other.args$b.lambdas)
      bmp.coef=coef(bmp,s = "lambda.min")
      # use penalized behavioral policy, for both penalty and value
      bnv = bmp.coef[2:(K+1),]
    }else{
      bm=glm(aa~sm-1,family='binomial') 
      bnv=bm$coefficients
    }
  
  if (other.args$absolute.spar.in.pen){
    # penalty is absolute
    bnp =  rep(0,K)
  }else{
    # set penalty to same as v
    bnp = bnv
  }

  
  #
    

  # to debug
  dont.est.beh=0
  if (dont.est.beh){
    b=b0
  }

  # set it again, since it is also used in variance calc
  dont.est.beh=0 
  
  #Nelder-Mead gives convergence 10 degeneracy sometimes with K=5
  #if (K>1){
  #method='Nelder-Mead' # maybe write gradient descent using l1 approx
  #}else{
  method='BFGS' #must use BFGS for one dimensional case
  #}
  
  maxit=1e6
  
  
  
  
  ######
  use.bowl=0
  if(use.bowl){
    fit.bowl = function(sm,rs,aa){
      
      sm.df = as.data.frame(sm)
      var.names = c(colnames(sm.df),-1)
      #var.names = c(colnames(sm.df))
      form = as.formula(paste("~", paste(var.names, collapse="+")))
      sm.df$R= rs
      sm.df$A = aa
      moPropen <- buildModelObj(model = form,
                                solver.method = 'glm',
                                solver.args = list('family'='binomial'),
                                predict.method = 'predict.glm',
                                predict.args = list(type='response'))
      
      fitSS <- bowl(moPropen = moPropen,
                    data = sm.df, reward = sm.df$R,  txName = 'A', 
                    regime = form,kernel="linear",lambdas=c(1))
      list(fitSS=fitSS,
           value.bowl=estimator(fitSS),
           regime=optTx(fitSS))
    }
    
    my.bowl = fit.bowl(so$sm,rs,aa)
    value.oracle.train = my.bowl$value.bowl
    
    
    #aa.test = unlist(lapply(eps.test,function(x){x$As}))
    #sm.test = so.test$sm.cs
    #eps.test=so.test$e
    #rs.test = unlist(lapply(eps.test,function(x){x$Rs}))
    
    
    
    #optTx(fitSS, sm.df)
    value.oracle.train = estimator(my.bowl$fitSS)
    mean.pred.oracle.train = mean(my.bowl$regime$decisionFunc)
  }else{
    # we know oracle is always treat
    oracle.beta = rep(1,K)
    oracle.beta[2] = 1e6
      value.oracle.train = vn(oracle.beta,bnv,eps,u)
    mean.pred.oracle.train = 1
  }
  
  
  
  
  #init at behavioral to minimize positivity viol at beginning

  # vn
  optim.res = optim(par=c(bnv),mn,
                    control=list(fnscale=-1,maxit=maxit),
                    method=method,
                    eps=eps,bv=bnv,
                    u=u,sigmas=sigmas
                    )
  

  bn.optim=NA
  betan=matrix(optim.res$par)
  
  # take min because sometimes inf
  rollout.betangamma = min(expit(sm%*%betan),1)


  rollout.bn = expit(sm%*%bnv)

  roll.diff.betangamma.bn.abs = sum(abs(rollout.betangamma-rollout.bn))/(n)
  vn.betagamma = vn(betan,bnv,eps,u)
  
  # jn
  optim.res.lamb = optim(par=c(bnv),wn,
                         control=list(fnscale=-1,maxit=maxit),
                         method=method,
                         bv=bnv,bp=bnp,eps=eps,
                         u=u,betahat=betan,
                         bhat=bnp,
                         lambda=lambda,sigmas=sigmas
                         )
  
  # convergence check
  print(c("converge?",optim.res.lamb$convergence))
  betanlambda=matrix(optim.res.lamb$par)
  
  # see what happens when use policy prospectively
  rollout.betangammalambda = expit(sm%*%betanlambda)
  roll.diff.betangammalambda.bn.ab = sum(abs(rollout.betangammalambda-rollout.bn))/(n)
  # could rep vn with mn here?
 
  # calc vn
  vn.betagammalambda=vn(betanlambda,bnv,eps,u)
  vn.betagammalambda.test = vn(betanlambda,bnv,eps.test,u)
  vn.bn = vn(bnp,bnv,eps,u)
  #calc jn
  wn.betagammalambda = wn(betanlambda,  bv=bnv,bp=bnp,eps,
                          u,betan,
                          bhat=bnp,
                          lambda,sigmas
                          )
  wn.betagammalambda.test = wn(betanlambda, bv=bnv,bp=bnp,eps.test,
                               u,betan,bhat=bnp,
                               lambda,sigmas
                               )
  
  tilde.optim.res = NA
  tilde.betan=NA
 
    var = matrix(NA)
    var.betanlambda= matrix(NA)
    var.b= matrix(NA)
  
  
  sigma.vn. = sigma.vn(betanlambda,bnv,eps,u)
  sigma.vn.test = sigma.vn(betanlambda,bnv,eps.test,u)
  diff = betanlambda-bnv
  diff.bn.v.betangamma = betan - bnv

  lamsel = abs(diff)>use.diff
  lp.rollout.betanlambda.sel = sm%*%(betanlambda*lamsel)
  lp.rollout.betanlambda.nosel = sm%*%(betanlambda*(1-lamsel))
  rollout.lps = expit(lp.rollout.betanlambda.sel+lp.rollout.betanlambda.nosel)[1:10]
  rollout.betangammalambda[1:10]
  comp.lp.sel.nosel = cbind(lp.rollout.betanlambda.sel,lp.rollout.betanlambda.nosel)


  list(var=var,var.betanlambda=var.betanlambda,
       var.b=var.b,
       betan=betan,
       betanlambda=betanlambda,
       bnv=bnv,
       bnp=bnp,
       bn.optim=bn.optim,
       tilde.betan=tilde.betan,
       vn.bn=vn.bn,
       vn.betagamma=vn.betagamma,
       vn.betagammalambda=vn.betagammalambda,
       vn.betagammalambda.test=vn.betagammalambda.test,
       wn.betagammalambda=wn.betagammalambda,
       wn.betagammalambda.test=wn.betagammalambda.test,
       sigam.vn. = sigma.vn.,
       sigam.vn.test = sigma.vn.test,
       split=split,
       n.train=n.train,
       n.test=n.test,
       rollout.betangamma =rollout.betangamma[1:100],
       rollout.betangammalambda = rollout.betangammalambda[1:100],
       rollout.betangamma.mean = mean(rollout.betangamma),
       rollout.betangammalambda.mean = mean(rollout.betangammalambda),
       rollout.lps=rollout.lps,
       rollout.bnmean = mean(rollout.bn),
       lambda=lambda,
       roll.diff.betangammalambda.bn.ab=roll.diff.betangammalambda.bn.ab,
       roll.diff.betangamma.bn.abs=roll.diff.betangamma.bn.abs,
       diff.bn.v.betangammalambda=diff,
       diff.bn.v.betangamma=diff.bn.v.betangamma,
       use.diff = use.diff,
       lamsel=lamsel*1,
       comp.lp.sel.nosel=head(comp.lp.sel.nosel,n = 10L),
       value.oracle.train = value.oracle.train,
       mean.pred.oracle.train = mean.pred.oracle.train
       )
}

lambda.path = function(eps.mimic,b0,
                       lambdas,names,
                       file,
                       split,scale.s,lambda.select.ix=FALSE,
                       gamma.sel=1,use.diff=1,other.args=NULL){
  
  # plots lambda vs coefficients
  
  n=length(eps.mimic) 
  K = dim(eps.mimic[[1]]$Ss)[1]
  ti=Sys.time()
  outer.res=list()
    dres = list()
    res = list()
    
    for (i in 1:length(lambdas)){
      print(c("lam",lambdas[i]))
      if ( (i==1)){plt.obj=1}else{plt.obj=0}
  
      res[[i]]=process.eps(eps.mimic,b0=b0,
                           u=1,
                           lambda=lambdas[i],plt.obj=0,scale.s=scale.s,
                           split=split,other.args=other.args)
    }

    dres = res
    # set this again; was because was previously loop here
    outer.res = dres
    
  te = Sys.time()-ti
  print(te)

  op = list(outer.res=outer.res,
            lambdas=lambdas,n=length(eps.mimic), 
            names=names,
            split=split,lambda.select.ix=lambda.select.ix,
            gamma.sel=gamma.sel,use.diff=use.diff,other.args=other.args)
  saveRDS(op,paste0(file,"op"))
  op
}

# get.possible.selects = function(K,nDiverg){
#   vecs = list()
#   numbernonzero = 1
#   combs=combn(K,nDiverg)
#   for (c in 1:dim(combs)[2]){
#     vec = rep(0,K) 
#     for (i in 1:dim(combs)[1]){
#       vec[combs[i,c]] = 1
#     }
#     vecs[[c]] = vec
#   }
#   vecs 
# }

plot.lambda.paths = function(outer.res,
                             max.n.div=1,plot.=TRUE,lambda.ns=NULL){
  # plots lambda against coefficients
  other.args = outer.res$other.args
  var.exists = !is.null(outer.res$or.shell.var)
  if (!var.exists){
    plot.emp.var = 0
  }else{
    outer.res.emp.var = outer.res$or.shell.var$outer.res
    plot.emp.var = other.args$plot.emp.var
    
  }
  #plot.adap.lasso.ci=other.args$plot.adap.lasso.ci
  #delta.sel = other.args$delta.sel
  #par(mfrow=c(2,1))
  split=outer.res$split
  gamma.sel=outer.res$gamma.sel
  print(split)
  # take c because sometimes it s a matrix if eg doing MC and averaging
  n.train = c(outer.res$outer.res[[1]]$n.train)
  n.test = c(outer.res$outer.res[[1]]$n.test)
  names = outer.res$names
  n = outer.res$n
  use.diff = outer.res$use.diff
  
  round.dig = 5
  if (split){
    n.train=n.train
    n.test=n.test
  }else{
    n.train=n
    n.test=n
  }
  #deltas = outer.res$deltas
  print(n)
  lambda.select.ix = outer.res$lambda.select.ix
  lambdas=outer.res$lambdas
  #gammas=outer.res$gammas
  outer.res= outer.res$outer.res
  nlam=length(lambdas)
  lambda.selects = list()
  jj=1
  #png("plts.png",width=500*8,height=500*2)
  value.ymin = NULL
  value.ymax = NULL
  #for (d in 1:length(deltas)){
  #delta = deltas[d]
  dres = outer.res#
  if (var.exists){
    dres.emp.var = outer.res.emp.var#
  }
  #for (j in 1:length(gammas)){
  pick.lam=1
  #gamma = gammas[j]
  res = dres#[[j]]
  betanlambdas.l=lapply(res,function(x){x[['betanlambda']]})
  betanlambdas=unlist.reshape(betanlambdas.l)
  #bs.l = lapply(res)
  vars = unlist.reshape(lapply(res, function(x){diag(x$var.betanlambda)}))
  
  vars.opt = unlist.reshape(lapply(res, function(x){diag(x$var.betan)}))
  
  
  #eps = gen.data(b=b0,n,init.state.mean,init.state.var,trans.var,T)
  #process.eps(eps,b0=NULL,gamma=1,u=1,delta=1)
  
  
  #must be vector for broadcasting later
  beh=c(res[[1]]$bnv)
  K = length(beh)
  all=c(beh,betanlambdas)
  
  
  zse = qnorm(0.975)*sqrt(vars/n.train)
  
  if(var.exists){
    
    res.emp.var = dres.emp.var
    
    M = length(res.emp.var)
    # so we take M*montecarloVariance(betan) \approx sigma^2
    
    #vars.emp =  M*unlist.reshape(lapply(res.emp.var, function(x){x$betanlambda}))
    # we then construct ci for beta0 as betan +- sqrt(sigma^2/n)*z
    #zse.emp = qnorm(0.975)*sqrt(vars.emp/n.train)
    
    #browser()
    
    vars.emp =  unlist.reshape(lapply(res.emp.var, function(x){x$betanlambda}))
    # we then construct ci for beta0 as betan +- sqrt(sigma^2/n)*z
    zse.emp = qnorm(0.975)*sqrt(vars.emp)
    
 
    print(zse.emp)
  }
  
  use.ci.for.ylim = 0
  if (use.ci.for.ylim){
    #if (plot.adap.lasso.ci){
    for (i in K:1){
      all = c(all,betanlambdas[,i]-zse[,i],betanlambdas[,i]+zse[,i])
    }
    #}
  }
  
  logx = other.args$logx
  if(logx){
    logxarg = 'x'
  }else{
    logxarg = ''
  }
  
  if (plot.){

    
    
    #if (other.args$DeltaOnlyTitle){  
    main=(bquote("n (train) =" ~ .(n.train) ~ "; " ~"n (test) =" ~ .(n.test) ~ "; " ~ Delta == .(use.diff))) 
    #}else{
    #  main=(bquote("n (train) =" ~ .(n.train) ~ "; " ~"n (test) =" ~ .(n.test) ~ "; " ~ ~ gamma == .(formatC(gamma, format = "e", digits = 2)) ~ '; ' ~ delta == .(delta)~ '; ' ~ Delta == .(use.diff)))
    #}
    
    # if we are not computing var, we don't have gamma
    #if (other.args$get.theor.var.of.ests==0){
    
    ylablambdas = TeX("Train $\\beta_{n,\\lambda}$ (dashed) or $b_{n}$ (solid)")
    #}else{  
    #ylablambdas = TeX("Train $\\beta_{n,\\gamma,\\lambda}$ (dashed) or $b_{n}$ (solid)")
    #}
    
    
    #############
    # bottom, left, top, right
    par(mar=c(1,5.1,2,2))
    plot(lambdas,rep(0,nlam),ylim=c(min(all),max(all)),lty=0,
         type='l',
         #main=(bquote("n=" ~ .(n.train) ~ "; " ~ gamma == .(formatC(gamma, format = "e", digits = 2)) ~ '; ' ~ delta == .(delta)~ '; ' ~ Delta == .(use.diff))),
         main=main,
         xlab=TeX("$\\lambda$"),
         ylab=ylablambdas,log=logxarg,xaxt='n', cex.lab=1.5,cex.axis=1.5,,cex.main=1.5,cex=1.5)
    
    
    colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                           "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                           c("#000000","#004949","#009292","#ff6db6","#ffb6db",
                             "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
                             "#920000","#924900","#db6d00","#24ff24","#ffff6d")[c(12,14)])
    
    
    if (plot.emp.var){
      for (i in K:1){
        polygon(c(rev(lambdas),lambdas),c(rev(betanlambdas[,i]-zse.emp[,i]),
                                          betanlambdas[,i]+zse.emp[,i]),
                col=adjustcolor(colorBlindBlack8[i],alpha.f=0.1),#rgb(0,0,0,alpha=0.1),
                lty=0)
      }
      
    }
    
    for (k in K:1){
      lines(lambdas,betanlambdas[,k],lty=2,col=colorBlindBlack8[k],lwd=2)
      abline(h=beh[k],lty=1,col=colorBlindBlack8[k],lwd=2)
    }
    
    
    legend("topright",TeX(names),lty=1,col=colorBlindBlack8[1:K],lwd=rep(2,K),cex=1.5)
    #}
    
    
    
    
  }
  
  
  
  vs = lapply(res,function(x){x$vn.betagammalambda})
  ws = lapply(res,function(x){x$wn.betagammalambda})
  vs = unlist(vs)
  sigma.vs = lapply(res,function(x){x$sigam.vn.})
  sigma.vs = unlist(sigma.vs)
  zse.v = qnorm(0.975)*sigma.vs/sqrt(n.train)
  
  
  ####### Lambda Sel ########
  
  # if we set lambda select, use it, but overridden below
  if(lambda.select.ix){
    star.ix=lambda.select.ix
  }
  
  
  if (is.null(other.args$num.se)){
    
    vmin = other.args$vmin
  }else{
    vmin = c(res[[length(res)]]$vn.betagammalambda + other.args$num.se*sqrt(res[[length(res)]]$sigam.vn./n.train))
  }
  
  
  
  select.lambda = function(betanlambdas,beh,zse,use.diff,vs,vmin){
    # instead let's find lambda that satisfies our criterion
    # of having a certain number diverging, but being the smallest lambda that does that
    if (use.diff){
      # diff big
      diff.from.behav = abs(t(betanlambdas)-beh)>use.diff
    }else{
      # diff big same as beh not in ci
      diff.from.behav = !(t(betanlambdas-zse)<=beh & t(betanlambdas+zse)>=beh)
    }
    select.lambda.adaptive = 1
    no.lambda.satisfies = 0
    if (select.lambda.adaptive){
      n.diff = apply(diff.from.behav,2,sum)
      # could possibly only take ones with positive here, so remove 000 selection
      # but in some ways 000 selection preferable "first do no harm
      
      dontLetEqualZero=0
      if (dontLetEqualZero){
        n.diff[n.diff==0] = 1e6
      }
      abs.diff.diff = abs(n.diff-max.n.div)
      
      #abs.min = which(abs.diff.diff==min(abs.diff.diff))
      
      
      gtvim  = vs>vmin
      # if it is not greater than minimum value, set it to inf
      abs.diff.diff[gtvim==FALSE] = Inf
      
      # take the min now of the ones, to get closest to target
      abs.min.gt.vmin = which((abs.diff.diff==min(abs.diff.diff)))
      
      # now take the max to get the most sparsity? 
      # or min to get most value. I think because
      # now we are over a value threshold, best to target sparsity
      # before, was no value threshold, so was taking largest value (min)
      # now we max. No now we still min
      
      # before, take with max value
      #star.ix = min(abs.min.gt.vmin)
      star.ix = max(abs.min.gt.vmin)
      
      exact = 0
      if (exact){
        lam.sat.crit = (n.diff == max.n.div)
        if (sum(lam.sat.crit)==0){
          # if this is the case, no lambda satisfies our criterion. so set to 1
          # then i guess we just don't use these parameter choices ..
          # need to check for this...
          star.ix = length(lambdas)
          no.lambda.satisfies = 1
        }else{
          star.ix = min(which(lam.sat.crit==TRUE))
        }
      }
    }
    selected.mask=(diff.from.behav*1)[,star.ix]
    list(star.ix=star.ix,no.lambda.satisfies=no.lambda.satisfies,selected.mask=selected.mask)
  }
  
  
  sel.res = select.lambda(betanlambdas,beh,zse,use.diff,vs,vmin)
  star.ix = sel.res$star.ix
  no.lambda.satisfies = sel.res$no.lambda.satisfies
  selected.mask = sel.res$selected.mask
  
  
  ####### Lambda Sel ########
  if (plot.){
    abline(v=lambdas[star.ix],lty=2)
  }
  
  
  
  
  ## test
  
  vst = lapply(res,function(x){x$vn.betagammalambda.test})
  wst =  lapply(res,function(x){x$wn.betagammalambda.test})
  vst = unlist(vst)
  sigma.vst = lapply(res,function(x){x$sigam.vn.test})
  sigma.vst = unlist(sigma.vst)
  zse.vt = qnorm(0.975)*sigma.vst/sqrt(n.test)
  
  
  #plt.train=1
  #if (plt.train){
  plot.w = 0
  if (plot.w){
    ylab = TeX('$V_{n}$ or $W_{n}$')
  }else{
    ylab = TeX('$V_{n}$')
  }
  use.first.val.for.ylim = 0
  if (use.first.val.for.ylim){
    #if ( d ==1){
    # take the value bounds from the first, and use them for every plot
    value.ymax = max(c(vs+zse.v,vst+zse.vt,vs))
    value.ymin = min(c(vs-zse.v,vst-zse.vt,vs))
    
    #}
  }else{
    value.ymax = max(c(vs+zse.v,vst+zse.vt,vs))
    value.ymin = min(c(vs-zse.v,vst-zse.vt,vs))
    
  }
  
  if (plot.){
    
    if (var.exists){
      plotx = 'n'
    }else{
      plotx = NULL
    }
    
    if (!var.exists){
      # bottom, left, top, right
      par(mar=c(4.4,5.1,1,2))
    }else{
      # bottom, left, top, right
      par(mar=c(1,5.1,.2,2))
    }
    
    plot(lambdas,vs,ylab=ylab,xlab=TeX("$\\lambda$"),type='l',lwd=3,
         ylim=c(value.ymin,value.ymax),
         log=logxarg,
         xaxt=plotx, cex.lab=1.5,cex.axis=1.5,,cex.main=1.5,cex=1.5
         #main=bquote("n=" ~ .(n.test))
         #main=(bquote("n=" ~ .(n.train) ~ "; " ~ gamma == .(formatC(gamma, format = "e", digits = 2)) ~ '; ' ~ delta == .(delta)~ '; ' ~ Delta == .(use.diff))),
    )
    if (var.exists){
      
      # rug, not used
      #rug(lambdas[lambda.ns],ticksize=.1,lty=2,col=1,lwd=1)
      #print(c("lambdans",lambda.ns))
      #print(table(lambdas[lambda.ns]))
      
      #for (lambda.n in lambda.ns){
      #abline(v=lambdas[lambda.n],lty=2,col=rgb(0,0,0,0.01))
      #}
      
    }
    if (plot.w){
      lines(lambdas,ws,lty=2,col=3)
    }
    #lines(lambdas,vs-zse.v,lty=2)
    polygon(c(rev(lambdas),lambdas),c(rev(vs-zse.v),
                                      vs+zse.v),
            col=adjustcolor(1,alpha.f=0.1),#rgb(0,0,0,alpha=0.1),
            lty=0)
    
    
    
    #lower.bound = vs-zse.v
    #star.ix = which(lower.bound==max(lower.bound))
    #points(lambdas[star.ix],lower.bound[star.ix],pch=8)
    
    #star.ix.v = which(vs==max(vs))
    #points(lambdas[star.ix],vs[star.ix],pch=8)
    
    
    #scale.plot(lambdas,sigma.vs,(vs[1]-vs),"train",plott=FALSE) #function(ls,x,y,label){
    
    lines(lambdas,vst,lty=3,col=2,lwd=3)#,ylab=TeX('$V_{n}(test)$'),xlab=TeX("$\\lambda$"),type='l',
    #ylim=c(min(c(vs-zse.v,vs)),max(c(vs+zse.v,vs))))
    if (plot.w){
      lines(lambdas,wst,lty=4,col=4)
    }
    #if ((d==deltaixleg)){
    
    #}
    #lines(lambdas,vs-zse.v,lty=2)
    polygon(c(rev(lambdas),lambdas),c(rev(vst-zse.vt),
                                      vst+zse.vt),
            col=adjustcolor(2,alpha.f=0.1),#rgb(0,0,0,alpha=0.1),
            lty=0)
    
    lepskii=0
    if(lepskii){
      lower.bound = vst-zse.vt
      star.ix = which(lower.bound==max(lower.bound))
      #points(lambdas[star.ix],lower.bound[star.ix],pch=8)
      #text(lambdas[star.ix],lower.bound[star.ix],"V-2sigma/sqrt(n)")
      
      star.ix.v = which(vst==max(vst))
    }
    
    if (plot.w){
      legend("topright",TeX(c("Train $V_n$","Test $V_n$","Train $W_n$", "Test $W_n$")),
             lty=c(1,3,2,4),col=c(1,2,3,4),lwd=c(3,3,1,1),cex=1.5)
    }else{
      
      if (var.exists){
        lamlab = TeX("$\\bar{\\lambda}_n$")
      }else{
        lamlab = TeX("$\\lambda_n$")
      }
      legend("topright",c("Train","Test",
                          lamlab
                          #TeX("$\\lambda_n$ (indiv.)")
      ),
      #lty=c(1,3,2,2),
      #col=c(1,2,1,4),
      #lwd=c(3,3,1,4)
      lty=c(1,3,2),
      col=c(1,2,1),
      lwd=c(3,3,1),cex=1.5
      )
    }
    
    
    
    if(lambda.select.ix){
      #pick.lam =0 here avoids scale plot plotting lambda
      pick.lam=0
    }
  }
  #lambda.star.ix=scale.plot(lambdas,sigma.vst,(vst[1]-vst),"test",plott=FALSE,pick.lam=pick.lam) #function(ls,x,y,label){
  #star.ix=lambda.star.ix$ix.star
  
  
  
  # for now just use second delta, assuming always will have one larger and smallr
  
  #if ((d==delta.sel)){  
  if (plot.){
    abline(v=lambdas[star.ix],lty=2)
  }
  #}
  #print(c("gamma",gammas[j],"delta",deltas[d]))
  #print(diff.from.behav)
  print(c("ix selected:",star.ix))
  print(c("selected mask:",selected.mask))
  
  lambda.selects[[jj]]=list(lambda.star.ix=star.ix,
                            #delta=delta,#gamma=gamma,
                            selected.mask=selected.mask, 
                            no.lambda.satisfies=no.lambda.satisfies)
  jj=jj+1
  #}
  #}
  lambda.selects
}

# this does selection/estimation for real data analysis (mimic.R)
select.and.est = function(eps,b0,
                          #gammas,
                          lambdas,
                          #deltas,
                          names,resfile,plotfile,scale.s,
                          #lambda.select.ix=FALSE,gamma.select.ix=2,
                          rdigits=3,use.diff=1e-3,other.args=NULL){

  K=dim(eps[[1]]$Ss)[1]
  
  n = length(eps)

  n.train = n
  #}
  n.test = n-n.train

  train.eps = eps[1:n.train]
  test.eps=eps[(n.train+1):n]

  # selection/estimation
  outer.res = lambda.path(train.eps,b0=b0,
                          lambdas=lambdas,
                          names=names,
                          file=paste0(resfile,"mimic.outer.res"),
                          split=1,scale.s=scale.s,
                          use.diff=use.diff,other.args=other.args)
  
  print(paste0("rollout of selected policy,
               outer res is outer.res[d(delta)][j(gamma)][i(lambda)]"))
  print("note that lambda does change here")
  print(outer.res$outer.res)
  
  nlam = length(lambdas)

  nplots = 2
  area.plots = other.args$area.plots
  res = 100*1.0
  png(plotfile,height=600*area.plots,width=1050,res=res)
  hist.la=0
  if(hist.la){
    w = 10
  }else{
    w = 20#15
  }
  ml=layout(matrix(c(1:(nplots)),ncol=1), 
            widths=rep(w,nplots), 
            heights=rep(c(12,5),1), TRUE)

  # : bottom, left, top, and right.
  par(mar=c(4.1,5.1,2,2))
  ixs=plot.lambda.paths(outer.res,
                        max.n.div=other.args$max.n.div)
  dev.off()

    outer.res.sel=NA
    n.sel=NA

  
  list(ixs=ixs,outer.res=outer.res)
}

check.pos = function(eps.mimic,cov.of.int,pen.b,train.test=1,plot.=TRUE,usesampler=FALSE,other.args){
  # checks positivity, etc for behavioral policy in real data analysis

  
  n = length(eps.mimic)

  
  # split data into test and training
 # if(split){
    # need to deal with this
  
    hf = as.integer(n/2)
    

    


   
    if (train.test){

      if (usesampler){  
      #train.ixs = sample(1:n,hf)  
    train.ix = as.logical(rbinom(n,1,prob=0.5))
    eps.test = eps.mimic[!train.ix]#[(hf+1):n]
    eps.train = eps.mimic[train.ix]#[1:hf]
      }else{
        eps.test = eps.mimic[(hf+1):n]
        eps.train = eps.mimic[1:hf]
        
      }
    
    }else{
      eps.test = eps.mimic
      eps.train = eps.mimic
    }
    

    n.test = length(eps.test)
    n.train = length(eps.train)
#  }else{
#    eps.test=eps
#  }
  
  # if(scale.s){
    so.train = center.scale.s(eps.train)
    eps.train = so.train$e
    
    so.test = center.scale.s(eps.test,means=so.train$means,sds=so.train$sds)
    eps.test = so.test$e

  # }else{
  #   # don t center and scale
  #   so = center.scale.s(eps,means=rep(0,K),sds = rep(1,K)) #not centering and scaling
  #   # but we still need the sds for the scaling factors.
  #   soc = center.scale.s(eps)
  #   sigmas = soc$sds
  #   # generally, we would do the above with scale factors. instead
  #   # override scale factors by setting sigmas to 1. we are not scaling for simulations. 
  #   # we just generate data with variance 1
  #   #we just scale for real data
  #   sigmas = rep(1,K)
  # }
  
  
  # first scale episodes
  #sc=center.scale.s(eps.mimic)
  
  get.matrices  = function(eps,sc,cov.of.int){  
    a=unlist(lapply(eps,function(x){x$As}))
    r=unlist(lapply(eps,function(x){x$Rs}))
    smm=sc$sm.cs
    colnames(smm)=cov.of.int
    list(a=a,r=r,smm=smm)
  }
  
  d.train = get.matrices(eps.train,so.train,cov.of.int)
  a.train = d.train$a
  r.train = d.train$r
  smm.train = d.train$smm 
  
  d.test = get.matrices(eps.test,so.test,cov.of.int)
  a.test = d.test$a
  r.test = d.test$r
  smm.test = d.test$smm 
  

  #hf = 
  
  if(pen.b){

    K = dim(smm.train)[2]
   #lambdas = seq(0,.01,length.out=100)
   
    beh.train=cv.glmnet(smm.train,a.train,intercept = FALSE,alpha=0,
                        family='binomial',
                        type.measure=other.args$type.measure,
                        lambda=other.args$b.lambdas)
    
   

    beh.coef=coef(beh.train,s = "lambda.min")[2:(K+1),]

    #beh.coef=coef(beh.train,s = "lambda.1se")[2:(K+1),]
    
    
    #print(beh.coef)
    pi.train = expit(smm.train%*%beh.coef)
    pi.test = expit(smm.test%*%beh.coef)


  }else{
  beh.train = glm(a.train~smm.train-1,family='binomial')
  
  pi.train = expit(smm.train%*%beh.train$coefficients)#beh.train$fitted.values
  pi.test = expit(smm.test%*%beh.train$coefficients)
  }
  


  #dev.off()
  
  
  
  get.treatment= function(r,a,smm,pi){

  n = length(eps.mimic)
  #na.states=apply(is.na(smm),1,sum)
  #knn(train_set=as.data.frame(cbind(a=a,smm),test_set=as.data.frame(smm),k=2,categorical_target = "a"))
  df.d=as.data.frame(cbind(r,a,smm))
  par.te=summary(lm(r~.,data=df.d))
  mu.1 = 1/n*sum(a*r/(pi)) 
  mu.0 = 1/n*sum((1-a)*r/(1-pi))
  treat.eff = mu.1 - mu.0
  df = as.data.frame(cbind(obs=a,pred=pi))
  #calibration_plot(data=df,obs="a",pred="pred")
  prop.r.less.than.0=mean(r<0)
  prop.a.1 = mean(a)

  # overlap plot
  #png("Overlap.png",width=1000,height=dim(msr)[2]*1000,res=150)
  #min.beh=min(c(1-pi,pi))
  #max.beh=max(c(1-pi,pi))
  #mean.beh=mean(c(1-pi,pi))
  min.beh=min(c(1-pi,pi))
  max.beh=max(c(1-pi,pi))
  mean.beh=mean(pi)
  

  list(par.te=par.te,treat.eff=treat.eff,prop.r.less.than.0=prop.r.less.than.0,
       prop.a.1=prop.a.1,min.beh=min.beh,max.beh=max.beh,mean.beh=mean.beh)
  }

  d.train=get.treatment(r.train,a.train,smm.train,pi.train)
  par.te.train = d.train$par.te
  treat.eff.train = d.train$treat.eff

  d.test=get.treatment(r.test,a.test,smm.test,pi.test)
  par.te.test = d.test$par.te
  treat.eff.test = d.test$treat.eff
  
  
  
  
  prop.r.less.than.0.train=d.train$prop.r.less.than.0
  prop.a.1.train = d.train$prop.a.1
  min.beh.train = d.train$min.beh
  max.beh.train = d.train$max.beh
  mean.beh.train = d.train$mean.beh
  
  

    #par(mfrow=c(2,1))
    #png(paste0("pen.b=",pen.b,"CalibrationTrain.png"),width=500,height=500,res=100)
    my.data.train = data.frame(a.train,pi.train)
    p.train=calibration_plot(data=my.data.train,obs="a.train",pred="pi.train",nTiles=5,title="Train Calibration",data_summary = TRUE)
    if (plot.){
    print(p.train)
    #dev.off()
    }
    #png(paste0("pen.b=",pen.b,"CalibrationTest.png"),width=500,height=500,res=100)
    my.data.test = data.frame(a.test,pi.test)
    p.test=calibration_plot(data=my.data.test,obs="a.test",pred="pi.test",nTiles=5,title="Test Calibration",data_summary = TRUE)
    if (plot.){
    print(p.test)
    }
    
    if (plot.){
    hist(pi.train[as.logical(1-a.train)],col=rgb(.9,.9,.9),xlab=TeX("$\\pi_{b_n}(A_0=1|s)$"),main="Overlap (Training set)")
    hist(pi.train[as.logical(a.train)],add=TRUE,col=rgb(.5,.5,.5))
    legend("topright",lty=c(1,1),lwd=c(4,4),c("Observed A=0","Observed A=1"),col=c(rgb(.9,.9,.9),rgb(.5,.5,.5)))
    #dev.off()
    }
  
  

  tb=c("n"=as.integer(n),
       #"T"=as.integer(T),
       "Prop. R<0"=round(prop.r.less.than.0.train,3),
       "Prop A=1"=round(prop.a.1.train,3),
       "min $\\pi_b$"=round(min.beh.train,3),
       "mean $\\pi_b$"=round(mean.beh.train,3),
       "max $\\pi_b$"=round(max.beh.train,3),
       "IPTW treat. eff. train"=round(treat.eff.train,3),
       "IPTW.treat.eff.test"=round(treat.eff.test,3),
       "Par treat eff train (pval)"=paste0(round(par.te.train$coefficients['a','Estimate'],3)," (",round(par.te.train$coefficients['a','Pr(>|t|)'],3),")"),
       "Par. treat. eff test (p-val)"=paste0(round(par.te.test$coefficients['a','Estimate'],3)," (",round(par.te.test$coefficients['a','Pr(>|t|)'],3),")")
  )

  dftb=as.data.frame(tb)
  #print(xtable(dftb),
  #      type="latex",sanitize.text.function = function(x){x},
  #      include.rownames=TRUE,digits=2)

  list(n=length(eps.mimic),
       #T=T,
       min.beh.train=min.beh.train,
       max.beh.train=max.beh.train,
       pi.train=pi.train,
       pi.test=pi.test,
       a.train=a.train,
       a.test=a.test,
       mean.beh.train=mean.beh.train,
       summ=summary(beh.train),
       IPTW.treat.eff.train=treat.eff.train,
       par.treat.eff.test=par.te.train$coefficients['a',],
       IPTW.treat.eff.test=treat.eff.test,
       par.treat.eff.test=par.te.test$coefficients['a',],
       prop.r.less.than.0.train=prop.r.less.than.0.train,
       prop.a.1.train=prop.a.1.train,
       dftb=dftb,beh.train=beh.train,
       p.train=p.train$data_summary,
       p.test=p.test$data_summary)
}


# runs a single selection
run.one.sel = function(param,m){
  
  # both sim and mimic
  
  exp.tag = param$exp.tag
  resfile = paste0(exp.tag,"/lambplots.outerress/",m)
  lambdas = param$lambdas
  nlam = length(lambdas)
  b0 = param$b0
  scale.s=param$scale.s
  res.dir.sel = param$res.dir.sel
  rdigits = param$rdigits
  use.diff = param$use.diff
  other.args= param$other.args
  
  # sim only
  n = param$n
  tau=param$tau
  R=param$R
  M= param$M
  K = param$K
  sd.epsi=param$sd.epsi
  mean.first=param$mean.first
  sd.first=param$sd.first
  u=param$u
  plt.obj=param$plt.obj
  mysel=param$mysel  
  names =  param$names  
  init.state.mean=param$init.state.mean
  init.state.var = param$init.state.var
  trans.var = param$trans.var
  
  #mimic only
  eps.mimic=param$eps.mimic
  renamed.cov.of.int=param$renamed.cov.of.int
  resfile = param$resfile
  plotfile=param$plotfile
  

  if (is.null(eps.mimic)){   # if it's simulation
  
  eps=gen.data.3arg(b0=b0,n=n,tau=tau,init.state.mean=init.state.mean,
                 init.state.var=init.state.var,trans.var=trans.var,R=R,seed=m,
                 other.args=other.args)
  

  n = length(eps)
  
  outer.res = lambda.path(eps,b=b0,
                          #gammas=gammas,
                          lambdas=lambdas,
                          names=names,
                          #deltas=deltas,
                          file="mc",split=1,scale.s=scale.s,
                          #sel=rep(1,K),
                          #lambda.select.ix=lambda.select.ix,
                          #gamma.sel=gamma.select.ix,
                          use.diff=use.diff,
                          other.args=other.args)
  
  
  ixs=plot.lambda.paths(outer.res,plot=FALSE)
  
  select.mask = ixs[[1]]$selected.mask
  #}
  print(ixs)
  nolam = ixs[[1]]$no.lambda.satisfies
  print(c("No lambda satisfies:",nolam))
  print("select mask")
  print(select.mask)

    n.sel=NULL
    gamma.tables=NULL
    outer.res.sel=NULL
    oplist = list(outer.res=outer.res, 
                  #outer.res.sels=outer.res.sels,
                  outer.res.sel=outer.res.sel,
                  gamma.tables=gamma.tables,
                  select.mask=select.mask,ixs=ixs)
  }else{ # it's mimic
    n=length(eps.mimic)
    set.seed(m)
    # just reorder indices

    reorder.ixs = sample(1:n,replace=FALSE) 
    # so when we take a test train split, we will now take different episodes each time
    eps.mimic = eps.mimic[reorder.ixs]
    
    ixs=select.and.est(eps.mimic,b0=rep(NA,dim(eps.mimic[[1]]$Ss)[1]),
                       lambdas=lambdas,
                       names=renamed.cov.of.int,
                       resfile=paste0(exp.tag,"res"),
                       plotfile=paste0(exp.tag,"grid.mimic.png"),
                       scale.s=scale.s,
                       use.diff=use.diff,other.args=other.args)
    
    n.sel=NULL
    gamma.tables=NULL
    outer.res.sel=NULL
    
    select.mask = ixs$ixs[[1]]$selected.mask
    
    oplist = list(outer.res=ixs$outer.res, 
                  #outer.res.sels=outer.res.sels,
                  outer.res.sel=outer.res.sel,
                  gamma.tables=gamma.tables,
                  select.mask=select.mask,ixs=ixs$ixs)
    
  }
  

  saveRDS(oplist,paste0(res.dir.sel,"/",m))
  
}




# this function not done, but started to do it
# check.pos.repeated = function(eps.mimic,cov.of.int,pen.b,train.test=1){
#   # checks positivity, etc for behavioral policy in real data analysis
#   
#   res = list()
#   
#   m = 3
#   for (i in 1:m){
#     n = length(eps.mimic)
#     
#     
#     # split data into test and training
#     # if(split){
#     # need to deal with this
#     
#     hf = as.integer(n/2)
#     
#     if (train.test){
#       usesampler=0
#       if (usesampler){  
#         #train.ixs = sample(1:n,hf)  
#         train.ix = as.logical(rbinom(n,1,prob=0.5))
#         eps.test = eps.mimic[!train.ix]#[(hf+1):n]
#         eps.train = eps.mimic[train.ix]#[1:hf]
#       }else{
#         eps.test = eps.mimic[(hf+1):n]
#         eps.train = eps.mimic[1:hf]
#         
#       }
#       
#     }else{
#       eps.test = eps.mimic
#       eps.train = eps.mimic
#     }
#     
#     
#     n.test = length(eps.test)
#     n.train = length(eps.train)
#     #  }else{
#     #    eps.test=eps
#     #  }
#     
#     # if(scale.s){
#     so.train = center.scale.s(eps.train)
#     eps.train = so.train$e
#     
#     so.test = center.scale.s(eps.test,means=so.train$means,sds=so.train$sds)
#     eps.test = so.test$e
#     
#     # }else{
#     #   # don t center and scale
#     #   so = center.scale.s(eps,means=rep(0,K),sds = rep(1,K)) #not centering and scaling
#     #   # but we still need the sds for the scaling factors.
#     #   soc = center.scale.s(eps)
#     #   sigmas = soc$sds
#     #   # generally, we would do the above with scale factors. instead
#     #   # override scale factors by setting sigmas to 1. we are not scaling for simulations. 
#     #   # we just generate data with variance 1
#     #   #we just scale for real data
#     #   sigmas = rep(1,K)
#     # }
#     
#     
#     # first scale episodes
#     #sc=center.scale.s(eps.mimic)
#     
#     get.matrices  = function(eps,sc,cov.of.int){  
#       a=unlist(lapply(eps,function(x){x$As}))
#       r=unlist(lapply(eps,function(x){x$Rs}))
#       smm=sc$sm.cs
#       colnames(smm)=cov.of.int
#       list(a=a,r=r,smm=smm)
#     }
#     
#     d.train = get.matrices(eps.train,so.train,cov.of.int)
#     a.train = d.train$a
#     r.train = d.train$r
#     smm.train = d.train$smm 
#     
#     d.test = get.matrices(eps.test,so.test,cov.of.int)
#     a.test = d.test$a
#     r.test = d.test$r
#     smm.test = d.test$smm 
#     
#     
#     #hf = 
#     
#     if(pen.b){
#       
#       K = dim(smm.train)[2]
#       #lambdas = seq(0,.01,length.out=100)
#       #lambdas = seq(0,.1,length.out=100)
#       beh.train=cv.glmnet(smm.train,a.train,intercept = FALSE,alpha=0,
#                           family='binomial')#,b.lambda=b.lambdas)
#       
#       
#       beh.coef=coef(beh.train,s = "lambda.min")[2:(K+1),]
#       #beh.coef=coef(beh,s = "lambda.1se")[2:(K+1),]
#       
#       
#       print(beh.coef)
#       pi.train = expit(smm.train%*%beh.coef)
#       pi.test = expit(smm.test%*%beh.coef)
#       
#       
#     }else{
#       beh.train = glm(a.train~smm.train-1,family='binomial')
#       
#       pi.train = expit(smm.train%*%beh.train$coefficients)#beh.train$fitted.values
#       pi.test = expit(smm.test%*%beh.train$coefficients)
#     }
#     
#     
#     
#     
#     
#     
#     
#     get.treatment= function(r,a,smm,pi){
#       
#       n = length(eps.mimic)
#       #na.states=apply(is.na(smm),1,sum)
#       #knn(train_set=as.data.frame(cbind(a=a,smm),test_set=as.data.frame(smm),k=2,categorical_target = "a"))
#       df.d=as.data.frame(cbind(r,a,smm))
#       par.te=summary(lm(r~.,data=df.d))
#       mu.1 = 1/n*sum(a*r/(pi)) 
#       mu.0 = 1/n*sum((1-a)*r/(1-pi))
#       treat.eff = mu.1 - mu.0
#       df = as.data.frame(cbind(obs=a,pred=pi))
#       #calibration_plot(data=df,obs="a",pred="pred")
#       prop.r.less.than.0=mean(r<0)
#       prop.a.1 = mean(a)
#       
#       # overlap plot
#       #png("Overlap.png",width=1000,height=dim(msr)[2]*1000,res=150)
#       #min.beh=min(c(1-pi,pi))
#       #max.beh=max(c(1-pi,pi))
#       #mean.beh=mean(c(1-pi,pi))
#       min.beh=min(c(1-pi,pi))
#       max.beh=max(c(1-pi,pi))
#       mean.beh=mean(pi)
#       
#       hist(pi[as.logical(1-a)],col=rgb(.9,.9,.9),xlab=TeX("$\\pi_{b_n}(A_0=1|s)$"),main="Overlap")
#       hist(pi[as.logical(a)],add=TRUE,col=rgb(.5,.5,.5))
#       legend("topright",lty=c(1,1),lwd=c(4,4),c("Observed A=0","Observed A=1"),col=c(rgb(.9,.9,.9),rgb(.5,.5,.5)))
#       #dev.off()
#       list(par.te=par.te,treat.eff=treat.eff,prop.r.less.than.0=prop.r.less.than.0,
#            prop.a.1=prop.a.1,min.beh=min.beh,max.beh=max.beh,mean.beh=mean.beh)
#     }
#     
#     d.train=get.treatment(r.train,a.train,smm.train,pi.train)
#     par.te.train = d.train$par.te
#     treat.eff.train = d.train$treat.eff
#     
#     d.test=get.treatment(r.test,a.test,smm.test,pi.test)
#     par.te.test = d.test$par.te
#     treat.eff.test = d.test$treat.eff
#     
#     
#     
#     
#     prop.r.less.than.0.train=d.train$prop.r.less.than.0
#     prop.a.1.train = d.train$prop.a.1
#     min.beh.train = d.train$min.beh
#     max.beh.train = d.train$max.beh
#     mean.beh.train = d.train$mean.beh
#     
#     res[[i]]=list(n=length(eps.mimic),
#                   #T=T,
#                   min.beh.train=min.beh.train,
#                   max.beh.train=max.beh.train,
#                   mean.beh.train=mean.beh.train,
#                   summ=summary(beh.train),
#                   IPTW.treat.eff.test=treat.eff.train,
#                   par.treat.eff.test=par.te.train$coefficients['a',],
#                   IPTW.treat.eff.test=treat.eff.test,
#                   par.treat.eff.test=par.te.test$coefficients['a',],
#                   prop.r.less.than.0.train=prop.r.less.than.0.train,
#                   prop.a.1.train=prop.a.1.train,
#                   beh.train=beh.train,
#                   pi.train=pi.train,
#                   pi.test=pi.test)
#     
#   }
#   
#   
#   browser()
#   lapply(res,function(x){x$})
#   
#   #par(mfrow=c(2,1))
#   #png(paste0("pen.b=",pen.b,"CalibrationTrain.png"),width=500,height=500,res=100)
#   my.data.train = data.frame(a.train,pi.train)
#   p=calibration_plot(data=my.data.train,obs="a.train",pred="pi.train",nTiles=5,title="Train Calibration")
#   print(p)
#   #dev.off()
#   
#   #png(paste0("pen.b=",pen.b,"CalibrationTest.png"),width=500,height=500,res=100)
#   my.data.test = data.frame(a.test,pi.test)
#   p=calibration_plot(data=my.data.test,obs="a.test",pred="pi.test",nTiles=5,title="Test Calibration")
#   print(p)
#   #dev.off()
#   
#   
#   
#   tb=c("n"=as.integer(n),
#        #"T"=as.integer(T),
#        "Prop. R<0"=round(prop.r.less.than.0.train,3),
#        "Prop A=1"=round(prop.a.1.train,3),
#        "min $\\pi_b$"=round(min.beh.train,3),
#        "mean $\\pi_b$"=round(mean.beh.train,3),
#        "max $\\pi_b$"=round(max.beh.train,3),
#        "IPTW treat. eff. train"=round(treat.eff.train,3),
#        "IPTW.treat.eff.test"=round(treat.eff.test,3),
#        "Par treat eff train (pval)"=paste0(round(par.te.train$coefficients['a','Estimate'],3)," (",round(par.te.train$coefficients['a','Pr(>|t|)'],3),")"),
#        "Par. treat. eff test (p-val)"=paste0(round(par.te.test$coefficients['a','Estimate'],3)," (",round(par.te.test$coefficients['a','Pr(>|t|)'],3),")")
#   )
#   
#   dftb=as.data.frame(tb)
#   
#   #print(xtable(dftb),
#   #      type="latex",sanitize.text.function = function(x){x},
#   #      include.rownames=TRUE,digits=2)
#   
#   
# }



# converts a list of MC experiment results into summary statistics (ie, mean, variance) 
# ie if we have 2 objects with results for 2 datasets, we will create an identical object
# that instead has the average of the 2 results in it (and also the variance as a sublist)
list.outer.res.to.summ = function(ors,mean=1){

  outer.res = ors[[1]]#$outer.res
  or.shell.av = outer.res
  or.shell.var = outer.res
  lambdas = outer.res$lambdas


  number.lam.gam.del.items = length(unlist(outer.res$outer.res))
  # for each MC iteration, make a column.  

  M =length(ors)
  r = matrix(nrow=number.lam.gam.del.items,ncol=M)
  for (m in 1:M){
    r[,m] = unlist(ors[[m]]$outer.res)
  }
  
  # then run stats. note that even taking averages 
  # of like n.train, which is just n.train
  #for (mean in c(0,1)){
  mean = 1
  if (mean){
    rav = apply(r,1,mean)
  }else{
    rav = apply(r,1,median)
  }
  rvar = apply(r,1,var)

  # one lambda
  # just get the names of one element of outer.res (one process.eps output)
  m.names=names(outer.res[[1]][[1]])#names(outer.res[[1]][[1]][[1]])
  k=1

      for (i in 1:length(lambdas)){
        for (na in m.names){
          
          # get the dimension of this named item (just use delta=gamma=lambda=1, since 
          # dimension is the same for all delta, gamma, lambda)
          # the last [[1]] is just to get the item itself
          m.dim = dim(as.matrix(or.shell.av$outer.res[[1]][na][[1]]))
          # cast as matrix so can get dim.  In the plot.lambda paths, we will uncast
          # to scalar again for some, such as sample size n
          if(is.null(m.dim)){m.dim=c(1,1)}
          # now because we unlisted this item, it should be rav dimensions 
          # k:k+(m.dim[1]*m.dim[2]-1), minus one is there because we start with k=1 not k=0
          av.res = matrix(rav[k:(k+(m.dim[1]*m.dim[2]-1))],
                          nrow=m.dim[1],ncol=m.dim[2])
          va.res = matrix(rvar[k:(k+(m.dim[1]*m.dim[2]-1))],
                          nrow=m.dim[1],ncol=m.dim[2])
          or.shell.av$outer.res[[i]][na][[1]] = av.res
          or.shell.var$outer.res[[i]][na][[1]] = va.res
          k=k+m.dim[1]*m.dim[2]
        }
      }

  # put variance in or.av, because so much code considers or.av
  
  print(c("outer.res av",or.shell.av))
  
  or.shell.av$or.shell.var = or.shell.var
  or.shell.av
}


# analyze results of a MC experiment, or of a summarized MC experiment
analyze.sel.res = function(res.dir,sel.res.dir,plot.dir,rdigits=3,other.args){
  hist.la = 1
  if(hist.la){
    w = 10
  }else{
    w = 8
  }
  print(sel.res.dir)
  print(plot.dir)
  fs = list.files(path=sel.res.dir,full.names=TRUE)
  print(fs)
  fs = mixedsort(fs)
  ors.and.tbs = lapply(fs,readRDS)

  
  masks = lapply(ors.and.tbs,function(x){x$select.mask})
  unique.masks = unique(masks)#masks[!duplicated(lapply(masks, sort))]
  
  print("Are all selected masks the same?")

  mm=unlist.reshape(masks)
  print(mm)
  print(apply(mm,2,var)==0)
  saveRDS(mm,paste0(res.dir,'/masks.matrix'))
  
  ms = matrix(apply(mm,2,mean),nrow=1)
  m.names=ors.and.tbs[[1]]$outer.res$names 
  colnames(ms) = m.names
  rownames(ms) = "$\\hat{P}(selected)$"
  print(xtable(ms),type="latex",sanitize.text.function = function(x){x})
  

  Cns = apply(mm,1,sum)
  if (length(Cns)>1){
  Cn.dist.tab=table(Cns)
  Cn.dist.tab = Cn.dist.tab/sum(Cn.dist.tab)
  Cns.dist.tab=data.frame(Cn.dist.tab)
  cns.dist.tab=t(Cns.dist.tab)
  rownames(cns.dist.tab) = c("$d_{n,\\lambda_n}$","$\\hat{P}(D_{n,\\lambda_n}=d_{n,\\lambda_n})$")
  print(xtable(cns.dist.tab,caption="Distribution of $D_{n,\\lambda_n},$ the number of empirically selected coefficients,
               over Monte-Carlo datasets."),
        include.colnames=FALSE,
        type="latex",sanitize.text.function = function(x){x})
  }
  #ms = apply(mm,2,mean)
  #print(xtable(ms))
  
  ors = lapply(ors.and.tbs,function(x){x$outer.res})
  ors.sel =  lapply(ors.and.tbs,function(x){x$outer.res.sel})
  lambda.ns = unlist(lapply(lapply(ors.and.tbs, function(x){x$ixs}),function(x){x[[1]]$lambda.star.ix}))
  
  vn.lambda.stars = rep(NA,length(ors))
  vns = rep(NA,length(ors))
  vn.bns =  rep(NA,length(ors))
  vn.bowls = rep(NA,length(ors))
  
  rollout.betanmeans = rep(NA,length(ors))
  rollout.betanlambdameans = rep(NA,length(ors.sel))
  rollout.bnmeans  = rep(NA,length(ors.sel))
  rollout.bowl.means  = rep(NA,length(ors.sel))
  
  
  spar.betangamma = rep(NA,length(ors))
  spar.betangammalambda = rep(NA,length(ors.sel))
  spar.bn  = rep(0,length(ors.sel))
                                 
  other.args = ors[[1]]$other.args
  
  
  

  for (i in 1:length(ors)){
   vn.lambda.stars[i] =  ors[[i]]$outer.res[[lambda.ns[i]]]$vn.betagammalambda
   vns[i] = ors[[i]]$outer.res[[lambda.ns[i]]]$vn.betagamma
   vn.bns[i] = ors[[i]]$outer.res[[lambda.ns[i]]]$vn.bn
   
   vn.bowls[i] = ors[[i]]$outer.res[[lambda.ns[i]]]$value.oracle.train
   
   rollout.betanmeans[i] = ors[[i]]$outer.res[[lambda.ns[i]]]$rollout.betangamma.mean
   rollout.betanlambdameans[i] = ors[[i]]$outer.res[[lambda.ns[i]]]$rollout.betangammalambda.mean
   rollout.bnmeans[i] = ors[[i]]$outer.res[[lambda.ns[i]]]$rollout.bnmean
   rollout.bowl.means[i] = ors[[i]]$outer.res[[lambda.ns[i]]]$mean.pred.oracle.train
   
   
   spar.betangammalambda[i] = sum(abs(ors[[i]]$outer.res[[lambda.ns[i]]]$diff.bn.v.betangammalambda)>other.args$use.diff)
   spar.betangamma[i] = sum(abs(ors[[i]]$outer.res[[lambda.ns[i]]]$diff.bn.v.betangamma)>other.args$use.diff)
   
  }
  
  
  

  # to get vars
  #unlist.reshape(lapply(res, function(x){x$vn.betagamma}))
  #unlist.reshape(lapply(res, function(x){x$vn.betagammalambda}))
  
  
  g.tbs =  lapply(ors.and.tbs,function(x){x$gamma.tables})
  
  M = length(ors)

  # get info to make plots
  outer.res = ors[[1]]
  
  gamma.select.ix=outer.res$gamma.sel
  K=dim(outer.res$outer.res[[1]]$betan)[1]
  #gammas = outer.res$gammas
  #ngam = length(gammas)
  #deltas = outer.res$deltas
  #ndel = length(deltas)
  lambdas = outer.res$lambdas
  nlam = length(lambdas)
  n = outer.res$n
  #T= outer.res$T
  use.mean = 1 # this is Ignored, and mean=1 in list.outer.res.to.summ
  ors.av = list.outer.res.to.summ(ors,mean=use.mean)
  res = 100*1.3
  #nplots=2

  tag = "noTag"
  other.args = outer.res$other.args
  area.plots = other.args$area.plots
  #plotfile = paste0(sel.res.dir,"/",tag,"mean",mean,
  # "_seed=",sd,"n=",n,"M=",M,"K=",K,"T=",T,"_MC.reps.png")
  plotfile = paste0(plot.dir,"/",tag,"use.mean",use.mean,
                    "n=",n,"M=",M,"K=",K,
                    #"T=",T,
                    "_MC.reps.png")


  print(plotfile)
  png(plotfile,height=600*area.plots,width=950,res=res)
  
  if (hist.la){
    nplots=3
    hts = rep(c(4.5,2.5,2))
  }else{
    nplots=2
    hts = rep(c(4.5,2.5))
  }
  
  ml=layout(matrix(c(1:(nplots)),ncol=1), 
            widths=rep(w,nplots), 
         heights=hts, TRUE)
  


  ixs=plot.lambda.paths(ors.av,
                        max.n.div=other.args$max.n.div,lambda.ns=lambda.ns)

  
  
  if (hist.la){
    # bottom, left, top, right
    par(mar=c(4,5.1,1,2))
    #par(mar=c(1,5.1,2,2))
    hist(lambdas[lambda.ns],xlab=TeX("$\\lambda$"), 
         xlim=c(min(lambdas),max(lambdas)),main="",ylab=TeX("Count $\\lambda_n$"),
         cex.lab=1.5,cex.axis=1.5,,cex.main=1.5,cex=1.5)
    box()
  }
  
  dev.off()
  saveRDS(ors.av,paste0(res.dir,'/ors.av'))


  all.v=c(vn.lambda.stars,vns,vn.bns)
  d = hist(all.v)
  xlim = range(d$breaks)
  
  
  all.vr=c(rollout.betanmeans,rollout.betanlambdameans,rollout.bnmeans)
  dr = hist(all.vr)
  xlimr = range(dr$breaks)
  
  plot.ora = 1-other.args$real.data
  plot.ora = 0
  
  if (plot.ora){
    height=1500
  }else{
    height=1000
  }
  
  png(paste0(plotfile,"VnDist.png"),res=300,width=2400,height=height)
  #bottom, left, top, and right. 
  par(mar=c(4,5,2,1))
  if (plot.ora){
  par(mfcol=c(3,3))
  }else{
    par(mfcol=c(2,3))
  }
  
  
  
  if (plot.ora){
  #hist(vn.bowls,xlim=xlim,xlab=TeX("Unconstrained (theor): $V_n(\\beta_0)$"),main='Value',breaks=10)
  #box()
  main=NULL

  
  hist(vns,main='Value',xlim=xlim,xlab=TeX("$V_n(\\beta_n)$ "),
       breaks=10,ylab=substitute(paste(bold("Unconstrained"))),cex.lab=1.2)
  box()
  
  }else{
    main = "Value" 
  }
  
  print("Value means:")
  print(c("mean value oracle",mean(vns)))
  print(c("mean value sugg",mean(vn.lambda.stars)))
  print(c("mean value beh",mean(vn.bns)))

  hist(vn.lambda.stars,xlim=xlim,main=main,xlab=TeX("$V_n( \\beta_{n, \\lambda_n})$"),
       breaks=10,ylab=substitute(paste(bold("Suggested"))),cex.lab=1.2)
  box()
  
  hist(vn.bns,xlim=xlim,main='',xlab=TeX("$V_n(b_n)$ "),
       breaks=10,ylab=substitute(paste(bold("Behavioral"))),cex.lab=1.2)
  box()

  
  
  #vn.bowls,rollout.bowl.means
  
  
  
  if (plot.ora){
  #hist(rollout.bowl.means,xlim=xlimr,main="Probability treatment",
  #     xlab=TeX("Unconstrained (theor): mean $\\pi(A_0=1|s_0;\\beta_0)$"),breaks=seq(0,1,by=0.01))
  #Mean $\\pi_{\\beta_{n, \\lambda = \\infty}}(A_0=1|S_0=s_0)$ (i.e.  $\\pi_{b_n}(A_0=1|S_0=s_0)$) "))
  #box()
  #main = NULL

  
  
  hist(rollout.betanmeans,xlim=xlimr,main='Prob(treatment)',
       xlab=TeX("Mean $\\pi(A_0=1|s_0;\\beta_n)$"),breaks=seq(0,1,by=0.01),ylab='')
                #$\\pi_{\\beta_{n, \\lambda = 0}}(A_0=1|S_0=s_0)$ (i.e. $\\pi_{\\beta_{n}}(A_0=1|S_0=s_0)$)"))
  box()
  }else{
    main="Prob(treatment)"
  }
  
  
  hist(rollout.betanlambdameans,xlim=xlimr,main=main,
       xlab=TeX("Mean $\\pi(A_0=1|s_0;\\beta_{n,\\lambda_n})$"),breaks=seq(0,1,by=0.01),ylab='')
                #Mean $\\pi_{\\beta_{n,\\lambda = \\lambda_n}}(A_0=1|S_0=s_0)$"),ylab='')
  box()
  hist(rollout.bnmeans,xlim=xlimr,main="",
       xlab=TeX("Mean $\\pi(A_0=1|s_0;b_n)$"),breaks=seq(0,1,by=0.01),ylab='')
                #Mean $\\pi_{\\beta_{n, \\lambda = \\infty}}(A_0=1|S_0=s_0)$ (i.e.  $\\pi_{b_n}(A_0=1|S_0=s_0)$) "))
  box()
  

  
  #bottom, left, top, and right. 
  #par(mar=c(4,5,1,2))
  
  if (plot.ora){
  #plot.new()
  
  
  hist(spar.betangamma,xlim=c(0,K),main="Relative sparsity",
       xlab=TeX("# coef. diff. from behavior ($d_{n,\\lambda=0}$)"),
       breaks=seq(0,K,by=.5),ylab='')
  #$\\pi_{\\beta_{n, \\lambda = 0}}(A_0=1|S_0=s_0)$ (i.e. $\\pi_{\\beta_{n}}(A_0=1|S_0=s_0)$)"))
  box()
  
  }else{
    main="Relative sparsity"
  }
  
  hist(spar.betangammalambda,xlim=c(0,K),main=main,
       xlab=TeX("# coef. diff. from behavior ($d_{n,\\lambda=\\lambda_n}$)"),
       breaks=seq(0,K,by=.5),ylab='')
  #Mean $\\pi_{\\beta_{n,\\lambda = \\lambda_n}}(A_0=1|S_0=s_0)$"))
  box()
  hist(spar.bn,xlim=c(0,K),main="",
       xlab=TeX("# coef. diff. from behavior ($d_{n,\\lambda=\\infty}$)"),breaks=seq(0,K,by=.5),ylab='')
  #Mean $\\pi_{\\beta_{n, \\lambda = \\infty}}(A_0=1|S_0=s_0)$ (i.e.  $\\pi_{b_n}(A_0=1|S_0=s_0)$) "))
  box()
  

  
  
  
  

  
  

  dev.off()
  
 

  
  # use the selection from the average of lambda plots

  nolam = NULL
  for (mysel in unique.masks){
    
  mask.indices = unlist(lapply(masks,function(x){all(x==mysel)}))

  # we only select the ones with the same mask indices
  ors.sel.av = list.outer.res.to.summ(ors.sel[mask.indices])
 
  # masks 
  n = ors.sel.av$n # or is it n.test? yes, it's 1/2 n

  }
  saveRDS(unique.masks,paste0(res.dir,'/unique.masks'))
  unique.masks
  
}

