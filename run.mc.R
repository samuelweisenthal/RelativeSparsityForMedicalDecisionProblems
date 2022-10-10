# main file to run Monte-Carlo experiments.  Generates a directory with results

library(latex2exp)
library(xtable)
library(gtools)

# useful functions
source('utils.R')
source('sim.R')
source('center.scale.R')
source('is.R') #vn
source('mc.utils.R')

# bh=1 means run on server
bh=0
debug =1

# whether to not estimate nuisance
dontestb=0

# other args that can be passed
other.args=list(logx=0,
                area.plots=1.3,
                plot.emp.var=1,
                use.diff=0.01,max.n.div=1,divbys=0)

# scale states
scale.s=1

# if use.diff is a number, then the difference between beta_nlambdak and bnk will be used
# to select, and that number will be the threshold
# ie if |beta-b|< use.diff, that coeff is not active
# if use.diff is NULL, then will see if the confidence interval of betangammalambdak captures bnk instead
use.diff = other.args$use.diff

# anonymousple size
n = 500
# number covariates
K=9

# number MC repetitions
plotsM = 500

# number lambdas in grid
nlam = 10
lambdas = seq(3,5,length.out=nlam)#seq(0.01,5,length.out=nlam)#seq(exp(-7),exp(0),length.out=3)#exp(seq(-7,0,length.out=5))#seq(0,.05,length.out=10)

# whether to censor observations, experimental
censor=0


if (debug){
  use.diff = 0.001#1e-1
  n = 22
  M=5
  plotsM = 2
  K=2
  nlam =2
  lambdas = seq(0.001,20,length.out=nlam)
}

# treatment effect
tau=matrix(rep(0,K))
tau[2,]=2

# name of experiment
str.tag = "Example_MC_Experiment"
other.arg.str= paste0(paste0(names(other.args),"=",other.args),collapse = "")
exp.tag = paste0(str.tag,"Res.Dir.n=",n,
                 #"T=",T,
                 "M=",M,"plotsM=",plotsM,"K=",K,
                 #"gammaselix=",gamma.select.ix,"lambda.select.ix=",lambda.select.ix,
                 "maxLam=",formatC(max(lambdas), format = "e", digits = 2),
                 "minLam=",formatC(min(lambdas),format="e",digits=2),
                 #"maxGamma=",formatC(max(gammas), format = "e", digits = 2),
                 #"minGamma=",formatC(min(gammas), format = "e", digits = 2),
                 "usediff=",use.diff,"tau2=",tau[2,])
# create directory specificially for this experiment
dir.create(exp.tag)

# create log files, etc
# main stdout
outfile = paste0(exp.tag,"/mainlog.txt")

#plot
plotfile = paste0(exp.tag,"/grid.mc.png")

#result object that is saved and can be read later
resfile = paste0(exp.tag,"/outer.res")

#parameters
param.path.sel = paste0(exp.tag,"/paramsel")

# selection tables
sel.tab.dir = paste0(exp.tag,"/sel.tables")

#selection result objects
res.dir.sel = paste0(exp.tag,"/res.sel")

#logs
mylogdir = paste0(exp.tag,"/logs")

# plot dir
res.dir.sel.plots =  paste0(exp.tag,"/plots")
dir.create(res.dir.sel)
dir.create(res.dir.sel.plots)
dir.create(mylogdir)

# save this parameter configuration in experiment directory as a record for reproducibility
if (!debug){
  sink(file=outfile)
  writeLines(readLines("run.mc.R"))
}


dir.create(sel.tab.dir)

# define reward
#R = function(s,a,sp,T){
R = function(s,a,sp){
  if (is.na(matrix(s)[1,])){
    r= 0
  }else{
    r = sp[2,]#(sp[2,]-s[2,])/s[2,]
  }
  if(length(r)>1){
    print("r/state dimension mismatch.")
    browser()
  }
  r
}


# true behavioral coefficients
b0=matrix(rep(0,K))
b0[1,]=-0.3
b0[2,]=0.2

init.state.mean = matrix(rep(3,K))

addcor = (1-diag(rep(1,K)))*0.01
if(K>2){
  addcor[,3:K] = 0 
  addcor[3:K,] = 0
}

init.state.var=diag(K) + addcor
trans.var = 1*diag(K) + addcor

sd.epsi=diag(rep(1,K)) # + addcor
mean.first=matrix(rep(0,K))
sd.first=diag(K) + addcor


# for tables, how many difits to use
rdigits = 2

# name of covariates in simulations
names = paste0("$S_{0,",1:K,"}$")

params.sel  =  list(
                    lambdas=lambdas,
                    b0=b0,n=n,tau=tau,
                    R=R,M=NULL,plotsM=plotsM,K=K,
                    sd.epsi=sd.epsi,mean.first=mean.first,
                    init.state.mean=init.state.mean,init.state.var=init.state.var,trans.var=trans.var,
                    sd.first=sd.first,
                    u=NULL,plt.obj=FALSE,scale.s=scale.s,
                    exp.tag=exp.tag,
                    res.dir.sel = res.dir.sel,rdigits=rdigits,names=names,use.diff=use.diff,bh=bh,
                    other.args=other.args)

saveRDS(params.sel,param.path.sel)


if (bh==1){
  # if on server, parallelize the creation and processing of individual datasets
  dir.create('../../sbatch_scripts')
  dir.create(paste0(exp.tag,"/logs.sel"))
  job.array.sbatch = "../../sbatch_scripts/run.sel.sbatch"
  custom.file.sel = paste0(job.array.sbatch,'.',exp.tag)
  unlink(custom.file.sel) 
  write(paste0("#!/bin/sh"),custom.file.sel,append=TRUE) #maybe needs newline?
  write(paste0("#SBATCH -n 1"),custom.file.sel,append=TRUE)
  write(paste0("#SBATCH --mail-type=begin"),custom.file.sel,append=TRUE)
  write(paste0("#SBATCH --mail-type=end"),custom.file.sel,append=TRUE)
  write(paste0("#SBATCH -a 1-",plotsM),custom.file.sel,append=TRUE)
  write(paste0("#SBATCH -o ", exp.tag,"/logs.sel/log%a.txt"),custom.file.sel,append=TRUE)
  write(paste0("#SBATCH -p interactive --time=0-10:00:00"),custom.file.sel,append=TRUE)
  write(paste0("#SBATCH --mem=10gb"),custom.file.sel,append=TRUE)
  write(paste0("module load r"),custom.file.sel,append=TRUE)
  write(paste0("cd ",getwd()),custom.file.sel,append=TRUE)
  write(paste0("Rscript ",getwd(),"/run.one.select.R ",
               param.path.sel," $SLURM_ARRAY_TASK_ID"),
        custom.file.sel,append=TRUE)
  system(paste0("sbatch ",custom.file.sel))
  
  analysis.sbatch = "../../sbatch_scripts/analyze.sel.sbatch"
  
  custom.file.sel.2 = paste0(analysis.sbatch,'.',exp.tag)
  unlink(custom.file.sel.2) 
  write(paste0("#!/bin/sh"),custom.file.sel.2,append=TRUE) #maybe needs newline?
  write(paste0("#SBATCH -p preempt --time=0-00:10:00"),custom.file.sel.2,append=TRUE)
  write(paste0("#SBATCH -n 1"),custom.file.sel.2,append=TRUE)
  write(paste0("#SBATCH --mem=10gb"),custom.file.sel.2,append = TRUE)
  write(paste0("#SBATCH --mail-type=begin"),custom.file.sel.2,append=TRUE)
  write(paste0("#SBATCH --mail-type=end"),custom.file.sel.2,append=TRUE)
  write(paste0("#SBATCH -o ", exp.tag,"/logs/sel.table.txt"),custom.file.sel.2,append=TRUE)
  write(paste0("module load r"),custom.file.sel.2,append=TRUE)
  write(paste0("cd ",getwd()),custom.file.sel.2,append=TRUE)
  write(paste0("Rscript ",getwd(),"/analyze.sel.R ", 
               exp.tag," ",
               res.dir.sel," ",res.dir.sel.plots),custom.file.sel.2,append=TRUE)
  
  mc.guide.sbatch = "../../sbatch_scripts/mc.guide.sbatch"
  
  
  custom.file.mc.guide = paste0(mc.guide.sbatch,'.',exp.tag)
  unlink(custom.file.mc.guide)
  write(paste0("#!/bin/sh"),custom.file.mc.guide,append=TRUE) #maybe needs newline?
  write(paste0("#SBATCH -n 1"),custom.file.mc.guide,append=TRUE)
  write(paste0("#SBATCH --mail-type=begin"),custom.file.mc.guide,append=TRUE)
  write(paste0("#SBATCH --mail-type=end"),custom.file.mc.guide,append=TRUE)
  write(paste0("#SBATCH -o ", exp.tag,"/log.mc.guide.txt"),custom.file.mc.guide,append=TRUE)
  write(paste0("#SBATCH -p interactive --time=0-10:00:00"),custom.file.mc.guide,append=TRUE)
  write(paste0("#SBATCH --mem=10gb"),custom.file.mc.guide,append=TRUE)
  write(paste0("module load r"),custom.file.mc.guide,append=TRUE)
  write(paste0("cd ",getwd()),custom.file.mc.guide,append=TRUE)
  write(paste0("Rscript ",getwd(),"/mc.guide.R ", 
               exp.tag," ", param.path.sel),custom.file.mc.guide,append=TRUE)
  
  pipeline =  paste0("../../sbatch_scripts/Infile.sel.",exp.tag)
  unlink(pipeline)
  write(custom.file.sel,pipeline,append=TRUE)
  write(custom.file.sel.2,pipeline,append=TRUE)
  write(custom.file.mc.guide,pipeline,append=TRUE)
  system(paste("sbatch-pipeline",pipeline))
  
}else{
  # if not on server just run dataset creation and processing in sequence
  for (m in 1:plotsM){
    run.one.sel(params.sel,m)
  }
  
  #analyze results
  analyze.sel.res(exp.tag,res.dir.sel,res.dir.sel.plots)
}



