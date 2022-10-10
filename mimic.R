# this file runs the real data analysis.
# there is also a mimic.prelim.analysis.Rmd file in this directory, 
# which runs more in depth
# and annotated preliminary analysis.  some preliminary analyses are 
# done here
# just as "unit tests" to make sure everything ok, but otherwise
# use this
# this file processes the outputs of np_load.py, in particular the stages
# multistage_rewards, multistage_actions, etc
# these are censored

library(latex2exp)
library(kableExtra)
library(xtable)

# some useful functions
source('utils.R')
source('mc.utils.R')
source('is.R')

# if debug, will just take eg first 300 episodes
# and will set fewer lambda, etc
debug = 1

# where are episodes, made by data_clean.py, np_load.py
stagesdirec = paste0("/stages")

# this tag just names output files (output will be an R object with summary data,
# a copy of this configuration file, and a plot)
str.tag = "mimic"

# miscellaneous arguments
other.args=list(logx=0,
                gammaixleg=1,deltaixleg=1,area.plots=1.3,
                plot.emp.var=1,
                max.n.div=2,use.diff=.5)

# threshold Delta for considering equal to behavior
use.diff = other.args$use.diff #1e-1

# lambda grid
min.lambda = 0.3
max.lambda = 1.4#.02
nlam = 10 #10#10
lambdas= seq(min.lambda,max.lambda,length.out=nlam)

if (debug){
  other.args=list(#get.theor.var.of.ests=0,plot.adap.lasso.ci=0,
    logx=0,
    gammaixleg=1,deltaixleg=1,area.plots=1.3,
    plot.emp.var=1,
    max.n.div=2,use.diff=.5)
    nlam = 2#5#3 #10#10
    lambdas= exp(seq(log(1e-5),log(2000),length.out=nlam))
  
}

nlam = length(lambdas)

ncov = 9 # 8 is max
# these are getting the correct stages, since we have some initial
# stages that are like 40 seconds just to collect other measurements
# from before episode starts

start.t=3 # beginning of traj
end.stage=3 # end of traj. we end early. so this is like first 30 min after hyp

other.arg.str= paste0(paste0(names(other.args),"=",other.args),collapse = "")
exp.tag = paste0("tag=",str.tag,"usediff=",use.diff,",start.t",start.t,"endStage",end.stage,"gammaselix=",
                 nlam,"minlam=",min(lambdas),"maxlam=",max(lambdas),
                 "ncov=",ncov,"use.dff=",use.diff)

print("Need to make sure in proper directory!!!")

# need to censor names
if(getwd()=="/Volumes/projects/Latent/anonymous/dir/relative_sparsity/code_betareg/relative_sparsity"){
  pth='/Users/anonymous/Box/MIMIC/mimic-iii-v1-4/hypotension-RL/model-data2'  
}else{
  print("stop")
  browser()
  pth='/scratch/sweisent/mimic'
  
  sink(file = paste0(exp.tag,".mimic.txt"), append = FALSE, 
       type = c("output", "message"),split = FALSE)
}

writeLines(readLines("mimic.R"))

# this is just an id for a unit test (checking that the data for this id is as expected).  
# This is done with much more description
# in the mimic.prelim.analysis.Rmd file, but we do it here just to be extra careful
id.check = 52

# read in stages
filenames <- list.files(paste0(pth,stagesdirec), pattern="*.csv", full.names=TRUE)
ldf <- lapply(filenames, read.csv)
ldf = lapply(ldf,as.matrix)
# note
state.dim = dim(ldf[[1]])[2] #K
nObs = dim(ldf[[1]])[1] #N

# rewards, actions, times (times aren't used, but useful to know them)
msr =  read.csv(paste0(pth,'/','multistage_rewards.csv'),header=TRUE, row.names = 1)
msa =  read.csv(paste0(pth,'/','multistage_actions.csv'),header=TRUE, row.names = 1)
mst =  read.csv(paste0(pth,'/','multistage_times.csv'),header=TRUE, row.names = 1)

# unit check should be 53,1,0.5067
msr[id.check,]
msa[id.check,]
mst[id.check,] # note that these times correspond to the rewards and actions, not states

# should be 48
for (i in start.t:end.stage){
print(ldf[[i]][id.check,]['map'])
}

# get dropout. should be 7
misa = rep(0,length(msa[,1]))

for (i in start.t:(end.stage)){
  misa = misa + is.na(msa[,i])
}

is.misa=misa>0

# don't take if dropout
take = !is.misa 

ixs = 1:dim(msa)[1]
newixs = ixs[take]


for (i in 1:length(ldf)){
  print(dim(ldf[[i]]))
}

# just take decisions betwwen start.t and end.stage
# starting at time start.t
# will be one stage

msr = msr[,start.t:(end.stage),drop=FALSE]
msa = msa[,start.t:(end.stage),drop=FALSE]
ldf = ldf[start.t:(end.stage)]
mst = mst[,start.t:(end.stage),drop=FALSE]

msr=msr[take,,drop=FALSE]
msa = msa[take,,drop=FALSE]
mst = mst[take,,drop=FALSE]
ldf=lapply(ldf,function(x){x[take,,drop=FALSE]})


# again, unit test, after doing all selection
id.check = which(newixs==id.check)
# should be 53
msr[id.check,]
# should be 1
msa[id.check,]

# should be 0.5067
mst[id.check,] # note that these times correspond to the rewards and actions, not states

# should be 48
for (i in 1:length(ldf)){
  print(ldf[[i]][id.check,]['map'])
}


# Joe’s paper  uses 9 variables, so we follow him. 
#He uses MAP, heart rate, urine output, lactate, Glasgow coma score, 
# serum creatinine, FiO2, total bilirubin, and platelet count. 

cov.of.int = c("map","hr","urine","lactate","GCS","creatinine","fio2",
               "bilirubin_total","platelets")#,"total_all_prev_fluids")#,"total_all_prev_vasos")#[1:ncov]

ldf2=lapply(ldf,function(x){x[,cov.of.int,drop=FALSE]})

# just make acronyms upper case, etc
renamed.cov.of.int  =  c("MAP","HR","urine","lactate","GCS","creatinine","Fio2",
                 "bilirubin","platelets")#,"past fluids")#,"past vasos")[1:ncov]

for (i in 1:length(ldf2)){
  colnames(ldf2[[i]])= cov.of.int
  
}

# put episodes into list rather than in arrays, list is processed by model
mseps = get.eps.ms(ldf2[[1]],msa,msr)

# do a positivity check. This is also done in prelim analysis,
# but do this here as well.
# note that it is pos over pi and 1-pi, so average is 0.5
check.pos(mseps,cov.of.int)

## just take 300 episodes if we are debugging
if (debug){
  eps.mimic = mseps[1:300]
}else{
  eps.mimic=mseps #mseps#[1:1000]#halfds
}

# scale states 
sc=center.scale.s(eps.mimic)
smm=sc$sm.cs

# tell the model to scale as well
scale.s = 1 # for real data, must scale

# for mimic, b0 doesn't matter, so just use NA. dimensions might be taken from this
ixs=select.and.est(eps.mimic,b0=rep(NA,dim(eps.mimic[[1]]$Ss)[1]),
                   lambdas=lambdas,
                   names=renamed.cov.of.int,
                   resfile=paste0(exp.tag,"res"),
                   plotfile=paste0(exp.tag,"grid.mimic.png"),
                   scale.s=scale.s,
                   use.diff=use.diff,other.args=other.args)



