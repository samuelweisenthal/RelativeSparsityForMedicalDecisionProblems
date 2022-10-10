# runs one Monte-Carlo result

library(xtable)
library(latex2exp)
source('mc.utils.R')
source('sim.R')

args = commandArgs(trailingOnly=TRUE)
print(args)

dire = args[1]
m = args[2]
params = readRDS(dire)
run.one.sel(params,m)


