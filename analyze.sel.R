# analyzes one Monte-Carlo result

library(gtools)
library(xtable)
library(latex2exp)
source('mc.utils.R')
source('utils.R')
# needs to go in sep file
args = commandArgs(trailingOnly=TRUE)
print("args begin")
print(args)
print("args end")
#if (!identical(args,character(0))){
res.dir = args[1]
sel.res.dir = args[2]
plot.dir = args[3]
#}
analyze.sel.res(res.dir,sel.res.dir,plot.dir,rdigits=3)
