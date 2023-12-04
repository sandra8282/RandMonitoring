library(psych); library(dplyr)
library(foreach); library(doParallel); library(tictoc)
load("D:/22ComparisonBaseline/simdat/sim5000_11.RData")
source('D:/22ComparisonBaseline/main_fct.R')
source('D:/22ComparisonBaseline/simAnalysisSOURCE.R')

ncores <- 8 #ncores <- detectCores()
niters <- length(sim1000)

################# ABTC and pvals
tic()
cl <- makeCluster(ncores)
registerDoParallel(cl)
bROC <- foreach(iter = 1:niters, .combine = rbind) %dopar% {
  simBROCanalysis(iter, sim1000 = sim1000)
}
stopCluster(cl)
toc()

tic()
cl <- makeCluster(ncores)
registerDoParallel(cl)
bROC2 <- foreach(iter = 1:niters, .combine = rbind) %dopar% {
  simBROCanalysis2(iter, sim1000 = sim1000)
}
stopCluster(cl)
toc()
tic()

write.csv(bROC, "D:/22ComparisonBaseline/simdat/results5000big11.csv")
write.csv(bROC2, "D:/22ComparisonBaseline/simdat/results5000simp11.csv")

################# curves
tic()
cl <- makeCluster(ncores)
registerDoParallel(cl)
curveROC <- foreach(iter = 1:niters) %dopar% {
  simBROC(iter, sim1000 = sim1000)
}
stopCluster(cl)
toc()

tic()
cl <- makeCluster(ncores)
registerDoParallel(cl)
curveROC2 <- foreach(iter = 1:niters) %dopar% {
  simBROC2(iter, sim1000 = sim1000)
}
stopCluster(cl)
toc()


save(curveROC, file = "D:/22ComparisonBaseline/simdat/curves5000big11.RData")
save(curveROC2, file = "D:/22ComparisonBaseline/simdat/curves5000simp11.RData")
tic()
cl <- makeCluster(ncores)
registerDoParallel(cl)
bysiteROC <- foreach(iter = 1:niters) %dopar% {
  simBROCanalysisBYSITE(iter)
}
stopCluster(cl)
toc()

save(bysiteROC, file = "D:/22ComparisonBaseline/simdat/results5000bysite11.RData")
