library(foreach); library(doParallel)
library(tictoc)
source('D:/22ComparisonBaseline/simNEWsource.R')

#Run simulations in parallel on <ncores> cores
ncores <- 8 #ncores <- detectCores()
seed <- 1979

########################################### SIM 1  ###########################################
Ns=5000
niters <- 1000
set.seed(seed);
seeds <- seed + sample(-2^25:2^25, niters)

tic()
cl <- makeCluster(ncores)
registerDoParallel(cl)
sim1000 <- foreach(iter = 1:niters) %dopar% {
  set.seed(seeds[iter])
  simfct1(seed=seeds[iter], iter=iter, Ns=Ns)
}
stopCluster(cl)
toc()

for (i in 1:niters) {
  tab <- data.frame(with(sim1000[[i]], table(block.id, block.size)))
  tab <- tab %>% filter(!(Freq %in% c(0, 2, 4, 6)))
}
(nrow(tab)/1000)*100

save(sim1000, file="D:/22ComparisonBaseline/simdat/sim5000_11.RData")
