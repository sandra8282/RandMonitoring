library(psych); library(dplyr); library(randomizr); library(tictoc)

#### function to get randomization schedule
rand_stratum <- function(stratum, n, levels, blk_sizes) {
  library(data.table)
  blk_sizes <- blk_sizes / 2
  dB <- data.table(blockrand::blockrand(
    n = n, num.levels = 2, levels = levels,
    id.prefix = stratum, block.prefix  = stratum,
    stratum = stratum, block.sizes = blk_sizes)
  )
  return(dB)
}

getrandsch <- function(nrand){
  schedule <- list()
  for (i in 1:nrow(nrand)){
    comb = nrand$comb[i]
    s = nrand$site[i]
    n = nrand$n[i]
    if (n==0){schedule[[i]]=NA
    } else {schedule[[i]] <- rand_stratum(stratum=i, n=n, levels=c(0, 1), 
                                            blk_sizes=c(2,4,6))} #block size 2,6,8
    names(schedule)[[i]] = i
  }
  return(schedule)
}

# blockrand::blockrand(n = max(nrand$n[i],2), 
#                                     num.levels = 2, # 2 treatments
#                                     levels = c(0, 1), # arm names
#                                     stratum = paste(i, comb, sep = "."), # stratum name
#                                     id.prefix = paste(i, ".", sep=""), # stratum abbrev
#                                     block.sizes = c(1,2,3), # 2,4,6
#                                     block.prefix = paste(i, ".", sep=""),
#                                     uneq.beg=TRUE) # stratum abbrev

#### function to create sample and then call randomizer function
simfct1 <- function(seed, iter, Ns){
  
  library(dplyr);
    set.seed(seed)
    site = data.frame(t(rmultinom(Ns, 1, c(0.4, 0.3, 0.3))))
    sampledat <- data.frame(site = ifelse(site$X1==1, 1, ifelse(site$X2==1,2,3)),
                    x1 = rnorm(Ns, 20 + 2.5*site$X1 + 0.5*site$X2, 2),
                    x2 = rbinom(Ns, 1, 0.8 - 0.3*site$X1 - 0.15*site$X2),
                    x3 = rbinom(Ns, 1, 0.6 - 0.1*site$X1 - 0.05*site$X2))
    sampledat = sampledat %>% mutate(randnum = 1:n()) 
    sampledat = sampledat %>% group_by(site, x3) %>% mutate(ordernum = 1:n())
    sampledat = data.frame(sampledat)
    sampledat$x1cat <- ifelse(sampledat$x1>21, 1, 0)
    sampledat$comb1 = ifelse(sampledat$x1cat ==1 & sampledat$x2==1 &sampledat$x3 == 1, 1, 
                             ifelse(sampledat$x1cat ==1 & sampledat$x2==1 &sampledat$x3 == 0, 2,
                                    ifelse(sampledat$x1cat ==1 & sampledat$x2==0 &sampledat$x3 == 1, 3,
                                           ifelse(sampledat$x1cat ==1 & sampledat$x2==0 &sampledat$x3 == 0, 4,
                                                  ifelse(sampledat$x1cat ==0 & sampledat$x2==1 &sampledat$x3 == 1, 5, 
                                                         ifelse(sampledat$x1cat ==0 & sampledat$x2==1 &sampledat$x3 == 0, 6,
                                                                ifelse(sampledat$x1cat ==0 & sampledat$x2==0 &sampledat$x3 == 1, 7, 8)
                                                         ))))))
    sampledat$schedule = ifelse(sampledat$site==1, sampledat$comb1, sampledat$comb1+8*(sampledat$site-1))
    nrand <- data.frame(table(sampledat$comb1, sampledat$site))
    colnames(nrand) = c("comb", "site", "n")
    nrand$comb = as.integer(as.character(nrand$comb))
    nrand$site = as.integer(as.character(nrand$site))
    #n4rand = (min(nrand$n))
  
  schedule = getrandsch(nrand=nrand)  ### Randomization schedules
  s_count = sort(unique(sampledat$schedule))
  newsdat = NULL
  for (schn in s_count){
    temp = sampledat %>% filter(schedule == schn)
    trtsch = schedule[[schn]]
    temp$trt = trtsch$treatment[1:nrow(temp)]
    temp$block.id = trtsch$block.id[1:nrow(temp)]
    temp$block.size = trtsch$block.size[1:nrow(temp)]
    temp$block.id = paste(schn, temp$block.id, sep=".")
    t2 = temp %>% select(!c(x1cat, comb1, schedule))
    newsdat = rbind(newsdat, t2)
  }
  sampledat <- newsdat
  sampledat$trt_error = sampledat$trt_error_mild = sampledat$trt
  sampledat <- sampledat %>% arrange(site, ordernum)
  inds = which(sampledat$site==2 & sampledat$ordernum >= 10)
  while(round(length(inds)/Ns,2)>0.1){
    out = c(1,length(inds))
    inds = inds[-sample(out,1)]
  }
  inds2=inds[1:(length(inds)/2)]
  while(round(length(inds2)/Ns,2)>0.05){
    inds2 = inds2[-length(inds2)]
  }
  sampledat$trt_error[inds] = 0
  sampledat$trt_error_mild[inds2] = 0
  return(sampledat)
}  