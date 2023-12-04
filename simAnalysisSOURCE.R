
########################
##### model with site
########################

#### ABCD
simBROCanalysis <- function(counti, sim1000){
  trialdat <- sim1000[[counti]]
  #Calculate propensity scores and logit
  m1 <- glm(trt ~ x1 + x2 + x3 + factor(site), family = binomial, data = trialdat)
  m2 <- glm(trt_error_mild ~ x1 + x2+ x3 + factor(site), family = binomial, data = trialdat)
  m3 <- glm(trt_error ~ x1 + x2+ x3 + factor(site), family = binomial, data = trialdat)
  trialdat$p <- predict(m1, newdata = trialdat, type = 'response')
  trialdat$p2 <- predict(m2, newdata = trialdat, type = 'response')
  trialdat$p3 <- predict(m3, newdata = trialdat, type = 'response')
  trialdat$xp <- predict(m1, newdata = trialdat)
  trialdat$xp2 <- predict(m2, newdata = trialdat)
  trialdat$xp3 <- predict(m3, newdata = trialdat)
  #obtain estimates for good randomization
  cd <- psych::cohen.d(trialdat$xp, trialdat$trt)
  r <- baseROC(t=trialdat$xp[trialdat$trt==1], 
               c=trialdat$xp[trialdat$trt==0], 
               AUC=FALSE, silent=TRUE)
  r1 <- c(cd1 = cd$cohen.d[2], 
          pval1 = t.test(trialdat$xp[trialdat$trt==0], 
                         trialdat$xp[trialdat$trt==1])$p.value, 
          ABTCD1 = r$ABTCD)
  #obtain estimates for randomization with a very mild error
  cd <- psych::cohen.d(trialdat$xp2, trialdat$trt_error_mild)
  r <- baseROC(t=trialdat$xp2[trialdat$trt_error_mild==1], 
               c=trialdat$xp2[trialdat$trt_error_mild==0], 
               AUC=FALSE, silent=TRUE)
  r2 <- c(cd2 = cd$cohen.d[2], 
          pval2 = t.test(trialdat$xp2[trialdat$trt_error_mild==0], 
                         trialdat$xp2[trialdat$trt_error_mild==1])$p.value, 
          ABTCD2 = r$ABTCD)
  #obtain estimates for randomization with a gross error
  cd <- psych::cohen.d(trialdat$xp3, trialdat$trt_error)
  r <- baseROC(t=trialdat$xp3[trialdat$trt_error==1], 
               c=trialdat$xp3[trialdat$trt_error==0], AUC=FALSE, silent=TRUE)
  r3 <- c(cd3 = cd$cohen.d[2], 
          pval3 = t.test(trialdat$xp3[trialdat$trt_error==0], 
                         trialdat$xp3[trialdat$trt_error==1])$p.value,
          ABTCD3 = r$ABTCD)
  return(r=c(r1, r2, r3))
}

#ROC
simBROC <- function(counti, sim1000){
  trialdat <- sim1000[[counti]]
  #Calculate propensity scores and logit
  m1 <- glm(trt ~ x1 + x2 + x3 + factor(site), family = binomial, data = trialdat)
  m2 <- glm(trt_error_mild ~ x1 + x2+ x3 + factor(site), family = binomial, data = trialdat)
  m3 <- glm(trt_error ~ x1 + x2+ x3 + factor(site), family = binomial, data = trialdat)
  trialdat$p <- predict(m1, newdata = trialdat, type = 'response')
  trialdat$p2 <- predict(m2, newdata = trialdat, type = 'response')
  trialdat$p3 <- predict(m3, newdata = trialdat, type = 'response')
  trialdat$xp <- predict(m1, newdata = trialdat)
  trialdat$xp2 <- predict(m2, newdata = trialdat)
  trialdat$xp3 <- predict(m3, newdata = trialdat)
  #obtain estimates for good randomization
  r <- baseROC(t=trialdat$xp[trialdat$trt==1], 
               c=trialdat$xp[trialdat$trt==0], AUC=FALSE, silent=TRUE)
  ROC1 <- r$ROC
  #obtain estimates for randomization with a very mild error
  r <- baseROC(t=trialdat$xp2[trialdat$trt_error_mild==1],
               c=trialdat$xp2[trialdat$trt_error_mild==0], AUC=FALSE, silent=TRUE)
  ROC2 <- r$ROC
  #obtain estimates for randomization with a gross error
  r <- baseROC(t=trialdat$xp3[trialdat$trt_error==1],
               c= trialdat$xp3[trialdat$trt_error==0], AUC=FALSE, silent=TRUE)
  ROC3 <- r$ROC
  return(ROC = list(ROC1, ROC2, ROC3))
}

########################
##### model without site
########################

#### ABCD
simBROCanalysis2 <- function(counti, sim1000){
  trialdat <- sim1000[[counti]]
  #Calculate propensity scores and logit
  m1 <- glm(trt ~ x1 + x2 + x3, family = binomial, data = trialdat)
  m2 <- glm(trt_error_mild ~ x1 + x2+ x3, family = binomial, data = trialdat)
  m3 <- glm(trt_error ~ x1 + x2+ x3, family = binomial, data = trialdat)
  trialdat$p <- predict(m1, newdata = trialdat, type = 'response')
  trialdat$p2 <- predict(m2, newdata = trialdat, type = 'response')
  trialdat$p3 <- predict(m3, newdata = trialdat, type = 'response')
  trialdat$xp <- predict(m1, newdata = trialdat)
  trialdat$xp2 <- predict(m2, newdata = trialdat)
  trialdat$xp3 <- predict(m3, newdata = trialdat)
  #obtain estimates for good randomization
  cd <- psych::cohen.d(trialdat$xp, trialdat$trt)
  r <- baseROC(t=trialdat$xp[trialdat$trt==1], 
               c=trialdat$xp[trialdat$trt==0], 
               AUC=FALSE, silent=TRUE)
  r1 <- c(cd1 = cd$cohen.d[2], 
          pval1 = t.test(trialdat$xp[trialdat$trt==0], 
                         trialdat$xp[trialdat$trt==1])$p.value, 
          ABTCD1 = r$ABTCD)
  #obtain estimates for randomization with a very mild error
  cd <- psych::cohen.d(trialdat$xp2, trialdat$trt_error_mild)
  r <- baseROC(t=trialdat$xp2[trialdat$trt_error_mild==1], 
               c=trialdat$xp2[trialdat$trt_error_mild==0], 
               AUC=FALSE, silent=TRUE)
  r2 <- c(cd2 = cd$cohen.d[2], 
          pval2 = t.test(trialdat$xp2[trialdat$trt_error_mild==0], 
                         trialdat$xp2[trialdat$trt_error_mild==1])$p.value, 
          ABTCD2 = r$ABTCD)
  #obtain estimates for randomization with a gross error
  cd <- psych::cohen.d(trialdat$xp3, trialdat$trt_error)
  r <- baseROC(t=trialdat$xp3[trialdat$trt_error==1], 
               c=trialdat$xp3[trialdat$trt_error==0], AUC=FALSE, silent=TRUE)
  r3 <- c(cd3 = cd$cohen.d[2], 
          pval3 = t.test(trialdat$xp3[trialdat$trt_error==0], 
                         trialdat$xp3[trialdat$trt_error==1])$p.value,
          ABTCD3 = r$ABTCD)
  return(r=c(r1, r2, r3))
}

#ROC
simBROC2 <- function(counti, sim1000){
  trialdat <- sim1000[[counti]]
  #Calculate propensity scores and logit
  m1 <- glm(trt ~ x1 + x2 + x3, family = binomial, data = trialdat)
  m2 <- glm(trt_error_mild ~ x1 + x2+ x3, family = binomial, data = trialdat)
  m3 <- glm(trt_error ~ x1 + x2+ x3, family = binomial, data = trialdat)
  trialdat$p <- predict(m1, newdata = trialdat, type = 'response')
  trialdat$p2 <- predict(m2, newdata = trialdat, type = 'response')
  trialdat$p3 <- predict(m3, newdata = trialdat, type = 'response')
  trialdat$xp <- predict(m1, newdata = trialdat)
  trialdat$xp2 <- predict(m2, newdata = trialdat)
  trialdat$xp3 <- predict(m3, newdata = trialdat)
  #obtain estimates for good randomization
  r <- baseROC(t=trialdat$xp[trialdat$trt==1], 
               c=trialdat$xp[trialdat$trt==0], AUC=FALSE, silent=TRUE)
  ROC1 <- r$ROC
  #obtain estimates for randomization with a very mild error
  r <- baseROC(t=trialdat$xp2[trialdat$trt_error_mild==1],
               c=trialdat$xp2[trialdat$trt_error_mild==0], AUC=FALSE, silent=TRUE)
  ROC2 <- r$ROC
  #obtain estimates for randomization with a gross error
  r <- baseROC(t=trialdat$xp3[trialdat$trt_error==1],
               c= trialdat$xp3[trialdat$trt_error==0], AUC=FALSE, silent=TRUE)
  ROC3 <- r$ROC
  return(ROC = list(ROC1, ROC2, ROC3))
}

##################################################################
#### Stratified by site
##################################################################

simBROCanalysisBYSITE <- function(counti){
  trialdat <- sim1000[[counti]]
  trialdat1 <- subset(trialdat, site==1)
  trialdat2 <- subset(trialdat, site==2)
  trialdat3 <- subset(trialdat, site==3)
  #No errors by site
  m11s <- glm(trt ~ x1 + x2 + x3, family = binomial, data = trialdat1)
  m12s <- glm(trt ~ x1 + x2 + x3, family = binomial, data = trialdat2)
  m13s <- glm(trt ~ x1 + x2 + x3, family = binomial, data = trialdat3)
  trialdat1$p1 <- predict(m11s, newdata = trialdat1, type = 'response')
  trialdat1$xp1 <- predict(m11s, newdata = trialdat1)
  trialdat2$p1 <- predict(m12s, newdata = trialdat2, type = 'response')
  trialdat2$xp1 <- predict(m12s, newdata = trialdat2)
  trialdat3$p1 <- predict(m13s, newdata = trialdat3, type = 'response')
  trialdat3$xp1 <- predict(m13s, newdata = trialdat3)
  r <- list(r11 = baseROC(t=trialdat1$xp1[trialdat1$trt==1], 
                          c=trialdat1$xp1[trialdat1$trt==0], 
                          AUC=FALSE, silent=TRUE),
            r12 = baseROC(t=trialdat2$xp1[trialdat2$trt==1], 
                          c=trialdat2$xp1[trialdat2$trt==0], 
                          AUC=FALSE, silent=TRUE),
            r13 = baseROC(t=trialdat3$xp1[trialdat3$trt==1], 
                          c=trialdat3$xp1[trialdat3$trt==0], 
                          AUC=FALSE, silent=TRUE))
  #obtain estimates for good randomization
  cds <- c(psych::cohen.d(trialdat1$xp1, trialdat1$trt)$cohen.d[2],
           psych::cohen.d(trialdat2$xp1, trialdat2$trt)$cohen.d[2],
           psych::cohen.d(trialdat3$xp1, trialdat3$trt)$cohen.d[2])
  pvals <- c(t.test(trialdat1$xp1[trialdat1$trt==0], 
                    trialdat1$xp1[trialdat1$trt==1])$p.value,
             t.test(trialdat2$xp1[trialdat2$trt==0], 
                    trialdat2$xp1[trialdat2$trt==1])$p.value,
             t.test(trialdat3$xp1[trialdat3$trt==0], 
                    trialdat3$xp1[trialdat3$trt==1])$p.value)
  result_noerror <- list(cds_ne = cds, pvals_ne = pvals, r_ne = r)
  
  #Mild error by site
  m21s <- glm(trt_error_mild ~ x1 + x2 + x3, family = binomial, data = trialdat1)
  m22s <- glm(trt_error_mild ~ x1 + x2 + x3, family = binomial, data = trialdat2)
  m23s <- glm(trt_error_mild ~ x1 + x2 + x3, family = binomial, data = trialdat3)
  trialdat1$p2 <- predict(m21s, newdata = trialdat1, type = 'response')
  trialdat1$xp2 <- predict(m21s, newdata = trialdat1)
  trialdat2$p2 <- predict(m22s, newdata = trialdat2, type = 'response')
  trialdat2$xp2 <- predict(m22s, newdata = trialdat2)
  trialdat3$p2 <- predict(m23s, newdata = trialdat3, type = 'response')
  trialdat3$xp2 <- predict(m23s, newdata = trialdat3)
  r <- list(r11 = baseROC(t=trialdat1$xp2[trialdat1$trt_error_mild==1], 
                          c=trialdat1$xp2[trialdat1$trt_error_mild==0], 
                          AUC=FALSE, silent=TRUE),
            r12 = baseROC(t=trialdat2$xp2[trialdat2$trt_error_mild==1], 
                          c=trialdat2$xp2[trialdat2$trt_error_mild==0], 
                          AUC=FALSE, silent=TRUE),
            r13 = baseROC(t=trialdat3$xp2[trialdat3$trt_error_mild==1], 
                          c=trialdat3$xp2[trialdat3$trt_error_mild==0], 
                          AUC=FALSE, silent=TRUE))
  #obtain estimates for good randomization
  cds <- c(psych::cohen.d(trialdat1$xp2, trialdat1$trt_error_mild)$cohen.d[2],
           psych::cohen.d(trialdat2$xp2, trialdat2$trt_error_mild)$cohen.d[2],
           psych::cohen.d(trialdat3$xp2, trialdat3$trt_error_mild)$cohen.d[2])
  pvals <- c(t.test(trialdat1$xp2[trialdat1$trt_error_mild==0], 
                    trialdat1$xp2[trialdat1$trt_error_mild==1])$p.value,
             t.test(trialdat2$xp2[trialdat2$trt_error_mild==0], 
                    trialdat2$xp2[trialdat2$trt_error_mild==1])$p.value,
             t.test(trialdat3$xp2[trialdat3$trt_error_mild==0], 
                    trialdat3$xp2[trialdat3$trt_error_mild==1])$p.value)
  result_milderror <- list(cds_me = cds, pvals_me = pvals, r_me = r)
  
  #Large error by site
  m31s <- glm(trt_error ~ x1 + x2 + x3, family = binomial, data = trialdat1)
  m32s <- glm(trt_error ~ x1 + x2 + x3, family = binomial, data = trialdat2)
  m33s <- glm(trt_error ~ x1 + x2 + x3, family = binomial, data = trialdat3)
  trialdat1$p3 <- predict(m31s, newdata = trialdat1, type = 'response')
  trialdat1$xp3 <- predict(m31s, newdata = trialdat1)
  trialdat2$p3 <- predict(m32s, newdata = trialdat2, type = 'response')
  trialdat2$xp3 <- predict(m32s, newdata = trialdat2)
  trialdat3$p3 <- predict(m33s, newdata = trialdat3, type = 'response')
  trialdat3$xp3 <- predict(m33s, newdata = trialdat3)
  r <- list(r11 = baseROC(t=trialdat1$xp3[trialdat1$trt_error==1], 
                          c=trialdat1$xp3[trialdat1$trt_error==0], AUC=FALSE, silent=TRUE),
            r12 = baseROC(t=trialdat2$xp3[trialdat2$trt_error==1], 
                          c=trialdat2$xp3[trialdat2$trt_error==0], AUC=FALSE, silent=TRUE),
            r13 = baseROC(t=trialdat3$xp3[trialdat3$trt_error==1], 
                          c=trialdat3$xp3[trialdat3$trt_error==0], AUC=FALSE, silent=TRUE))
  #obtain estimates for good randomization
  cds <- c(psych::cohen.d(trialdat1$xp3, trialdat1$trt_error)$cohen.d[2],
           psych::cohen.d(trialdat2$xp3, trialdat2$trt_error)$cohen.d[2],
           psych::cohen.d(trialdat3$xp3, trialdat3$trt_error)$cohen.d[2])
  pvals <- c(t.test(trialdat1$xp3[trialdat1$trt_error==0], 
                    trialdat1$xp3[trialdat1$trt_error==1])$p.value,
             t.test(trialdat2$xp3[trialdat2$trt_error==0], 
                    trialdat2$xp3[trialdat2$trt_error==1])$p.value,
             t.test(trialdat3$xp3[trialdat3$trt_error==0], 
                    trialdat3$xp3[trialdat3$trt_error==1])$p.value)
  result_lerror <- list(cds_le = cds, pvals_le = pvals, r_le = r)
  return(r=list(result_noerror, result_milderror, result_lerror))
}


 
# 
# simBROCanalysis3 <- function(counti){
#   trialdat <- sim1000[[counti]]
#   #Calculate propensity scores and logit
#   m1 <- glm(trt ~ (x1 + x2 + x3)*factor(site), family = binomial, data = trialdat)
#   m2 <- glm(trt_error_mild ~ (x1 + x2 + x3)*factor(site), family = binomial, data = trialdat)
#   m3 <- glm(trt_error ~ (x1 + x2 + x3)*factor(site), family = binomial, data = trialdat)
#   trialdat$p <- predict(m1, newdata = trialdat, type = 'response')
#   trialdat$p2 <- predict(m2, newdata = trialdat, type = 'response')
#   trialdat$p3 <- predict(m3, newdata = trialdat, type = 'response')
#   trialdat$xp <- predict(m1, newdata = trialdat)
#   trialdat$xp2 <- predict(m2, newdata = trialdat)
#   trialdat$xp3 <- predict(m3, newdata = trialdat)
#   #obtain estimates for good randomization
#   cd <- psych::cohen.d(trialdat$xp, trialdat$trt)
#   r <- baseROC(t=trialdat$xp[trialdat$trt==1], 
#                c=trialdat$xp[trialdat$trt==0], 
#                AUC=FALSE, silent=TRUE)
#   r1 <- c(cd1 = cd$cohen.d[2], 
#           pval1 = t.test(trialdat$xp[trialdat$trt==0], 
#                          trialdat$xp[trialdat$trt==1])$p.value, 
#           ABTCD1 = r$ABTCD)
#   #obtain estimates for randomization with a very mild error
#   cd <- psych::cohen.d(trialdat$xp2, trialdat$trt_error_mild)
#   r <- baseROC(t=trialdat$xp2[trialdat$trt_error_mild==1], 
#                c=trialdat$xp2[trialdat$trt_error_mild==0], 
#                AUC=FALSE, silent=TRUE)
#   r2 <- c(cd2 = cd$cohen.d[2], 
#           pval2 = t.test(trialdat$xp2[trialdat$trt_error_mild==0], 
#                          trialdat$xp2[trialdat$trt_error_mild==1])$p.value, 
#           ABTCD2 = r$ABTCD)
#   #obtain estimates for randomization with a gross error
#   cd <- psych::cohen.d(trialdat$xp3, trialdat$trt_error)
#   r <- baseROC(t=trialdat$xp3[trialdat$trt_error==1], 
#                c=trialdat$xp3[trialdat$trt_error==0], AUC=FALSE, silent=TRUE)
#   r3 <- c(cd3 = cd$cohen.d[2], 
#           pval3 = t.test(trialdat$xp3[trialdat$trt_error==0], 
#                          trialdat$xp3[trialdat$trt_error==1])$p.value,
#           ABTCD3 = r$ABTCD)
#   return(r=c(r1, r2, r3))
# }
# 
# simBROC3 <- function(counti){
#   trialdat <- sim1000[[counti]]
#   #Calculate propensity scores and logit
#   m1 <- glm(trt ~ (x1 + x2 + x3)*factor(site), family = binomial, data = trialdat)
#   m2 <- glm(trt_error_mild ~ (x1 + x2 + x3)*factor(site), family = binomial, data = trialdat)
#   m3 <- glm(trt_error ~ (x1 + x2 + x3)*factor(site), family = binomial, data = trialdat)
#   trialdat$p <- predict(m1, newdata = trialdat, type = 'response')
#   trialdat$p2 <- predict(m2, newdata = trialdat, type = 'response')
#   trialdat$p3 <- predict(m3, newdata = trialdat, type = 'response')
#   trialdat$xp <- predict(m1, newdata = trialdat)
#   trialdat$xp2 <- predict(m2, newdata = trialdat)
#   trialdat$xp3 <- predict(m3, newdata = trialdat)
#   #obtain estimates for good randomization
#   r <- baseROC(t=trialdat$xp[trialdat$trt==1], 
#                c=trialdat$xp[trialdat$trt==0], AUC=FALSE, silent=TRUE)
#   ROC1 <- r$ROC
#   #obtain estimates for randomization with a very mild error
#   r <- baseROC(t=trialdat$xp2[trialdat$trt_error_mild==1],
#                c=trialdat$xp2[trialdat$trt_error_mild==0], AUC=FALSE, silent=TRUE)
#   ROC2 <- r$ROC
#   #obtain estimates for randomization with a gross error
#   r <- baseROC(t=trialdat$xp3[trialdat$trt_error==1],
#                c= trialdat$xp3[trialdat$trt_error==0], AUC=FALSE, silent=TRUE)
#   ROC3 <- r$ROC
#   return(ROC = list(ROC1, ROC2, ROC3))
# }