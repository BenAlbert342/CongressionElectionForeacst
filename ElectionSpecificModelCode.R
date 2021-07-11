rm(list = ls())

suppressWarnings(library(tidyverse))
suppressWarnings(library(verification))
suppressWarnings(library(car))

Calc.Mcfadden.R.squared<- function(y,mod,x=NULL) {
  nullmod <- glm(y~1, family="binomial")
  LL <- logLik(mod)
  R2 <- 1-LL/logLik(nullmod)
  
  return(R2)
}

Con.Mat <- function(y,y_hat,threshold=0.5) {
  y_hat <- ifelse(y_hat>threshold,1,0)
  Y <- data.frame(y,y_hat)
  Classification <- Y %>%
    mutate(n=1) %>%
    group_by(y,y_hat) %>%
    summarise(n=sum(n)) %>%
    pivot_wider(names_from = y, values_from = n)
  
  return(Classification)
}

set.seed(1234)

HOR <- readRDS('ElectionSpecificVarsHOR.rds')
Senate <- readRDS('ElectionSpecificVarsSenate.rds')

nrow(HOR)
summary(HOR)

nrow(Senate)
summary(Senate)

inds <- sample(1:nrow(HOR),size = round(0.75*nrow(HOR)))
test <- HOR[-inds,]
HOR <- HOR[inds,]

mod.HOR <- glm(Dem~Pres+Mid+Pres.Mid+VS+Incum+PA,HOR,family = binomial())
summary(mod.HOR)

tmp <- data.frame(rep(1,2),c(T,F),c(1,0),rep(mean(HOR$VS),2),rep(0,2),rep(mean(HOR$PA[HOR$Pres==1],2),))
names(tmp) <- c("Pres","Mid","Pres.Mid","VS","Incum","PA")
preds <- predict(mod.HOR,tmp,type = 'response')
preds[1] - preds[2]

tmp <- data.frame(rep(1,2),c(T,T),c(1,1),rep(mean(HOR$VS),2),c(1,0),rep(mean(HOR$PA[HOR$Pres==1],2),))
names(tmp) <- c("Pres","Mid","Pres.Mid","VS","Incum","PA")
preds <- predict(mod.HOR,tmp,type = 'response')
preds[1] - preds[2]

Calc.Mcfadden.R.squared(HOR$Dem,mod.HOR)

Con.Mat(HOR$Dem,predict(mod.HOR,type = 'response'))

roc.area(HOR$Dem,predict(mod.HOR,type = 'response'))
roc.plot(HOR$Dem,predict(mod.HOR,type = 'response'))

y_hat <- predict(mod.HOR,test,type = "response")
Con.Mat(test$Dem,y_hat)

roc.area(test$Dem,y_hat)
roc.plot(test$Dem,y_hat) 

inds <- sample(1:nrow(Senate),size = round(0.75*nrow(Senate)))
test <- Senate[-inds,]
Senate <- Senate[inds,]

mod.Senate <- glm(Dem~Pres+Mid+Pres.Mid+VS+Incum+PA,Senate,family = binomial())
summary(mod.Senate)

linearHypothesis(mod.Senate,c("MidTRUE = 0","Pres.Mid = 0"))

tmp <- data.frame(rep(1,2),c(T,F),c(1,0),rep(mean(Senate$VS),2),rep(0,2),rep(mean(Senate$PA[Senate$Pres==1],2),))
names(tmp) <- c("Pres","Mid","Pres.Mid","VS","Incum","PA")
preds <- predict(mod.Senate,tmp,type = 'response')
preds[1] - preds[2]

tmp <- data.frame(rep(1,2),c(T,T),c(1,1),rep(mean(Senate$VS),2),c(1,0),rep(mean(Senate$PA[Senate$Pres==1],2),))
names(tmp) <- c("Pres","Mid","Pres.Mid","VS","Incum","PA")
preds <- predict(mod.Senate,tmp,type = 'response')
preds[1] - preds[2]

Calc.Mcfadden.R.squared(Senate$Dem,mod.Senate)

Con.Mat(Senate$Dem,predict(mod.Senate,type = 'response'))

roc.area(Senate$Dem,predict(mod.Senate,type = 'response'))
roc.plot(Senate$Dem,predict(mod.Senate,type = 'response')) 

y_hat <- predict(mod.Senate,test,type = "response")
Con.Mat(test$Dem,y_hat)

roc.area(test$Dem,y_hat)
roc.plot(test$Dem,y_hat) 




