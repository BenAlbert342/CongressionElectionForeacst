rm(list=ls())

library(tidyverse)
library(xts)
library(zoo)
library(forecast)
library(mgcv)
library(lubridate)

setwd('C:/Users/15307/Desktop/CongressionalElelctionSimulator/Data')

fit.gam.arima <- function(df) {
  GAM <- gam(Approval~s(Years.In.Office,bs = "cr",fx=T)+s(rgdp.pc,bs = "cr",fx=T)+s(unemp,bs = "cr",fx=T)+s(cpi,bs = "cr",fx=T)+
               PopularVote+Elected,data=df,family = gaussian)
  df$res <- residuals(GAM)
  arima <- auto.arima(df$res)
  df$arima.fit <- as.numeric(fitted(arima))
  GAM <- gam(Approval~s(Years.In.Office,bs = "cr",fx=T)+s(rgdp.pc,bs = "cr",fx=T)+s(unemp,bs = "cr",fx=T)+s(cpi,bs = "cr",fx=T)+
               PopularVote+Elected+arima.fit,data=df,family =gaussian)
  mod <- list(GAM,arima)
  names(mod) <- c('GAM','ARIMA')
  
  return(mod)
}

Plot.Trend <- function(df,mod) {
  plot(df$Approval)
  df$trend.line <- as.numeric(predict(mod))
  lines(df$trend.line,col='red')
  addLegend(legend.loc = "bottomleft", legend.names = c("Approval","Forecast"), 
            pch=15,col = c("black","red"))
}

predict.gam.arima <- function(df,mod,k=1,sd.cap=F,lb,ub,se=F) {
  arima.forecast <- as.numeric(forecast(mod$ARIMA,k)[["mean"]])
  df$arima.fit <- as.numeric(arima.forecast)
  preds <- predict(mod$GAM,df)
  if (se) {
    df$arima.fit <- as.numeric(forecast(mod$ARIMA, k)[["lower"]][,"95%"])
    lower <- predict(mod$GAM,df,se.fit=T)
    lower <- lower$fit - 1.96*lower$se.fit
    
    df$arima.fit <- as.numeric(forecast(mod$ARIMA, k)[["upper"]][,"95%"])
    upper <- predict(mod$GAM,df,se.fit=T)
    upper <- upper$fit+1.96*upper$se.fit
    
    preds <- list(preds,lower,upper)
    names(preds) <- c('fit','lower','upper')
    if (sd.cap) {
      preds$fit <- pmax(preds$fit,lb)
      preds$fit <- pmin(preds$fit,ub)
    }
  }
  if (sd.cap&!se) {
    preds <- pmax(preds,lb)
    preds <- pmin(preds,ub)
  }
  
  return(preds)
}

Pres.Approval.BT <- function(df,test,p=12,Constrain.Prediction,m=1.5) {
  Mod <- fit.gam.arima(df)
  dates <- index(test)
  wf.predictions <- NULL
  wf.df <- df[,c("Approval","Years.In.Office","rgdp.pc","cpi","unemp","Elected","PopularVote")]
  
  for (i in 2:length(dates)) {
    if (length(dates)-i>=(p-1)) {
      n <- p-1
      pred.dat <- as.data.frame(test[i:c(i+n),c('PopularVote','Years.In.Office','rgdp.pc','cpi','unemp','Elected')])
    } else {
      n <- length(dates) - i
      pred.dat <- as.data.frame(test[i:c(i+n),c('PopularVote','Years.In.Office','rgdp.pc','cpi','unemp','Elected')])
    }
    pred.dat$rgdp.pc <- pred.dat$rgdp.pc[1]# on date i will not know what gdp was for date i + 12, only gdp for days 1:i
    pred.dat$cpi <- pred.dat$cpi[1]
    pred.dat$unemp <- pred.dat$unemp[1]
    pred.dat$PopularVote <- pred.dat$PopularVote[1] # will not know popular vote results before it happend 
    if (Constrain.Prediction) {
      ub <- as.numeric(wf.df$Approval[nrow(wf.df)])+m*sd(wf.df$Approval)
      lb <- as.numeric(wf.df$Approval[nrow(wf.df)])-m*sd(wf.df$Approval)
      preds <- predict.gam.arima(pred.dat,Mod,k=c(1+n),sd.cap = T,
                                 lb=max(lb,0),ub=max(ub,100))
    } else {
      preds <- predict.gam.arima(pred.dat,Mod,k=c(1+n),sd.cap = F)
    }
    if (length(preds)<p) {
      preds <- c(preds,rep(NA,p-c(n+1)))
    }
    names(preds) <- paste0('fc.',seq(1:p))
    preds <- as.xts(t(preds),dates[i])
    preds$wf.mean <- mean(wf.df$Approval)
    wf.predictions <- rbind(wf.predictions,preds)
    wf.df <- rbind(wf.df,test[i,c("Approval","Years.In.Office","rgdp.pc","cpi","unemp","Elected","PopularVote")])
    Mod <- fit.gam.arima(wf.df)
    if (i%%25==0|i==length(dates)) {
      cat(paste0(round(100*i/length(dates)),'% Complete. \n'))
    }
  }
  
  return(wf.predictions)
}

Pres.Aproval.BT.Res <- function(BT.Res,test) {
  OOSR2 <- NULL
  RMSE <- NULL
  for (i in 1:c(ncol(BT.Res)-1)) {
    pred <- na.omit(as.numeric(BT.Res[,paste0('fc.',i)]))
    pred <- as.xts(pred,index(test)[c(i+1):nrow(test)])
    tmp <- test$Approval
    tmp <- na.omit(cbind(tmp,pred,BT.Res$wf.mean))
    plot(tmp,main=paste(i,'Month Forecast'),col=c("black","red","green"))
    print(addLegend(legend.loc = "bottomleft", legend.names = c("Approval","Forecast","WF.Mean"), 
                    pch=15,col = c("black","red","green")))
    tmp2 <- sqrt(mean((tmp$Approval-tmp$pred)^2))
    tmp <- 1 - (sum((tmp$Approval-tmp$pred)^2)/sum((tmp$Approval-tmp$wf.mean)^2))
    OOSR2 <- c(OOSR2,tmp)
    RMSE <- c(RMSE,tmp2)
  }
  
  print(plot(1:c(ncol(BT.Res)-1),OOSR2))
  print(abline(0,0))
  print(plot(1:c(ncol(BT.Res)-1),RMSE))
  
  out <- data.frame(1:c(ncol(BT.Res)-1),OOSR2,RMSE)
  names(out) <- c('FC.Period','OOSR2','RMSE')
  return(out)
}


df <- readRDS('PresidentApprovalData.rds')
head(df)

df <- as.xts(df[,c("Approval","PopularVote",'Years.In.Office','rgdp.pc','cpi','unemp','Elected')],as.Date(df$Date))
plot(df$Approval)

df <- readRDS('PresidentApprovalData.rds')
df$Years.In.Office <- round(df$Years.In.Office,4)
Biden <- df[df$President=="Joseph Biden",c('Years.In.Office','Approval')]
names(Biden) <- c('Years.In.Office','Biden')
Trump <- df[df$President=="Donald Trump",c('Years.In.Office','Approval')]
names(Trump) <- c('Years.In.Office','Trump')
Obama <- df[df$President=="Barack Obama",c('Years.In.Office','Approval')]
names(Obama) <- c('Years.In.Office','Obama')
names <- c('Years.In.Office','Obama')
Average.PA <- df %>%
  dplyr::select(Approval,Years.In.Office) %>%
  dplyr::filter(Years.In.Office<=8) %>%
  group_by(Years.In.Office) %>%
  dplyr::summarise(Approval.mean=mean(Approval),
                   Approval.sd=sd(Approval)) %>%
  ungroup()
Average.PA <- left_join(Average.PA,Biden,by='Years.In.Office')
Average.PA <- left_join(Average.PA,Trump,by='Years.In.Office')
Average.PA <- left_join(Average.PA,Obama,by='Years.In.Office')

head(Average.PA)
Average.PA %>%
  dplyr::select(Years.In.Office,'Mean'=Approval.mean,Biden,Trump,Obama) %>%
  pivot_longer(-Years.In.Office) %>%
  ggplot(aes(Years.In.Office,value,color=name)) + geom_point() + geom_line() + theme_classic() + labs(title = "Average Presidential Approval and Time in Office")

df <- as.xts(df[,c("Approval","PopularVote",'Years.In.Office','rgdp.pc','cpi','unemp','Elected')],as.Date(df$Date))
head(df)

test <- df[round(0.75*nrow(df)):nrow(df),]
df <- df[1:round(0.75*nrow(df)),]

Mod <- fit.gam.arima(df)
summary(Mod$ARIMA)
summary(Mod$GAM)

par(mfrow=c(2,2))
plot.gam(Mod$GAM)

par(mfrow=c(2,2))
gam.check(Mod$GAM)
par(mfrow=c(1,1))
Plot.Trend(df,Mod$GAM)

wf.predictions <- Pres.Approval.BT(df,test,p=12,Constrain.Prediction = F)

Pres.Aproval.BT.Res(wf.predictions,test)

wf.predictions <- Pres.Approval.BT(df,test,p=12,Constrain.Prediction=T,m=1.5)
Pres.Aproval.BT.Res(wf.predictions,test)

df <- readRDS('PresidentApprovalData.rds')
Biden <- df[df$President=="Joseph Biden",]

df <- as.xts(df[,c("Approval","PopularVote",'Years.In.Office','rgdp.pc','cpi','unemp','Elected')],as.Date(df$Date))

Years.In.Office <- seq(as.numeric(Biden$Years.In.Office[nrow(Biden)])+0.0834,by=0.0834,length.out = 8)
rgdp.pc <- rep(as.numeric(Biden$rgdp.pc[nrow(Biden)]),8)
cpi <- rep(as.numeric(Biden$cpi[nrow(Biden)]),8)
unemp <- rep(as.numeric(Biden$unemp[nrow(Biden)]),8)
Elected <- rep(as.numeric(Biden$Elected[nrow(Biden)]),8)
PopularVote <- rep(as.numeric(Biden$PopularVote[nrow(Biden)]),8)
pred <- data.frame(Years.In.Office,rgdp.pc,unemp,cpi,PopularVote,Elected)
Mod <- fit.gam.arima(df)
ub <- as.numeric(df$Approval[nrow(df)])+1.5*sd(df$Approval)
lb <- as.numeric(df$Approval[nrow(df)])-1.5*sd(df$Approval)

tmp <- predict.gam.arima(pred,Mod,k=8,sd.cap=T,lb=lb,ub=ub,se=T)
Biden <- Biden[,c('Date','Approval')]
Biden$lower <- Biden$Approval
Biden$upper <- Biden$Approval
Biden <- as.xts(Biden[,-1],Biden$Date)

pred$Approval <- tmp$fit
pred$lower <- tmp$lower
pred$upper <- tmp$upper
pred <- pred[,c('Approval','lower','upper')]
pred <- as.xts(pred,seq.Date(as.Date('2021-06-01'),length.out = nrow(pred)+1,by='month')[2:9])

Biden <- rbind(Biden,pred)
plot(Biden,col = c("black","red","red"),lty=c(1,2,2))
Biden
