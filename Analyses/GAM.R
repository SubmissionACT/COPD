# this file provides the codes for GAM main analyses(O3-adjusted model)
# For single-factor model, just delete the O3 covariant (delete "+cb.coplt3" from line 22, 73, 123) 
# 99.9%: 3.29; 99%: 2.58
library(tidyverse); library(dlnm);
library(lubridate); library(tsModel);library(ggsci)
library(mgcv); library(scales); library(splines)

################################## Load data ###################################
df <- read.csv("COPD0516.csv")
df$t <- 1:length(df$date)
data=df

################################## Pollutant ###################################

FTjust <- function(dis, plt, coplt1, coplt2,df_t, coplt3){
  
  cb.plt = crossbasis(data[[plt]], lag=6, argvar=list(fun="lin"), arglag=list(fun="integer"))
  cb.coplt1 = crossbasis(data[[coplt1]], lag=3, argvar=list(fun="ns",df = 3), arglag=list(fun="integer"))
  cb.coplt2 = crossbasis(data[[coplt2]], lag=3, argvar=list(fun="ns",df = 3), arglag=list(fun="integer"))
  cb.coplt3 = crossbasis(data[[coplt3]], lag=6, argvar=list(fun="lin"), arglag=list(fun="integer"))
  
  plist <- ifelse(plt == coplt3, "cb.plt+cb.coplt1+cb.coplt2", "cb.plt+cb.coplt1+cb.coplt2+cb.coplt3")
  
  eval(parse(text=paste("model <- gam(", dis," ~ ", plist, "+
        s(t,k=",df_t,"*12,fx=T)+
        as.factor(dow)+as.factor(holiday),
        family=quasipoisson(link = 'log'), data=data)", sep='')))
  
  iqr<-10
  eval(parse(text=paste("pred <- crosspred(cb.plt, model,
                        at=iqr, cumul=TRUE)",sep='')))
  matRRfit <- data.frame(fit = "matfit", dis = dis, plt = plt, iqr=iqr,
                         coplt = coplt3,df_t=df_t,
                         pred$matfit, stringsAsFactors = FALSE)
  matRRhigh <- data.frame(fit = "mathigh", dis = dis, plt = plt, iqr=iqr,
                          coplt = coplt3,df_t=df_t,
                          pred$matfit+3.29*pred$matse, stringsAsFactors = FALSE)
  matRRlow <- data.frame(fit = "matlow", dis = dis, plt = plt,iqr=iqr,
                         coplt = coplt3,df_t=df_t,
                         pred$matfit-3.29*pred$matse, stringsAsFactors = FALSE)
  out <- bind_rows(matRRfit, matRRhigh, matRRlow)
  out
}

plist <- list(dis = c('copd'),
              plt = c("co","fsp","no2","o3h8max"),
              coplt1 = c('temp'),
              coplt2 = c("rh"),
              df_t=c(7),
              coplt3 = c("o3h8max"))

bb <- plist %>% cross_df()  %>%  pmap_df(FTjust)

bbA <- bb %>%
  gather(lagA, value, starts_with("lag")) %>%
  spread(fit, value) %>%
  mutate(lag = parse_number(lagA)) %>%
  mutate(lag=factor(lag))

output <- bbA %>% mutate(
  sig=as.factor(ifelse(mathigh<=0,1,ifelse(matlow>=0,1,0))))

output1 <- output

################################## temperature ###################################

FTjust <- function(dis, plt, coplt1, coplt2,df_t, coplt3){
  
  cb.plt  = crossbasis(data[[coplt1]], lag=6, argvar=list(fun="ns",df=3), arglag=list(fun="integer"))
  cb.coplt2 = crossbasis(data[[coplt2]], lag=3, argvar=list(fun="ns",df=3), arglag=list(fun="integer"))
  cb.coplt3 = crossbasis(data[[coplt3]], lag=6, argvar=list(fun="lin"), arglag=list(fun="integer"))
  
  plist <- ifelse(plt==coplt3,"cb.plt+cb.coplt2","cb.plt+cb.coplt2+cb.coplt3")
  
  eval(parse(text=paste("model <- gam(", dis," ~ ", plist, "+
        s(t,k=",df_t,"*12,fx=T)+
        as.factor(dow)+as.factor(holiday),
        family=quasipoisson(link = 'log'),  data=data)", sep='')))
  
  iqr<-mean(data[["temp"]],na.rm=TRUE)
  eval(parse(text=paste("pred <- crosspred(cb.plt, model,
                        at=iqr+1, cumul=TRUE)",sep='')))
  matRRfit <- data.frame(fit = "matfit", dis = dis, plt = "temp", iqr=iqr,
                         coplt = coplt3,df_t=df_t,
                         pred$matfit, stringsAsFactors = FALSE)
  matRRhigh <- data.frame(fit = "mathigh", dis = dis, plt = "temp", iqr=iqr,
                          coplt = coplt3,df_t=df_t,
                          pred$matfit+3.29*pred$matse, stringsAsFactors = FALSE)
  matRRlow <- data.frame(fit = "matlow", dis = dis, plt = "temp",iqr=iqr,
                         coplt = coplt3,df_t=df_t,
                         pred$matfit-3.29*pred$matse, stringsAsFactors = FALSE)
  out <- bind_rows(matRRfit, matRRhigh, matRRlow)
  out
}

plist <- list(dis = c('copd'),
              plt = c("temp"),
              coplt1 = c('temp'),
              coplt2 = c("rh"),
              df_t=c(7),
              coplt3 = c("o3h8max"))

bb <- plist %>% cross_df()  %>%  pmap_df(FTjust)

bbA <- bb %>%
  gather(lagA, value, starts_with("lag")) %>%
  spread(fit, value) %>%
  mutate(lag = parse_number(lagA)) %>%
  mutate(lag=factor(lag))

output <- bbA %>% mutate(
  sig=as.factor(ifelse(mathigh<=0,1,ifelse(matlow>=0,1,0))))
output2 <- output

################################## Relative humidity ###################################

FTjust <- function(dis, plt, coplt1, coplt2,df_t, coplt3){
  
  cb.coplt2 = crossbasis(data[[coplt1]], lag=3, argvar=list(fun="ns",df=3), arglag=list(fun="integer"))
  cb.plt = crossbasis(data[[coplt2]], lag=6, argvar=list(fun="ns",df=3), arglag=list(fun="integer"))
  cb.coplt3 = crossbasis(data[[coplt3]], lag=6, argvar=list(fun="lin"), arglag=list(fun="integer"))
  
  plist <- ifelse(plt==coplt3,"cb.plt+cb.coplt2","cb.plt+cb.coplt2+cb.coplt3")
  
  eval(parse(text=paste("model <- gam(", dis," ~ ", plist, "+
        s(t,k=",df_t,"*12,fx=T)+
        as.factor(dow)+as.factor(holiday),
        family=quasipoisson(link = 'log'), data=data)", sep='')))
  
  iqr<-mean(data[["rh"]],na.rm=TRUE)
  
  eval(parse(text=paste("pred <- crosspred(cb.plt, model,
                        at=10+iqr, cumul=TRUE)",sep='')))
  matRRfit <- data.frame(fit = "matfit", dis = dis, plt = "rh", iqr=iqr,
                         coplt = coplt3,df_t=df_t,
                         pred$matfit, stringsAsFactors = FALSE)
  matRRhigh <- data.frame(fit = "mathigh", dis = dis, plt = "rh", iqr=iqr,
                          coplt = coplt3,df_t=df_t,
                          pred$matfit+3.29*pred$matse, stringsAsFactors = FALSE)
  matRRlow <- data.frame(fit = "matlow", dis = dis, plt = "rh",iqr=iqr,
                         coplt = coplt3,df_t=df_t,
                         pred$matfit-3.29*pred$matse, stringsAsFactors = FALSE)
  out <- bind_rows(matRRfit, matRRhigh, matRRlow)
  out
}

plist <- list(dis = c('copd'),
              plt = c("rh"),
              coplt1 = c('temp'),
              coplt2 = c("rh"),
              df_t=c(7),
              coplt3 = c("o3h8max"))

bb <- plist %>% cross_df()  %>%  pmap_df(FTjust)

bbA <- bb %>%
  gather(lagA, value, starts_with("lag")) %>%
  spread(fit, value) %>%
  mutate(lag = parse_number(lagA)) %>%
  mutate(lag=factor(lag))

output <- bbA %>% mutate(
  sig=as.factor(ifelse(mathigh<=0,1,ifelse(matlow>=0,1,0))))
output3 <- output

################################## Plotting ###################################

output <- rbind(output1, output2, output3)
names(output) <- c("dis","plt","iqr","coplt","df_t","lagA","beta","h","l","lag","Significance")

output$lag <- factor(output$lag)

output$Significance <- factor(output$Significance,levels = c(0,1),labels = c("Non-sig","Sig"))

output$Coplt <- output$coplt
output$plt <- factor(output$plt,levels = c("o3h8max","co","fsp","no2","temp","rh"),labels = c("O3","CO","PM2.5","NO2","Temp","RH"))
output$df_t <- factor(output$df_t)

P2_2 <- ggplot(output, aes(lag,beta,ymin = l, ymax = h, group=control)) +
  geom_hline(yintercept = 0,linetype='dashed') +
  geom_errorbar(aes(group = plt, col = plt, width=0), 
                size=1, position = position_dodge((width=0.5))) +
  scale_color_jama()+
  geom_point(aes(group = plt, col = plt,shape=Significance),size = 2, 
             position = position_dodge(width=0.5),,fontface = "bold") +
  scale_shape_manual(values = c(16,8)) +
  facet_wrap(~plt,scales='free',labeller = "label_parsed",ncol = 3 ) +
  ylab(expression(paste(beta,"(Risk estimates air pollutants → COPD)", sep = "")))+
  xlab("Lag")+  
  theme_bw() +
  theme(axis.text.x=element_text(size = 12, color= "black"),
        axis.text.y=element_text(size = 12, color= "black"),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.ticks.length=unit(0.2,'cm'),
        axis.line = element_line(colour = "black"),
        legend.position=c("bottom"),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        legend.box="vertical",
        legend.margin=margin(-15,unit="pt"),
        legend.box.spacing = margin(15.5),
        legend.background = element_rect(fill="transparent"),
        strip.background = element_rect(
          color = "white", fill = "white"),
        panel.border = element_blank(),
        legend.key = element_rect(fill = "transparent"),
        strip.text = element_text(size = 14),
        panel.grid = element_blank(),
        plot.title = element_text(size = 18, hjust=0.5))
P2_2        

  
