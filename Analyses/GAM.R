# this file provides the codes for GAM main analyses(O3-adjusted model)
# For single-factor model, just delete the O3 covariant (delete "+cb.coplt3" from line 22, 76, 128) 
# 99.9%: 3.29; 99%: 2.58
library(tidyverse); library(dlnm);library(scales)
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
  
  iqr <- data %>%
    pull(!!sym(plt)) %>%
    IQR(na.rm = TRUE)
  
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
  
  iqr <- data %>%
    pull(!!sym(coplt1)) %>%
    IQR(na.rm = TRUE)
  eval(parse(text=paste("pred <- crosspred(cb.plt, model,
                        at=iqr+17.5, cumul=TRUE)",sep='')))
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
  
  iqr <- data %>%
    pull(!!sym(coplt2)) %>%
    IQR(na.rm = TRUE)
  
  eval(parse(text=paste("pred <- crosspred(cb.plt, model,
                        at=60+iqr, cumul=TRUE)",sep='')))
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

output <- output %>%
  mutate(
    beta = if_else(plt == "temp", beta / 5, beta),
    h = if_else(plt == "temp", h / 5, h),
    l = if_else(plt == "temp", l / 5, l)
  )

output$Coplt <- output$coplt
output$plt <- factor(output$plt, levels=c("o3h8max","co","fsp","no2","temp","rh"),
                     labels = c(expression(paste("O"[3],sep = " ")),
                                expression(paste("CO",sep = " ")),
                                expression(paste("PM"[2.5],sep = " ")),
                                expression(paste("NO"[2],sep = " ")),
                                expression(paste("Temp",sep = " ")),
                                expression(paste("Humid",sep = " "))))


color_values <- c('gray20', "#d62728") 
output$lag <- as.numeric(output$lag)

P2_2 <- ggplot(output, aes(lag, beta, ymin = l, ymax = h)) +  
  geom_hline(yintercept = 0, linetype = 'dashed', color = "black", size = 0.8) +
  geom_ribbon( fill = "gray70", color = NA, alpha = 0.6) +  
  geom_point(aes(shape = factor(Significance),col = Significance, shape = Significance), size = 2, stroke = 1.5, alpha=0.7) +  
  geom_line(col="darkgray") +
  geom_errorbar(data = subset(output, Significance == "Sig"), 
                aes(ymin = l, ymax = h), 
                width = 0.2, 
                color = "#d62728", 
                size = 0.8) +  
  scale_color_manual(values = color_values) +
  scale_shape_manual(values = c(16, 8)) +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  facet_wrap(~plt, scales = 'fixed', labeller = label_parsed, ncol = 6) +
  ylab(expression(paste("Risk estimates β (env → COPD)", sep = ""))) +
  scale_x_continuous(breaks = unique(output$lag), labels = c(0:6)) +
  xlab("Lag") +
  theme_minimal(base_size = 14) + 
  theme(
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    axis.ticks.length = unit(0.2, 'cm'),
    axis.line = element_line(colour = "black"),
    legend.position = "None",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.box = "vertical",
    legend.margin = margin(-15, unit = "pt"),
    legend.box.spacing = unit(10, "pt"),
    legend.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent"),
    strip.background = element_rect(color = "white", fill = "white"),
    strip.text = element_text(size = 14),
    panel.grid.major = element_line(color = "gray90", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold")
  )

P2_2
    