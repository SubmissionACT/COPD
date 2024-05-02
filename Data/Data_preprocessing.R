library(tidyverse); library(dlnm);
library(lubridate); library(tsModel);library(ggsci)
library(mgcv); library(scales); library(splines)

################################## Load data ###################################
df <- read.csv("COPD0516.csv")
df1 <- df %>% select(date, copd, temp, rh, co, fsp, no2, o3h8max)

### data standardization --- dtrend, dseason, and normalizaiton
nomz <- function(x, normalization=T, dseason=T, season_sd=F, sea=365, dtrend=T, dTtype="linear"){
  x <- as.numeric(x)
  xt <- x
  # Detrend
  if(dtrend==T & dTtype=="first"){xt <- diff(xt)} else if (dtrend==T & dTtype=="linear"){
    lm.t <- lm(xt~c(1:length(xt)))
    xt <- xt-(lm.t$coefficients[1]+lm.t$coefficients[2]*c(1:length(xt)))}
  # Deseason
  if(dseason==T){
    xs <- as.numeric(apply(matrix(xt[1:(sea*length(xt)%/%sea)],ncol=sea,byrow=T),2,mean,na.rm=T))
    xsd <- as.numeric(apply(matrix(xt[1:(sea*length(xt)%/%sea)],ncol=sea,byrow=T),2,sd,na.rm=T))
    xt <- xt-c(rep(xs,1+length(xt)%/%sea))[1:length(xt)]
    if(season_sd==T){xt <- xt/(c(rep(xsd,1+length(xt)%/%sea))[1:length(xt)])}}
  # Normalization (zero mean & unity variance)
  if(normalization==T){xt <- (xt-mean(xt,na.rm=T))/sd(xt,na.rm=T)}
  return(xt)
}

df1 <- df1 %>% as.data.frame() %>% mutate_at(2:ncol(.),nomz) 

write.csv(df1, file = "COPD0516_preprocessed.csv")
