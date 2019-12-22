##############################################################
#Code for estimating one step ahead errors

#mydir <- "C:/Users/Felix/Google Drive/Sync drive/Research/NationalHomicide/Data and Analysis" #set to your personal machine
mydir <- "C:\\Users\\apwhe\\Dropbox\\School_Projects\\National_Homicide_ErrorBars\\Analysis"
setwd(mydir)

library(forecast)
library(ggplot2)
#load the data
safety <- read.csv(file="CrimeRate_to2018.csv", header=TRUE)

#The rates per 100,000 are truncated, so replacing with floating point values
safety$HomRate <- (safety$Murder/safety$Population)*100000
safety$AggRate <- (safety$Aggravated_assault/safety$Population)*100000
safety$RobRate <- (safety$Robbery/safety$Population)*100000
safety$RapRate <- (safety$Legacy_rape/safety$Population)*100000

#################################################
#FUNCTIONS

#Returns low,cover,or high
lch <- function(val,low,high){
  if(val < low){res <- "low"}
  else if (val > high){res <- "high"}
  else {res <- "cover"}
  return(res)
}

#Function to do the fitting and step ahead intervals
MyForecast <- function(series, years, train, intervals){
  #Creating original time series
  my_ts <- ts(series, start=years[1], end = years[2])
  #Data frame to stuff the results in
  end_df <- data.frame(years=years[1]:years[2],series=series) 
  #maybe append one more row
  #to do one step ahead 
  #And add Pred/SE?s
  for (i in intervals){ 
    nm_low <- paste0("lower_",i)
    nm_high <- paste0("higher_",i)
    nm_cov <- paste0("cov_",i)
    end_df[,c(nm_low)] <- NA
    end_df[,c(nm_high)] <- NA	
    end_df[,c(nm_cov)] <- ""
  }
  #Splitting up training and testing data
  mark_train <- train - years[1] + 1
  after_train <- mark_train + 1
  end <- years[2] - years[1] + 1
  trainingdata <- subset(my_ts, end=mark_train)
  #Fitting the model on the test data to select ARIMA
  mod_fit <- auto.arima(trainingdata, lambda=0, stepwise=FALSE, approximation=FALSE)
  #Applying the predictions to the future data, one step ahead
  for (ii in mark_train:(end-1)){
    curr_dat <- subset(my_ts, end=ii)
    ar_model_update <- Arima(curr_dat, lambda=0, model=mod_fit)
    for (j in intervals){
      nm_low <- paste0("lower_",j)
      nm_high <- paste0("higher_",j)
      nm_cov <- paste0("cov_",j)
      fcst <- forecast(ar_model_update,h=1,level=j/100)
      end_df[ii+1,c(nm_low)] <- fcst$lower
      end_df[ii+1,c(nm_high)] <- fcst$upper
      end_df[ii+1,c(nm_cov)] <- lch(end_df[ii+1,2],fcst$lower,fcst$upper)
    }
  }
  #Creating a nice cover table
  cov_vars <- paste0("cov_",intervals)
  lt <- reshape(data=end_df[after_train:end,c("years",cov_vars)],times=cov_vars,
                   timevar="Int",idvar="years",direction="long",
                   varying=cov_vars,v.names="Cover")
  lt$Cover <- factor(lt$Cover,levels=c("low","cover","high"))
  cov_tab <- table(lt$Int,lt$Cover) 
  #Returning the data and the model fit
  return(list(mod_fit,end_df,cov_tab))
}

#Function to make a nice ggplot
my_bgg <- function(data,title,ylab){
  my_plot <- ggplot() +
    geom_ribbon(data=data, aes(x = years, ymin=lower_98, ymax=higher_98, colour = "98% PI"), fill = "#CCCCCC") + 
    geom_ribbon(data=data, aes(x = years, ymin=lower_90, ymax=higher_90, colour = "90% PI"), fill = "#B2B2B2") + 
    geom_ribbon(data=data, aes(x = years, ymin=lower_80, ymax=higher_80, colour = "80% PI"), fill = "#929292") + 
    geom_ribbon(data=data, aes(x = years, ymin=lower_50, ymax=higher_50, colour = "50% PI"), fill = "#666666") +
    geom_line(data=data, aes(x= years, y=series), color = "#000000") +
    geom_point(data=data, aes(x= years, y=series, colour = "Observed Crime"), fill = "#000000") + 
    scale_x_continuous(breaks=seq(1960,2018,by=5)) + 
    ggtitle(title,"adjusted for population") +
    labs(x = "Year", y= ylab, caption = "UCR data compiled by FBI")  + 
    scale_colour_manual(name = "PI", labels = c("50% PI", "80% PI","90% PI","98% PI", "Observed Crime"), 
                        values=c("Observed Crime" = "#000000", "98% PI" = "#CCCCCC", "90% PI" = "#B2B2B2", 
                                 "80% PI" = "#929292", "50% PI" = "#666666")) + 
    theme(legend.key=element_rect(fill="white")) + theme_bw() + theme_linedraw() + 
    theme(panel.grid = element_line(color="black", linetype = "dotted")) + 
    guides(colour = guide_legend(override.aes = list(size=3, stroke=3))) + theme(text = element_text(size=8)) +
    theme(legend.position=c(.5,.2))
  return(my_plot)
}
#################################################

#################################################
#Now doing the analysis for the individual series

my_years <- range(safety$Year)
train_end <- 1990
pints <- c(50,80,90,98)

##################
#Homicide Rate Analysis
hr <- MyForecast(series=safety$HomRate,years=my_years,train=train_end,intervals=pints)
summary(hr[[1]])

#######################
#Table 1, fit stats for different ARIMA structure
my_lag <- 6 #Taken from checkresiduals function for this length

orders <- list(c(1,1,0),c(0,1,0),c(0,1,1),c(1,1,1),c(2,1,0),c(2,1,1))
tot_mods <- length(orders)
AIC_Tab <- data.frame(rows=1:tot_mods)
AIC_Tab$Ord <- ""
AIC_Tab$Xsq <- NA
AIC_Tab$p <- NA
AIC_Tab$Aicc <- NA
for (i in 1:tot_mods){
  arima_mod <- Arima(hr[[1]]$x, lambda=0, order = orders[[i]])
  LB <- Box.test(zoo::na.approx(residuals(arima_mod)), fitdf = length(coef(a1)), 
                 lag = my_lag, type = "Ljung")
  AIC_Tab[i,c('Ord')] <- paste0("ARIMA (",paste0(orders[[i]],collapse=","),"):")
  AIC_Tab[i,c('Xsq')] <- LB$statistic
  AIC_Tab[i,c('p')] <- LB$p.value
  AIC_Tab[i,c('Aicc')] <- arima_mod$aicc
}

#Table 1
AIC_Tab
######################

#Figure 1
png('HomResids.png',height=5,width=7,units="in",res=1000,type="cairo")
checkresiduals(hr[[1]])
dev.off()

#Figure 2
png('Homi_forecast.png', height=5, width=7, units="in", res=1000, type="cairo") 
my_bgg(data=hr[[2]],title="Homicide Rates in the United States 1960-2018",
       ylab="Homicide Rate per 100,000")
dev.off()

#Cover Table
hr[[3]]
##################

##################
#Aggravated Assault Rate Analysis
ar <- MyForecast(series=safety$AggRate,years=my_years,train=train_end,intervals=pints)
summary(ar[[1]])

#Fig 3
png('AA_forecast.png', height=5, width=7, units="in", res=1000, type="cairo") 
my_bgg(data=ar[[2]],title="Aggravated Assault Rates in the United States 1960-2018",
       ylab="Aggravated Assault Rate per 100,000")
dev.off()

#Cover Table
ar[[3]]
##################

##################
#Rape
ra <- MyForecast(series=safety$RapRate,years=my_years,train=train_end,intervals=pints)
summary(ra[[1]])

#Fig 4
png('RP_forecast.png', height=5, width=7, units="in", res=1000, type="cairo") 
my_bgg(data=ra[[2]],title="Legacy Rape Rates in the United States 1960-2018",
       ylab="Rape Rate per 100,000")
dev.off()

#Cover Table
ra[[3]]
##################

##################
#Robbery
ro <- MyForecast(series=safety$RobRate,years=my_years,train=train_end,intervals=pints)
summary(ro[[1]])

#Fig 5
png('RB_forecast.png', height=5, width=7, units="in", res=1000, type="cairo") 
my_bgg(data=ro[[2]],title="Robbery Rates in the United States 1960-2018",
       ylab="Robber Rate per 100,000")
dev.off()
#Maybe combine 3,4,5 figures into one via patchwork?

#Cover Table
ro[[3]]
##################

#Can add cover tables together to get overall coverage
ar[[3]] + hr[[3]] + ra[[3]] + ro[[3]] #Not too shabby
#################################################


#################################################
#Supplemental analysis when varying the start year
#For the model fitting process, calculate coverage
#For each crime type and cumulative

test_years <- 1980:2000
fill_set <- expand.grid(paste0("cov_",pints),c("low","cover","high"),years=test_years)
fill_set$Hom <- NA
fill_set$Agg <- NA
fill_set$Rob <- NA
fill_set$Rap <- NA

#Range of years to check 
for (y in test_years){
  print(y) #To check if it messes up inside loop
  hr_y <- MyForecast(series=safety$HomRate,years=my_years,train=y,intervals=pints)
  ar_y <- MyForecast(series=safety$AggRate,years=my_years,train=y,intervals=pints)
  ro_y <- MyForecast(series=safety$RobRate,years=my_years,train=y,intervals=pints)
  ra_y <- MyForecast(series=safety$RapRate,years=my_years,train=y,intervals=pints)
  hr_cdf <- as.data.frame(hr_y[[3]])
  ar_cdf <- as.data.frame(ar_y[[3]])
  ro_cdf <- as.data.frame(ro_y[[3]])
  ra_cdf <- as.data.frame(ra_y[[3]])
  fill_set[fill_set$years == y,c('Hom')] <- hr_cdf$Freq
  fill_set[fill_set$years == y,c('Agg')] <- ar_cdf$Freq
  fill_set[fill_set$years == y,c('Rob')] <- ro_cdf$Freq
  fill_set[fill_set$years == y,c('Rap')] <- ra_cdf$Freq
}

fill_set$Tot <- fill_set$Hom + fill_set$Agg + fill_set$Rob + fill_set$Rap
fill_set$Tot_Cov <- fill_set$Tot * (fill_set$Var2 == "cover")

agg_cov <- aggregate(cbind(Tot_Cov,Tot)~Var1+years,data=fill_set,FUN=sum)
agg_cov$Perc <- (agg_cov$Tot_Cov / agg_cov$Tot)*100

cov_stats <- reshape(data=agg_cov,timevar='Var1',idvar=c('Tot','years'),direction="wide",drop=c("Tot_Cov"))
cov_stats


#Save results
save.image(file="HomErrorBar_Results.RData")
#################################################

