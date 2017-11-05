
########################################################################
## Functions to perform martingale residual-based imputation          ##
## "Martingale residual-based method to control for confounders       ##
##  measured only in a validation sample in time-to-event analysis",  ##
##  Burne R.M., Abrahamowicz M., Stat-in-Med 2016                     ##
##                                                                    ##
## Code: R.M. Burne                                                   ##
## Date edited: 22/11/2016  BY  R.M. Burne                            ##
########################################################################

# This code contains the functions MRBasedImputation() and ImputeFunction()
# MRBasedImputation() calls ImputeFUnction(), so it's necessary to load both functions before calling

###############################
## MARTINGALE MODEL FUNCTION ##
###############################
## Description:
# Uses MR-based imputation to impute unmeasured confounders, then fits a final CoxPH model dependent 
# on exposure, measured and unmeasured confounders. Returns model output as well as data with imputed U values
# Works for both internal and external VS - but input is still a single data frame (just specify external.vs term)
# Works for time-fixed OR time-varying data

## Usage:
# MRBasedImputation(dat=, var.obs=, var.unmeas=, exposure=, time=, event=, external.vs=TRUE)

## Arguments:
# dat: data frame with combined main data and VS (must be in data.frame form) 
#       NB - there can be no missingness other than unmeasured confounders for the main data
# var.obs: character vector with names of measured confounders eg c("C1","C2")
# var.unmeas: character vector with names of unmeasured confounders eg c("U1","U2")
# exposure: character vector length 1 with name of exposure eg "X"
# time: character vector length 1 with name of time variable eg "y". 
#         For right-censored (non-time-varying) data this is the follow-up time
#         For interval data, this is the starting time for the interval
# time2: (optional) character vector length 1 with name of STOP variable eg "stop"
#         Specify if using interval notation (i.e. Surv(start,stop,event)) for time-varying data
# event: character vector length 1 with name of censoring indicator variable eg "e"
# external.vs: TRUE/FALSE - If validation sample is external to the main data (default=TRUE)
# MI: TRUE/FALSE - If false, perform single imputation. If true, multiple imputation (default=FALSE)
# m: Specify if MI=TRUE: number of imputations
# return.data: TRUE/FALSE - whether or not to return the data with imputed values (default=TRUE)


MRBasedImputation <- function(dat, var.obs, var.unmeas, exposure, time, time2, event, external.vs=TRUE, MI=FALSE, m, return.data=TRUE){
  require(survival)
  
  # Specify "Surv" object dependent on whether time-varying data - (start,stop) notation - or time-fixed data is present
  if(missing(time2)){
    surv.obj <- paste("Surv(",time,",",event,")")
  } else {
    surv.obj <- paste("Surv(",time,",",time2,",",event,")")
  }
  
  # Get martingale residual based on whether validation sample is external or not
  if(external.vs==TRUE){
    # Split data into main and VS, since models are run separately for External VS
    datamain <- dat[!(complete.cases(dat)),]
    dataVS <- dat[complete.cases(dat),]
    rm(dat);gc()
    
    # Get martingale residuals in data
    formula <- paste0(surv.obj,"~",exposure)
    for (i in 1:length(var.obs)){
      formula<-paste0(formula,"+",var.obs[i])
    }
    
    naivemodelVS <- coxph(eval(parse(text=formula)),data=dataVS)
    naivemodelmain <- coxph(eval(parse(text=formula)),data=datamain)
    
    # Obtain martingale residuals #
    dataVS$MR <- residuals(naivemodelVS,type="martingale")
    datamain$MR <- residuals(naivemodelmain,type="martingale")
    
  } else { # (If external.vs=FALSE)
    # Get martingale residuals in (full) data
    formula <- paste0(surv.obj,"~",exposure)
    for (i in 1:length(var.obs)){
      formula<-paste0(formula,"+",var.obs[i])
    }
    
    naivemodel <- coxph(eval(parse(text=formula)),data=dat)
    
    # Obtain martingale residuals #
    dat$MR <- residuals(naivemodel,type="martingale")
    
    # Get main and VS
    datamain <- dat[!(complete.cases(dat)),]
    dataVS <- dat[complete.cases(dat),]
    rm(dat);gc()
  }
  
  # Imputation steps
  temp <- datamain[,!names(datamain) %in% var.unmeas]
  if(MI==FALSE){ # if not multiple imputation, data.MR will be a single data frame with imputed values for var.unmeas
    U.MR <- sapply(var.unmeas,ImputeFunction,var.obs=var.obs,exposure=exposure,VS=dataVS,main=datamain,event=event)
    data.MR <- rbind(cbind(temp,U.MR),dataVS)
  } else { # if multiple imputation, data.MR will be a list of data.frames with imputed values for var.unmeas
    data.MR <- list()
    for(i in 1:m){
      U.MR <- sapply(var.unmeas,ImputeFunction,var.obs=var.obs,exposure=exposure,VS=dataVS,main=datamain,event=event)
      data.MR[[i]] <- rbind(cbind(temp,U.MR),dataVS)
      rm(U.MR)
    }
  }
  
  
  ## FINAL MODELS AND OUTPUT
  to.return <- list()
  to.return$models <- list()
  if(return.data==TRUE){
    to.return$data <- data.MR
    }
  
  
  formula.full<-paste(surv.obj,"~",exposure,sep="")
  for (i in 1:length(var.obs)){
    formula.full<-paste(formula.full,"+",var.obs[i],sep="")
  }
  for (i in 1:length(var.unmeas)){
    formula.full<-paste(formula.full,"+",var.unmeas[i],sep="")
  }
  
  # Imputation with MR:
  # Depending on whether MI=FALSE or TRUE, perform imputation once or m (specified) times
  if(MI==FALSE){
    # single imputation:
    model.MRImp <- coxph(eval(parse(text=formula.full)),data=data.MR)
    to.return$models$MRImp <- summary(model.MRImp)
  } else { 
    # multiple imputation:
    # Get model results for each of the m imputed datasets
    templistofmodels <- lapply(1:m,function(i){ summary(coxph(eval(parse(text=formula.full)), data=data.MR[[i]])) })
    # Use Rubin's rules to combine these results
    MI.est <- colMeans(do.call("rbind",lapply(lapply(templistofmodels,"coef"),"[",,"coef")))  
    w.sd <- do.call("rbind",lapply(lapply(templistofmodels,"coef"),"[",,"se(coef)"))
    w.var <- colMeans(w.sd^2)
    b.var <- unlist(lapply(1:(sum(length(var.obs)+length(var.unmeas))+1),
                           function(j){(1 + 1/m)*sum((do.call("rbind",lapply(lapply(templistofmodels,"coef"),"[",,"coef"))[,j] - MI.est[j])^2)}))
    MI.sd <- sqrt(w.var + b.var)
    
    MI.results <- cbind(MI.est,MI.sd); colnames(MI.results) <- c("coef","se(coef)")
    model.MRimp <- list()
    to.return$models$MRimp <- MI.results
  }
  
  return(to.return)
}


#######################
##IMPUTATION FUNCTION##
#######################
## Description:
# Performs the MR-based imputation step. Called within the MRBasedImputation function above.
#
## Arguments:
# U.name: character vector of length 1 with name of the unmeasured confounder to impute
# var.obs: character vector with observed variables (to be used in imputation model)
# exposure: character vector length 1 with name of exposure eg "X"
# VS: validation sample (data.frame)
# main: main data (data.frame)
# timevar: If time variable is used in imputation model (used for LogT imputation)
# event: character vector length 1 with name of censoring indicator variable eg "e"
# imputation: options "martingale" or "timevar" - for either MR-based imputation or LogT imputation


ImputeFunction <- function(U.name,var.obs,exposure,VS,main,event){
  
  ## Create variable "type" - binomial / gaussian (continuous) ##
  # Note that currently this does NOT consider any sort of adjustment for normality, it just considers 
  # continuous confounders as normal
  type <- if(all(VS[,names(VS)==U.name] %in% c(0,1))){"binomial"} else {"gaussian"}
  
  ## Fit model dependent on type ## 
  # Define the formula: #
  formula <- paste0(exposure, "+MR")
  
  for (i in 1:length(var.obs)){
    formula<-paste(formula,"+",var.obs[i],sep="")
  }
  
  imp.model <- switch(type,
                      gaussian = lm(eval(parse(text=paste(U.name,"~",formula))),data=VS),
                      binomial = glm(eval(parse(text=paste(U.name,"~",formula))),data=VS,family=binomial))
  
  coef <- imp.model$coef
  
  temp <- names(coef)[2:length(coef)]
  XmatVS <- as.matrix(cbind(1,subset(VS,select=temp)))
  Xmatmain <- as.matrix(cbind(1,subset(main,select=temp)))
  
  ## Impute step. ##
  u.imp <- switch(type,
                  binomial = apply(Xmatmain, 1,
                                   function(x){
                                     mu <- sum(coef*x)
                                     p <- exp(mu)/(1+exp(mu))
                                     new.u <- rbinom(1,1,p)
                                   }),
                  gaussian = apply(Xmatmain, 1,
                                   function(x){
                                     mu <- sum(coef*x)
                                     var <- summary(imp.model)$sigma^2*
                                       (1+t(as.matrix(x))%*%solve(t(XmatVS)%*%XmatVS)%*%as.matrix(x))    
                                     # (This is the prediction variance)
                                     new.u <- rnorm(1,mean=mu,sd=sqrt(var))}))
  
  rm(Xmatmain);rm(XmatVS) 
  gc()
  return(u.imp)  
}
