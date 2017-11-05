
############################################################################
## How to use functions to perform martingale residual-based imputation   ##
##                                                                        ##
## For version of code: MR-based-imputation.R                             ##
##                                                                        ##
## "Martingale residual-based method to control for confounders           ##
##  measured only in a validation sample in time-to-event analysis",      ##
##  Burne R.M., Abrahamowicz M., Stat-in-Med 2016                         ##
##                                                                        ##
## Code: R.M. Burne                                                       ##
## Date edited: 29/11/2016  BY  R.M. Burne                                ##
############################################################################


## DESCRIPTION OF THE FUNCTION:
# Uses MR-based imputation to impute unmeasured confounders, then fits a final CoxPH model dependent 
# on exposure, measured and unmeasured confounders. Returns model output as well as data with imputed U values
# Works for both internal and external VS - but input is still a single data frame (just specify external.vs term)
# Works for time-fixed OR time-varying data (see examples below)

## Usage:
# MRBasedImputation(dat=, var.obs=, var.unmeas=, exposure=, time=, event=, external.vs=TRUE)

## Arguments:
# dat: data frame with *combined* main data and VS (must be in data.frame form) 
#       ***NB*** - there can be no missingness (other than missing unmeasured confounders in the main data)
# var.obs: character vector with names of measured confounders eg c("C1","C2")
# var.unmeas: character vector with names of unmeasured confounders eg c("U1","U2")
#       These may be continuous (normal) or binary (0/1) ONLY. No groups of dummy variables.
# exposure: character vector length 1 with name of exposure eg "X"
# time: character vector length 1 with name of time variable eg "y". 
#         For right-censored (non-time-varying) data this is the follow-up time
#         For interval data, this is the starting time for the interval
# time2: (optional) character vector length 1 with name of STOP variable eg "stop"
#         Specify only if using interval notation for time-varying data
# event: character vector length 1 with name of censoring indicator variable eg "e"
# external.vs: TRUE/FALSE - If validation sample is external to the main data (default=TRUE)
# MI: TRUE/FALSE - If false, perform single imputation. If true, multiple imputation (default=FALSE)
# m: Specify if MI=TRUE: number of imputations
# return.data: TRUE/FALSE - whether or not to return the data with imputed values (default=TRUE)


### EXAMPLES ###

# Source the functions:
setwd("C://path//to//your//directory")
source("MR-based-imputation.R")


### 1. Time-fixed data (single imputation) ###
# Load simulated data
load("simdata_TF.RData")

# To see the set-up of data:
str(data)
head(data)
tail(data)

nrow(data[!complete.cases(data),])    # 9000: the main database
nrow(data[complete.cases(data),])     # 1000: the validation sample

# NOTE: U1 & U2 are missing (NA) for ALL subjects in the main data, and available for ALL subjects
# in the VS. This is how the VS and the main data are distinguished in the function.
# There is NO other missingness (this is important).
# In this case (non-time-varying data) each line denotes a subject.
# U1 is continuous, and U2 is binary
# As long as all values of an unmeasured variable are 0/1 in the VS, the function will treat it as a binary variable,
# and if not, it will be treated as Normal.

# Use the function to obtain the output:
output <- MRBasedImputation(dat=data, var.obs=c("C1","C2"), var.unmeas=c("U1","U2"), exposure="X", time="time", event="event", external.vs=FALSE)
# (In this example, we treated the validation sample as an *internal* sample, 
# from the same source population. Could just as easily specify external.vs=TRUE if not.)

names(output)
# [1] "models" "data" 

# Get output from model with imputed U1 & U2 values.
output$models
# $MRImp
# Call:
#  coxph(formula = eval(parse(text = formula.full)), data = data.MR)
# ... (model output) ...

str(output$data)
# 'data.frame':	10000 obs. of  9 variables: ...
# Note that this will be the data with imputed values of the confounders U1 & U2

head(output$data)




### 2. Time-varying data (single imputation) ###
# Load simulated data:
load("simdata_TD.RData")

# To see the set-up of data:
str(data.td)
head(data.td)
tail(data.td)
length(unique(data.td$Id[complete.cases(data.td)]))   # 1000 subjects in the VS (complete cases)
length(unique(data.td$Id[!complete.cases(data.td)]))   # 10000 subjects in the main data (non-complete cases)

# NOTE: Again, U1 & U2 are missing (NA) for ALL subjects in the main data, and available for ALL subjects
# in the VS. This is how the VS and the main data are distinguished in the function.
# There is NO other missingness (this is important).
# In this case (time-varying data) each line denotes an interval for a subject.
# The interval beginning and end variables here are called "Start" and "Stop"

# Use the function to obtain the output:
output <- MRBasedImputation(dat=data.td, var.obs=c("C1","C2"), var.unmeas=c("U1","U2"), exposure="X", time="Start", time2="Stop", event="Event", external.vs=FALSE)

names(output)
# [1] "models" "data" 

# Get output from model with imputed U1 & U2 values.
output$models
# $MRImp
# Call:
#  coxph(formula = eval(parse(text = formula.full)), data = data.MR)
# ... (model output) ...

str(output$data)
# 'data.frame':	104508 obs. of  10 variables: ...
# this is the data with imputed values of the confounders U1 & U2

head(output$data)



### Multiple imputation ###
output <- MRBasedImputation(dat=data, var.obs=c("C1","C2"), var.unmeas=c("U1","U2"), exposure="X", time="y", event="e", external.vs=FALSE, MI=TRUE, m=5)

names(output)
# [1] "models" "data" 

output$models
# $MRImp
#coef   se(coef)
# ...

# Note that this output returns only coefficient and se(coef) at present - not all of the output returned from e.g. a summary() of a model

length(output$data)
# 5
# (data is returned in a list, where each element is a data.frame corresponding to a single imputation)

str(output$data[[1]])
# 'data.frame':	10000 obs. of  9 variables: ...

head(output$data[[1]])
