## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(CGNM)

## ----define the model function------------------------------------------------
library(RxODE)

model_text="
d/dt(X_1)=-ka*X_1
d/dt(C_2)=(ka*X_1-CL_2*C_2)/V1"

# here the model defined as above is compiled
model=RxODE(model_text)

# here we define the model function where takes in the parameter vector x and return the model simulation
model_function=function(x){

observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)

theta <- c(ka=x[1],V1=x[2],CL_2=x[3])
ev <- eventTable()
ev$add.dosing(dose = 1000, start.time =0)
ev$add.sampling(observation_time)
odeSol=model$solve(theta, ev)
log10(odeSol[,"C_2"])
}

## -----------------------------------------------------------------------------
observation=log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))

## ---- warning = FALSE---------------------------------------------------------
CGNM_result=Cluster_Gauss_Newton_method(nonlinearFunction=model_function,
targetVector = observation,
initial_lowerRange = c(0.1,0.1,0.1),initial_upperRange =  c(10,10,10))

## -----------------------------------------------------------------------------
head(acceptedApproximateMinimizers(CGNM_result))

## -----------------------------------------------------------------------------
CGNM_bootstrap=Cluster_Gauss_Newton_Bootstrap_method(CGNM_result, nonlinearFunction=model_function)

## -----------------------------------------------------------------------------
library(ggplot2)

## ---- fig.width=6, fig.height=3.5---------------------------------------------
plot_Rank_SSR(CGNM_result)

## ---- fig.width=6, fig.height=3.5---------------------------------------------
plot_paraDistribution_byHistogram(CGNM_bootstrap,  ParameterNames=c("Ka","V1","CL_2"), ReparameterizationDef=c("x1","x2","x3"), bins = 50)

## -----------------------------------------------------------------------------
plot_goodnessOfFit(CGNM_result, plotType = 1, independentVariableVector = c(0.1,0.2,0.4,0.6,1,2,3,6,12), plotRank = seq(1,50))

## -----------------------------------------------------------------------------
plot_goodnessOfFit(CGNM_bootstrap, plotType = 1, independentVariableVector = c(0.1,0.2,0.4,0.6,1,2,3,6,12))

