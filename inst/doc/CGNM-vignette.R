## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(CGNM)

## ----define the model function------------------------------------------------
model_function=function(x){

  observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
  Dose=1000
  F=1

  ka=10^x[1]
  V1=10^x[2]
  CL_2=10^x[3]
  t=observation_time

  Cp=ka*F*Dose/(V1*(ka-CL_2/V1))*(exp(-CL_2/V1*t)-exp(-ka*t))

  log10(Cp)
}


## -----------------------------------------------------------------------------
observation=log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))

## ---- warning = FALSE---------------------------------------------------------
CGNM_result=Cluster_Gauss_Newton_method(nonlinearFunction=model_function,
targetVector = observation,
initial_lowerRange = c(-2,-2,-2),initial_upperRange =  c(2,2,2), saveLog=TRUE, num_minimizersToFind = 500)

## -----------------------------------------------------------------------------
head(acceptedApproximateMinimizers(CGNM_result))

## -----------------------------------------------------------------------------
CGNM_bootstrap=Cluster_Gauss_Newton_Bootstrap_method(CGNM_result, nonlinearFunction=model_function)

## -----------------------------------------------------------------------------
library(ggplot2)

## ---- fig.width=6, fig.height=3.5---------------------------------------------
plot_Rank_SSR(CGNM_result)

## ---- fig.width=6, fig.height=3.5---------------------------------------------
plot_paraDistribution_byHistogram(CGNM_bootstrap,  ParameterNames=c("Ka","V1","CL_2"), ReparameterizationDef=c("10^x1","10^x2","10^x3"), bins = 50)

## ---- fig.width = 7-----------------------------------------------------------
plot_goodnessOfFit(CGNM_result, plotType = 1, independentVariableVector = c(0.1,0.2,0.4,0.6,1,2,3,6,12), plotRank = seq(1,50))

## ---- fig.width = 7-----------------------------------------------------------
plot_goodnessOfFit(CGNM_bootstrap, plotType = 1, independentVariableVector = c(0.1,0.2,0.4,0.6,1,2,3,6,12))

## ---- fig.width = 7-----------------------------------------------------------
plot_profileLikelihood(c("CGNM_log","CGNM_log_bootstrap"))

## ---- fig.width = 7, fig.height = 6-------------------------------------------
plot_2DprofileLikelihood("CGNM_log", numBins = 50, ParameterNames=c("log10(Ka)","log10(V1)","log10(CL_2)"), ReparameterizationDef = c("x1","x2","x3"))

