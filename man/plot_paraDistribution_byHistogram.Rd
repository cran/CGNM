% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PostProcess.R
\name{plot_paraDistribution_byHistogram}
\alias{plot_paraDistribution_byHistogram}
\title{plot_paraDistribution_byHistogram}
\usage{
plot_paraDistribution_byHistogram(
  CGNM_result,
  indicesToInclude = NA,
  ParameterNames = NA,
  ReparameterizationDef = NA,
  bins = 30
)
}
\arguments{
\item{CGNM_result}{(required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.}

\item{indicesToInclude}{(default: NA) \emph{A vector of integers} indices to include in the plot (if NA, use indices chosen by the acceptedIndices() function with default setting).}

\item{ParameterNames}{(default: NA) \emph{A vector of strings} the user can supply so that these names are used when making the plot. (Note if it set as NA or vector of incorrect length then the parameters are named as theta1, theta2, ... or as in ReparameterizationDef)}

\item{ReparameterizationDef}{(default: NA) \emph{A vector of strings} the user can supply definition of reparameterization where each string follows R syntax}

\item{bins}{(default: 30) \emph{A natural number} Number of bins used for plotting histogram.}
}
\value{
\emph{A ggplot object} including the violin plot, interquartile range and median, minimum and maximum.
}
\description{
Make histograms to visualize the initial distribution and distribition of the accepted approximate minimizers found by the CGNM.
}
\examples{


model_analytic_function=function(x){

 observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
 Dose=1000
 F=1

 ka=x[1]
 V1=x[2]
 CL_2=x[3]
 t=observation_time

 Cp=ka*F*Dose/(V1*(ka-CL_2/V1))*(exp(-CL_2/V1*t)-exp(-ka*t))

 log10(Cp)
}

observation=log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))



CGNM_result=Cluster_Gauss_Newton_method(
nonlinearFunction=model_analytic_function,
targetVector = observation,
initial_lowerRange = rep(0.01,3), initial_upperRange =  rep(100,3),
lowerBound=rep(0,3), ParameterNames = c("Ka","V1","CL"),
num_iter = 10, num_minimizersToFind = 100, saveLog = FALSE)

plot_paraDistribution_byHistogram(CGNM_result)
plot_paraDistribution_byHistogram(CGNM_result,
     ReparameterizationDef=c("log10(Ka)","log10(V1)","log10(CL)"))

}
