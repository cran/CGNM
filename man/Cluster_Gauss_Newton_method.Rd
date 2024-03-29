% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Cluster_Gauss_Newton_method.R
\name{Cluster_Gauss_Newton_method}
\alias{Cluster_Gauss_Newton_method}
\title{Cluster_Gauss_Newton_method}
\usage{
Cluster_Gauss_Newton_method(
  nonlinearFunction,
  ...,
  targetVector,
  initial_lowerRange,
  initial_upperRange,
  lowerBound = NA,
  upperBound = NA,
  ParameterNames = NA,
  stayIn_initialRange = FALSE,
  num_minimizersToFind = 250,
  num_iteration = 25,
  saveLog = TRUE,
  runName = "",
  textMemo = "",
  algorithmParameter_initialLambda = 1,
  algorithmParameter_gamma = 2,
  algorithmVersion = 3,
  initialIterateMatrix = NA,
  targetMatrix = NA,
  keepInitialDistribution = NA
)
}
\arguments{
\item{nonlinearFunction}{(required input) \emph{A function with input of a vector x of real number of length n and output a vector y of real number of length m.} In the context of model fitting the nonlinearFunction is \strong{the model}.  Given the CGNM does not assume the uniqueness of the minimizer, m can be less than n.  Also CGNM does not assume any particular form of the nonlinear function and also does not require the function to be continuously differentiable (see Appendix D of our publication for an example when this function is discontinuous). Also this function can be matrix to matrix equation.  This can be used for parallerization, see vignettes for examples.}

\item{...}{Further arguments to be supplied to nonlinearFunction}

\item{targetVector}{(required input) \emph{A vector of real number of length m} where we minimize the Euclidean distance between the nonlinearFuncition and targetVector.  In the context of curve fitting targetVector can be though as \strong{the observational data}.}

\item{initial_lowerRange}{(required input) \emph{A vector of real number of length n} where each element represents  \strong{the lower range of the initial iterate}. Similarly to regular Gauss-Newton method, CGNM iteratively reduce the residual to find minimizers.  Essential differences is that CGNM start from the initial RANGE and not an initial point.}

\item{initial_upperRange}{(required input) \emph{A vector of real number of length n} where each element represents  \strong{the upper range of the initial iterate}.}

\item{lowerBound}{(default: NA) \emph{A vector of real number or NA of length n} where each element represents  \strong{the lower bound of the parameter search}.  If no lower bound set that element NA. Note that CGNM is an unconstraint optimization method so the final minimizer can be anywhere.  In the parameter estimation problem, there often is a constraints to the parameters (e.g., parameters cannot be negative). So when the upper or lower bound is set using this option, parameter transformation is conducted internally (e.g., if either the upper or lower bound is given parameters are log transformed, if the upper and lower bounds are given logit transform is used.)}

\item{upperBound}{(default: NA) \emph{A vector of real number or NA of length n} where each element represents  \strong{the upper bound of the parameter search}.  If no upper bound set that element NA.}

\item{ParameterNames}{(default: NA) \emph{A vector of string} of length n User can specify names of the parameters that will be used for the plots.}

\item{stayIn_initialRange}{(default: FALSE) \emph{TRUE or FALSE} if set TRUE, the parameter search will conducted strictly within the range specified by initial_lowerRange and initial_upperRange.}

\item{num_minimizersToFind}{(default: 250) \emph{A positive integer} defining number of approximate minimizers CGNM will find. We usually \strong{use 250 when testing the model and 1000 for the final analysis}.  The computational cost increase proportionally to this number; however, larger number algorithm becomes more stable and increase the chance of finding more better minimizers. See Appendix C of our paper for detail.}

\item{num_iteration}{(default: 25)  \emph{A positive integer} defining maximum number of iterations. We usually \strong{set 25 while model building and 100 for final analysis}.  Given each point terminates the computation when the convergence criterion is met the computation cost does not grow proportionally to the number of iterations (hence safe to increase this without significant increase in the computational cost).}

\item{saveLog}{(default: TRUE) \emph{TRUE or FALSE} indicating either or not to save computation result from each iteration in CGNM_log folder. It requires disk write access right in the current working directory. \strong{Recommended to set TRUE if the computation is expected to take long time} as user can retrieve intrim computation result even if the computation is terminated prematurely (or even during the computation).}

\item{runName}{(default: "") \emph{string} that user can ue to identify the CGNM runs. The run history will be saved in the folder name CGNM_log_<runName>.  If this is set to "TIME" then runName is automatically set by the run start time.}

\item{textMemo}{(default: "") \emph{string} that user can write an arbitrary text (without influencing computation). This text is stored with the computation result so that can be used for example to describe model so that the user can recognize the computation result.}

\item{algorithmParameter_initialLambda}{(default: 1) \emph{A positive number} for initial value for the regularization coefficient lambda see Appendix B of of our paper for detail.}

\item{algorithmParameter_gamma}{(default: 2) \emph{A positive number} a positive scalar value for adjusting the strength of the weighting for the linear approximation see Appendix A of our paper for detail.}

\item{algorithmVersion}{(default: 3.0) \emph{A positive number} user can choose different version of CGNM algorithm currently 1.0 and 3.0 are available.  If number chosen other than 1.0 or 3.0 it will choose 1.0.}

\item{initialIterateMatrix}{(default: NA) \emph{A matrix} with dimension num_minimizersToFind x n.  User can provide initial iterate as a matrix  This input is used when the user wishes not to generate initial iterate randomly from the initial range.  The user is responsible for ensuring all function evaluation at each initial iterate does not produce NaN.}

\item{targetMatrix}{(default: NA) \emph{A matrix} with dimension num_minimizersToFind x m  User can define multiple target vectors in the matrix form.  This input is mainly used when running bootstrap method and not intended to be used for other purposes.}

\item{keepInitialDistribution}{(default: NA) \emph{A vector of TRUE or FALSE} of length n User can specify if the initial distribution of one of the input variable (e.g. parameter) to be kept as the initial iterate throughout CGNM iterations.}
}
\value{
list of a matrix X, Y,residual_history and initialX, as well as a list runSetting
\enumerate{\item X: \emph{a num_minimizersToFind by n matrix} which stores the approximate minimizers of the nonlinear least squares in each row. In the context of model fitting they are \strong{the estimated parameter sets}.
\item Y: \emph{a num_minimizersToFind by m matrix} which stores the nonlinearFunction evaluated at the corresponding approximate minimizers in matrix X above. In the context of model fitting each row corresponds to \strong{the model simulations}.
\item residual_history: \emph{a num_iteration by num_minimizersToFind matrix} storing sum of squares residual for all iterations.
\item initialX: \emph{a num_minimizersToFind by n matrix} which stores the set of initial iterates.
\item runSetting: a list containing all the input variables to Cluster_Gauss_Newton_method (i.e., nonlinearFunction, targetVector, initial_lowerRange, initial_upperRange ,algorithmParameter_initialLambda, algorithmParameter_gamma, num_minimizersToFind, num_iteration, saveLog, runName, textMemo).}
}
\description{
Find multiple minimisers of the nonlinear least squares problem.
\deqn{argmin_x ||f(x)-y*||}
where
\enumerate{\item f: nonlinear function (e.g., mathematical model)
\item y*: target vector (e.g., observed data to fit the mathematical model)
\item x: variable of the nonlinear function that we aim to find the values that minimize (minimizers) the differences between the nonlinear function and target vector (e.g., model parameter)
}
Parameter estimation problems of mathematical models can often be formulated as nonlinear least squares problems.  In this context f can be thought at a model, x is the parameter, and y* is the observation.
CGNM iteratively estimates the minimizer of the nonlinear least squares problem from various initial estimates hence finds multiple minimizers.
Full detail of the algorithm and comparison with conventional method is available in the following publication, also please cite this publication when this algorithm is used in your research: Aoki et al. (2020) <doi.org/10.1007/s11081-020-09571-2>. Cluster Gauss–Newton method. Optimization and Engineering, 1-31.  As illustrated in this paper, CGNM is faster and more robust compared to repeatedly applying the conventional optimization/nonlinear least squares algorithm from various initial estimates. In addition, CGNM can realize this speed assuming the nonlinear function to be a black-box function (e.g. does not use things like adjoint equation of a system of ODE as the function does not have to be based on a system of ODEs.).
}
\examples{
##lip-flop kinetics (an example known to have two distinct solutions)

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
targetVector = observation, num_iteration = 10, num_minimizersToFind = 100,
initial_lowerRange = c(0.1,0.1,0.1), initial_upperRange =  c(10,10,10),
saveLog = FALSE)

acceptedApproximateMinimizers(CGNM_result)

\dontrun{
library(RxODE)

model_text="
d/dt(X_1)=-ka*X_1
d/dt(C_2)=(ka*X_1-CL_2*C_2)/V1"

model=RxODE(model_text)
#define nonlinearFunction
model_function=function(x){

observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)

theta <- c(ka=x[1],V1=x[2],CL_2=x[3])
ev <- eventTable()
ev$add.dosing(dose = 1000, start.time =0)
ev$add.sampling(observation_time)
odeSol=model$solve(theta, ev)
log10(odeSol[,"C_2"])

}

observation=log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))

CGNM_result=Cluster_Gauss_Newton_method(nonlinearFunction=model_function,
targetVector = observation, saveLog = FALSE,
initial_lowerRange = c(0.1,0.1,0.1),initial_upperRange =  c(10,10,10))}

}
