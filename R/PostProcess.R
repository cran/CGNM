#' @title acceptedMaxSSR
#' @description
#' CGNM find multiple sets of minimizers of the nonlinear least squares (nls) problem by solving nls from various initial iterates.  Although CGNM is shown to be robust compared to other conventional multi-start algorithms, not all initial iterates minimizes successfully.  By assuming sum of squares residual (SSR) follows the chai-square distribution we first reject the approximated minimiser who SSR is statistically significantly worse than the minimum SSR found by the CGNM.  Then use elbow-method (a heuristic often used in mathematical optimisation to balance the quality and the quantity of the solution found) to find the "acceptable" maximum SSR.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param cutoff_pvalue (default: 0.05) \emph{A number} defines the rejection p-value for the first stage of acceptable computational result screening.
#' @param numParametersIncluded (default: NA) \emph{A natural number} defines the number of parameter sets to be included in the assessment of the acceptable parameters.  If set NA then use all the parameters found by the CGNM.
#' @return \emph{A positive real number} that is the maximum sum of squares residual (SSR) the algorithm has selected to accept.
#' @examples
#' library(parallel) #if parallel library is loaded CGNM paralleizes the algorithm automatically
#'
#'model_analytic_function=function(x){
#'
#'  observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
#'  Dose=1000
#'  F=1
#'
#'  ka=x[1]
#'  V1=x[2]
#'  CL_2=x[3]
#'  t=observation_time
#'
#'  Cp=ka*F*Dose/(V1*(ka-CL_2/V1))*(exp(-CL_2/V1*t)-exp(-ka*t))
#'
#'  log10(Cp)
#'}
#'
#' observation=log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))
#'
#' CGNM_result=Cluster_Gauss_Newton_method(
#' nonlinearFunction=model_analytic_function,
#' targetVector = observation,
#' initial_lowerRange = c(0.1,0.1,0.1), initial_upperRange =  c(10,10,10),
#' num_iter = 10, num_minimizersToFind = 100)
#'
#' acceptedMaxSSR(CGNM_result)
#' @export
acceptedMaxSSR=function(CGNM_result, cutoff_pvalue=0.05, numParametersIncluded=NA){

SSR_vec=CGNM_result$residual_history[,dim(CGNM_result$residual_history)[2]]
min_R=min(SSR_vec)
minIndex=which(SSR_vec==min_R)[1]
residual_vec=CGNM_result$Y[minIndex,]-CGNM_result$runSetting$targetVector

acceptMaxSSR=qchisq(1-cutoff_pvalue, df=length(residual_vec)) *(sd(residual_vec))^2+min_R

accept_index=which(SSR_vec<acceptMaxSSR)

accept_vec=as.vector(SSR_vec<acceptMaxSSR)

numAccept=sum(accept_vec, na.rm = TRUE)



  if(!is.na(numParametersIncluded)){
    sortedAcceptedSSR_value=sort(SSR_vec[accept_vec])[seq(1,numParametersIncluded)]
  }else{
    sortedAcceptedSSR=sort(SSR_vec[accept_vec])
  }


  trapizoido_area=c()
  for(i in seq(1,length(sortedAcceptedSSR)-1))
    trapizoido_area=c(trapizoido_area, (sortedAcceptedSSR[1]+sortedAcceptedSSR[i])*i+(sortedAcceptedSSR[i+1]+sortedAcceptedSSR[length(sortedAcceptedSSR)])*(length(sortedAcceptedSSR)-i))

  strictMaxAcceptSSR=sortedAcceptedSSR[which(trapizoido_area==min(trapizoido_area))]

  accept_index=which(SSR_vec<strictMaxAcceptSSR)
  accept_vec=as.vector(SSR_vec<strictMaxAcceptSSR)
  numAccept=sum(accept_vec, na.rm = TRUE)
  acceptMaxSSR=strictMaxAcceptSSR

  return(acceptMaxSSR)
}

#' @title acceptedApproximateMinimizers
#' @description
#' CGNM find multiple sets of minimizers of the nonlinear least squares (nls) problem by solving nls from various initial iterates.  Although CGNM is shown to be robust compared to other conventional multi-start algorithms, not all initial iterates minimizes successfully.  By assuming sum of squares residual (SSR) follows the chai-square distribution we first reject the approximated minimiser who SSR is statistically significantly worse than the minimum SSR found by the CGNM.  Then use elbow-method (a heuristic often used in mathematical optimisation to balance the quality and the quantity of the solution found) to find the "acceptable" maximum SSR. This function outputs the acceptable approximate minimizers of the nonlinear least squares problem found by the CGNM.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param cutoff_pvalue (default: 0.05) \emph{A number} defines the rejection p-value for the first stage of acceptable computational result screening.
#' @param numParametersIncluded (default: NA) \emph{A natural number} defines the number of parameter sets to be included in the assessment of the acceptable parameters.  If set NA then use all the parameters found by the CGNM.
#' @return \emph{A matrix} that each row stores the accepted approximate minimizers found by CGNM.
#' @examples
#' library(parallel) #if parallel library is loaded CGNM paralleizes the algorithm automatically
#'
#'model_analytic_function=function(x){
#'
#'  observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
#'  Dose=1000
#'  F=1
#'
#'  ka=x[1]
#'  V1=x[2]
#'  CL_2=x[3]
#'  t=observation_time
#'
#'  Cp=ka*F*Dose/(V1*(ka-CL_2/V1))*(exp(-CL_2/V1*t)-exp(-ka*t))
#'
#'  log10(Cp)
#'}
#'
#' observation=log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))
#'
#' CGNM_result=Cluster_Gauss_Newton_method(
#' nonlinearFunction=model_analytic_function,
#' targetVector = observation,
#' initial_lowerRange = c(0.1,0.1,0.1), initial_upperRange =  c(10,10,10),
#' num_iter = 10, num_minimizersToFind = 100)
#'
#' acceptedApproximateMinimizers(CGNM_result)
#' @export
acceptedApproximateMinimizers=function(CGNM_result, cutoff_pvalue=0.05, numParametersIncluded=NA){
  CGNM_result$X[acceptedIndecies(CGNM_result, cutoff_pvalue, numParametersIncluded),]
}

#' @title acceptedIndecies
#' @description
#' CGNM find multiple sets of minimizers of the nonlinear least squares (nls) problem by solving nls from various initial iterates.  Although CGNM is shown to be robust compared to other conventional multi-start algorithms, not all initial iterates minimizes successfully.  By assuming sum of squares residual (SSR) follows the chai-square distribution we first reject the approximated minimiser who SSR is statistically significantly worse than the minimum SSR found by the CGNM.  Then use elbow-method (a heuristic often used in mathematical optimisation to balance the quality and the quantity of the solution found) to find the "acceptable" maximum SSR. This function outputs the indices of acceptable approximate minimizers of the nonlinear least squares problem found by the CGNM.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param cutoff_pvalue (default: 0.05) \emph{A number} defines the rejection p-value for the first stage of acceptable computational result screening.
#' @param numParametersIncluded (default: NA) \emph{A natural number} defines the number of parameter sets to be included in the assessment of the acceptable parameters.  If set NA then use all the parameters found by the CGNM.
#' @return \emph{A vector of natural number} that contains the indices of accepted approximate minimizers found by CGNM.
#' @examples
#' library(parallel) #if parallel library is loaded CGNM paralleizes the algorithm automatically
#'
#'model_analytic_function=function(x){
#'
#'  observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
#'  Dose=1000
#'  F=1
#'
#'  ka=x[1]
#'  V1=x[2]
#'  CL_2=x[3]
#'  t=observation_time
#'
#'  Cp=ka*F*Dose/(V1*(ka-CL_2/V1))*(exp(-CL_2/V1*t)-exp(-ka*t))
#'
#'  log10(Cp)
#'}
#'
#' observation=log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))
#'
#' CGNM_result=Cluster_Gauss_Newton_method(
#' nonlinearFunction=model_analytic_function,
#' targetVector = observation,
#' initial_lowerRange = c(0.1,0.1,0.1), initial_upperRange =  c(10,10,10),
#' num_iter = 10, num_minimizersToFind = 100)
#'
#' acceptedIndecies(CGNM_result)
#' @export
acceptedIndecies=function(CGNM_result, cutoff_pvalue=0.05, numParametersIncluded=NA){
  acceptMaxSSR_value=acceptedMaxSSR(CGNM_result, cutoff_pvalue, numParametersIncluded)

  SSR_vec=CGNM_result$residual_history[,dim(CGNM_result$residual_history)[2]]

  return(which(SSR_vec<acceptMaxSSR_value))
}

#' @title acceptedIndecies_binary
#' @description
#' CGNM find multiple sets of minimizers of the nonlinear least squares (nls) problem by solving nls from various initial iterates.  Although CGNM is shown to be robust compared to other conventional multi-start algorithms, not all initial iterates minimizes successfully.  By assuming sum of squares residual (SSR) follows the chai-square distribution we first reject the approximated minimiser who SSR is statistically significantly worse than the minimum SSR found by the CGNM.  Then use elbow-method (a heuristic often used in mathematical optimisation to balance the quality and the quantity of the solution found) to find the "acceptable" maximum SSR. This function outputs the indices of acceptable approximate minimizers of the nonlinear least squares problem found by the CGNM. (note that acceptedIndecies(CGNM_result) is equal to seq(1,length(acceptedIndecies_binary(CGNM_result)))[acceptedIndecies_binary(CGNM_result)])
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param cutoff_pvalue (default: 0.05) \emph{A number} defines the rejection p-value for the first stage of acceptable computational result screening.
#' @param numParametersIncluded (default: NA) \emph{A natural number} defines the number of parameter sets to be included in the assessment of the acceptable parameters.  If set NA then use all the parameters found by the CGNM.
#' @return \emph{A vector of TRUE and FALSE} that indicate if the each of the approximate minimizer found by CGNM is acceptable or not.
#' @examples
#' library(parallel) #if parallel library is loaded CGNM paralleizes the algorithm automatically
#'
#'model_analytic_function=function(x){
#'
#'  observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
#'  Dose=1000
#'  F=1
#'
#'  ka=x[1]
#'  V1=x[2]
#'  CL_2=x[3]
#'  t=observation_time
#'
#'  Cp=ka*F*Dose/(V1*(ka-CL_2/V1))*(exp(-CL_2/V1*t)-exp(-ka*t))
#'
#'  log10(Cp)
#'}
#'
#' observation=log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))
#'
#' CGNM_result=Cluster_Gauss_Newton_method(
#' nonlinearFunction=model_analytic_function,
#' targetVector = observation,
#' initial_lowerRange = c(0.1,0.1,0.1), initial_upperRange =  c(10,10,10),
#' num_iter = 10, num_minimizersToFind = 100)
#'
#' acceptedIndecies_binary(CGNM_result)
#' @export
acceptedIndecies_binary=function(CGNM_result, cutoff_pvalue=0.05, numParametersIncluded=NA){
  acceptMaxSSR_value=acceptedMaxSSR(CGNM_result, cutoff_pvalue, numParametersIncluded)

  SSR_vec=CGNM_result$residual_history[,dim(CGNM_result$residual_history)[2]]

  return(as.vector(SSR_vec<acceptMaxSSR_value))
}


#' @title plot_Rank_SSR
#' @description
#' Make SSR v.s. rank plot. This plot is often used to visualize the maximum accepted SSR.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param cutoff_pvalue (default: 0.05) \emph{A number} defines the rejection p-value for the first stage of acceptable computational result screening.
#' @param numParametersIncluded (default: NA) \emph{A natural number} defines the number of parameter sets to be included in the assessment of the acceptable parameters.  If set NA then use all the parameters found by the CGNM.
#' @return \emph{A ggplot object} of SSR v.s. rank.
#' @examples
#' library(parallel) #if parallel library is loaded CGNM paralleizes the algorithm automatically
#'
#'model_analytic_function=function(x){
#'
#'  observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
#'  Dose=1000
#'  F=1
#'
#'  ka=x[1]
#'  V1=x[2]
#'  CL_2=x[3]
#'  t=observation_time
#'
#'  Cp=ka*F*Dose/(V1*(ka-CL_2/V1))*(exp(-CL_2/V1*t)-exp(-ka*t))
#'
#'  log10(Cp)
#'}
#'
#' observation=log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))
#'
#' CGNM_result=Cluster_Gauss_Newton_method(
#' nonlinearFunction=model_analytic_function,
#' targetVector = observation,
#' initial_lowerRange = c(0.1,0.1,0.1), initial_upperRange =  c(10,10,10),
#' num_iter = 10, num_minimizersToFind = 100)
#'
#' plot_Rank_SSR(CGNM_result)
#' @export
#' @import ggplot2

plot_Rank_SSR=function(CGNM_result, cutoff_pvalue=0.05, numParametersIncluded=NA){
  SSR_value=NULL
  is_accepted=NULL
  SSR_vec=CGNM_result$residual_history[,dim(CGNM_result$residual_history)[2]]

  acceptedIndecies_b=acceptedIndecies_binary(CGNM_result,cutoff_pvalue, numParametersIncluded)

  numAccept=sum(acceptedIndecies_b)
  acceptMaxSSR=acceptedMaxSSR(CGNM_result,cutoff_pvalue, numParametersIncluded)
  min_R=min(SSR_vec)

  plot_df=data.frame(SSR_value=SSR_vec, rank=rank(SSR_vec), is_accepted=acceptedIndecies_b)
  if("ggplot2" %in% tolower((.packages()))){
  ggplot(plot_df, aes(x=rank,y=SSR_value, colour=is_accepted))+geom_point()+coord_cartesian(ylim=c(0,acceptMaxSSR*2))+geom_vline(xintercept = numAccept, color="grey")+ annotate(geom="text", x=numAccept, y=acceptMaxSSR*0.5, label=paste("Accepted: ",numAccept,"\n Accepted max SSR: ",formatC(acceptMaxSSR, format = "g", digits = 3)),angle = 90,
                                                                                                                                                                       color="black")+ annotate(geom="text", x=length(SSR_vec)*0.1, y=min_R*1.1, label=paste("min SSR: ",formatC(min_R, format = "g", digits = 3)),
                                                                                                                                                                                                color="black")


  }else{
    print("ggplot2 package needs to be loaded before using this function")
  }
}


#' @title plot_paraDistribution_byViolinPlots
#' @description
#' Make violin plot to compare the initial distribution and distribition of the accepted approximate minimizers found by the CGNM. Bars in the violin plots indicates the interquartile range. The solid line connects the interquartile ranges of the initial distribution and the distribution of the accepted approximate minimizer at the final iterate.  The blacklines connets the minimums and maximums of the initial distribution and the distribution of the accepted approximate minimizer at the final iterate. The black dots indicate the median.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param cutoff_pvalue (default: 0.05) \emph{A number} defines the rejection p-value for the first stage of acceptable computational result screening.
#' @param numParametersIncluded (default: NA) \emph{A natural number} defines the number of parameter sets to be included in the assessment of the acceptable parameters.  If set NA then use all the parameters found by the CGNM.
#' @param ParameterNames (default: NA) \emph{A vector of string} the user can supply so that these names are used when making the plot. (Note if it set as NA or vector of incorrect length then the parameters are named as x_1, x_2, ...)
#' @return \emph{A ggplot object} including the violin plot, interquartile range and median, minimum and maximum.
#' @examples
#' library(parallel) #if parallel library is loaded CGNM paralleizes the algorithm automatically
#'
#'model_analytic_function=function(x){
#'
#'  observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
#'  Dose=1000
#'  F=1
#'
#'  ka=x[1]
#'  V1=x[2]
#'  CL_2=x[3]
#'  t=observation_time
#'
#'  Cp=ka*F*Dose/(V1*(ka-CL_2/V1))*(exp(-CL_2/V1*t)-exp(-ka*t))
#'
#'  log10(Cp)
#'}
#'
#' observation=log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))
#'
#' CGNM_result=Cluster_Gauss_Newton_method(
#' nonlinearFunction=model_analytic_function,
#' targetVector = observation,
#' initial_lowerRange = c(0.1,0.1,0.1), initial_upperRange =  c(10,10,10),
#' num_iter = 10, num_minimizersToFind = 100)
#'
#' plot_paraDistribution_byViolinPlots(CGNM_result)
#' @export
#' @import ggplot2

plot_paraDistribution_byViolinPlots=function(CGNM_result, cutoff_pvalue=0.05, numParametersIncluded=NA, ParameterNames=NA){
  Kind_iter=NULL
  X_value=NULL
  freeParaValues=data.frame()
  acceptedIndecies_b=acceptedIndecies_binary(CGNM_result,cutoff_pvalue, numParametersIncluded)

  if(length(ParameterNames)==dim(CGNM_result$initialX)[2]){

  }else{
    ParameterNames=paste0("x_",seq(1,dim(CGNM_result$initialX)[2]))
  }

    for(i in seq(1,dim(CGNM_result$initialX)[2])){
      freeParaValues=rbind(freeParaValues, data.frame(Name=ParameterNames[i],X_value=CGNM_result$initialX[,i], Kind_iter="Initial"))
    }

    for(i in seq(1,dim(CGNM_result$X)[2])){
      freeParaValues=rbind(freeParaValues, data.frame(Name=ParameterNames[i],X_value=CGNM_result$X[acceptedIndecies_b,i], Kind_iter="Final Accepted"))
    }

  freeParaValues$Kind_iter=factor(freeParaValues$Kind_iter, levels = c("Initial", "Final Accepted"))
  freeParaValues$Name=factor(freeParaValues$Name, levels = ParameterNames)

  if("ggplot2" %in% tolower((.packages()))){

    p<-ggplot(freeParaValues,aes(x=Kind_iter,y=X_value))+facet_wrap(Name~., scales = "free")+scale_y_continuous(trans = 'log10')

    p+geom_violin(trim=T,fill="#999999",linetype="blank",alpha=I(1/2))+
      stat_summary(geom="pointrange",fun = median, fun.min = function(x) quantile(x,probs=0.25), fun.max = function(x) quantile(x,probs=0.75), size=0.5,alpha=.5)+
      stat_summary(geom="line",fun =  function(x) quantile(x,probs=0), aes(group=1),size=0.5,alpha=.3,linetype=2)+
      stat_summary(geom="line",fun =  function(x) quantile(x,probs=1), aes(group=1),size=0.5,alpha=.3,linetype=2)+
      stat_summary(geom="line",fun =  function(x) quantile(x,probs=0.25), aes(group=1),size=0.5,alpha=.3)+
      stat_summary(geom="line",fun =  function(x) quantile(x,probs=0.75), aes(group=1),size=0.5,alpha=.3)+
      theme(legend.position="none")+xlab("")

  }else{
    print("ggplot2 package needs to be loaded before using this function")
  }
}



