## TODO make autodiagnosis function, is Lambda all big?  Is the minimum SSR parameter set within the initial range?
## TODO implement best parameter
## TODO implement worst accepted parameter

#' @title acceptedMaxSSR
#' @description
#' CGNM find multiple sets of minimizers of the nonlinear least squares (nls) problem by solving nls from various initial iterates.  Although CGNM is shown to be robust compared to other conventional multi-start algorithms, not all initial iterates minimizes successfully.  By assuming sum of squares residual (SSR) follows the chai-square distribution we first reject the approximated minimiser who SSR is statistically significantly worse than the minimum SSR found by the CGNM.  Then use elbow-method (a heuristic often used in mathematical optimisation to balance the quality and the quantity of the solution found) to find the "acceptable" maximum SSR.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param cutoff_pvalue (default: 0.05) \emph{A number} defines the rejection p-value for the first stage of acceptable computational result screening.
#' @param numParametersIncluded (default: NA) \emph{A natural number} defines the number of parameter sets to be included in the assessment of the acceptable parameters.  If set NA then use all the parameters found by the CGNM.
#' @param useAcceptedApproximateMinimizers (default: TRUE) \emph{TRUE or FALSE} If true then use chai-square and elbow method to choose maximum accepted SSR.  If false returnsnumParametersIncluded-th smallest SSR (or if numParametersIncluded=NA then returns the largest SSR).
#' @return \emph{A positive real number} that is the maximum sum of squares residual (SSR) the algorithm has selected to accept.
#' @examples
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
acceptedMaxSSR=function(CGNM_result, cutoff_pvalue=0.05, numParametersIncluded=NA, useAcceptedApproximateMinimizers=TRUE){

SSR_vec=CGNM_result$residual_history[,dim(CGNM_result$residual_history)[2]]

if(useAcceptedApproximateMinimizers){
  min_R=min(SSR_vec)
  minIndex=which(SSR_vec==min_R)[1]

  targetVector=as.numeric(CGNM_result$runSetting$targetVector)
  targetVector[is.na(targetVector)]=0

  residual_vec=CGNM_result$Y[minIndex,]-targetVector

  acceptMaxSSR=max(qchisq(1-cutoff_pvalue, df=length(residual_vec)) *(sd(residual_vec))^2+min_R, sqrt(sqrt(.Machine$double.eps)))

  accept_index=which(SSR_vec<acceptMaxSSR)

  accept_vec=as.vector(SSR_vec<acceptMaxSSR)

  numAccept=sum(accept_vec, na.rm = TRUE)



  if(!is.na(numParametersIncluded)&(sum(accept_vec)>numParametersIncluded)){
    sortedAcceptedSSR=sort(SSR_vec[accept_vec])[seq(1,numParametersIncluded)]
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

}else if(is.na(numParametersIncluded)|numParametersIncluded>length(SSR_vec)){
  acceptMaxSSR=max(SSR_vec, na.rm = TRUE)
}else{
  acceptMaxSSR=sort(SSR_vec)[numParametersIncluded]
}

  return(acceptMaxSSR)
}

#' @title acceptedApproximateMinimizers
#' @description
#' CGNM find multiple sets of minimizers of the nonlinear least squares (nls) problem by solving nls from various initial iterates.  Although CGNM is shown to be robust compared to other conventional multi-start algorithms, not all initial iterates minimizes successfully.  By assuming sum of squares residual (SSR) follows the chai-square distribution we first reject the approximated minimiser who SSR is statistically significantly worse than the minimum SSR found by the CGNM.  Then use elbow-method (a heuristic often used in mathematical optimisation to balance the quality and the quantity of the solution found) to find the "acceptable" maximum SSR. This function outputs the acceptable approximate minimizers of the nonlinear least squares problem found by the CGNM.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param cutoff_pvalue (default: 0.05) \emph{A number} defines the rejection p-value for the first stage of acceptable computational result screening.
#' @param numParametersIncluded (default: NA) \emph{A natural number} defines the number of parameter sets to be included in the assessment of the acceptable parameters.  If set NA then use all the parameters found by the CGNM.
#' @param useAcceptedApproximateMinimizers (default: TRUE) \emph{TRUE or FALSE} If true then use chai-square and elbow method to choose maximum accepted SSR.  If false returns the parameters upto numParametersIncluded-th smallest SSR (or if numParametersIncluded=NA then use all the parameters found by the CGNM).
#' @return \emph{A matrix} that each row stores the accepted approximate minimizers found by CGNM.
#' @examples
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
acceptedApproximateMinimizers=function(CGNM_result, cutoff_pvalue=0.05, numParametersIncluded=NA, useAcceptedApproximateMinimizers=TRUE){
  CGNM_result$X[acceptedIndecies(CGNM_result, cutoff_pvalue, numParametersIncluded, useAcceptedApproximateMinimizers),]
}

#' @title acceptedIndecies
#' @description
#' CGNM find multiple sets of minimizers of the nonlinear least squares (nls) problem by solving nls from various initial iterates.  Although CGNM is shown to be robust compared to other conventional multi-start algorithms, not all initial iterates minimizes successfully.  By assuming sum of squares residual (SSR) follows the chai-square distribution we first reject the approximated minimiser who SSR is statistically significantly worse than the minimum SSR found by the CGNM.  Then use elbow-method (a heuristic often used in mathematical optimisation to balance the quality and the quantity of the solution found) to find the "acceptable" maximum SSR. This function outputs the indices of acceptable approximate minimizers of the nonlinear least squares problem found by the CGNM.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param cutoff_pvalue (default: 0.05) \emph{A number} defines the rejection p-value for the first stage of acceptable computational result screening.
#' @param numParametersIncluded (default: NA) \emph{A natural number} defines the number of parameter sets to be included in the assessment of the acceptable parameters.  If set NA then use all the parameters found by the CGNM.
#' @param useAcceptedApproximateMinimizers (default: TRUE) \emph{TRUE or FALSE} If true then use chai-square and elbow method to choose maximum accepted SSR.  If false returns the parameters upto numParametersIncluded-th smallest SSR (or if numParametersIncluded=NA then use all the parameters found by the CGNM).
#' @return \emph{A vector of natural number} that contains the indices of accepted approximate minimizers found by CGNM.
#' @examples
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
acceptedIndecies=function(CGNM_result, cutoff_pvalue=0.05, numParametersIncluded=NA, useAcceptedApproximateMinimizers=TRUE){
  acceptMaxSSR_value=acceptedMaxSSR(CGNM_result, cutoff_pvalue, numParametersIncluded, useAcceptedApproximateMinimizers)

  SSR_vec=CGNM_result$residual_history[,dim(CGNM_result$residual_history)[2]]

  return(which(SSR_vec<=acceptMaxSSR_value))
}

#' @title acceptedIndecies_binary
#' @description
#' CGNM find multiple sets of minimizers of the nonlinear least squares (nls) problem by solving nls from various initial iterates.  Although CGNM is shown to be robust compared to other conventional multi-start algorithms, not all initial iterates minimizes successfully.  By assuming sum of squares residual (SSR) follows the chai-square distribution we first reject the approximated minimiser who SSR is statistically significantly worse than the minimum SSR found by the CGNM.  Then use elbow-method (a heuristic often used in mathematical optimisation to balance the quality and the quantity of the solution found) to find the "acceptable" maximum SSR. This function outputs the indices of acceptable approximate minimizers of the nonlinear least squares problem found by the CGNM. (note that acceptedIndecies(CGNM_result) is equal to seq(1,length(acceptedIndecies_binary(CGNM_result)))[acceptedIndecies_binary(CGNM_result)])
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param cutoff_pvalue (default: 0.05) \emph{A number} defines the rejection p-value for the first stage of acceptable computational result screening.
#' @param numParametersIncluded (default: NA) \emph{A natural number} defines the number of parameter sets to be included in the assessment of the acceptable parameters.  If set NA then use all the parameters found by the CGNM.
#' @param useAcceptedApproximateMinimizers (default: TRUE) \emph{TRUE or FALSE} If true then use chai-square and elbow method to choose maximum accepted SSR.  If false returns the indicies upto numParametersIncluded-th smallest SSR (or if numParametersIncluded=NA then use all the parameters found by the CGNM).
#' @return \emph{A vector of TRUE and FALSE} that indicate if the each of the approximate minimizer found by CGNM is acceptable or not.
#' @examples
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
acceptedIndecies_binary=function(CGNM_result, cutoff_pvalue=0.05, numParametersIncluded=NA, useAcceptedApproximateMinimizers=TRUE){
  acceptMaxSSR_value=acceptedMaxSSR(CGNM_result, cutoff_pvalue, numParametersIncluded, useAcceptedApproximateMinimizers)

  SSR_vec=CGNM_result$residual_history[,dim(CGNM_result$residual_history)[2]]

  return(as.vector(SSR_vec<=acceptMaxSSR_value))
}


#' @title plot_Rank_SSR
#' @description
#' Make SSR v.s. rank plot. This plot is often used to visualize the maximum accepted SSR.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param cutoff_pvalue (default: 0.05) \emph{A number} defines the rejection p-value for the first stage of acceptable computational result screening.
#' @param numParametersIncluded (default: NA) \emph{A natural number} defines the number of parameter sets to be included in the assessment of the acceptable parameters.  If set NA then use all the parameters found by the CGNM.
#' @param useAcceptedApproximateMinimizers (default: TRUE) \emph{TRUE or FALSE} If true then use chai-square and elbow method to choose maximum accepted SSR.  If false use parameters upto numParametersIncluded-th smallest SSR (or if numParametersIncluded=NA then use all the parameters found by the CGNM).
#' @return \emph{A ggplot object} of SSR v.s. rank.
#' @examples
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
plot_Rank_SSR=function(CGNM_result, cutoff_pvalue=0.05, numParametersIncluded=NA, useAcceptedApproximateMinimizers=TRUE){
  SSR_value=NULL
  is_accepted=NULL
  SSR_vec=CGNM_result$residual_history[,dim(CGNM_result$residual_history)[2]]

  acceptedIndecies_b=acceptedIndecies_binary(CGNM_result,cutoff_pvalue, numParametersIncluded, useAcceptedApproximateMinimizers)

  numAccept=sum(acceptedIndecies_b)
  acceptMaxSSR=acceptedMaxSSR(CGNM_result,cutoff_pvalue, numParametersIncluded, useAcceptedApproximateMinimizers)
  min_R=min(SSR_vec)

  plot_df=data.frame(SSR_value=SSR_vec, rank=rank(SSR_vec), is_accepted=acceptedIndecies_b)

  ggplot2::ggplot(plot_df, ggplot2::aes(x=rank,y=SSR_value, colour=is_accepted))+ggplot2::geom_point()+ggplot2::coord_cartesian(ylim=c(0,acceptMaxSSR*2))+ggplot2::geom_vline(xintercept = numAccept, color="grey")+ ggplot2::annotate(geom="text", x=numAccept, y=acceptMaxSSR*0.5, label=paste("Accepted: ",numAccept,"\n Accepted max SSR: ",formatC(acceptMaxSSR, format = "g", digits = 3)),angle = 90,
                                                                                                                                                                       color="black")+ ggplot2::annotate(geom="text", x=length(SSR_vec)*0.1, y=min_R*1.1, label=paste("min SSR: ",formatC(min_R, format = "g", digits = 3)),
                                                                                                                                                                                                color="black")



}


#' @title plot_paraDistribution_byViolinPlots
#' @description
#' Make violin plot to compare the initial distribution and distribition of the accepted approximate minimizers found by the CGNM. Bars in the violin plots indicates the interquartile range. The solid line connects the interquartile ranges of the initial distribution and the distribution of the accepted approximate minimizer at the final iterate.  The blacklines connets the minimums and maximums of the initial distribution and the distribution of the accepted approximate minimizer at the final iterate. The black dots indicate the median.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param cutoff_pvalue (default: 0.05) \emph{A number} defines the rejection p-value for the first stage of acceptable computational result screening.
#' @param numParametersIncluded (default: NA) \emph{A natural number} defines the number of parameter sets to be included in the assessment of the acceptable parameters.  If set NA then use all the parameters found by the CGNM.
#' @param ParameterNames (default: NA) \emph{A vector of string} the user can supply so that these names are used when making the plot. (Note if it set as NA or vector of incorrect length then the parameters are named as x_1, x_2, ...)
#' @param useAcceptedApproximateMinimizers (default: TRUE) \emph{TRUE or FALSE} If true then use chai-square and elbow method to choose maximum accepted SSR.  If false use parameters upto numParametersIncluded-th smallest SSR (or if numParametersIncluded=NA then use all the parameters found by the CGNM).
#' @return \emph{A ggplot object} including the violin plot, interquartile range and median, minimum and maximum.
#' @examples
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

plot_paraDistribution_byViolinPlots=function(CGNM_result, cutoff_pvalue=0.05, numParametersIncluded=NA, ParameterNames=NA, useAcceptedApproximateMinimizers=TRUE){
  Kind_iter=NULL
  X_value=NULL
  cluster=NULL
  freeParaValues=data.frame()

  if(useAcceptedApproximateMinimizers){
    useIndecies_b=acceptedIndecies_binary(CGNM_result,cutoff_pvalue, numParametersIncluded, useAcceptedApproximateMinimizers)
  }else{
    if(is.na(numParametersIncluded)|numParametersIncluded>dim(CGNM_result$X)[1]){
      useIndecies_b=rep(TRUE, dim(CGNM_result$X)[1])
    }else{
      SSR_vec=CGNM_result$residual_history[,dim(CGNM_result$residual_history)[2]]
      useIndecies_b=(SSR_vec<=sort(SSR_vec)[numParametersIncluded])
    }
  }

  if(length(ParameterNames)==dim(CGNM_result$initialX)[2]){

  }else{
    ParameterNames=paste0("x_",seq(1,dim(CGNM_result$initialX)[2]))
  }

    for(i in seq(1,dim(CGNM_result$initialX)[2])){
      freeParaValues=rbind(freeParaValues, data.frame(Name=ParameterNames[i],X_value=CGNM_result$initialX[,i], Kind_iter="Initial"))
    }

    for(i in seq(1,dim(CGNM_result$X)[2])){
      freeParaValues=rbind(freeParaValues, data.frame(Name=ParameterNames[i],X_value=CGNM_result$X[useIndecies_b,i], Kind_iter="Final Accepted"))
    }

  freeParaValues$Kind_iter=factor(freeParaValues$Kind_iter, levels = c("Initial", "Final Accepted"))
  freeParaValues$Name=factor(freeParaValues$Name, levels = ParameterNames)


    p<-ggplot2::ggplot(freeParaValues,ggplot2::aes(x=Kind_iter,y=X_value))+ggplot2::facet_wrap(Name~., scales = "free")

    p+ggplot2::geom_violin(trim=T,fill="#999999",linetype="blank",alpha=I(1/2))+
      ggplot2::stat_summary(geom="pointrange",fun = median, fun.min = function(x) quantile(x,probs=0.25), fun.max = function(x) quantile(x,probs=0.75), size=0.5,alpha=.5)+
      ggplot2::stat_summary(geom="line",fun =  function(x) quantile(x,probs=0), ggplot2::aes(group=1),size=0.5,alpha=.3,linetype=2)+
      ggplot2::stat_summary(geom="line",fun =  function(x) quantile(x,probs=1), ggplot2::aes(group=1),size=0.5,alpha=.3,linetype=2)+
      ggplot2::stat_summary(geom="line",fun =  function(x) quantile(x,probs=0.25), ggplot2::aes(group=1),size=0.5,alpha=.3)+
      ggplot2::stat_summary(geom="line",fun =  function(x) quantile(x,probs=0.75), ggplot2::aes(group=1),size=0.5,alpha=.3)+
      ggplot2::theme(legend.position="none")+xlab("")


}


#' @title plot_SSR_parameterValue
#' @description
#' Make SSR v.s. parameterValue plot of the accepted approximate minimizers found by the CGNM. Bars in the violin plots indicates the interquartile range.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param cutoff_pvalue (default: 0.05) \emph{A number} defines the rejection p-value for the first stage of acceptable computational result screening.
#' @param numParametersIncluded (default: NA) \emph{A natural number} defines the number of parameter sets to be included in the assessment of the acceptable parameters.  If set NA then use all the parameters found by the CGNM.
#' @param ParameterNames (default: NA) \emph{A vector of string} the user can supply so that these names are used when making the plot. (Note if it set as NA or vector of incorrect length then the parameters are named as x_1, x_2, ...)
#' @param showInitialRange (default: TRUE) \emph{TRUE or FALSE} if TRUE then the initial range appears in the plot.
#' @param useAcceptedApproximateMinimizers (default: TRUE) \emph{TRUE or FALSE} If true then use chai-square and elbow method to choose maximum accepted SSR.  If false use parameters upto numParametersIncluded-th smallest SSR (or if numParametersIncluded=NA then use all the parameters found by the CGNM).
#' @return \emph{A ggplot object} including the violin plot, interquartile range and median, minimum and maximum.
#' @examples
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
#' plot_SSR_parameterValue(CGNM_result)
#' @export
#' @import ggplot2
plot_SSR_parameterValue=function(CGNM_result, cutoff_pvalue=0.05, numParametersIncluded=NA, ParameterNames=NA, showInitialRange=TRUE, useAcceptedApproximateMinimizers=TRUE){
  SSR=NULL
  X_value=NULL
  cluster=NULL
  boundValue=NULL

  lastIter=dim(CGNM_result$residual_history)[2]

 # if(useAcceptedApproximateMinimizers){
    useIndecies_b=acceptedIndecies_binary(CGNM_result,cutoff_pvalue, numParametersIncluded, useAcceptedApproximateMinimizers)
  # }else{
  #   if(is.na(numParametersIncluded)|numParametersIncluded>dim(CGNM_result$X)[1]){
  #     useIndecies_b=rep(TRUE, dim(CGNM_result$X)[1])
  #   }else{
  #     SSR_vec=CGNM_result$residual_history[,dim(CGNM_result$residual_history)[2]]
  #     useIndecies_b=(SSR_vec<=sort(SSR_vec)[numParametersIncluded])
  #   }
  # }

  if(length(ParameterNames)==dim(CGNM_result$initialX)[2]){

  }else{
    ParameterNames=paste0("x_",seq(1,dim(CGNM_result$initialX)[2]))
  }

  freeParaValues=data.frame()

  for(i in seq(1,dim(CGNM_result$X)[2])){
    freeParaValues=rbind(freeParaValues, data.frame(Name=ParameterNames[i],X_value=CGNM_result$X[useIndecies_b,i], SSR=CGNM_result$residual_history[useIndecies_b,lastIter]))
  }

  freeParaValues$Name=factor(freeParaValues$Name, levels = ParameterNames)

  #kmeans_result=optimal_kmeans(cbind(CGNM_result$X[useIndecies_b,], freeParaValues$SSR))
  kmeans_result=optimal_kmeans(CGNM_result$X)

  freeParaValues$cluster=as.factor(kmeans_result$cluster[useIndecies_b])

  g=ggplot2::ggplot(freeParaValues, ggplot2::aes(y=SSR, x=X_value, colour=cluster))+ggplot2::geom_point(alpha=0.3)

  if(showInitialRange){
    g=g+ggplot2::geom_vline(data=data.frame(boundValue=CGNM_result$runSetting$initial_lowerRange, Name=ParameterNames), ggplot2::aes(xintercept=boundValue), colour="gray")
    g=g+ggplot2::geom_vline(data=data.frame(boundValue=CGNM_result$runSetting$initial_upperRange, Name=ParameterNames), ggplot2::aes(xintercept=boundValue), colour="gray")
  }

  g+ggplot2::facet_wrap(.~Name,scales="free")
}



#' @title plot_parameterValue_scatterPlots
#' @description
#' Make scatter plots of the accepted approximate minimizers found by the CGNM. Bars in the violin plots indicates the interquartile range.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param cutoff_pvalue (default: 0.05) \emph{A number} defines the rejection p-value for the first stage of acceptable computational result screening.
#' @param numParametersIncluded (default: NA) \emph{A natural number} defines the number of parameter sets to be included in the assessment of the acceptable parameters.  If set NA then use all the parameters found by the CGNM.
#' @param ParameterNames (default: NA) \emph{A vector of string} the user can supply so that these names are used when making the plot. (Note if it set as NA or vector of incorrect length then the parameters are named as x_1, x_2, ...)
#' @param useAcceptedApproximateMinimizers (default: TRUE) \emph{TRUE or FALSE} If true then use chai-square and elbow method to choose maximum accepted SSR.  If false use parameters upto numParametersIncluded-th smallest SSR (or if numParametersIncluded=NA then use all the parameters found by the CGNM).
#' @return \emph{A ggplot object} including the violin plot, interquartile range and median, minimum and maximum.
#' @examples
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
#' plot_parameterValue_scatterPlots(CGNM_result)
#' @export
#' @import ggplot2
plot_parameterValue_scatterPlots=function(CGNM_result, cutoff_pvalue=0.05, numParametersIncluded=NA, ParameterNames=NA, useAcceptedApproximateMinimizers=TRUE){
  SSR=NULL
  X_value=NULL
  cluster=NULL
  Y_value=NULL
  cluster=NULL

  lastIter=dim(CGNM_result$residual_history)[2]

  plot_df=data.frame()

 # if(useAcceptedApproximateMinimizers){
    useIndecies_b=acceptedIndecies_binary(CGNM_result,cutoff_pvalue, numParametersIncluded, useAcceptedApproximateMinimizers)
  # }else{
  #   if(is.na(numParametersIncluded)|numParametersIncluded>dim(CGNM_result$X)[1]){
  #     useIndecies_b=rep(TRUE, dim(CGNM_result$X)[1])
  #   }else{
  #     SSR_vec=CGNM_result$residual_history[,dim(CGNM_result$residual_history)[2]]
  #     useIndecies_b=(SSR_vec<=sort(SSR_vec)[numParametersIncluded])
  #   }
  # }

  if(length(ParameterNames)==dim(CGNM_result$initialX)[2]){

  }else{
    ParameterNames=paste0("x_",seq(1,dim(CGNM_result$initialX)[2]))
  }


  kmeans_result=optimal_kmeans(cbind(CGNM_result$X[useIndecies_b,], plot_df$SSR))
  cluster_vec=as.factor(kmeans_result$cluster)

  for(i in seq(1,dim(CGNM_result$X)[2])){
    for(j in seq(1,dim(CGNM_result$X)[2])){
      plot_df=rbind(plot_df, data.frame(xName=ParameterNames[i],X_value=CGNM_result$X[useIndecies_b,i],yName=ParameterNames[j],Y_value=CGNM_result$X[useIndecies_b,j],cluster=cluster_vec))
    }
  }



  ggplot2::ggplot(plot_df, ggplot2::aes(y=Y_value, x=X_value, colour=cluster))+ggplot2::geom_point(alpha=0.3)+ggplot2::facet_wrap(xName~yName,scales="free")+
    ggplot2::theme(
      strip.background = element_blank(),
      strip.text.x = element_blank()
    )

}


#' @title bestApproximateMinimizers
#' @description
#' Returns the approximate minimizers with minimum SSR found by CGNM.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param numParameterSet (default 1) \emph{A natural number} number of parameter sets to output (chosen from the smallest SSR to numParameterSet-th smallest SSR) .
#' @return \emph{A vector} a vector of accepted approximate minimizers with minimum SSR found by CGNM.
#' @examples
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
#' bestApproximateMinimizers(CGNM_result,10)
#' @export
bestApproximateMinimizers=function(CGNM_result, numParameterSet=1){
  SSRorder=order(CGNM_result$residual_history[,dim(CGNM_result$residual_history)[2]])

  CGNM_result$X[SSRorder[1:numParameterSet],]
}
