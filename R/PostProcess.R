## TODO make autodiagnosis function, is Lambda all big?  Is the minimum SSR parameter set within the initial range?
## TODO implement best parameter
## TODO implement worst accepted parameter


makeParaDistributionPlotDataFrame=function(CGNM_result, indicesToInclude=NA, cutoff_pvalue=0.05, numParametersIncluded=NA, ParameterNames=NA, ReparameterizationDef=NA, useAcceptedApproximateMinimizers=TRUE){
  Kind_iter=NULL
  X_value=NULL
  cluster=NULL
  freeParaValues=data.frame()
  SSR_vec=CGNM_result$residual_history[,dim(CGNM_result$residual_history)[2]]


  if(is.na(ParameterNames[1])|length(ParameterNames)!=length(ReparameterizationDef)){
    ParameterNames=ReparameterizationDef
  }

  if(is.na(ParameterNames[1])){
    ParameterNames=paste0("x",seq(1,dim(CGNM_result$initialX)[2]))
    ReparameterizationDef=ParameterNames
  }

  if(!is.na(indicesToInclude[1])){
    useIndecies_b=seq(1,dim(CGNM_result$X)[1])%in%indicesToInclude
  }else if(useAcceptedApproximateMinimizers){
    useIndecies_b=acceptedIndices_binary(CGNM_result,cutoff_pvalue, numParametersIncluded, useAcceptedApproximateMinimizers)
  }else{
    if(is.na(numParametersIncluded)|numParametersIncluded>dim(CGNM_result$X)[1]){
      useIndecies_b=rep(TRUE, dim(CGNM_result$X)[1])
    }else{
      useIndecies_b=(SSR_vec<=sort(SSR_vec)[numParametersIncluded])
    }
  }

  initialX_df=CGNM_result$initialX
  colnames(initialX_df)=paste0("x",seq(1,dim(CGNM_result$initialX)[2]))
  initialX_df=data.frame(initialX_df)

  repara_initialX_df=data.frame(row.names = seq(1,dim(initialX_df)[1]))
  for(i in seq(1,length(ParameterNames))){
    repara_initialX_df[,ParameterNames[i]]=with(initialX_df, eval(parse(text=ReparameterizationDef[i])))
  }


  finalX_df=CGNM_result$X
  colnames(finalX_df)=paste0("x",seq(1,dim(CGNM_result$X)[2]))
  finalX_df=data.frame(finalX_df)

  repara_finalX_df=data.frame(row.names = seq(1,dim(finalX_df)[1]))
  for(i in seq(1,length(ParameterNames))){
    repara_finalX_df[,ParameterNames[i]]=with(finalX_df, eval(parse(text=ReparameterizationDef[i])))
  }


  for(i in seq(1,length(ParameterNames))){
    freeParaValues=rbind(freeParaValues, data.frame(Name=ParameterNames[i],X_value=repara_initialX_df[,i], Kind_iter="Initial", SSR=NA))
  }

  for(i in seq(1,length(ParameterNames))){
    freeParaValues=rbind(freeParaValues, data.frame(Name=ParameterNames[i],X_value=repara_finalX_df[useIndecies_b,i], Kind_iter="Final Accepted", SSR=SSR_vec[useIndecies_b]))
  }

  if(!is.null(CGNM_result$bootstrapX)){

    bootstrapX_df=CGNM_result$bootstrapX
    colnames(bootstrapX_df)=paste0("x",seq(1,dim(CGNM_result$bootstrapX)[2]))
    bootstrapX_df=data.frame(bootstrapX_df)

    repara_bootstrapX_df=data.frame(row.names = seq(1,dim(bootstrapX_df)[1]))
    for(i in seq(1,length(ParameterNames))){
      repara_bootstrapX_df[,ParameterNames[i]]=with(bootstrapX_df, eval(parse(text=ReparameterizationDef[i])))
    }

    for(i in seq(1,length(ParameterNames))){
      freeParaValues=rbind(freeParaValues, data.frame(Name=ParameterNames[i],X_value=repara_bootstrapX_df[,i], Kind_iter="Bootstrap",SSR=NA))
    }

    freeParaValues$Kind_iter=factor(freeParaValues$Kind_iter, levels = c("Initial", "Final Accepted", "Bootstrap"))

  }else{
    freeParaValues$Kind_iter=factor(freeParaValues$Kind_iter, levels = c("Initial", "Final Accepted"))

  }

  freeParaValues$Name=factor(freeParaValues$Name, levels = ParameterNames)

  return(freeParaValues)
}

#' @title acceptedMaxSSR
#' @description
#' CGNM find multiple sets of minimizers of the nonlinear least squares (nls) problem by solving nls from various initial iterates.  Although CGNM is shown to be robust compared to other conventional multi-start algorithms, not all initial iterates minimizes successfully.  By assuming sum of squares residual (SSR) follows the chai-square distribution we first reject the approximated minimiser who SSR is statistically significantly worse than the minimum SSR found by the CGNM.  Then use elbow-method (a heuristic often used in mathematical optimisation to balance the quality and the quantity of the solution found) to find the "acceptable" maximum SSR.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param cutoff_pvalue (default: 0.05) \emph{A number} defines the rejection p-value for the first stage of acceptable computational result screening.
#' @param numParametersIncluded (default: NA) \emph{A natural number} defines the number of parameter sets to be included in the assessment of the acceptable parameters.  If set NA then use all the parameters found by the CGNM.
#' @param useAcceptedApproximateMinimizers (default: TRUE) \emph{TRUE or FALSE} If true then use chai-square and elbow method to choose maximum accepted SSR.  If false returnsnumParametersIncluded-th smallest SSR (or if numParametersIncluded=NA then returns the largest SSR).
#' @param algorithm (default: 2) \emph{1 or 2} specify the algorithm used for obtain accepted approximate minimizers. (Algorithm 1 uses elbow method, Algorithm 2 uses Grubbs' Test for Outliers.)
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
acceptedMaxSSR=function(CGNM_result, cutoff_pvalue=0.05, numParametersIncluded=NA, useAcceptedApproximateMinimizers=TRUE, algorithm=2){

SSR_vec=CGNM_result$residual_history[,dim(CGNM_result$residual_history)[2]]

if(useAcceptedApproximateMinimizers){

  if(algorithm==2){

    numInInitialSet=dim(CGNM_result$residual_history)[1]

    orderedIndex=order(SSR_vec)

    bestIndex=topIndices(CGNM_result,1)

    ptest_vec=c()

    cumResidualFromBestFit=c()

    pvalue_vec=c()
    for(i in seq(2,numInInitialSet)){

      numInSample=i-1
      cumResidualFromBestFit=c(cumResidualFromBestFit, sum((CGNM_result$Y[orderedIndex[i-1],]-CGNM_result$Y[orderedIndex[i],])))

      if(numInSample>3){

        tStat=qt(1-cutoff_pvalue/(2*numInSample), df=(numInSample-2))

        #test statistics from Grubbs' Test for Outliers
        testStatistics=(numInSample-1)/sqrt(numInSample)*sqrt(tStat^2/(numInSample-2+tStat^2))

        #test hypothesis where the null hypothesis is that the last added sum of residual is outlier compared to the ones that are already accepted
        ptest_vec=c(ptest_vec,abs(mean(cumResidualFromBestFit)-cumResidualFromBestFit[length(cumResidualFromBestFit)])/sd(cumResidualFromBestFit)>testStatistics)

      }else{
        ptest_vec=c(ptest_vec,FALSE)
      }
    }

    acceptMaxSSR=sort(SSR_vec)[which(ptest_vec)[1]+1]
  }else{
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
  }


}else if(is.na(numParametersIncluded)|numParametersIncluded>length(SSR_vec)){
  acceptMaxSSR=max(SSR_vec, na.rm = TRUE)
}else{
  acceptMaxSSR=sort(SSR_vec)[numParametersIncluded]
}

acceptMaxSSR=max(acceptMaxSSR, min(SSR_vec)+sqrt(.Machine$double.eps))
  return(acceptMaxSSR)
}

#' @title acceptedApproximateMinimizers
#' @description
#' CGNM find multiple sets of minimizers of the nonlinear least squares (nls) problem by solving nls from various initial iterates.  Although CGNM is shown to be robust compared to other conventional multi-start algorithms, not all initial iterates minimizes successfully.  By assuming sum of squares residual (SSR) follows the chai-square distribution we first reject the approximated minimiser who SSR is statistically significantly worse than the minimum SSR found by the CGNM.  Then use elbow-method (a heuristic often used in mathematical optimisation to balance the quality and the quantity of the solution found) to find the "acceptable" maximum SSR. This function outputs the acceptable approximate minimizers of the nonlinear least squares problem found by the CGNM.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param cutoff_pvalue (default: 0.05) \emph{A number} defines the rejection p-value for the first stage of acceptable computational result screening.
#' @param numParametersIncluded (default: NA) \emph{A natural number} defines the number of parameter sets to be included in the assessment of the acceptable parameters.  If set NA then use all the parameters found by the CGNM.
#' @param useAcceptedApproximateMinimizers (default: TRUE) \emph{TRUE or FALSE} If true then use chai-square and elbow method to choose maximum accepted SSR.  If false returns the parameters upto numParametersIncluded-th smallest SSR (or if numParametersIncluded=NA then use all the parameters found by the CGNM).
#' @param algorithm (default: 2) \emph{1 or 2} specify the algorithm used for obtain accepted approximate minimizers. (Algorithm 1 uses elbow method, Algorithm 2 uses Grubbs' Test for Outliers.)
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
acceptedApproximateMinimizers=function(CGNM_result, cutoff_pvalue=0.05, numParametersIncluded=NA, useAcceptedApproximateMinimizers=TRUE, algorithm=2){
  CGNM_result$X[acceptedIndices(CGNM_result, cutoff_pvalue, numParametersIncluded, useAcceptedApproximateMinimizers, algorithm=algorithm),]
}

#' @title acceptedIndices
#' @description
#' CGNM find multiple sets of minimizers of the nonlinear least squares (nls) problem by solving nls from various initial iterates.  Although CGNM is shown to be robust compared to other conventional multi-start algorithms, not all initial iterates minimizes successfully.  By assuming sum of squares residual (SSR) follows the chai-square distribution we first reject the approximated minimiser who SSR is statistically significantly worse than the minimum SSR found by the CGNM.  Then use elbow-method (a heuristic often used in mathematical optimisation to balance the quality and the quantity of the solution found) to find the "acceptable" maximum SSR. This function outputs the indices of acceptable approximate minimizers of the nonlinear least squares problem found by the CGNM.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param cutoff_pvalue (default: 0.05) \emph{A number} defines the rejection p-value for the first stage of acceptable computational result screening.
#' @param numParametersIncluded (default: NA) \emph{A natural number} defines the number of parameter sets to be included in the assessment of the acceptable parameters.  If set NA then use all the parameters found by the CGNM.
#' @param useAcceptedApproximateMinimizers (default: TRUE) \emph{TRUE or FALSE} If true then use chai-square and elbow method to choose maximum accepted SSR.  If false returns the parameters upto numParametersIncluded-th smallest SSR (or if numParametersIncluded=NA then use all the parameters found by the CGNM).
#' @param algorithm (default: 2) \emph{1 or 2} specify the algorithm used for obtain accepted approximate minimizers. (Algorithm 1 uses elbow method, Algorithm 2 uses Grubbs' Test for Outliers.)
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
#' acceptedIndices(CGNM_result)
#' @export
acceptedIndices=function(CGNM_result, cutoff_pvalue=0.05, numParametersIncluded=NA, useAcceptedApproximateMinimizers=TRUE, algorithm=2){
  acceptMaxSSR_value=acceptedMaxSSR(CGNM_result, cutoff_pvalue, numParametersIncluded, useAcceptedApproximateMinimizers, algorithm=algorithm)

  SSR_vec=CGNM_result$residual_history[,dim(CGNM_result$residual_history)[2]]

  return(which(SSR_vec<=acceptMaxSSR_value))
}


#' @title topIndices
#' @description
#' CGNM find multiple sets of minimizers of the nonlinear least squares (nls) problem by solving nls from various initial iterates.  Although CGNM is shown to be robust compared to other conventional multi-start algorithms, not all initial iterates minimizes successfully.  One can visually inspect rank v.s. SSR plot and manually choose number of best fit acceptable parameters.  By using this function "topIndices", we can obtain the indices of the "numTopIndices" best fit parameter combinations.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param numTopIndices (required input) \emph{An integer} .
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
#' topInd=topIndices(CGNM_result, 10)
#'
#' ## This gives top 10 approximate minimizers
#' CGNM_result$X[topInd,]
#' @export
topIndices=function(CGNM_result, numTopIndices){

  SSR_vec=CGNM_result$residual_history[,dim(CGNM_result$residual_history)[2]]
  sortedSSR=sort(SSR_vec)
  outIndices=which(SSR_vec<=sortedSSR[numTopIndices])

  return(outIndices)
}



#' @title acceptedIndices_binary
#' @description
#' CGNM find multiple sets of minimizers of the nonlinear least squares (nls) problem by solving nls from various initial iterates.  Although CGNM is shown to be robust compared to other conventional multi-start algorithms, not all initial iterates minimizes successfully.  By assuming sum of squares residual (SSR) follows the chai-square distribution we first reject the approximated minimiser who SSR is statistically significantly worse than the minimum SSR found by the CGNM.  Then use elbow-method (a heuristic often used in mathematical optimisation to balance the quality and the quantity of the solution found) to find the "acceptable" maximum SSR. This function outputs the indices of acceptable approximate minimizers of the nonlinear least squares problem found by the CGNM. (note that acceptedIndices(CGNM_result) is equal to seq(1,length(acceptedIndices_binary(CGNM_result)))[acceptedIndices_binary(CGNM_result)])
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param cutoff_pvalue (default: 0.05) \emph{A number} defines the rejection p-value for the first stage of acceptable computational result screening.
#' @param numParametersIncluded (default: NA) \emph{A natural number} defines the number of parameter sets to be included in the assessment of the acceptable parameters.  If set NA then use all the parameters found by the CGNM.
#' @param useAcceptedApproximateMinimizers (default: TRUE) \emph{TRUE or FALSE} If true then use chai-square and elbow method to choose maximum accepted SSR.  If false returns the indicies upto numParametersIncluded-th smallest SSR (or if numParametersIncluded=NA then use all the parameters found by the CGNM).
#' @param algorithm (default: 2) \emph{1 or 2} specify the algorithm used for obtain accepted approximate minimizers. (Algorithm 1 uses elbow method, Algorithm 2 uses Grubbs' Test for Outliers.)
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
#' acceptedIndices_binary(CGNM_result)
#' @export
acceptedIndices_binary=function(CGNM_result, cutoff_pvalue=0.05, numParametersIncluded=NA, useAcceptedApproximateMinimizers=TRUE, algorithm=2){
  acceptMaxSSR_value=acceptedMaxSSR(CGNM_result, cutoff_pvalue, numParametersIncluded, useAcceptedApproximateMinimizers, algorithm= algorithm)

  SSR_vec=CGNM_result$residual_history[,dim(CGNM_result$residual_history)[2]]

  return(as.vector(SSR_vec<=acceptMaxSSR_value))
}


#' @title plot_Rank_SSR
#' @description
#' Make SSR v.s. rank plot. This plot is often used to visualize the maximum accepted SSR.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param indicesToInclude (default: NA) \emph{A vector of integers} indices to include in the plot (if NA, use indices chosen by the acceptedIndices() function with default setting).
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
plot_Rank_SSR=function(CGNM_result, indicesToInclude=NA){
  SSR_value=NULL
  cutoff_pvalue=0.05
  numParametersIncluded=NA
  useAcceptedApproximateMinimizers=TRUE

  is_accepted=NULL
  SSR_vec=CGNM_result$residual_history[,dim(CGNM_result$residual_history)[2]]

  if(!is.na(indicesToInclude[1])){
    acceptedIndices_b=seq(1,dim(CGNM_result$X)[1])%in%indicesToInclude
  }else{
    acceptedIndices_b=acceptedIndices_binary(CGNM_result,cutoff_pvalue, numParametersIncluded, useAcceptedApproximateMinimizers)
  }

  numAccept=sum(acceptedIndices_b)
  acceptMaxSSR=acceptedMaxSSR(CGNM_result,cutoff_pvalue, numParametersIncluded, useAcceptedApproximateMinimizers)
  min_R=min(SSR_vec)

  plot_df=data.frame(SSR_value=SSR_vec, rank=rank(SSR_vec), is_accepted=acceptedIndices_b)

  ggplot2::ggplot(plot_df, ggplot2::aes(x=rank,y=SSR_value, colour=is_accepted))+ggplot2::geom_point()+ggplot2::coord_cartesian(ylim=c(0,acceptMaxSSR*2))+ggplot2::geom_vline(xintercept = numAccept, color="grey")+ ggplot2::annotate(geom="text", x=numAccept, y=acceptMaxSSR*0.5, label=paste("Accepted: ",numAccept,"\n Accepted max SSR: ",formatC(acceptMaxSSR, format = "g", digits = 3)),angle = 90,
                                                                                                                                                                       color="black")+ ggplot2::annotate(geom="text", x=length(SSR_vec)*0.1, y=min_R*1.1, label=paste("min SSR: ",formatC(min_R, format = "g", digits = 3)),
                                                                                                                                                                                                color="black")+ylab("SSR")



}


#' @title plot_paraDistribution_byViolinPlots
#' @description
#' Make violin plot to compare the initial distribution and distribition of the accepted approximate minimizers found by the CGNM. Bars in the violin plots indicates the interquartile range. The solid line connects the interquartile ranges of the initial distribution and the distribution of the accepted approximate minimizer at the final iterate.  The blacklines connets the minimums and maximums of the initial distribution and the distribution of the accepted approximate minimizer at the final iterate. The black dots indicate the median.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param indicesToInclude (default: NA) \emph{A vector of integers} indices to include in the plot (if NA, use indices chosen by the acceptedIndices() function with default setting).
#' @param ParameterNames (default: NA) \emph{A vector of strings} the user can supply so that these names are used when making the plot. (Note if it set as NA or vector of incorrect length then the parameters are named as x1, x2, ...)
#' @param ReparameterizationDef (default: NA) \emph{A vector of strings} the user can supply definition of reparameterization where each string follows R syntax and also refer to the ith element of the x vector (the input variable to the nonlinear function) as xi (e.g., if the first input variable to the nonlinear function is defined as x1=log10(Ka), then by setting "10^x1" as one of the strings in this vector, you can plot the violin plot of Ka)
#' @return \emph{A ggplot object} including the violin plot, interquartile range and median, minimum and maximum.
#' @examples
#'
#'model_analytic_function=function(x){
#'
#'  observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
#'  Dose=1000
#'  F=1
#'
#'  ka=10^x[1]
#'  V1=10^x[2]
#'  CL_2=10^x[3]
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
#' initial_lowerRange = c(-2,-2,-2), initial_upperRange =  c(1,2,2),
#' num_iter = 10, num_minimizersToFind = 100)
#'
#' a_Indices=acceptedIndices(CGNM_result)
#' plot_paraDistribution_byViolinPlots(CGNM_result, indicesToInclude=a_Indices)
#' plot_paraDistribution_byViolinPlots(CGNM_result, indicesToInclude=a_Indices,
#'      ParameterNames=c("Ka","V1","CL_2"),
#'      ReparameterizationDef=c("10^x1","10^x2","10^x3"))
#' @export
#' @import ggplot2

plot_paraDistribution_byViolinPlots=function(CGNM_result, indicesToInclude=NA, ParameterNames=NA, ReparameterizationDef=NA){

  Kind_iter=NULL
  X_value=NULL

  cutoff_pvalue=0.05
  numParametersIncluded=NA
  useAcceptedApproximateMinimizers=TRUE

  freeParaValues=makeParaDistributionPlotDataFrame(CGNM_result, indicesToInclude, cutoff_pvalue, numParametersIncluded, ParameterNames, ReparameterizationDef, useAcceptedApproximateMinimizers)

  p<-ggplot2::ggplot(freeParaValues,ggplot2::aes(x=Kind_iter,y=X_value))+ggplot2::facet_wrap(Name~., scales = "free")

  p+ggplot2::geom_violin(trim=T,fill="#999999",linetype="blank",alpha=I(1/2))+
    ggplot2::stat_summary(geom="pointrange",fun = median, fun.min = function(x) quantile(x,probs=0.25), fun.max = function(x) quantile(x,probs=0.75), size=0.5,alpha=.5)+
    ggplot2::stat_summary(geom="line",fun =  function(x) quantile(x,probs=0), ggplot2::aes(group=1),size=0.5,alpha=.3,linetype=2)+
    ggplot2::stat_summary(geom="line",fun =  function(x) quantile(x,probs=1), ggplot2::aes(group=1),size=0.5,alpha=.3,linetype=2)+
    ggplot2::stat_summary(geom="line",fun =  function(x) quantile(x,probs=0.25), ggplot2::aes(group=1),size=0.5,alpha=.3)+
    ggplot2::stat_summary(geom="line",fun =  function(x) quantile(x,probs=0.75), ggplot2::aes(group=1),size=0.5,alpha=.3)+
    ggplot2::theme(legend.position="none",axis.text.x = element_text(angle = 25, hjust = 1))+xlab("")+ylab("Value")


}


#' @title plot_paraDistribution_byHistogram
#' @description
#' Make histograms to visualize the initial distribution and distribition of the accepted approximate minimizers found by the CGNM.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param indicesToInclude (default: NA) \emph{A vector of integers} indices to include in the plot (if NA, use indices chosen by the acceptedIndices() function with default setting).
#' @param ParameterNames (default: NA) \emph{A vector of strings} the user can supply so that these names are used when making the plot. (Note if it set as NA or vector of incorrect length then the parameters are named as x1, x2, ...)
#' @param ReparameterizationDef (default: NA) \emph{A vector of strings} the user can supply definition of reparameterization where each string follows R syntax and also refer to the ith element of the x vector (the input variable to the nonlinear function) as xi (e.g., if the first input variable to the nonlinear function is defined as x1=log10(Ka), then by setting "10^x1" as one of the strings in this vector, you can plot the violin plot of Ka)
#' @param bins (default: 30) \emph{A natural number} Number of bins used for plotting histogram.
#' @return \emph{A ggplot object} including the violin plot, interquartile range and median, minimum and maximum.
#' @examples
#'
#'model_analytic_function=function(x){
#'
#'  observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
#'  Dose=1000
#'  F=1
#'
#'  ka=10^x[1]
#'  V1=10^x[2]
#'  CL_2=10^x[3]
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
#' initial_lowerRange = c(-2,-2,-2), initial_upperRange =  c(1,2,2),
#' num_iter = 10, num_minimizersToFind = 100)
#'
#' plot_paraDistribution_byHistogram(CGNM_result)
#' plot_paraDistribution_byHistogram(CGNM_result,
#'      ParameterNames=c("Ka","V1","CL_2"),
#'      ReparameterizationDef=c("10^x1","10^x2","10^x3"))
#' @export
#' @import ggplot2

plot_paraDistribution_byHistogram=function(CGNM_result, indicesToInclude=NA, ParameterNames=NA, ReparameterizationDef=NA,  bins=30){

  X_value=NULL
  Kind_iter=NULL

  cutoff_pvalue=0.05
  numParametersIncluded=NA
  useAcceptedApproximateMinimizers=TRUE
  freeParaValues=makeParaDistributionPlotDataFrame(CGNM_result, indicesToInclude, cutoff_pvalue, numParametersIncluded, ParameterNames, ReparameterizationDef, useAcceptedApproximateMinimizers)

  p<-ggplot2::ggplot(freeParaValues,ggplot2::aes(X_value))+ggplot2::geom_histogram(bins = bins)+ggplot2::facet_grid(Kind_iter~Name, scales = "free")+xlab("")

  p
}



#' @title plot_profileLikelihood
#' @description
#' Draw profile likelihood surface using the function evaluations conducted during CGNM computation. Note plot_SSRsurface can only be used when log is saved by setting saveLog=TRUE option when running Cluster_Gauss_Newton_method().  The grey horizontal line is the threshold for 95% pointwise confidence interval.
#' @param logLocation (required input) \emph{A string} of folder directory where CGNM computation log files exist.
#' @param numBins (default: NA) \emph{A positive integer} SSR surface is plotted by finding the minimum SSR given one of the parameters is fixed and then repeat this for various values.  numBins specifies the number of different parameter values to fix for each parameter. (if set NA the number of bins are set as num_minimizersToFind/10)
#' @param ParameterNames (default: NA) \emph{A vector of strings} the user can supply so that these names are used when making the plot. (Note if it set as NA or vector of incorrect length then the parameters are named as x1, x2, ...)
#' @param ReparameterizationDef (default: NA) \emph{A vector of strings} the user can supply definition of reparameterization where each string follows R syntax and also refer to the ith element of the x vector (the input variable to the nonlinear function) as xi (e.g., if the first input variable to the nonlinear function is defined as x1=log10(Ka), then by setting "10^x1" as one of the strings in this vector, you can plot the violin plot of Ka)
#' @param showInitialRange (default: FALSE) \emph{TRUE or FALSE} if TRUE then the initial range appears in the plot.
#' @return \emph{A ggplot object} including the violin plot, interquartile range and median, minimum and maximum.
#' @examples
#'\dontrun{
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
#' num_iter = 10, num_minimizersToFind = 100, saveLog=TRUE)
#'
#' plot_profileLikelihood("CGNM_log")
#' }
#' @export
#' @import ggplot2
plot_profileLikelihood=function(logLocation, numBins=NA,  ParameterNames=NA, ReparameterizationDef=NA, showInitialRange=FALSE){
  plot_SSRsurface(logLocation, profile_likelihood=TRUE, numBins=numBins, ParameterNames=ParameterNames, ReparameterizationDef=ReparameterizationDef, showInitialRange=showInitialRange)
  }

#' @title plot_SSRsurface
#' @description
#' Make minimum SSR v.s. parameterValue plot using the function evaluations used during CGNM computation. Note plot_SSRsurface can only be used when log is saved by setting saveLog=TRUE option when running Cluster_Gauss_Newton_method().
#' @param logLocation (required input) \emph{A string or a list of strings} of folder directory where CGNM computation log files exist.
#' @param profile_likelihood (default: FALSE) \emph{TRUE or FALSE} If set TRUE plot profile likelihood (assuming normal distribution of residual) instead of SSR surface.
#' @param numBins (default: NA) \emph{A positive integer} SSR surface is plotted by finding the minimum SSR given one of the parameters is fixed and then repeat this for various values.  numBins specifies the number of different parameter values to fix for each parameter. (if set NA the number of bins are set as num_minimizersToFind/10)
#' @param maxSSR (default: NA) \emph{A positive number} the maximum SSR that will be plotted on SSR surface plot.  This option is used to zoom into the SSR surface near the minimum SSR.
#' @param ParameterNames (default: NA) \emph{A vector of strings} the user can supply so that these names are used when making the plot. (Note if it set as NA or vector of incorrect length then the parameters are named as x1, x2, ...)
#' @param ReparameterizationDef (default: NA) \emph{A vector of strings} the user can supply definition of reparameterization where each string follows R syntax and also refer to the ith element of the x vector (the input variable to the nonlinear function) as xi (e.g., if the first input variable to the nonlinear function is defined as x1=log10(Ka), then by setting "10^x1" as one of the strings in this vector, you can plot the violin plot of Ka)
#' @param showInitialRange (default: FALSE) \emph{TRUE or FALSE} if TRUE then the initial range appears in the plot.
#' @return \emph{A ggplot object} including the violin plot, interquartile range and median, minimum and maximum.
#' @examples
#'\dontrun{
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
#' num_iter = 10, num_minimizersToFind = 100, saveLog=TRUE)
#'
#' plot_SSRsurface("CGNM_log") + scale_y_continuous(trans='log10')
#' }
#' @export
#' @import ggplot2
plot_SSRsurface=function(logLocation, profile_likelihood=FALSE, numBins=NA, maxSSR=NA, ParameterNames=NA, ReparameterizationDef=NA, showInitialRange=FALSE){

  CGNM_result=NULL
  minSSR=NULL
  negative2LogLikelihood=NULL
  value=NULL
  boundValue=NULL


  if(typeof(logLocation)=="character"){
    DirectoryName_vec=c(logLocation)

  }else if(typeof(logLocation)=="list"){
      if(logLocation$runSetting$runName==""){
        DirectoryName_vec=c("CGNM_log","CGNM_log_bootstrap")
      }else{
        DirectoryName_vec=c(paste0("CGNM_log_",logLocation$runSetting$runName),paste0("CGNM_log_",logLocation$runSetting$runName,"bootstrap"))
      }



  }else{
    warning("DirectoryName need to be either the CGNM_result object or a string")
  }


  load(paste0(DirectoryName_vec[1],"/iteration_1.RDATA"))


  if(is.na(ParameterNames[1])|length(ParameterNames)!=length(ReparameterizationDef)){
    ParameterNames=ReparameterizationDef
  }

  if(is.na(ParameterNames[1])){
    ParameterNames=paste0("x",seq(1,dim(CGNM_result$initialX)[2]))
    ReparameterizationDef=ParameterNames
  }




  if(is.na(numBins)){
    numBins=round(CGNM_result$runSetting$num_minimizersToFind/10)
  }


  rawX_nu=CGNM_result$initialX
  R_nu=CGNM_result$residual_history[,1]

  initiLNLfunc=CGNM_result$runSetting$nonlinearFunction

  residual_variance_vec=c()

  for(DirectoryName in DirectoryName_vec){

    checkNL=TRUE
    for(i in seq(1,CGNM_result$runSetting$num_iteration)){
      fileName=paste0(DirectoryName,"/iteration_",i,".RDATA")
      if (file.exists(fileName)){
        load(fileName)
        if(identical(initiLNLfunc, CGNM_result$runSetting$nonlinearFunction)){
          rawX_nu=rbind(rawX_nu,CGNM_result$X)
          R_nu=c(R_nu, CGNM_result$residual_history[,dim(CGNM_result$residual_history)[2]])
          if(checkNL){
            print(paste0("log saved in ",getwd(),"/",DirectoryName," is used to draw SSR/likelihood surface"))
            checkNL=FALSE
          }
        }else if(checkNL){
          print(paste0("log saved in ",getwd(),"/",DirectoryName," will NOT be used as the nonlinear function used in this log is not the same as ",getwd(),"/",DirectoryName[1]))
          checkNL=FALSE
        }
      }
    }

    residual_variance_vec=c(residual_variance_vec, min(CGNM_result$residual_history,na.rm = TRUE)/sum(!is.na(CGNM_result$runSetting$targetVector)))

  }

  residual_variance=min(residual_variance_vec)

  rawXinit=CGNM_result$initialX
  reparaXinit=data.frame(row.names = seq(1,dim(rawXinit)[1]))

  colnames(rawXinit)=paste0("x",seq(1,dim(rawXinit)[2]))
  rawXinit=data.frame(rawXinit)

  for(i in seq(1,length(ParameterNames))){
    reparaXinit[,ParameterNames[i]]=with(rawXinit, eval(parse(text=ReparameterizationDef[i])))
  }


#
#   for(i in seq(1,CGNM_result$runSetting$num_iteration)){
#     fileName=paste0(DirectoryName,"bootstrap/iteration_",i,".RDATA")
#
#
#     if (file.exists(fileName)){
#       load(fileName)
#       rawX_nu=rbind(rawX_nu,CGNM_result$X)
#
#       R_temp=c()
#       for(j in seq(1, dim(CGNM_result$Y)[1])){
#         R_temp=c(R_temp,sum((CGNM_result$Y[j,]-CGNM_result$runSetting$targetVector)^2,na.rm = TRUE))
#       }
#
#       R_nu=c(R_nu, R_temp)
#     }
#   }



  tempMatrix=unique(cbind(rawX_nu,R_nu))
  rawX_nu=tempMatrix[,seq(1,dim(rawX_nu)[2])]
  R_nu=tempMatrix[,dim(tempMatrix)[2]]

  colnames(rawX_nu)=paste0("x",seq(1,dim(rawX_nu)[2]))
  rawX_nu=data.frame(rawX_nu)

  X_nu=data.frame(row.names = seq(1,dim(rawX_nu)[1]))
  for(i in seq(1,length(ParameterNames))){
    X_nu[,ParameterNames[i]]=with(rawX_nu, eval(parse(text=ReparameterizationDef[i])))
  }


  likelihoodSurfacePlot_df=data.frame()

  for(j in seq(1,length(ParameterNames))){
    for(i in seq(1,numBins)){
      binUpper=quantile(X_nu[,j], probs = c((i-1)*(1/numBins),(i)*(1/numBins)))[2]
      binLower=quantile(X_nu[,j], probs = c((i-1)*(1/numBins),(i)*(1/numBins)))[1]
      indexToUse=((X_nu[,j]<=binUpper)&(X_nu[,j]>binLower))
      likelihoodSurfacePlot_df=rbind(likelihoodSurfacePlot_df, data.frame(value=X_nu[which(R_nu==min(R_nu[indexToUse],na.rm = TRUE)),j], minSSR=min(R_nu[indexToUse]), parameterName=ParameterNames[j]))

    }
  }

#  ggplot(likelihoodSurfacePlot_df, aes(x=value,y=minSSR))+geom_point()+facet_wrap(.~parameterName,scales = "free")+geom_smooth()+ylim(0.5,1)

  if(!is.na(maxSSR)){
    likelihoodSurfacePlot_df=subset(likelihoodSurfacePlot_df, minSSR<=maxSSR)
  }

  likelihoodSurfacePlot_df$negative2LogLikelihood=likelihoodSurfacePlot_df$minSSR/residual_variance


  if(profile_likelihood){
    likelihoodSurfacePlot_df=subset(likelihoodSurfacePlot_df, negative2LogLikelihood<=(2*(3.84+min(likelihoodSurfacePlot_df$negative2LogLikelihood))))

    g=ggplot2::ggplot(likelihoodSurfacePlot_df, ggplot2::aes(x=value,y=negative2LogLikelihood))+ggplot2::geom_point()+ggplot2::geom_line()

  }else{
    g=ggplot2::ggplot(likelihoodSurfacePlot_df, ggplot2::aes(x=value,y=minSSR))+ggplot2::geom_point()+ggplot2::geom_line()

  }

#  g=ggplot2::ggplot(likelihoodSurfacePlot_df, ggplot2::aes(x=value,y=minSSR))+ggplot2::geom_point()+ggplot2::geom_line()



  if(showInitialRange){
    Xinit_min=c()
    Xinit_max=c()
    for(i in seq(1,dim(reparaXinit)[2])){
      Xinit_min=c(Xinit_min,min(reparaXinit[,i]))
      Xinit_max=c(Xinit_max,max(reparaXinit[,i]))
    }

    g=g+ggplot2::geom_vline(data=data.frame(boundValue=Xinit_min, parameterName=ParameterNames), ggplot2::aes(xintercept=boundValue), colour="gray")
    g=g+ggplot2::geom_vline(data=data.frame(boundValue=Xinit_max, parameterName=ParameterNames), ggplot2::aes(xintercept=boundValue), colour="gray")
  }

  if(profile_likelihood){
    g=g+ggplot2::facet_wrap(.~parameterName,scales = "free_x")+ggplot2::ylab("-2log likelihood")+ggplot2::xlab("Parameter Value")
    g=g+ggplot2::geom_hline(yintercept = 3.84+min(likelihoodSurfacePlot_df$negative2LogLikelihood), colour="grey")

  }else{
    g=g+ggplot2::facet_wrap(.~parameterName,scales = "free_x")+ggplot2::ylab("SSR")+ggplot2::xlab("Parameter Value")

  }

  return(g)
}

#' @title plot_SSR_parameterValue
#' @description
#' Make SSR v.s. parameterValue plot of the accepted approximate minimizers found by the CGNM. Bars in the violin plots indicates the interquartile range.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param indicesToInclude (default: NA) \emph{A vector of integers} indices to include in the plot (if NA, use indices chosen by the acceptedIndices() function with default setting).
#' @param ParameterNames (default: NA) \emph{A vector of strings} the user can supply so that these names are used when making the plot. (Note if it set as NA or vector of incorrect length then the parameters are named as x1, x2, ...)
#' @param ReparameterizationDef (default: NA) \emph{A vector of strings} the user can supply definition of reparameterization where each string follows R syntax and also refer to the ith element of the x vector (the input variable to the nonlinear function) as xi (e.g., if the first input variable to the nonlinear function is defined as x1=log10(Ka), then by setting "10^x1" as one of the strings in this vector, you can plot the violin plot of Ka)
#' @param showInitialRange (default: TRUE) \emph{TRUE or FALSE} if TRUE then the initial range appears in the plot.
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
plot_SSR_parameterValue=function(CGNM_result, indicesToInclude=NA, ParameterNames=NA, ReparameterizationDef=NA, showInitialRange=TRUE){

  Kind_iter=NULL
  SSR=NULL
  X_value=NULL
  cluster=NULL
  boundValue=NULL

  cutoff_pvalue=0.05
  numParametersIncluded=NA
  useAcceptedApproximateMinimizers=TRUE

  freeParaValues=makeParaDistributionPlotDataFrame(CGNM_result, indicesToInclude, cutoff_pvalue, numParametersIncluded, ParameterNames, ReparameterizationDef, useAcceptedApproximateMinimizers)

  plotdata_df=subset(freeParaValues, Kind_iter=="Final Accepted")

  kmeans_result=optimal_kmeans(matrix(plotdata_df$X_value, ncol=dim(CGNM_result$X)[2]))

  plotdata_df$cluster=as.factor(kmeans_result$cluster)

  g=ggplot2::ggplot(plotdata_df, ggplot2::aes(y=SSR, x=X_value, colour=cluster))+ggplot2::geom_point(alpha=0.3)

  if(showInitialRange&&is.na(ParameterNames)){
    ParameterNames=paste0("x",seq(1,dim(CGNM_result$X)[2]))
    g=g+ggplot2::geom_vline(data=data.frame(boundValue=CGNM_result$runSetting$initial_lowerRange, Name=ParameterNames), ggplot2::aes(xintercept=boundValue), colour="gray")
    g=g+ggplot2::geom_vline(data=data.frame(boundValue=CGNM_result$runSetting$initial_upperRange, Name=ParameterNames), ggplot2::aes(xintercept=boundValue), colour="gray")
  }

  g+ggplot2::facet_wrap(.~Name,scales="free")
}



#' @title plot_parameterValue_scatterPlots
#' @description
#' Make scatter plots of the accepted approximate minimizers found by the CGNM. Bars in the violin plots indicates the interquartile range.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param indicesToInclude (default: NA) \emph{A vector of integers} indices to include in the plot (if NA, use indices chosen by the acceptedIndices() function with default setting).
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
plot_parameterValue_scatterPlots=function(CGNM_result, indicesToInclude=NA){

  cutoff_pvalue=0.05
  numParametersIncluded=NA
  ParameterNames=NA
  useAcceptedApproximateMinimizers=TRUE

  SSR=NULL
  X_value=NULL
  cluster=NULL
  Y_value=NULL
  cluster=NULL

  lastIter=dim(CGNM_result$residual_history)[2]

  plot_df=data.frame()

  if(!is.na(indicesToInclude[1])){
    useIndecies_b=seq(1,dim(CGNM_result$X)[1])%in%indicesToInclude
  }else{
    useIndecies_b=acceptedIndices_binary(CGNM_result,cutoff_pvalue, numParametersIncluded, useAcceptedApproximateMinimizers)
  }

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


#' @title plot_goodnessOfFit
#' @description
#' Make goodness of fit plots to assess the model-fit and bias in residual distribution.\cr\cr
#' Explanation of the terminologies in terms of PBPK model fitting to the time-course drug concentration measurements:
#' \enumerate{\item "independent variable" is time
#' \item"dependent variable" is the concentration.
#' \item "Residual" is the difference between the measured concentration and the model simulation with the parameter fond by the CGNM.
#' \item "m" is number of observations}
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param plotType (default: 1) \emph{1,2 or 3}\cr specify the kind of goodness of fit plot to create\enumerate{\item dependent variable v.s. independent variable with overlay of the target as red dots (e.g.,plots of time-course concentration profile with overlay of observed concentration in red dots). When CGNM_result include bootstrap analysis result, then the model simulation with median, 5 percentile and 95 percentile will be plotted.
#'  \item residual v.s. dependent variable (used to check to make sure one has chosen the right "shape" of residual distribution, i.e., additive, proportional etc., check to make sure there is no noticeable trend.)
#'  \item residual v.s. independent variable (e.g., use to check if the model-fit is equally good throughout different phases of time-course profile.)}
#' @param plotRank (default: c(1)) \emph{an integer of a vector of integers}\cr Specify which rank of the parameter to use for the goodness of fit plots. (e.g., if one wishes to use rank 1 to 100 then set it to be seq(1,100), or if one wish to use 88th rank parameters then set this as 88.)
#' @param independentVariableVector (default: NA) \emph{a vector of numerics of length m} \cr set independent variables that target values are associated with (e.g., time of the drug concentration measurement one is fitting PBPK model to) \cr(when this variable is set to NA, seq(1,m) will be used as independent variable when appropriate).
#' @param dependentVariableTypeVector (default: NA) \emph{a vector of text of length m} \cr when this variable is set (i.e., not NA) then the goodness of fit analyses is done for each variable type.  For example, if we are fitting the PBPK model to data with multiple dose arms, one can see the goodness of fit for each dose arm by specifying which dose group the observations are from.
#' @return \emph{A ggplot object} of the goodness of fit plot.
#' @examples
#'
#'model_analytic_function=function(x){
#'
#'  observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
#'  Dose=1000
#'  F=1
#'
#'  ka=10^x[1]
#'  V1=10^x[2]
#'  CL_2=10^x[3]
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
#' initial_lowerRange = c(-2,-2,-2), initial_upperRange =  c(1,2,2),
#' num_iter = 10, num_minimizersToFind = 100)
#'
#' plot_paraDistribution_byHistogram(CGNM_result)
#' plot_paraDistribution_byHistogram(CGNM_result,
#'      ParameterNames=c("Ka","V1","CL_2"),
#'      ReparameterizationDef=c("10^x1","10^x2","10^x3"))
#' @export
#' @import ggplot2

plot_goodnessOfFit=function(CGNM_result, plotType=1, plotRank=c(1), independentVariableVector=NA, dependentVariableTypeVector=NA){

  independent_variable=NULL
  lower=NULL
  upper=NULL
  target=NULL
  model_fit=NULL
  ind=NULL
  residual=NULL
  SSR=NULL


  independentVariableVector_in=independentVariableVector

  SSR_vec=CGNM_result$residual_history[,dim(CGNM_result$residual_history)[2]]
  indexToInclude=seq(1,length(SSR_vec))[rank(SSR_vec,na.last = TRUE)%in%plotRank]

  if(is.na(independentVariableVector[1])||length(independentVariableVector)!=length(CGNM_result$runSetting$targetVector)){
    independentVariableVector=seq(1, length(CGNM_result$runSetting$targetVector))
  }else{
    independentVariableVector=as.numeric(independentVariableVector)
  }

  if(is.na(dependentVariableTypeVector[1])||length(dependentVariableTypeVector)!=length(CGNM_result$runSetting$targetVector)){
    dependentVariableTypeVector=NA
  }else{
    dependentVariableTypeVector=as.factor(dependentVariableTypeVector)
  }

  residualPlot_df=data.frame()

  for(index in indexToInclude){
    residualPlot_df=rbind(residualPlot_df, data.frame(independent_variable=independentVariableVector,
                                                      residual=CGNM_result$Y[index,]-CGNM_result$runSetting$targetVector,
                                                      target=CGNM_result$runSetting$targetVector,
                                                      model_fit=CGNM_result$Y[index,],
                                                      dependent_variable_type=dependentVariableTypeVector,
                                                      square_residual=(CGNM_result$Y[index,]-CGNM_result$runSetting$targetVector)^2,
                                                      ind=index))
  }
  residualPlot_df$ind=as.factor(residualPlot_df$ind)


  withBootstrap=!is.null(CGNM_result$bootstrapY)

  median_vec=c()
  percentile5_vec=c()
  percentile95_vec=c()
  uncertaintyBound_df=data.frame()

  if(withBootstrap){
    for(i in seq(1, dim(CGNM_result$bootstrapY)[2])){
      tempQ=quantile(CGNM_result$bootstrapY[,i], prob=c(0.05,0.5,0.95), na.rm=TRUE)
      median_vec=c(median_vec, tempQ[2])
      percentile5_vec=c(percentile5_vec, tempQ[1])
      percentile95_vec=c(percentile95_vec, tempQ[3])
    }
    uncertaintyBound_df=data.frame(independent_variable=independentVariableVector, dependent_variable_type=dependentVariableTypeVector, median=median_vec, lower=percentile5_vec, upper=percentile95_vec, target=CGNM_result$runSetting$targetVector)
  }


  if(plotType==1){ #dependent variable v.s. independent variable (e.g., concentration-time profile)
    if(withBootstrap){
      p<-ggplot2::ggplot(uncertaintyBound_df,ggplot2::aes(x=independent_variable, y=median))+ggplot2::geom_line(colour="blue")+ggplot2::geom_ribbon(ggplot2::aes(ymin =lower, ymax =upper, x=independent_variable), alpha=0.2,fill="blue")+ggplot2::geom_point(ggplot2::aes(x=independent_variable, y=target), colour="red")

    }else if(length(unique(residualPlot_df$ind))==1){
      p<-ggplot2::ggplot(residualPlot_df,ggplot2::aes(x=independent_variable, y=model_fit))+ggplot2::geom_line()+ggplot2::geom_point(ggplot2::aes(x=independent_variable, y=target), colour="red")

    }else{
      p<-ggplot2::ggplot(residualPlot_df,ggplot2::aes(x=independent_variable, y=model_fit, group=ind))+ggplot2::geom_line()+ggplot2::geom_point(ggplot2::aes(x=independent_variable, y=target), colour="red")

    }

    p=p+ggplot2::xlab("Independent Variable")
    p=p+ggplot2::ylab("Dependent Variable")

    if(!is.na(dependentVariableTypeVector[1])){
      p=p+ggplot2::facet_grid(.~dependent_variable_type, scales = "free")
    }
  }else if(plotType==2){ #residual v.s. dependent variable (e.g., residual-fitted model profile)
    p<-ggplot2::ggplot(residualPlot_df,ggplot2::aes(x=model_fit, y=residual))+ggplot2::geom_point()
    p=p+ggplot2::geom_smooth()
    p=p+ggplot2::xlab("Dependent Variable")
    p=p+ggplot2::ylab("Residual")

    if(!is.na(dependentVariableTypeVector[1])){
      p=p+ggplot2::facet_grid(.~dependent_variable_type, scales = "free")
    }

    p=p+geom_hline(yintercept = 0)
  }else if(plotType==3){ #residual v.s. independent variable (e.g., residual-time profile)

    if(is.na(independentVariableVector_in[1])){
      p<-ggplot2::ggplot(residualPlot_df,ggplot2::aes(x=residual))+ggplot2::geom_histogram()
      p=p+ggplot2::xlab("Independent Variable")
      p=p+ggplot2::ylab("Residual")

      if(!is.na(dependentVariableTypeVector[1])){
        p=p+ggplot2::facet_grid(.~dependent_variable_type, scales = "free")
      }
    }else{
      p<-ggplot2::ggplot(residualPlot_df,ggplot2::aes(x=independent_variable, y=residual))+ggplot2::geom_point()
      p=p+ggplot2::geom_smooth()
      p=p+ggplot2::xlab("Independent Variable")
      p=p+ggplot2::ylab("Residual")

      if(!is.na(dependentVariableTypeVector[1])){
        p=p+ggplot2::facet_grid(.~dependent_variable_type, scales = "free")
      }

      p=p+geom_hline(yintercept = 0)
    }

  }else{
    p=ggplot2::ggplot(residualPlot_df)
  }


  if(plotType!=4&&(length(unique(dependentVariableTypeVector ))>1&length(plotRank)==1)){
    SSR_byVariable=aggregate(residualPlot_df$square_residual, by=list(dependent_variable_type=residualPlot_df$dependent_variable_type), FUN=sum)
    SSR_byVariable$SSR=formatC(SSR_byVariable$x, format = "g", digits = 3)

    p=p+
      geom_text(
        size    = 5,
        data    = SSR_byVariable,
        mapping = aes(x = Inf, y = Inf, label = paste0("SSR=",SSR)),
        hjust   = 1.05,
        vjust   = 1.5,
        inherit.aes = FALSE
      )
  }

  p
}




#' @title table_parameterSummary
#' @description
#' Make summary table of the approximate local minimizers found by CGNM.  If bootstrap analysis result is available, relative standard error (RSE: standard deviation/mean) will also be included in the table.
#' @param CGNM_result (required input) \emph{A list} stores the computational result from Cluster_Gauss_Newton_method() function in CGNM package.
#' @param indicesToInclude (default: NA) \emph{A vector of integers} indices to include in the plot (if NA, use indices chosen by the acceptedIndices() function with default setting).
#' @param ParameterNames (default: NA) \emph{A vector of strings} the user can supply so that these names are used when making the plot. (Note if it set as NA or vector of incorrect length then the parameters are named as x1, x2, ...)
#' @param ReparameterizationDef (default: NA) \emph{A vector of strings} the user can supply definition of reparameterization where each string follows R syntax and also refer to the ith element of the x vector (the input variable to the nonlinear function) as xi (e.g., if the first input variable to the nonlinear function is defined as x1=log10(Ka), then by setting "10^x1" as one of the strings in this vector, you can plot the violin plot of Ka)
#' @return \emph{A ggplot object} including the violin plot, interquartile range and median, minimum and maximum.
#' @examples
#'
#'model_analytic_function=function(x){
#'
#'  observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
#'  Dose=1000
#'  F=1
#'
#'  ka=10^x[1]
#'  V1=10^x[2]
#'  CL_2=10^x[3]
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
#' initial_lowerRange = c(-2,-2,-2), initial_upperRange =  c(1,2,2),
#' num_iter = 10, num_minimizersToFind = 100)
#'
#' table_parameterSummary(CGNM_result)
#' table_parameterSummary(CGNM_result, ParameterNames=c("Ka","V1","CL_2"),
#'      ReparameterizationDef=c("10^x1","10^x2","10^x3"))
#' @export
#' @import ggplot2

table_parameterSummary=function(CGNM_result, indicesToInclude=NA, ParameterNames=NA, ReparameterizationDef=NA){

  Name=NULL
  Kind_iter=NULL

  cutoff_pvalue=0.05
  numParametersIncluded=NA
  useAcceptedApproximateMinimizers=TRUE
  freeParaValues=makeParaDistributionPlotDataFrame(CGNM_result, indicesToInclude, cutoff_pvalue, numParametersIncluded, ParameterNames, ReparameterizationDef, useAcceptedApproximateMinimizers)

  summaryData_df=data.frame()

  variableNames=unique(freeParaValues$Name)
  for(vName in variableNames){

    tempVec=quantile(subset(freeParaValues, Name==vName&Kind_iter=="Final Accepted")$X_value, probs = c(0,0.25,0.5,0.75,1), na.rm = TRUE)

    if(!is.null(CGNM_result$bootstrapX)){
      bootstrapDS=subset(freeParaValues, Name==vName&Kind_iter=="Bootstrap")$X_value
      tempVec=c(tempVec, sd(bootstrapDS,na.rm = TRUE)/mean(bootstrapDS,na.rm = TRUE)*100)
    }

    summaryData_df=rbind(summaryData_df,as.numeric(tempVec))
  }
  rownames(summaryData_df)=variableNames

  if(!is.null(CGNM_result$bootstrapX)){

    colnames(summaryData_df)=c("Minimum","25 percentile", "Median","75 percentile", "Maximum", "Bootstrap RSE (%)")
  }else{
    colnames(summaryData_df)=c("Minimum","25 percentile", "Median","75 percentile", "Maximum")

  }

  return(summaryData_df)
}
