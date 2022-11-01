#source("R/Cluster_Gauss_Newton_method.R")


plot_simulation=function(simulationResult_matrix, independentVariableVector=NA, dependentVariableTypeVector=NA, plotIndividual=5){

  independent_variable=NULL
  lower=NULL
  upper=NULL
  target=NULL
  model_fit=NULL
  ind=NULL
  residual=NULL
  SSR=NULL

  individual=NULL
  label=NULL
  newvalue=NULL

  independentVariableVector_in=independentVariableVector


  if(is.na(independentVariableVector[1])||length(independentVariableVector)!=dim(simulationResult_matrix)[2]){
    independentVariableVector=seq(1, dim(simulationResult_matrix)[2])
  }else{
    independentVariableVector=as.numeric(independentVariableVector)
  }

  if(is.na(dependentVariableTypeVector[1])||length(dependentVariableTypeVector)!=dim(simulationResult_matrix)[2]){
    dependentVariableTypeVector=NA
  }else{
    dependentVariableTypeVector=as.factor(dependentVariableTypeVector)
  }

  median_vec=c()
  percentile5_vec=c()
  percentile95_vec=c()
  uncertaintyBound_df=data.frame()

  for(i in seq(1, dim(simulationResult_matrix)[2])){
    tempQ=quantile(simulationResult_matrix[,i], prob=c(0.05,0.5,0.95), na.rm=TRUE)
    median_vec=c(median_vec, tempQ[2])
    percentile5_vec=c(percentile5_vec, tempQ[1])
    percentile95_vec=c(percentile95_vec, tempQ[3])
  }

  uncertaintyBound_df=data.frame(independent_variable=independentVariableVector, dependent_variable_type=dependentVariableTypeVector, median=median_vec, lower=percentile5_vec, upper=percentile95_vec)

  uncertaintyBound_df$label=with(uncertaintyBound_df, paste0(independent_variable, dependent_variable_type))


  uncertaintyBound_mean_df=data.frame()
  for(nowLabel in unique(uncertaintyBound_df$label)){

    tempDF=subset(uncertaintyBound_df, label==nowLabel)
    tempDF$median=mean(tempDF$median, na.rm=TRUE)
    tempDF$percentile5_vec=mean(tempDF$percentile5_vec, na.rm=TRUE)
    tempDF$percentile95_vec=mean(tempDF$percentile95_vec, na.rm=TRUE)
    uncertaintyBound_mean_df=rbind(uncertaintyBound_mean_df, tempDF[1,])

  }


  p<-ggplot2::ggplot(uncertaintyBound_mean_df,ggplot2::aes(x=independent_variable, y=median))+ggplot2::geom_line(colour="blue")+ggplot2::geom_ribbon(ggplot2::aes(ymin =lower, ymax =upper, x=independent_variable), alpha=0.2,fill="blue")

  if(plotIndividual>0){

    if(plotIndividual<dim(simulationResult_matrix)[1]){
      colSD_simY=col_sd(simulationResult_matrix)
      maxSD_ind=which(colSD_simY==max(colSD_simY))

      quantileAtMaxSD=quantile(simulationResult_matrix[,maxSD_ind], prob=(seq(0,1,1/(plotIndividual-1))), type=3)

      indexToPloe_vec=which(simulationResult_matrix[,maxSD_ind]%in%as.numeric(quantileAtMaxSD))
    }else{
      indexToPloe_vec=seq(1,dim(simulationResult_matrix)[1])
    }

    print(indexToPloe_vec)


    individualPlot_df=data.frame()
    for(i in indexToPloe_vec){
      individualPlot_df=rbind(individualPlot_df, data.frame(individual=i,independent_variable=independentVariableVector, dependent_variable_type=dependentVariableTypeVector, value=as.numeric(simulationResult_matrix[i,])))
    }

    individualPlot_df$label=with(individualPlot_df, paste0(independent_variable, dependent_variable_type, "ind:", individual))

    individualPlot_new_df=data.frame()
    for(nowLabel in unique(individualPlot_df$label)){

      tempDF=subset(individualPlot_df, label==nowLabel)
      tempDF$newvalue=mean(tempDF$value, na.rm=TRUE)
      individualPlot_new_df=rbind(individualPlot_new_df, tempDF[1,])

    }


    p=p+ggplot2::geom_line(data=individualPlot_new_df,ggplot2::aes(x=independent_variable, y=newvalue, group=individual) )
  }


  p=p+ggplot2::facet_wrap(.~dependent_variable_type, scales = "free")

  p
}



chooseXof=function(choice_vec, numToChoose){
  out=matrix(choice_vec,nrow=length(choice_vec), ncol=1)
  for(j in seq(1,numToChoose-1)){
    oldOut=out

    out=matrix(NA,nrow = 1,ncol = dim(oldOut)[2]+1)
    for(i in seq(1,dim(oldOut)[1])){
      for(toAdd in choice_vec[(choice_vec>oldOut[i,dim(oldOut)[2]])]){
        out=rbind(out, c(oldOut[i,],toAdd))
      }

    }

    out=out[seq(2,dim(out)[1]),]
  }
  return(out)
}



miniumDesign=function(jacobian){
  jac_abs=abs(jacobian)

  chosenIndex=c()

  jac_temp=-jac_abs
  for(i in seq(1, dim(jac_temp)[2])){
    chosenIndex=c(chosenIndex,which(rank(jac_temp[,i])==1))
    jac_temp[chosenIndex,]=0
  }

  return(chosenIndex)
}

detFIM=function(jacobian_in, design, obsVariance_vec=NA){

  jacobian_ToUse=jacobian_in[design,]


  if(!is.na(obsVariance_vec)){
    if(length(obsVariance_vec)==dim(jacobian_ToUse)[1]){
      jacobian_ToUse=jacobian_ToUse*diag(sqrt(obsVariance_vec[design]))
    }
  }
  #jacobian_svd=svd((jacobian))

  #print(jacobian_svd$d)
  #return(sum(log(jacobian_svd$d)))

#  out=determinant(t(jacobian_ToUse)%*%jacobian_ToUse,logarithm=TRUE)$modulus

  FIM_chol=chol(t(jacobian_ToUse)%*%jacobian_ToUse)
  diagFIMinv=diag(chol2inv(FIM_chol))
  detFIM=sum(log(diag(FIM_chol)))

#  FIM_chol_update=cholupdate(FIM_chol, jacobian_in[1,])
#  out=sum(log(diag(FIM_chol_update)))


 # out=sum(log(svd(jacobian_ToUse)$d))
  return(list(detFIM=detFIM, diagFIMinv=diagFIMinv))
}

designEvaluations=function(jacobian, design_list, obsVariance_vec=NA){
  detFIM_vec=c()
  diagInvFIM_matrix=matrix(NA, nrow = length(design_list), ncol=dim(jacobian)[2])
  for(i in seq(1,length(design_list))){
    detFIM_result=detFIM(jacobian, design_list[[i]], obsVariance_vec)

    detFIM_vec=c(detFIM_vec, detFIM_result$detFIM)
    diagInvFIM_matrix[i,]=detFIM_result$diagFIMinv
  }
  return(list(detFIM_vec=detFIM_vec,diagInvFIM_matrix=diagInvFIM_matrix))
}

optimizeDesign=function(X_in, Y_in){
  jac_list=createApproxJacobianlist(X_in, Y_in)
  min_design_list=list()

  for(i in seq(1,length(jac_list))){
    min_design=miniumDesign(jac_list[[i]])
    min_design_list[[i]]=min_design
  }

  return(((min_design_list)))

 # return(designEvaluations(jac_list[[i]],min_design_list[[i]]))


}




createApproxJacobianlist <-  function(X_in, Y_in,  algorithmParameter_gamma=2){

  #// aliveIndex: index of the parameter set in the cluster whose SSR have decreased in the previous iteration.
  #// X_in: all set of parameters in the cluster from the previous itaration  X_in[i][j]: jth parameter value in the ith parameter set in the cluster
  #// Y_in: all set of the function value (solution of the forward problem, f(parameter)) from the previous iteration
  # global algorithmParameter_gamma


  ##//STEP 1: Linear Approximation

  YisNotNaN <- matrix(1,dim(X_in)[1],1)

  for (k in 1:dim(Y_in)[1]){
    for (i in 1:dim(Y_in)[2]){
      if (is.nan(Y_in[k,i])|is.infinite(Y_in[k,i])){
        YisNotNaN[k] <- 0
      }
    }
  }

  #// Linear approximation will be A x + y_o \approx f(x) = y for all y in Y matrix (stored in each row of Y matrix)
  #// This can be written as A_with_y0 X_append_one_column^T \approx Y^T


  kmean_result=optimal_kmeans(X_in)

  while(min(table(kmean_result$cluster))<dim(X_in)[2]|min(kmean_result$withinss)==0){
    kmean_result=kmeans(col_normalize(X_in),centers = length(unique(kmean_result$cluster))-1)
  }

  notAcceptableClusters=as.numeric(which(table(kmean_result$cluster)==1,arr.ind = TRUE))

  relative_distance <- matrix(0,dim(X_in)[1],dim(X_in)[1])
  rec_relative_distance <- matrix(0,dim(X_in)[1],dim(X_in)[1])


  numClusters=length(unique(kmean_result$cluster))
  sdMatrix=matrix(0,  dim(X_in)[2],dim(X_in)[1])

  for(i in seq(1, numClusters)){
    for(j in seq(1, dim(X_in)[2])){
      sdMatrix[j,kmean_result$cluster==i]=sd(X_in[kmean_result$cluster==i,j])
    }
  }


  tX_in=t(X_in)

  for (i in seq(1,dim(X_in)[1])){

    useIndex=(kmean_result$cluster==kmean_result$cluster[i])
    relative_distance[i,useIndex] = colSums(((tX_in[,useIndex]-tX_in[,i])/(sdMatrix[,useIndex]))^2)

  }


  rec_relative_distance=1/relative_distance
  rec_relative_distance[relative_distance==0]=0



  # Linear approximation will be A_with_y0  X_append_one_column^T \approx Y^T now solve it for A_with_y0
  # We solve for \min || X_append_one_column  A_with_y0^T - Y ||_F
  # use CGNR and solve each row of A_with_y0 (i.e., each column of A_with_y0^T)
  A <- matrix(0, dim(Y_in)[2],dim(X_in)[2])
  deltaX <- matrix(0, dim(X_in)[1],dim(X_in)[2])
  delta_X <- matrix(0, dim(X_in)[1],dim(X_in)[2])
  tempOnes <- matrix(1, dim(Y_in)[1],1)

  A_list=list()

  for (k in seq(1,dim(Y_in)[1])){
    if((!(kmean_result$cluster[k] %in% notAcceptableClusters))){
      ## 2-1): Construct weighted linear approximation of the nonlinear function
      for (i in seq(1,dim(X_in)[1])){
        deltaX[i,] <- X_in[i,]-X_in[k,]
      }

      # use the following code if one wishes to use MASS's matrix
      # inverse function (ginv) instead of CGNR



      rec_relative_distance_tempvec=(rec_relative_distance[ ,k])^algorithmParameter_gamma
      rec_relative_distance_tempvec[is.infinite(rec_relative_distance_tempvec)]=max(c(10^10,rec_relative_distance_tempvec[!is.infinite(rec_relative_distance_tempvec)]))

      AT=t(diag(rec_relative_distance_tempvec)%*%deltaX)

      ATAinv=MASS::ginv(AT%*%t(AT))

      delta_Y=Y_in

      for(i in seq(1, dim(Y_in)[2])){
        delta_Y[,i]=Y_in[,i]-Y_in[k,i]
      }

      A=t(ATAinv%*%(AT%*%diag((rec_relative_distance[ ,k])^algorithmParameter_gamma)%*%(delta_Y)));


    }
    A_list[[k]]=A

  }

  return(A_list)
}
