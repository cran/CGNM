repmat=function(subM,numrows,numcol){

  matrix_out=subM

  if(numrows>1){
    for(i in seq(1,numrows-1)){
      matrix_out=rbind(matrix_out,subM)
    }
  }

  tempM=matrix_out

  if(numcol>1){
    for(i in seq(1,numcol-1)){
      matrix_out=cbind(matrix_out,tempM)
    }
  }


  return(matrix_out)
}

matlabSum <-function(matrix,direction){
  if (direction==1){
    return(t(as.matrix(colSums(matrix))))
  }else{
    return(as.matrix(rowSums(matrix)))
  }
}

matlabMedian<-function(matrix,direction){
  out=c()
  if (direction==1){
    for(i in seq(1,dim(matrix)[2])){
      out=c(out,median(matrix[,i], na.rm =TRUE))
    }
    return(t(as.matrix(out)))
  }else{
    for(i in seq(1,dim(matrix)[1])){
      out=c(out,median(matrix[i,], na.rm =TRUE))
    }
    return((as.matrix(out)))
  }
}

dot=function(vec1,vec2){
  return(as.numeric(t(as.matrix(vec1))%*%as.matrix(vec2)))
}

matlabRand<-function(numRow,numCol){
  outMat=runif(numCol)
  for(i in seq(2,numRow)){
    outMat=rbind(outMat,runif(numCol))
  }
  return(outMat)
}

tryCatch_nonlinearFunction=function(x,num_observations, nonlinearFunction){
  out <- tryCatch(
    {
      nonlinearFunction(x)

    },
    error=function(cond) {
      return(rep(NaN,num_observations))
    },
    warning=function(cond) {
      nonlinearFunction(x)
    },
    finally={

    }
  )

  if(length(out)!=num_observations){
    out=rep(NaN,num_observations)
  }else if(is.na(sum(out))){
    out=rep(NaN,num_observations)
  }else if(is.infinite(sum(abs(out)))){
    out=rep(NaN,num_observations)
  }

  return(out)
}

#' @title Cluster_Gauss_Newton_method
#' @description Find multiple minimisers of the nonlinear least squares problem.
#' \deqn{argmin_x ||f(x)-y*||}
#' where
#' \enumerate{\item f: nonlinear function (e.g., mathematical model)
#' \item y*: target vector (e.g., observed data to fit the mathematical model)
#' \item x: variable of the nonlinear function that we aim to find the values that minimize (minimizers) the differences between the nonlinear function and target vector (e.g., model parameter)
#' }
#' Parameter estimation problems of mathematical models can often be formulated as nonlinear least squares problems.  In this context f can be thought at a model, x is the parameter, and y* is the observation.
#' CGNM iteratively estimates the minimizer of the nonlinear least squares problem from various initial estimates hence finds multiple minimizers.
#' Full detail of the algorithm and comparison with conventional method is available in the following publication, also please cite this publication when this algorithm is used in your research: Aoki et al. (2020) <doi.org/10.1007/s11081-020-09571-2>. Cluster Gauss–Newton method. Optimization and Engineering, 1-31.  As illustrated in this paper, CGNM is faster and more robust compared to repeatedly applying the conventional optimization/nonlinear least squares algorithm from various initial estimates. In addition, CGNM can realize this speed assuming the nonlinear function to be a black-box function (e.g. does not use things like adjoint equation of a system of ODE as the function does not have to be based on a system of ODEs.).
#' @param nonlinearFunction (required input) \emph{A function with input of a vector x of real number of length n and output a vector y of real number of length m.} In the context of model fitting the nonlinearFunction is \strong{the model}.  Given the CGNM does not assume the uniqueness of the minimizer, m can be less than n.  Also CGNM does not assume any particular form of the nonlinear function and also does not require the function to be continuously differentiable (see Appendix D of our publication for an example when this function is discontinuous).
#' @param targetVector (required input) \emph{A vector of real number of length m} where we minimize the Euclidean distance between the nonlinearFuncition and targetVector.  In the context of curve fitting targetVector can be though as \strong{the observational data}.
#' @param initial_lowerRange (required input) \emph{A vector of real number of length n} where each element represents  \strong{the lower range of the initial iterate}. Similarly to regular Gauss-Newton method, CGNM iteratively reduce the residual to find minimizers.  Essential differences is that CGNM start from the initial RANGE and not an initial point. Note that CGNM is an unconstraint optimization method so the final minimizer can be anywhere (and outside of this specified range).  In the parameter estimation problem, there often is a constraints to the parameters (e.g., parameters cannot be negative).. If you wish to constraint the parameter domain do so via parameter transformation (e.g., if parameter needs to be positive do log transform, if there is upper and lower bounds consider using logit transform.)
#' @param initial_upperRange (required input) \emph{A vector of real number of length n} where each element represents  \strong{the upper range of the initial iterate}.
#' @param num_minimizersToFind  (default: 250) \emph{A positive integer} defining number of approximate minimizers CGNM will find. We usually \strong{use 250 when testing the model and 1000 for the final analysis}.  The computational cost increase proportionally to this number; however, larger number algorithm becomes more stable and increase the chance of finding more better minimizers. See Appendix C of our paper for detail.
#' @param num_iteration (default: 25)  \emph{A positive integer} defining maximum number of iterations. We usually \strong{set 25 while model building and 100 for final analysis}.  Given each point terminates the computation when the convergence criterion is met the computation cost does not grow proportionally to the number of iterations (hence safe to increase this without significant increase in the computational cost).
#' @param saveLog (default: FALSE) \emph{TRUE or FALSE} indicating either or not to save computation result from each iteration in CGNM_log folder. It requires disk write access right in the current working directory. \strong{Recommended to set TRUE if the computation is expected to take long time} as user can retrieve intrim computation result even if the computation is terminated prematurely (or even during the computation).
#' @param runName (default: "") \emph{string} that user can ue to identify the CGNM runs. The run history will be saved in the folder name CGNM_log_<runName>.  If this is set to "TIME" then runName is automatically set by the run start time.
#' @param textMemo (default: "") \emph{string} that user can write an arbitrary text (without influencing computation). This text is stored with the computation result so that can be used for example to describe model so that the user can recognize the computation result.
#' @param algorithmParameter_initialLambda (default: 1) \emph{A positive number} for initial value for the regularization coefficient lambda see Appendix B of of our paper for detail.
#' @param algorithmParameter_gamma_in (default: 2) \emph{A positive number} a positive scalar value for adjusting the strength of the weighting for the linear approximation see Appendix A of our paper for detail.
#' @return list of a matrix X, Y,residual_history and initialX, as well as a list runSetting
#' \enumerate{\item X: \emph{a num_minimizersToFind by n matrix} which stores the approximate minimizers of the nonlinear least squares in each row. In the context of model fitting they are \strong{the estimated parameter sets}.
#' \item Y: \emph{a num_minimizersToFind by m matrix} which stores the nonlinearFunction evaluated at the corresponding approximate minimizers in matrix X above. In the context of model fitting each row corresponds to \strong{the model simulations}.
#' \item residual_history: \emph{a num_iteration by num_minimizersToFind matrix} storing sum of squares residual for all iterations.
#' \item initialX: \emph{a num_minimizersToFind by n matrix} which stores the set of initial iterates.
#' \item runSetting: a list containing all the input variables to Cluster_Gauss_Newton_method (i.e., nonlinearFunction, targetVector, initial_lowerRange, initial_upperRange ,algorithmParameter_initialLambda, algorithmParameter_gamma_in, num_minimizersToFind, num_iteration, saveLog, runName, textMemo).}
#' @examples
#' ##lip-flop kinetics (an example known to have two distinct solutions)
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
#' targetVector = observation, num_iteration = 10, num_minimizersToFind = 100,
#' initial_lowerRange = c(0.1,0.1,0.1), initial_upperRange =  c(10,10,10))
#'
#' acceptedApproximateMinimizers(CGNM_result)
#'
#  ## flip-flop kinetics using RxODE (an example known to have two distinct solutions)
#' \dontrun{
#' library(RxODE)
#' library(parallel) #if parallel library is loaded CGNM paralleizes the algorithm automatically
#'
#' model_text="
#' d/dt(X_1)=-ka*X_1
#' d/dt(C_2)=(ka*X_1-CL_2*C_2)/V1"
#'
#' model=RxODE(model_text)
#' #define nonlinearFunction
#' model_function=function(x){
  #'
  #' observation_time=c(0.1,0.2,0.4,0.6,1,2,3,6,12)
  #'
  #' theta <- c(ka=x[1],V1=x[2],CL_2=x[3])
  #' ev <- eventTable()
  #' ev$add.dosing(dose = 1000, start.time =0)
  #' ev$add.sampling(observation_time)
  #' odeSol=model$solve(theta, ev)
  #' log10(odeSol[,"C_2"])
  #'
  #' }
#'
#' observation=log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238))
#'
#' CGNM_result=Cluster_Gauss_Newton_method(nonlinearFunction=model_function,
#' targetVector = observation,
#' initial_lowerRange = c(0.1,0.1,0.1),initial_upperRange =  c(10,10,10))}
#'
#' @export
#' @import stats parallel


Cluster_Gauss_Newton_method <- function(nonlinearFunction, targetVector, initial_lowerRange, initial_upperRange , num_minimizersToFind=250, num_iteration=25, saveLog=FALSE, runName="", textMemo="",algorithmParameter_initialLambda=1, algorithmParameter_gamma_in=2){

  runSetting=list(nonlinearFunction=nonlinearFunction, targetVector=targetVector, initial_lowerRange=initial_lowerRange, initial_upperRange=initial_upperRange ,algorithmParameter_initialLambda=algorithmParameter_initialLambda, algorithmParameter_gamma_in=algorithmParameter_gamma_in, num_minimizersToFind=num_minimizersToFind, num_iteration=num_iteration, saveLog=saveLog, runName=runName, textMemo=textMemo)

  showIntermetiateResults=FALSE

  X_ul_in=t(matrix(c(initial_lowerRange,initial_upperRange), nrow=length(initial_upperRange)))

  saveFolderName="CGNM_log"
  if(runName=="TIME"){
    saveFolderName=paste0(saveFolderName,"_",Sys.time())
  }else if(runName!=""){
    saveFolderName=paste0(saveFolderName,"_",runName)
  }

  if(saveLog){
    dir.create(saveFolderName)
  }

  X_history=c()
  Y_history=c()

    algorithmParameter_gamma=algorithmParameter_gamma_in

    method <- 'CGNM'# set algorithm name for the log files)
    descriptionText <- paste0(toString(num_minimizersToFind),'samples_initLambda',toString(algorithmParameter_initialLambda),'_distanceOrder',toString(algorithmParameter_gamma))

    # setting the name depending on if we do restart or not.  Restart is an
    # experimental feature so I suggest not to use at this point.

    num_parameters <- dim(X_ul_in)[2]
    num_observations = length(targetVector)

    timeOneParaTook <- c()

  #  tic

    ############################
    # 1) Pre-iteration process #
    ############################
    X <- matrix(0,num_minimizersToFind,num_parameters)# initialize X (We keep the same naming for the variables as in the manuscript)
    Y <- matrix(1,num_minimizersToFind,num_observations)# initialize Y
    Y_new <- matrix(1,num_minimizersToFind,num_observations)# initialize Y
    lambda_vec <- matrix(1,1,num_minimizersToFind)*algorithmParameter_initialLambda# initialise regularisation parameter lambda

    residual_history <- c()# initialise a matrix that stores the SSR for all iterations


    ## Generate initial cluster
    is_alive <- matrix(0,1,num_minimizersToFind)# initialise the vector to keep track of if the nonlinear function was able to be evaluated at the randomly generated x.


    # repeat the randomsampling of x until num_minimizersToFind of 'valid' x are sampled. Where we consider x to be valid if f(x) can be evaluated.
    while(sum(is_alive)<num_minimizersToFind){

      print("Generating initial cluster.")

        # random sampling of X
        X_temp <- matlabRand(num_minimizersToFind,length(X_ul_in[1,]))*(repmat(X_ul_in[2,],num_minimizersToFind,1)-repmat(X_ul_in[1,],num_minimizersToFind,1))+repmat(X_ul_in[1,],num_minimizersToFind,1)

        # replace the rows of X matrix with randomly sampled x if x was
        # determined to be not valid
        X[is_alive==0,]=X_temp[is_alive==0,]

        # if("parallel" %in% tolower((.packages()))){
        #
        #   Y_list=mclapply(split(X, rep(seq(1:nrow(X)),ncol(X))), tryCatch_nonlinearFunction, num_observations=num_observations, nonlinearFunction=nonlinearFunction, mc.cores = 4)
        #
        # }else{
        Y_list=lapply(split(X, rep(seq(1:nrow(X)),ncol(X))), tryCatch_nonlinearFunction, num_observations=num_observations, nonlinearFunction=nonlinearFunction)
        #}

        Y=t(matrix(unlist(Y_list),ncol=length(Y_list)))


        # compute SSR
        residual <- matlabSum((Y-repmat(t(as.matrix(targetVector)),num_minimizersToFind,1))^2,2)

        if(showIntermetiateResults){
          print(residual)
        }

        # determine valid x (in this case if the residual is not a number
        # we consider x to be valid)
        is_alive=is.nan(residual)==0

    }

    # store the residual vector for the initial cluster into the
    # residual_history matrix
    residual_history <- rbind(residual_history,residual)

    # store the matrix X and Y to X_history and Y_history, respectively
    X_history[[1]] <- X
    Y_history[[1]] <- Y

    lambda_history <- t(as.matrix(lambda_vec))

    prev_residual <- residual

    for (k in seq(1,num_iteration)){

        out_temp <- main_iteration(X, Y, lambda_vec, X_ul_in[1,], X_ul_in[2,], targetVector, algorithmParameter_gamma)

        X_new=X+out_temp

        # if("parallel" %in% tolower((.packages()))){
        #   Y_list=mclapply(split(X_new, rep(seq(1:nrow(X_new)),ncol(X_new))), tryCatch_nonlinearFunction, num_observations=num_observations, nonlinearFunction=nonlinearFunction, mc.cores = 4)
        #
        #   }else{
        Y_list=lapply(split(X_new, rep(seq(1:nrow(X_new)),ncol(X_new))), tryCatch_nonlinearFunction, num_observations=num_observations, nonlinearFunction=nonlinearFunction)
        #}
        Y_new=t(matrix(unlist(Y_list),ncol=length(Y_list)))

        if(showIntermetiateResults){
          print(Y_new)
        }

        residual_new <- t(matlabSum((Y_new-repmat(t(targetVector),num_minimizersToFind,1))^2,2))


        residual <- prev_residual
        for (i in seq(1,num_minimizersToFind)){
            if (!is.nan(residual_new[i])&&(prev_residual[i]>residual_new[i])){
                X[i,] <- X[i,]+out_temp[i,]
                residual[i] <- residual_new[i]
                Y[i,] <- Y_new[i,]
                lambda_vec[i] <- lambda_vec[i]/10
            } else {
                lambda_vec[i] <- lambda_vec[i]*10
            }
        }


        X_history[[k+1]] <- X
        Y_history[[k+1]] <- Y
        lambda_history <- cbind(lambda_history,t(lambda_vec))

        residual_history <- cbind(residual_history,residual)

        prev_residual <- residual

       # print(matlabMedian(residual_history,1))
          print(paste0("Iteration:",toString(k),"  Median sum of squares residual=", toString(median(prev_residual))))


        if(saveLog){

          result_list=list(X=X,Y=Y,residual_history=residual_history, initialX=X_history[[1]], runSetting=runSetting )
          save(file=paste0(saveFolderName,'/iteration_',toString(k),'.RDATA'), result_list)

          # write.csv(file=paste0('CGNM_log/',method,'_',descriptionText,'_xHistory_',toString(k),'.csv'),X_history[[k]])
          # write.csv(file=paste0('CGNM_log/',method,'_',descriptionText,'_yHistory_',toString(k),'.csv'),Y_history[[k]])
          #
          # write.csv(file=paste0('CGNM_log/',method,'_',descriptionText,'_residualHistory.csv'),t(residual_history))
          # write.csv(file=paste0('CGNM_log/',method,'_',descriptionText,'_compTime.csv'),timeOneParaTook)
          # write.csv(file=paste0('CGNM_log/',method,'_',descriptionText,'_lambdaHistory.csv'),t(lambda_history))
        }




        if(median(prev_residual)==0){
          break
        }
    }

    # write.csv(file=paste0(saveFolderName,'/',method,'_',descriptionText,'_xHistory_',toString(k),'.csv'),X_history[[k]])
    # write.csv(file=paste0(saveFolderName,'/',method,'_',descriptionText,'_yHistory_',toString(k),'.csv'),Y_history[[k]])
    #
    # write.csv(file=paste0(saveFolderName,'/',method,'_',descriptionText,'_residualHistory.csv'),t(residual_history))
    # write.csv(file=paste0(saveFolderName,'/',method,'_',descriptionText,'_compTime.csv'),timeOneParaTook)
    # write.csv(file=paste0(saveFolderName,'/',method,'_',descriptionText,'_lambdaHistory.csv'),t(lambda_history))

    result_list=list(X=X,Y=Y,residual_history=residual_history, runSetting=runSetting, initialX=X_history[[1]] )

    if(saveLog){
      save(file=paste0(saveFolderName,'/iteration_final.RDATA'), result_list)
    }

  return(result_list)

}


main_iteration <-  function(X_in, Y_in, lambdaV, minX, maxX, targetVector, algorithmParameter_gamma){

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


    x_width <- maxX-minX

    relative_distance <- matrix(0,dim(X_in)[1],dim(X_in)[1])
    rec_relative_distance <- matrix(0,dim(X_in)[1],dim(X_in)[1])

    for (i in seq(1,dim(X_in)[1])){
        relative_distance[i,i] <- 0
        rec_relative_distance[i,i] <- 0

        if(i>1){
          for (j in seq(1,(i-1))){
            relative_distance[i,j] <- 0
            relative_distance[j,i] <- 0

            for (k in seq(1,dim(X_in)[2])){
              relative_distance[i,j] <- relative_distance[i,j]+((X_in[i,k]-X_in[j,k])/x_width[k])^2
            }

            if (YisNotNaN[i]==1&YisNotNaN[j]==1){

              relative_distance[i,j] <- (relative_distance[i,j])

              relative_distance[j,i] <- relative_distance[i,j]
              rec_relative_distance[i,j] <- 1/relative_distance[i,j]
              rec_relative_distance[j,i] <- rec_relative_distance[i,j]
            } else {
              relative_distance[i,j] <- 0
              rec_relative_distance[i,j] <- 0
              rec_relative_distance[j,i] <- 0
            }

          }

        }

    }




    # Linear approximation will be A_with_y0  X_append_one_column^T \approx Y^T now solve it for A_with_y0
    # We solve for \min || X_append_one_column  A_with_y0^T - Y ||_F
    # use CGNR and solve each row of A_with_y0 (i.e., each column of A_with_y0^T)
    A <- matrix(0, dim(Y_in)[2],dim(X_in)[2])
    deltaX <- matrix(0, dim(X_in)[1],dim(X_in)[2])
    delta_X <- matrix(0, dim(X_in)[1],dim(X_in)[2])
    tempOnes <- matrix(1, dim(Y_in)[1],1)

    for (k in seq(1,dim(Y_in)[1])){

      if(lambdaV[k]<10^10){
        ## 2-1): Construct weighted linear approximation of the nonlinear function
        for (i in seq(1,dim(X_in)[1])){
          deltaX[i,] <- X_in[i,]-X_in[k,]
        }

        # use the following code if one wishes to use matlab's matrix
        # inverse function instead of CGNR
        #        AT=(diag((rec_relative_distance[ ,k]).^algorithmParameter_gamma)*deltaX)';
        #        ATAinv=inv(AT*AT');
        #        A=(ATAinv*(AT*diag((rec_relative_distance[ ,k]).^algorithmParameter_gamma)*(Y_in-repmat(Y_in[k,],dim(Y_in)[1],1))))';

        for (i in seq(1,dim(Y_in)[2])){
          A[i,] <- CGNR_ATAx_ATb_with_weight(deltaX,Y_in[ ,i]-Y_in[k,i]*tempOnes,rec_relative_distance[ ,k]^algorithmParameter_gamma)
        }

        ## 2-2): Solve for x that minimizes the resodual using the weighted linear approximation
        delta_X[k,] <- CGNR_ATAx_ATb_with_reg(A, t(targetVector-Y_in[k,]),lambdaV[k])
      }

    }

    return(delta_X)
}

CGNR_ATAx_ATb <- function(A, b){#solves A^TA x = A^Tb

    num_max_iter <- dim(A)[2]


    x <- matrix(0,dim(A)[2],1)

    if (dim(A)[1]==length(b)){

        r <- t(A)%*%b
        p <- r

        innerProd_rr_old <- dot(r,r)

        for (i in seq(1,num_max_iter)){
            tempVector <- A%*%p
            tempD <- dot(tempVector,tempVector)
            if (is.na(tempD)||tempD==0){
                break
            }
            alpha <- innerProd_rr_old/tempD
            x <- x+alpha*p
            r <- r-alpha*t(A)%*%tempVector
            innerProd_rr_new <- dot(r,r)
            beta <- innerProd_rr_new/innerProd_rr_old
            innerProd_rr_old <- innerProd_rr_new
            if (innerProd_rr_old==0){
                break
            }
            p <- r+beta*p
        }
    } else {
        print('CGNR error: matrix A and vetor b sizes do not make sense.')
        print(paste('Size of A: ', toString(dim(A)[1])))
        print('Size of b: ', toString(length(b)))
    }
  return(x)
}

CGNR_ATAx_ATb_with_weight=function(A, b, weight){

    for (i in seq(1,min(length(weight),dim(A)[1]))){
        if (weight[i]==0){
            for (j in 1:dim(A)[2]){
                A[i,j] <- 0
            }
            b[i] <- 0
        } else {
            A[i,] <- weight[i]*A[i,]
            b[i] <- weight[i]*b[i]
        }
    }
    return(CGNR_ATAx_ATb(A,b))
}


CGNR_ATAx_ATb_with_reg=function(A, b, lambda){#   //solves (A^TA+lambda I) x = A^Tb

    num_max_iter <- dim(A)[1]

    x <- matrix(0,dim(A)[2],1)

    if (dim(A)[1]==length(b)){

        r <- t(A)%*%t(b)
        p <- r

        innerProd_rr_old <- dot(r,r)
        AtA <- t(A)%*%A
        AtA_lambda <- AtA

        for (i in seq(1,dim(AtA)[1])){
            AtA_lambda[i,i] <- AtA_lambda[i,i]+lambda
        }


        for (i in seq(1,num_max_iter)){
            tempD <- dot(p,(AtA_lambda%*%p))
            if (is.na(tempD)||tempD==0){
                break
            }
            alpha <- innerProd_rr_old/tempD
            x <- x+alpha*p
            r <- r-alpha*(AtA_lambda%*%p)
            innerProd_rr_new <- dot(r,r)
            beta <- innerProd_rr_new/innerProd_rr_old
            innerProd_rr_old <- innerProd_rr_new
            if (innerProd_rr_old==0){
                break
            }
            p <- r+beta*p
        }
    } else {
        print('CGNR with reg error: matrix A and vetor b sizes do not make sense.')
        print(paste('Size of A: ', toString(dim(A)[1])))
        print(paste('Size of b: ', toString(length(b))))
    }
    return(x)
}

