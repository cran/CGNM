
#' @title make_ShinyCGNM_observationData
#' @description
#' A helper function to write out the csv file that can be read in as the observation file in shinyCGNM
#' @param ID (required input) \emph{string or vector}
#' @param time (required input) \emph{number or numeric vector}
#' @param Observation_expression (required input) \emph{string or string vector}
#' @param Observed_value (required input) \emph{number or numeric vector}
#' @param ResidualError_model (required input) \emph{0 or 1} 0: additive residual model, 1: relative residual model
#' @param Memo (default: NA)  \emph{NA, string, or string vector} If TRUE plot absolute values of the residual.
#' @param fileName (default: NA) \emph{NA or string}
#' @return \emph{data.frame} if fileName is NA \emph{Null} if fileName is not NA but instead write out the csv file
#' @examples
#'
#' make_ShinyCGNM_observationData(
#' ID=1,
#' time=c(1,2,3,6,12,24),
#' Observation_expression="C_central",
#' Observed_value=c(0.1, 0.3, 0.6, 0.1, 0.05, 0.01),
#' ResidualError_model=1
#' )
#'
#' @export
#' @import utils

make_ShinyCGNM_observationData=function(ID, time, Observation_expression, Observed_value, ResidualError_model, Memo=NA, fileName=NA){

  out_df=data.frame(ID=ID,
                    time=time,
                    Observation_expression=Observation_expression,
                    Observed_value=Observed_value,
                    ResidualError_model=ResidualError_model,
                    Memo=Memo)

  if (is.na(fileName)){
    return(out_df)
  }else{
    write.csv(file = fileName, out_df, row.names = FALSE)
  }
}


#' @title make_ShinyCGNM_doseData
#' @description
#' A helper function to write out the csv file that can be read in as the dose file in shinyCGNM
#' @param ID (required input) \emph{string or vector}
#' @param dose (required input) \emph{number or numeric vector}
#' @param dosing.to (required input) \emph{string or string vector}
#' @param start.time (required input) \emph{number or numeric vector}
#' @param rate (default: NA) \emph{NA, number or numeric vector} infusion rate; set to NA for a bolus dose.
#' @param nbr.doses (default: 1) \emph{number or numeric vector} number of doses to administer starting at start.time, spaced dosing.interval apart.
#' @param dosing.interval (default: NA) \emph{NA, number or numeric vector} time between repeated doses; only relevant when nbr.doses is greater than 1.
#' @param fileName (default: NA) \emph{NA or string}
#' @return \emph{data.frame} if fileName is NA \emph{Null} if fileName is not NA but instead write out the csv file
#' @examples
#'
#' make_ShinyCGNM_doseData(
#' ID=seq(1,5),
#' dose=10,
#' dosing.to="A_admin",
#' start.time=0,
#' rate=NA
#' )
#'
#' @export
#' @import utils

make_ShinyCGNM_doseData=function(ID, dose, dosing.to, start.time, rate=NA, nbr.doses=1, dosing.interval=NA, fileName=NA){

  out_df=data.frame(ID=ID,
                    dose=dose,
                    dosing.to=dosing.to,
                    start.time=start.time,
                    rate=rate,
                    nbr.doses=nbr.doses,
                    dosing.interval=dosing.interval)

  if (is.na(fileName)){
    return(out_df)
  }else{
    write.csv(file = fileName, out_df, row.names = FALSE)
  }
}


#' @title make_ShinyCGNM_parameterInfo
#' @description
#' A helper function to write out the csv file that can be read in as the parameter info file in shinyCGNM
#' @param ParameterName (required input) \emph{string or vector} Name of the parameters. Must cover every parameter used in the ODE_text, every free symbol used in the dose data's dose, start.time, and rate columns, and every free symbol used in the initial condition data's value column.
#' @param Initial_lower_range (required input) \emph{number or numeric vector} lower range of the initial guess used in \code{\link{Cluster_Gauss_Newton_method}}'s \code{initial_lowerRange}.
#' @param Initial_upper_range (required input) \emph{number or numeric vector} upper range of the initial guess used in \code{\link{Cluster_Gauss_Newton_method}}'s \code{initial_upperRange}.
#' @param Lower_bound (default: 0) \emph{number, NA, or numeric vector} used in \code{\link{Cluster_Gauss_Newton_method}}'s \code{lowerBound}. Set NA if there is no lower bound.
#' @param Upper_bound (default: NA) \emph{number, NA, or numeric vector} used in \code{\link{Cluster_Gauss_Newton_method}}'s \code{upperBound}. Set NA if there is no upper bound.
#' @param VaryByID (default: 0) \emph{0, 3, or another positive integer, can be a vector} 0: the parameter is shared across all IDs. 3: the parameter varies by the first 3 characters of ID (e.g., to represent an occasion/period). Any other nonzero value: the parameter varies by the full ID.
#' @param MO_weight (default: 0) \emph{number or numeric vector} used in \code{\link{Cluster_Gauss_Newton_method}}'s \code{MO_weights} for the middle-out method. Set 0 if the parameter is not used in a middle-out constraint.
#' @param MO_value (default: NA) \emph{number, NA, or numeric vector} used in \code{\link{Cluster_Gauss_Newton_method}}'s \code{MO_values}. Required (non-NA) when the corresponding \code{MO_weight} is not 0.
#' @param Unit (default: NA) \emph{NA, string, or string vector} unit of the parameter, only used for documentation purposes in the generated script.
#' @param fileName (default: NA) \emph{NA or string}
#' @return \emph{data.frame} if fileName is NA \emph{Null} if fileName is not NA but instead write out the csv file
#' @examples
#'
#' make_ShinyCGNM_parameterInfo(
#' ParameterName=c("ka","CL","V1"),
#' Initial_lower_range=c(0.01,0.01,0.01),
#' Initial_upper_range=c(100,100,100)
#' )
#'
#' @export
#' @import utils

make_ShinyCGNM_parameterInfo=function(ParameterName, Initial_lower_range, Initial_upper_range, Lower_bound=0, Upper_bound=NA, VaryByID=0, MO_weight=0, MO_value=NA, Unit=NA, fileName=NA){

  out_df=data.frame(ParameterName=ParameterName,
                    Initial_lower_range=Initial_lower_range,
                    Initial_upper_range=Initial_upper_range,
                    Lower_bound=Lower_bound,
                    Upper_bound=Upper_bound,
                    VaryByID=VaryByID,
                    MO_weight=MO_weight,
                    MO_value=MO_value,
                    Unit=Unit)

  if (is.na(fileName)){
    return(out_df)
  }else{
    write.csv(file = fileName, out_df, row.names = FALSE)
  }
}


#' @title make_ShinyCGNM_initialCondition
#' @description
#' A helper function to build the (optional) initial condition table accepted by
#' \code{\link{generateCGNM_script}}'s \code{initialConditionData_table}. This lets you set
#' the value of an ODE compartment at time 0 explicitly, independently of (or together with)
#' any dosing into that compartment from a dose table.
#' @param ID (required input) \emph{string or vector} which individual(s) this initial condition applies to.
#' @param state (required input) \emph{string or vector} name of the ODE compartment (state variable) whose initial condition is being set.
#' @param value (required input) \emph{number, string, or vector} the initial condition. May be a plain number, or a string expression that references a parameter name (mirroring \code{\link{make_ShinyCGNM_doseData}}'s \code{dose} argument), e.g. \code{"V1*C0"}.
#' @param fileName (default: NA) \emph{NA or string}
#' @return \emph{data.frame} if fileName is NA \emph{Null} if fileName is not NA but instead write out the csv file
#' @examples
#'
#' make_ShinyCGNM_initialCondition(
#' ID="1",
#' state="central",
#' value=5
#' )
#'
#' @export
#' @import utils

make_ShinyCGNM_initialCondition=function(ID, state, value, fileName=NA){

  out_df=data.frame(ID=ID,
                    state=state,
                    value=value)

  if (is.na(fileName)){
    return(out_df)
  }else{
    write.csv(file = fileName, out_df, row.names = FALSE)
  }
}


#' @title make_ShinyCGNM_simulationTimepoints
#' @description
#' A helper function to build the simulation time point table accepted by
#' \code{\link{generateCGNM_script}}'s \code{simulationTimepoints_table} (used when
#' \code{includeSimulationCode = TRUE}). Unlike \code{\link{make_ShinyCGNM_observationData}},
#' there is no \code{Observed_value}/\code{ResidualError_model}, since nothing is being fit here
#' - these are just the time points and expressions to simulate and plot.
#' @param ID (required input) \emph{string or vector} which individual(s)/scenario(s) to simulate.
#' @param time (required input) \emph{number or numeric vector} time points to simulate at.
#' @param Observation_expression (required input) \emph{string or string vector} the name of a compartment (state variable) or derived (LHS) variable from the ODE, exactly as it appears in ODE_text - unlike \code{\link{make_ShinyCGNM_observationData}}'s Observation_expression, this is matched by name, not evaluated as an expression, so e.g. \code{"C_central"} works but \code{"log10(C_central)"} does not.
#' @param Memo (default: NA) \emph{NA, string, or string vector}
#' @param fileName (default: NA) \emph{NA or string}
#' @return \emph{data.frame} if fileName is NA \emph{Null} if fileName is not NA but instead write out the csv file
#' @examples
#'
#' make_ShinyCGNM_simulationTimepoints(
#' ID=1,
#' time=seq(0,24,by=0.5),
#' Observation_expression="C_central"
#' )
#'
#' @export
#' @import utils

make_ShinyCGNM_simulationTimepoints=function(ID, time, Observation_expression, Memo=NA, fileName=NA){

  out_df=data.frame(ID=ID,
                    time=time,
                    Observation_expression=Observation_expression,
                    Memo=Memo)

  if (is.na(fileName)){
    return(out_df)
  }else{
    write.csv(file = fileName, out_df, row.names = FALSE)
  }
}
