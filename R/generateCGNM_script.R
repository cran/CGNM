# ---------------------------------------------------------------------------
# Internal code-generation helpers.
#
# pasteWithCollapse_WithApprox80charLimits(), puttogether_model_code(),
# makeCGNM_runCode(), and simulation_code_text() below are ported unchanged
# from inst/shinyCGNM/CodeGenerationRelatedFunctions.R, which shinyCGNM() uses
# to turn ODE_text plus its data-entry tables into a runnable CGNM script.
# They are plain functions (no Shiny dependency), so generateCGNM_script()
# below reuses them directly instead of re-implementing the code generation.
# .makeCodeForPosthoc_middleout() is adapted (not verbatim; see its own
# comment) from inst/shinyCGNM/server.R's makeCodeForPosthoc_middleout().
# ---------------------------------------------------------------------------

pasteWithCollapse_WithApprox80charLimits=function(paste_vec, collapse){
  numTotalChar=nchar(paste(paste_vec, collapse = collapse))
  numToDivide=ceiling(numTotalChar/80)
  numPerDivide=ceiling(length(paste_vec)/numToDivide)

  splitted_list=split(paste_vec, ceiling(seq_along(paste_vec)/numPerDivide))
  outText_vec=c()
  for(splitNameNu in names(splitted_list)){
    outText_vec=c(outText_vec,  paste(splitted_list[[splitNameNu]], collapse = collapse))
  }

  return(paste(outText_vec, collapse = gsub(",",",\n",collapse)))
}


puttogether_model_code = function(input, ll, rv, for_simulation=FALSE) {
  ParameterName_vec = ll$parameterInfo_dat$ParameterName
  Initial_upper_range_vec = ll$parameterInfo_dat$Initial_upper_range
  Initial_lower_range_vec = ll$parameterInfo_dat$Initial_lower_range
  Upper_bound_vec = ll$parameterInfo_dat$Upper_bound
  Lower_bound_vec = ll$parameterInfo_dat$Lower_bound
  MO_weight_vec = ll$parameterInfo_dat$MO_weight
  MO_value_vec = ll$parameterInfo_dat$MO_value
  runName=input$runNameText
  runName=gsub("[[:punct:]]","", runName)
  runName=gsub(" ","_", runName)

  unique_IDs = sort(unique(c(ll$ObservedData_dat$ID, ll$Dose_dat$ID)))

  unique_3_IDs=unique(substr(unique_IDs,1,3))

  UseIndVar = FALSE
  IndividualParameter_text = ""
  if (!is.null(ll$parameterInfo_dat$VaryByID)) {
    if (sum(ll$parameterInfo_dat$VaryByID != 0) > 0) {
      UseIndVar = TRUE

      ParameterName_vec = subset(ll$parameterInfo_dat, VaryByID == 0)$ParameterName
      Initial_upper_range_vec = subset(ll$parameterInfo_dat, VaryByID ==
                                         0)$Initial_upper_range
      Initial_lower_range_vec = subset(ll$parameterInfo_dat, VaryByID ==
                                         0)$Initial_lower_range
      Upper_bound_vec = subset(ll$parameterInfo_dat, VaryByID == 0)$Upper_bound
      Lower_bound_vec = subset(ll$parameterInfo_dat, VaryByID == 0)$Lower_bound

      MO_weight_vec = subset(ll$parameterInfo_dat, VaryByID == 0)$MO_weight
      MO_value_vec = subset(ll$parameterInfo_dat, VaryByID == 0)$MO_value


      temp_paraInfo = subset(ll$parameterInfo_dat, VaryByID != 0&VaryByID != 3)
      indParameterNames = temp_paraInfo$ParameterName


      temp_paraInfo_3 = subset(ll$parameterInfo_dat, VaryByID == 3)
      indParameterNames_3 = temp_paraInfo_3$ParameterName


      Parameter_text = paste(paste0(
        subset(ll$parameterInfo_dat, VaryByID == 0)$ParameterName,
        "=x[",
        seq(1, dim(
          subset(ll$parameterInfo_dat, VaryByID == 0)
        )[1]),
        "]"
      ),
      collapse = "\n")

      indParaNames = c()

      if(dim(temp_paraInfo)[1]>0){
        for (i in seq(1, dim(temp_paraInfo)[1])) {
          indParaNames = c(indParaNames,
                           paste0(temp_paraInfo$ParameterName[i], "_ID", unique_IDs))
          ParameterName_vec = c(
            ParameterName_vec,
            paste0(temp_paraInfo$ParameterName[i], "_ID", unique_IDs)
          )
          Initial_upper_range_vec = c(
            Initial_upper_range_vec,
            rep(
              temp_paraInfo$Initial_upper_range[i],
              length(unique_IDs)
            )
          )
          Initial_lower_range_vec = c(
            Initial_lower_range_vec,
            rep(
              temp_paraInfo$Initial_lower_range[i],
              length(unique_IDs)
            )
          )
          Upper_bound_vec = c(Upper_bound_vec,
                              rep(temp_paraInfo$Upper_bound[i], length(unique_IDs)))
          Lower_bound_vec = c(Lower_bound_vec,
                              rep(temp_paraInfo$Lower_bound[i], length(unique_IDs)))
          MO_weight_vec = c(MO_weight_vec,
                            rep(temp_paraInfo$MO_weight[i], length(unique_IDs)))
          MO_value_vec = c(MO_value_vec, rep(temp_paraInfo$MO_value[i], length(unique_IDs)))
        }
      }

      if(dim(temp_paraInfo_3)[1]>0){
        for (i in seq(1, dim(temp_paraInfo_3)[1])) {
          indParaNames = c(indParaNames,
                           paste0(temp_paraInfo_3$ParameterName[i], "_ID", unique_3_IDs))
          ParameterName_vec = c(
            ParameterName_vec,
            paste0(temp_paraInfo_3$ParameterName[i], "_ID", unique_3_IDs)
          )
          Initial_upper_range_vec = c(
            Initial_upper_range_vec,
            rep(
              temp_paraInfo_3$Initial_upper_range[i],
              length(unique_3_IDs)
            )
          )
          Initial_lower_range_vec = c(
            Initial_lower_range_vec,
            rep(
              temp_paraInfo_3$Initial_lower_range[i],
              length(unique_3_IDs)
            )
          )
          Upper_bound_vec = c(Upper_bound_vec,
                              rep(temp_paraInfo_3$Upper_bound[i], length(unique_3_IDs)))
          Lower_bound_vec = c(Lower_bound_vec,
                              rep(temp_paraInfo_3$Lower_bound[i], length(unique_3_IDs)))
          MO_weight_vec = c(MO_weight_vec,
                            rep(temp_paraInfo_3$MO_weight[i], length(unique_3_IDs)))
          MO_value_vec = c(MO_value_vec, rep(temp_paraInfo_3$MO_value[i], length(unique_3_IDs)))
        }
      }


      Parameter_text = paste0(Parameter_text, "\n", paste0(paste0(
        indParaNames,
        "=x[",
        seq(
          dim(subset(ll$parameterInfo_dat, VaryByID == 0))[1] + 1,
          dim(subset(ll$parameterInfo_dat, VaryByID == 0))[1] + length(indParaNames)
        ),
        "]"
      ),
      collapse = "\n"))

    } else{
      Parameter_text = paste(paste0(
        ll$parameterInfo_dat$ParameterName,
        "=x[",
        seq(1, dim(ll$parameterInfo_dat)[1]),
        "]"
      ),
      collapse = "\n")
    }
  } else{
    Parameter_text = paste(paste0(
      ll$parameterInfo_dat$ParameterName,
      "=x[",
      seq(1, dim(ll$parameterInfo_dat)[1]),
      "]"
    ),
    collapse = "\n")
  }

  odeCodeText = paste0("ODE_text=\"\n",
                       input$ODE_text,
                       "\"\n compiledModel=RxODE(ODE_text)")


  parameterText = paste0(ifelse(for_simulation,paste0("\nsimulation_function_",runName,"=function(x){\n"),paste0("\n\n\nmodel_function_",runName,"=function(x){\n")), Parameter_text)

  if (!UseIndVar) {
    parameterText = paste0(
      "\n\n",
      parameterText,
      "\n\nmodelingParameter=c(",
      pasteWithCollapse_WithApprox80charLimits(paste0(
        rv$ODEparameter, "=", rv$ODEparameter
      ), collapse = ","),
      ")\n\n"
    )
  }

  dataSet_text = paste0(
    "\n\n",ifelse(for_simulation,"simulationDataSkelton", "dataSet"),"<<-subset(data.frame(seq=seq(1,",

    dim(ll$ObservedData_dat)[1],
    "),
      time=c(",
    pasteWithCollapse_WithApprox80charLimits(ll$ObservedData_dat$time, collapse = ","),
    "),
      ID=c(\"",
    pasteWithCollapse_WithApprox80charLimits(ll$ObservedData_dat$ID, collapse = "\",\""),
    "\"),
    ",ifelse(for_simulation,"",
      paste0("Observed_value=c(",
             pasteWithCollapse_WithApprox80charLimits(ll$ObservedData_dat$Observed_value, collapse = ","),
    "),")),"
      Observation_expression=c(\"",
    pasteWithCollapse_WithApprox80charLimits(ll$ObservedData_dat$Observation_expression, collapse = "\",\""),
    "\"))",ifelse(!for_simulation,
",!is.na(Observed_value))",")\n\n
uniqueTime=unique(simulationDataSkelton$time)\n"
  ))

  dose_obs_Text = paste0("\n\nsimResult_df=data.frame()\n")

  for (ID_nu in sort(unique_IDs)) {
    dose_df = subset(ll$Dose_dat, ID == ID_nu)
    obs_df = subset(ll$ObservedData_dat, ID == ID_nu)
    init_df = subset(ll$InitialCondition_dat, ID == ID_nu)

    if (UseIndVar) {
      dose_obs_Text = paste0(dose_obs_Text, "\n\n## ID: ", ID_nu, "     ", paste( unique(obs_df$Memo), collapse = ", "),
                             "\n\n")


      dose_obs_Text = paste0(
        dose_obs_Text,
        ifelse(length(indParameterNames)>0,
        paste(
          paste0(indParameterNames, "=", indParameterNames, "_ID", ID_nu),
          collapse = "\n"
        ), ""),"\n",
        ifelse(length(indParameterNames_3)>0, paste(
          paste0(indParameterNames_3, "=", indParameterNames_3, "_ID", substr(ID_nu,1,3)),
          collapse = "\n"
        ),""),
        "\n\nmodelingParameter=c(",
        pasteWithCollapse_WithApprox80charLimits(
          paste0(rv$ODEparameter, "=", rv$ODEparameter),
          collapse = ","
        ),
        ")\n\n ev <- eventTable()"
      )

    } else{
      dose_obs_Text = paste0(dose_obs_Text, "\n\n ## ID: ", ID_nu, "     ", paste( unique(obs_df$Memo), collapse = ", "),"
      ev <- eventTable()")
    }

    if (dim(dose_df)[1] > 0) {
      for (i in seq(1, dim(dose_df)[1])) {
        if (is.na(dose_df$rate[i])) {
          dose_obs_Text = paste0(
            dose_obs_Text,
            "
      ev$add.dosing(dose = ",
            dose_df$dose[i],
            ", start.time=",
            dose_df$start.time[i],
            ", dosing.to=\"",
            dose_df$dosing.to[i],
            "\", nbr.doses=",
            dose_df$nbr.doses[i],
            ", dosing.interval=",
            dose_df$dosing.interval[i],
            ")"
          )
        } else{
          dose_obs_Text = paste0(
            dose_obs_Text,
            "
      ev$add.dosing(dose = ",
            dose_df$dose[i],
            ", rate = ",
            dose_df$rate[i],
            ", start.time=",
            dose_df$start.time[i],
            ", dosing.to=\"",
            dose_df$dosing.to[i],
            "\", nbr.doses=",
            dose_df$nbr.doses[i],
            ", dosing.interval=",
            dose_df$dosing.interval[i],
            ")"


          )
        }
      }
    }

    if (dim(obs_df)[1] > 0) {
      initsArgText=""
      if (dim(init_df)[1] > 0) {
        dose_obs_Text = paste0(
          dose_obs_Text,
          "\n      initCondition_vec=c(",
          paste(paste0(init_df$state, "=", init_df$value), collapse = ","),
          ")"
        )
        initsArgText=", inits=initCondition_vec"
      }

      dose_obs_Text =paste0(
        dose_obs_Text,
        "\n      ev$add.sampling(",ifelse(for_simulation,"uniqueTime",paste0("c(",
        paste(sort(as.numeric(
          unique(obs_df$time)
        )), collapse = ", "),
        ")")),
        ")\n      odeSol=compiledModel$solve(modelingParameter, ev", initsArgText, ")"
      )

      uniqueObsVariables = unique(obs_df$Observation_expression)

      if(!for_simulation){
        for (obsVariable_nu in uniqueObsVariables) {
          dose_obs_Text = paste0(
            dose_obs_Text,
            "\n      simResult_df=rbind(simResult_df,data.frame(value=with(data.frame(odeSol), ",
            obsVariable_nu ,
            "), Observation_expression=\"",
            obsVariable_nu,
            "\", time=odeSol[,\"time\"], ID=\"",
            ID_nu,
            "\"))"
          )
        }
      }

      if(for_simulation){
        dose_obs_Text=paste0(dose_obs_Text, paste0("\n\n      for(i in seq(2,dim(odeSol)[2])){
        simResult_df=rbind(simResult_df,data.frame(value=odeSol[,i], Observation_expression = colnames(odeSol)[i], time=odeSol[,1], ID=\"",ID_nu,"\"))
      }")
        )
      }
    }
  }
  dose_obs_Text = paste0(
    dose_obs_Text,
    "

      mergedData=merge(",ifelse(for_simulation,"simulationDataSkelton", "dataSet"),",simResult_df, all.x = TRUE)
      mergedData=mergedData[order(mergedData$seq),]
      return(mergedData$value)
}"
  )

  if(for_simulation){
    dose_obs_Text = paste0(
      dose_obs_Text,
      "\n\n
      simulationModel=simulation_function_",runName
    )
  }


  CGNM_run_text = paste0(
    paste0(
      "
ParaNames=c(\"",
      pasteWithCollapse_WithApprox80charLimits(ParameterName_vec, collapse = "\",\""),
      "\")
"
    ),

    paste0("UR=c(",
           paste(Initial_upper_range_vec, collapse = ","),
           ")
"),
    paste0("LR=c(",
           paste(Initial_lower_range_vec, collapse = ","),
           ")
"),
    paste0("U_bound=c(",
           paste(
             as.numeric(Upper_bound_vec), collapse = ","
           ),
           ")
"),
    paste0("L_bound=c(",
           paste(
             as.numeric(Lower_bound_vec), collapse = ","
           ),
           ")
"),
    paste0("observation=dataSet$Observed_value
")
  )

  if (ll$UseMO) {
    CGNM_run_text = paste0(
      CGNM_run_text,
      "MO_weights_vec=c(",
      paste(MO_weight_vec, collapse = ",") ,
      ")
MO_values_vec=c(",
      paste(MO_value_vec, collapse = ",") ,
      ")"
    )

  }

  testCode_text <-
    paste0(
      "library(CGNM)\nlibrary(rxode2)\n\n",ifelse(for_simulation,"","#If you wish to rerun CGNM regardless of if there is a log file available set it to TRUE\nForceReEstimation=FALSE\n\n"),
      odeCodeText,
      dataSet_text,
      parameterText,
      dose_obs_Text,
      ifelse(for_simulation,"",CGNM_run_text),
      "\n\n"
    )

  return(testCode_text)
}


makeCGNM_runCode = function(parallel = "none",
                            ll,
                            numIter = 25,
                            numMinimizersTofind = 250,
                            bootstrap = TRUE,
                            runName="") {

  runName=gsub("[[:punct:]]","", runName)
  runName=gsub(" ","_", runName)

  useResidualFunction = sum(ll$ObservedData_dat$ResidualError_model != 0) >
    0

  if (useResidualFunction) {
    CGNM_runOptions = "targetVector = rep(0,length(observation)), "
  } else{
    CGNM_runOptions = "targetVector = observation, "
  }

  if (ll$UseMO) {
    CGNM_runOptions = paste0(
      CGNM_runOptions,
      "initial_lowerRange = LR, initial_upperRange = UR, lowerBound = L_bound, upperBound = U_bound,ParameterNames = ParaNames, MO_weights=MO_weights_vec, MO_values=MO_values_vec, num_minimizersToFind = ",
      numMinimizersTofind,
      ", num_iteration=",
      numIter
    )

  } else{
    CGNM_runOptions = paste0(
      CGNM_runOptions,
      "initial_lowerRange = LR, initial_upperRange = UR, lowerBound = L_bound, upperBound = U_bound,ParameterNames = ParaNames, num_minimizersToFind = ",
      numMinimizersTofind,
      ", num_iteration=",
      numIter
    )
  }

  CGNM_runOptions = paste0(
    CGNM_runOptions,
    ",runName=\"",runName,"\""
  )


  out = ""
  if (useResidualFunction) {
    if (sum(ll$ObservedData_dat$ResidualError_model != 1) == 0) {
      out = paste0(
        out,
        "
residual_model=function(y_sim){
  out=(y_sim-dataSet$Observed_value)/y_sim

  return(out)
}"
      )
    } else{
      out = paste0(
        out,
        "
residual_model=function(y_sim){
  out=(y_sim-dataSet$Observed_value)/y_sim
  out[c(",
        paste(
          which(ll$ObservedData_dat$ResidualError_model == 0),
          collapse = ","
        )

        ,
        ")]=y_sim[c(",
        paste(
          which(ll$ObservedData_dat$ResidualError_model == 0),
          collapse = ","
        )

        ,
        ")]

  return(out)
}"
      )
    }
  }

  if (ll$UseMO) {
    out = paste0(
      out,
      "\n\n## CGNM R package above or equal to version 0.8.1 is necessary to run middle out method using MO_weights, MO_values options as implemented below."
    )
  }
  if(bootstrap){
    out=paste0(out, "\n\n", "if(!file.exists(\"",runName,"_CGNM_log_bootstrap/CGNM_bootstrapResult.RDATA\")|ForceReEstimation){"
    )
  }else{
    out=paste0(out, "\n\n", "if(!file.exists(\"",runName,"_CGNM_log/CGNM_result.RDATA\")|ForceReEstimation){"
    )
  }

  if (parallel == "none") {
    if (useResidualFunction) {
      out = paste0(
        out,
        "

model_function_withResidualmodel=function(x){
  return(residual_model(model_function_",runName,"(x)))
}

CGNM_result=Cluster_Gauss_Newton_method(model_function_withResidualmodel, ", CGNM_runOptions, ")
", ifelse(
  bootstrap,
  "CGNM_result=Cluster_Gauss_Newton_Bootstrap_method(CGNM_result,model_function_withResidualmodel)",
  ""
))

    } else{
      out = paste0(
        out,
        "

CGNM_result=Cluster_Gauss_Newton_method(model_function_",runName,", ",
        CGNM_runOptions,
        ")
",
        ifelse(
          bootstrap,
          paste0("CGNM_result=Cluster_Gauss_Newton_Bootstrap_method(CGNM_result,model_function_",runName,")"),
          ""
        )
      )
    }

  } else if (parallel == "win") {
    out = paste0(
      out,
      "

library(foreach)
library(doParallel)

numCoretoUse=detectCores()-1
registerDoParallel(numCoretoUse)
cluster=makeCluster(numCoretoUse, type =\"PSOCK\")
registerDoParallel(cl=cluster)

obsLength=length(observation)

# Given CGNM searches through wide range of parameter combination, it can encounter
# parameter combinations that is not feasible to evaluate. This try catch function
# is implemented within CGNM for regular functions but for the matrix functions
# user needs to implement outside of CGNM

modelFunction_tryCatch=function(x){
 out=tryCatch({",
      ifelse(
        useResidualFunction,
        paste0("residual_model(model_function_",runName,"(x))"),
        paste0("model_function_",runName,"(x)")
      ),
      "},
              error=function(cond) {rep(NA, obsLength)}
 )
 return(out)
}

model_matrix_function=function(x){
  X=as.matrix(x)

  if(is.matrix(X)){
    Y_list=foreach(i=1:dim(X)[1], .export = c(\"model_function_",runName,"\",",
ifelse(useResidualFunction, "\"residual_model\",", ""),
" \"modelFunction_tryCatch\", \"dataSet\", \"obsLength\", \"compiledModel\"), .packages = c(\"rxode2\"))%dopar%{
      modelFunction_tryCatch(as.numeric(X[i,]))
    }

    Y=t(matrix(unlist(Y_list),ncol=length(Y_list)))

  }else{

   Y= modelFunction_tryCatch(X)
  }

  return(Y)

}

CGNM_result=Cluster_Gauss_Newton_method(model_matrix_function, ",
CGNM_runOptions,
")
"
,
ifelse(
  bootstrap,
  "CGNM_result=Cluster_Gauss_Newton_Bootstrap_method(CGNM_result,model_matrix_function)",
  ""
),
"

stopCluster(cl=cluster)")

  } else if (parallel == "mac") {
    out = paste0(
      out,
      "

library(parallel)

    obsLength=length(observation)

    ## Given CGNM searches through wide range of parameter combination, it can encounter
    ## parameter combinations that is not feasible to evaluate. This try catch function
    ## is implemented within CGNM for regular functions but for the matrix functions
    ## user needs to implement outside of CGNM

    modelFunction_tryCatch=function(x){
      out=tryCatch({",
      ifelse(
        useResidualFunction,
        paste0("residual_model(model_function_",runName,"(x))"),
        paste0("model_function_",runName,"(x)")
      ),
      "},
                   error=function(cond) {rep(NA, obsLength)}
      )
      return(out)
    }

    model_matrix_function=function(x){
      Y_list=mclapply(split(x, rep(seq(1:nrow(x)),ncol(x))), modelFunction_tryCatch,mc.cores = (parallel::detectCores()-1), mc.preschedule = FALSE)

      Y=t(matrix(unlist(Y_list),ncol=length(Y_list)))

      return(Y)
    }

CGNM_result=Cluster_Gauss_Newton_method(model_matrix_function, ",
      CGNM_runOptions,
      ")
"
      ,
      ifelse(
        bootstrap,
        "CGNM_result=Cluster_Gauss_Newton_Bootstrap_method(CGNM_result,model_matrix_function)",
        ""
      )
    )

  }

if(bootstrap){
  out=paste0(out, "\n\n}else{\nload( file=\"",runName,"_CGNM_log_bootstrap/CGNM_bootstrapResult.RDATA\")\n}\n")
}else{
  out=paste0(out, "\n\n}else{\nload( file=\"",runName,"_CGNM_log/iteration_final.RDATA\")\n}\n")
}


if (ll$UseMO) {
  out = paste0(
    out,
    "


MO_para_names=CGNM_result$runSetting$ParameterNames[CGNM_result$runSetting$MO_weights!=0]
MO_values=CGNM_result$runSetting$MO_values[CGNM_result$runSetting$MO_weights!=0]

plot_goodnessOfFit(CGNM_result, independentVariableVector = c(dataSet$time, rep(0,length(MO_para_names))) ,dependentVariableTypeVector = c(paste(\"ID:\",dataSet$ID, dataSet$Observation_expression), MO_para_names) )+ggplot2::geom_point(colour=\"blue\")+ggplot2::labs(caption = \"Note the middleout values shown here are after transformation\")
plot_profileLikelihood(CGNM_result)+scale_x_continuous(trans=\"log10\")+ggplot2::geom_vline(data=data.frame(value=MO_values, parameterName=MO_para_names), aes(xintercept=value), colour=\"darkgrey\")
"
  )
} else{
  out = paste0(
    out,
    "\n\nplot_goodnessOfFit(CGNM_result, independentVariableVector = dataSet$time ,dependentVariableTypeVector = paste(\"ID:\",dataSet$ID, dataSet$Observation_expression))
plot_profileLikelihood(CGNM_result)+scale_x_continuous(trans=\"log10\")
"
  )
}


if (sum(ll$ObservedData_dat$ResidualError_model != 0) > 0) {
  out = paste0(
    out,
    "
plot_simulationWithCI(model_function_",runName,",parameter_matrix = CGNM_result$bootstrapParameterCombinations, independentVariableVector = dataSet$time, dependentVariableTypeVector =  paste( dataSet$Observation_expression,\"ID:\",dataSet$ID),
                      observationVector = dataSet$Observed_value, observationIndpendentVariableVector = dataSet$time, observationDependentVariableTypeVector =  paste( dataSet$Observation_expression,\"ID:\",dataSet$ID))+scale_y_continuous(trans=\"log10\")
"
  )

}
return(out)
}


simulation_code_text=function(runName){
  runName=gsub("[[:punct:]]","", runName)
  runName=gsub(" ","_", runName)

  paste0("
load(\"",runName,"_CGNM_log_bootstrap/CGNM_bootstrapResult.RDATA\")
plot_simulationWithCI(simulation_function",ifelse(runName=="","","_"),runName,",parameter_matrix = CGNM_result$bootstrapParameterCombinations, independentVariableVector = simulationDataSkelton$time, dependentVariableTypeVector =  paste(simulationDataSkelton$Observation_expression, \"ID:\",simulationDataSkelton$ID))+scale_y_continuous(trans=\"log10\")
")
}


# .makeCodeForPosthoc_middleout() is adapted from shinyCGNM's makeCodeForPosthoc_middleout()
# (inst/shinyCGNM/server.R), parameterized on parameterInfo_table instead of closing over the
# app's reactive ll$parameterInfo_dat. Two changes from the original, both needed to make the
# generated postHoc_likelihood() actually runnable without hand-editing first:
#  1. the "middle out value" variables are named middleOutValue_<ParameterName> from the start
#     (the original GUI code defined plain <ParameterName>=... at the top level but referenced
#     middleOutValue_<ParameterName> inside postHoc_likelihood(), which only worked once the
#     user had manually renamed the top-level variables to match);
#  2. postHoc_likelihood() now binds <ParameterName>=middleOutValue_<ParameterName> locally
#     before evaluating CGNM_result$runSetting$ParameterTransFormationDef[i] (a string like
#     "log(ka-0)"), since that expression is written in terms of the natural parameter name and
#     otherwise has nothing to resolve "ka" to (the original GUI code left this undefined,
#     which errors with "object 'ka' not found" the first time postHoc_likelihood() is actually
#     called, i.e. as soon as plot_profileLikelihood() is run).
.makeCodeForPosthoc_middleout=function(parameterInfo_table){
  ParameterName=parameterInfo_table$ParameterName
  nPara=length(ParameterName)

  paste0("## Use the code below to draw profile likelihood where the likelihood is defined post hoc. Works on CGNM version 0.7.0 or above.


## Change the values below to be >0 if wish to do middleout.  Set to 0 for the parameters not wishing to do the middle out.
",
         paste(paste0("weight_",ParameterName,"=0"), collapse = "\n"),"

## Change the values below to be the middle out values. (The values used below initially are the mean of the lower and upper range.)
",
         paste(paste0("middleOutValue_",ParameterName,"=", (as.numeric(parameterInfo_table$Initial_lower_range)+as.numeric(parameterInfo_table$Initial_upper_range))/2), collapse = "\n"),"

postHoc_likelihood=function(CGNM_result, initial=FALSE){

if(initial){
    out=CGNM_result$initialY[,!is.na(CGNM_result$runSetting$targetVector)]-t(matrix(rep(CGNM_result$runSetting$targetVector[!is.na(CGNM_result$runSetting$targetVector)],dim(CGNM_result$initialY)[1]),nrow=dim(CGNM_result$initialY)[2]))

    ",paste(paste0("x",seq(1,nPara),"=CGNM_result$initialX[,",seq(1,nPara),"]"),  collapse="\n    "),"

  }else{
    out=CGNM_result$Y[,!is.na(CGNM_result$runSetting$targetVector)]-t(matrix(rep(CGNM_result$runSetting$targetVector[!is.na(CGNM_result$runSetting$targetVector)],dim(CGNM_result$Y)[1]),nrow=dim(CGNM_result$Y)[2]))

    ",paste(paste0("x",seq(1,nPara),"=CGNM_result$X[,",seq(1,nPara),"]"),  collapse="\n    "),"

}

## these are the values x1, x2, ... are compared against (in transformed space) below;
## edit the middleOutValue_<ParameterName> variables above, not these bindings
",paste(paste0("    ",ParameterName,"=middleOutValue_",ParameterName), collapse="\n"),"

## If one wishes to do middle out in normal scale make changes to the following code by removing log10
out=cbind(out,
",
         paste0(paste0("    weight_",ParameterName,"*(x",seq(1,nPara),"-", "eval(parse(text=CGNM_result$runSetting$ParameterTransFormationDef[",seq(1,nPara),"])))"), collapse = ",\n"),
         "
)

return(rowSums((out)^2,na.rm = TRUE))
}
middleOutValue_df=data.frame(value=c(",paste(paste0("middleOutValue_", ParameterName), collapse=","),"),
                           parameterName=c(\"",paste(ParameterName, collapse="\",\""),"\"))

middleOutValue_df=middleOutValue_df[c(",paste(paste0("weight_", ParameterName), collapse="!=0,"),
"!=0),]

plot_profileLikelihood(CGNM_result,Likelihood_function = postHoc_likelihood)+geom_vline(data=middleOutValue_df,aes(xintercept=value))")
}


# ---------------------------------------------------------------------------
# Internal validation helpers for generateCGNM_script().
#
# These reproduce, outside of Shiny, the checks that inst/shinyCGNM/server.R
# performs on cell edits and CSV uploads for the ODE text and the three
# tables (parameterInfo, doseData, observedData).
# ---------------------------------------------------------------------------

.CGNM_compileODE=function(ODE_text){
  if(!is.character(ODE_text)||length(ODE_text)!=1||is.na(ODE_text)||trimws(ODE_text)==""){
    stop("ODE_text must be a single non-empty character string.")
  }
  if(!requireNamespace("rxode2", quietly = TRUE)){
    stop("The 'rxode2' package is required to compile ODE_text and is not installed.\nInstall it with: install.packages(\"rxode2\")")
  }

  compiledModel=tryCatch(rxode2::RxODE(ODE_text), error=function(e) e, warning=function(w) w)
  if(inherits(compiledModel, "condition")){
    stop(paste0("ODE_text could not be compiled by rxode2::RxODE(); check the ODE syntax.\nUnderlying message:\n", conditionMessage(compiledModel)))
  }

  modelVars=compiledModel$get.modelVars()
  list(stateVariables=modelVars$state, ODEparameters=modelVars$params, LHSvariables=modelVars$lhs)
}

.CGNM_noQuotes=function(x, label){
  bad=grepl('"', x, fixed=TRUE)
  if(any(bad)){
    stop(paste0(label, " cannot contain the character '\"'. Offending value(s): ", paste(unique(x[bad]), collapse=", ")))
  }
}

.CGNM_validate_doseTable=function(doseData_table, stateVariables, label="doseData_table"){
  requiredCols=c("ID","dose","dosing.to","start.time")
  optionalDefaults=list(rate=NA, nbr.doses=1, dosing.interval=NA)

  if(is.data.frame(doseData_table)){
    missingReq=setdiff(requiredCols, names(doseData_table))
    if(length(missingReq)>0){
      stop(paste0(label, " is missing required column(s): ", paste(missingReq, collapse=", "),
                   ". Required columns are ID, dose, dosing.to, start.time (rate, nbr.doses, and dosing.interval are optional)."))
    }
    out=doseData_table
  } else if(length(doseData_table)==1 && is.na(doseData_table)){
    out=data.frame(ID=character(0), dose=character(0), dosing.to=character(0), start.time=character(0))
  } else{
    stop(paste0(label, " must be a data.frame, or NA if the model has no dosing events."))
  }

  for(colNu in names(optionalDefaults)){
    if(!colNu %in% names(out)) out[[colNu]]=rep(optionalDefaults[[colNu]], dim(out)[1])
  }

  out$ID=as.character(out$ID)
  out$dosing.to=as.character(out$dosing.to)
  out$dose=as.character(out$dose)
  out$start.time=as.character(out$start.time)
  out$rate=as.character(out$rate)
  suppressWarnings(out$nbr.doses<-as.numeric(out$nbr.doses))
  suppressWarnings(out$dosing.interval<-as.numeric(out$dosing.interval))

  if(dim(out)[1]>0){
    if(any(is.na(out$dose)|out$dose=="")){
      stop(paste0(label, "$dose cannot be NA/blank; every dosing row needs a dose amount or expression."))
    }

    badTo = is.na(out$dosing.to) | out$dosing.to=="" | !(out$dosing.to %in% stateVariables)
    if(any(badTo)){
      stop(paste0(label, "$dosing.to must name a compartment (state variable) defined by ODE_text. Invalid entries: ",
                   paste(unique(out$dosing.to[badTo]), collapse=", "),
                   ". Compartments defined by ODE_text: ", paste(stateVariables, collapse=", ")))
    }

    .CGNM_noQuotes(out$ID, paste0(label, "$ID"))
    .CGNM_noQuotes(out$dosing.to, paste0(label, "$dosing.to"))
  }

  out[,c("ID","dose","dosing.to","start.time","rate","nbr.doses","dosing.interval")]
}

.CGNM_extractDoseParameters=function(doseData_table){
  validEntries=unique(c(doseData_table$dose, doseData_table$start.time, doseData_table$rate))
  validEntries=validEntries[!is.na(validEntries)]
  validEntries=validEntries[!(validEntries %in% c("", " ", "\n"))]
  if(length(validEntries)==0) return(character(0))

  parsedExpr=tryCatch(parse(text=paste(validEntries, collapse="+")), error=function(e) e)
  if(inherits(parsedExpr, "condition")){
    stop(paste0("Could not parse the dose/start.time/rate entries in doseData_table as R expressions.\nUnderlying message:\n", conditionMessage(parsedExpr)))
  }
  all.vars(parsedExpr)
}

.CGNM_validate_initialConditionTable=function(initialConditionData_table, stateVariables){
  requiredCols=c("ID","state","value")

  if(is.data.frame(initialConditionData_table)){
    missingReq=setdiff(requiredCols, names(initialConditionData_table))
    if(length(missingReq)>0){
      stop(paste0("initialConditionData_table is missing required column(s): ", paste(missingReq, collapse=", "),
                   ". Required columns are ID, state, value."))
    }
    out=initialConditionData_table[,requiredCols]
  } else if(length(initialConditionData_table)==1 && is.na(initialConditionData_table)){
    out=data.frame(ID=character(0), state=character(0), value=character(0))
  } else{
    stop("initialConditionData_table must be a data.frame, or NA if no initial conditions are set explicitly (all compartments then start at 0, subject to any dosing).")
  }

  out$ID=as.character(out$ID)
  out$state=as.character(out$state)
  out$value=as.character(out$value)

  if(dim(out)[1]>0){
    if(any(is.na(out$value)|out$value=="")){
      stop("initialConditionData_table$value cannot be NA/blank; every row needs an initial-condition amount or expression.")
    }

    badState = is.na(out$state) | out$state=="" | !(out$state %in% stateVariables)
    if(any(badState)){
      stop(paste0("initialConditionData_table$state must name a compartment (state variable) defined by ODE_text. Invalid entries: ",
                   paste(unique(out$state[badState]), collapse=", "),
                   ". Compartments defined by ODE_text: ", paste(stateVariables, collapse=", ")))
    }

    dupKey=paste(out$ID, out$state)
    if(any(duplicated(dupKey))){
      stop(paste0("initialConditionData_table has more than one row for the same ID/state combination: ",
                   paste(unique(dupKey[duplicated(dupKey)]), collapse=", ")))
    }

    .CGNM_noQuotes(out$ID, "initialConditionData_table$ID")
    .CGNM_noQuotes(out$state, "initialConditionData_table$state")
  }

  out
}

.CGNM_extractInitialConditionParameters=function(initialConditionData_table){
  validEntries=unique(initialConditionData_table$value)
  validEntries=validEntries[!is.na(validEntries)]
  validEntries=validEntries[!(validEntries %in% c("", " ", "\n"))]
  if(length(validEntries)==0) return(character(0))

  parsedExpr=tryCatch(parse(text=paste(validEntries, collapse="+")), error=function(e) e)
  if(inherits(parsedExpr, "condition")){
    stop(paste0("Could not parse the value entries in initialConditionData_table as R expressions.\nUnderlying message:\n", conditionMessage(parsedExpr)))
  }
  all.vars(parsedExpr)
}

.CGNM_validate_parameterInfoTable=function(parameterInfo_table, requiredParameterNames){
  if(!is.data.frame(parameterInfo_table)){
    stop("parameterInfo_table must be a data.frame.")
  }
  requiredCols=c("ParameterName","Initial_lower_range","Initial_upper_range")
  missingReq=setdiff(requiredCols, names(parameterInfo_table))
  if(length(missingReq)>0){
    stop(paste0("parameterInfo_table is missing required column(s): ", paste(missingReq, collapse=", ")))
  }

  out=parameterInfo_table
  out$ParameterName=as.character(out$ParameterName)

  if(any(is.na(out$ParameterName))||any(out$ParameterName=="")){
    stop("parameterInfo_table$ParameterName cannot contain NA or empty entries.")
  }
  .CGNM_noQuotes(out$ParameterName, "parameterInfo_table$ParameterName")
  if(any(duplicated(out$ParameterName))){
    stop(paste0("parameterInfo_table$ParameterName contains duplicated entries: ",
                 paste(unique(out$ParameterName[duplicated(out$ParameterName)]), collapse=", ")))
  }

  nParam=dim(out)[1]
  if(!"Lower_bound" %in% names(out)) out$Lower_bound=rep(0, nParam)
  if(!"Upper_bound" %in% names(out)) out$Upper_bound=rep(NA, nParam)
  if(!"VaryByID" %in% names(out)) out$VaryByID=rep(0, nParam)
  if(!"MO_weight" %in% names(out)) out$MO_weight=rep(0, nParam)
  if(!"MO_value" %in% names(out)) out$MO_value=rep(NA, nParam)
  if(!"Unit" %in% names(out)) out$Unit=rep(NA, nParam)

  requireNumeric=c("Initial_lower_range","Initial_upper_range","VaryByID","MO_weight")
  allowNANumeric=c("Lower_bound","Upper_bound","MO_value")

  for(colNu in requireNumeric){
    suppressWarnings(numericVersion<-as.numeric(out[[colNu]]))
    if(any(is.na(numericVersion))){
      stop(paste0("parameterInfo_table$", colNu, " must be numeric (NA not allowed) for every parameter. Problem entries for: ",
                   paste(out$ParameterName[is.na(numericVersion)], collapse=", ")))
    }
    out[[colNu]]=numericVersion
  }
  for(colNu in allowNANumeric){
    suppressWarnings(numericVersion<-as.numeric(out[[colNu]]))
    badIndex=is.na(numericVersion)&!is.na(out[[colNu]])
    if(any(badIndex)){
      stop(paste0("parameterInfo_table$", colNu, " must be numeric or NA. Problem entries for: ",
                   paste(out$ParameterName[badIndex], collapse=", ")))
    }
    out[[colNu]]=numericVersion
  }

  if(any(out$Initial_lower_range>out$Initial_upper_range)){
    bad=out$Initial_lower_range>out$Initial_upper_range
    stop(paste0("parameterInfo_table$Initial_lower_range must be <= Initial_upper_range. Problem entries for: ",
                 paste(out$ParameterName[bad], collapse=", ")))
  }
  lbBad=!is.na(out$Lower_bound)&(out$Lower_bound>out$Initial_lower_range)
  if(any(lbBad)){
    stop(paste0("parameterInfo_table$Lower_bound must be <= Initial_lower_range. Problem entries for: ",
                 paste(out$ParameterName[lbBad], collapse=", ")))
  }
  ubBad=!is.na(out$Upper_bound)&(out$Upper_bound<out$Initial_upper_range)
  if(any(ubBad)){
    stop(paste0("parameterInfo_table$Upper_bound must be >= Initial_upper_range. Problem entries for: ",
                 paste(out$ParameterName[ubBad], collapse=", ")))
  }
  moBad=out$MO_weight!=0&is.na(out$MO_value)
  if(any(moBad)){
    stop(paste0("parameterInfo_table$MO_value must be a non-NA number whenever the corresponding MO_weight is not 0. Problem entries for: ",
                 paste(out$ParameterName[moBad], collapse=", ")))
  }

  missingParams=setdiff(requiredParameterNames, out$ParameterName)
  if(length(missingParams)>0){
    stop(paste0("parameterInfo_table is missing row(s) for parameter(s) required by ODE_text, doseData_table, and/or initialConditionData_table: ",
                 paste(missingParams, collapse=", ")))
  }
  extraParams=setdiff(out$ParameterName, requiredParameterNames)
  if(length(extraParams)>0){
    message(paste0("Note: parameterInfo_table contains parameter(s) not referenced by ODE_text, doseData_table, or initialConditionData_table (kept as-is, unused in the generated model function): ",
                    paste(extraParams, collapse=", ")))
  }

  out[,c("ParameterName","Initial_lower_range","Initial_upper_range","Lower_bound","Upper_bound","VaryByID","MO_weight","MO_value","Unit")]
}

.CGNM_validate_observedDataTable=function(observedData_table){
  if(!is.data.frame(observedData_table)){
    stop("observedData_table must be a data.frame.")
  }
  requiredCols=c("ID","time","Observation_expression","Observed_value")
  missingReq=setdiff(requiredCols, names(observedData_table))
  if(length(missingReq)>0){
    stop(paste0("observedData_table is missing required column(s): ", paste(missingReq, collapse=", ")))
  }
  if(dim(observedData_table)[1]==0){
    stop("observedData_table must have at least one row.")
  }

  out=observedData_table
  if(!"ResidualError_model" %in% names(out)) out$ResidualError_model=0
  if(!"Memo" %in% names(out)) out$Memo=NA

  out$ID=as.character(out$ID)
  .CGNM_noQuotes(out$ID, "observedData_table$ID")

  suppressWarnings(out$time<-as.numeric(out$time))
  if(any(is.na(out$time))){
    stop("observedData_table$time must be numeric for every row.")
  }

  out$Observation_expression=as.character(out$Observation_expression)
  if(any(is.na(out$Observation_expression)|out$Observation_expression=="")){
    stop("observedData_table$Observation_expression cannot be NA/blank.")
  }
  .CGNM_noQuotes(out$Observation_expression, "observedData_table$Observation_expression")
  parseErrors=vapply(out$Observation_expression, function(x){
    inherits(tryCatch(parse(text=x), error=function(e) e), "error")
  }, logical(1))
  if(any(parseErrors)){
    stop(paste0("observedData_table$Observation_expression contains invalid R expression(s): ",
                 paste(unique(out$Observation_expression[parseErrors]), collapse=", ")))
  }

  badResidual=!(out$ResidualError_model %in% c(0,1))
  if(any(badResidual)){
    stop("observedData_table$ResidualError_model may only contain 0 (additive) or 1 (relative); currently only these two residual error models are implemented.")
  }

  evaluatedValue=suppressWarnings(vapply(out$Observed_value, function(v){
    if(is.numeric(v)) return(as.numeric(v))
    val=tryCatch(eval(parse(text=as.character(v))), error=function(e) NA_real_)
    if(!is.numeric(val)||length(val)!=1) return(NA_real_)
    as.numeric(val)
  }, numeric(1)))
  badValue=is.na(evaluatedValue)|!is.finite(evaluatedValue)
  if(any(badValue)){
    stop(paste0("observedData_table$Observed_value must be numeric, or an expression that evaluates to a single finite numeric value. Problem entries: ",
                 paste(out$Observed_value[badValue], collapse=", ")))
  }

  out[,c("ID","time","Observation_expression","Observed_value","ResidualError_model","Memo")]
}

.CGNM_validate_simulationTimepointsTable=function(simulationTimepoints_table, validExpressions){
  if(!is.data.frame(simulationTimepoints_table)){
    stop("simulationTimepoints_table must be a data.frame.")
  }
  requiredCols=c("ID","time","Observation_expression")
  missingReq=setdiff(requiredCols, names(simulationTimepoints_table))
  if(length(missingReq)>0){
    stop(paste0("simulationTimepoints_table is missing required column(s): ", paste(missingReq, collapse=", ")))
  }
  if(dim(simulationTimepoints_table)[1]==0){
    stop("simulationTimepoints_table must have at least one row.")
  }

  out=simulationTimepoints_table
  if(!"Memo" %in% names(out)) out$Memo=NA

  out$ID=as.character(out$ID)
  .CGNM_noQuotes(out$ID, "simulationTimepoints_table$ID")

  suppressWarnings(out$time<-as.numeric(out$time))
  if(any(is.na(out$time))){
    stop("simulationTimepoints_table$time must be numeric for every row.")
  }

  out$Observation_expression=as.character(out$Observation_expression)
  if(any(is.na(out$Observation_expression)|out$Observation_expression=="")){
    stop("simulationTimepoints_table$Observation_expression cannot be NA/blank.")
  }
  .CGNM_noQuotes(out$Observation_expression, "simulationTimepoints_table$Observation_expression")

  # Unlike observedData_table$Observation_expression (evaluated via with(data.frame(odeSol), .)
  # in the fit's model function), the simulation code dumps every compartment/LHS variable from
  # rxode2's solve() output as-is and merges by exact name -- so an arbitrary expression like
  # "log10(C_central)" silently merges to NA here instead of being evaluated.
  badExpr=!(out$Observation_expression %in% validExpressions)
  if(any(badExpr)){
    stop(paste0("simulationTimepoints_table$Observation_expression must exactly match one of ODE_text's compartment or derived (LHS) variable names -- it is not evaluated as an R expression the way observedData_table$Observation_expression is. Invalid entries: ",
                 paste(unique(out$Observation_expression[badExpr]), collapse=", "),
                 ". Valid names: ", paste(validExpressions, collapse=", ")))
  }

  out[,c("ID","time","Observation_expression","Memo")]
}


#' @title generateCGNM_script
#' @description
#' Programmatic (non-GUI) equivalent of the shinyCGNM app's "make CGNM script" workflow.
#' Takes the ODE model text and the same three tables the shinyCGNM app uses
#' (parameter info, observed data, and dose data), runs the same validation
#' checks the app performs on those tables (ODE compiles, dosing targets a
#' real compartment, every referenced parameter has a range, observed values
#' are numeric, residual error model is 0/1, ...), and returns the generated
#' R script (an rxode2-based model function plus a
#' \code{\link{Cluster_Gauss_Newton_method}}/\code{\link{Cluster_Gauss_Newton_Bootstrap_method}}
#' call) as a character string, ready to be written to a file, \code{cat()}'d,
#' or \code{eval(parse(text = .))}'d directly.
#' @param ODE_text (required input) \emph{string} the ODE model, written in rxode2 syntax (see \code{rxode2::RxODE()}). Its state variables become the compartments doseData_table$dosing.to and initialConditionData_table$state can target, and its parameters (together with any free symbols used in doseData_table's dose/start.time/rate columns and initialConditionData_table's value column) are the parameters parameterInfo_table must provide ranges for.
#' @param parameterInfo_table (required input) \emph{data.frame} one row per parameter. Required columns: ParameterName, Initial_lower_range, Initial_upper_range. Optional columns (defaulted if absent): Lower_bound (default 0), Upper_bound (default NA), VaryByID (default 0; 0=shared across IDs, 3=vary by the first 3 characters of ID, any other nonzero value=vary by the full ID), MO_weight (default 0), MO_value (default NA), Unit (default NA). See \code{\link{make_ShinyCGNM_parameterInfo}}.
#' @param observedData_table (required input) \emph{data.frame} one row per observation. Required columns: ID, time, Observation_expression, Observed_value. Optional columns (defaulted if absent): ResidualError_model (default 0; 0=additive, 1=relative), Memo (default NA). See \code{\link{make_ShinyCGNM_observationData}}.
#' @param doseData_table (default: NA) \emph{data.frame or NA} one row per dosing event. Required columns if supplied: ID, dose, dosing.to, start.time. Optional columns (defaulted if absent): rate (default NA, i.e. bolus dose), nbr.doses (default 1), dosing.interval (default NA). Set to NA (the default) if the model has no dosing events. See \code{\link{make_ShinyCGNM_doseData}}.
#' @param initialConditionData_table (default: NA) \emph{data.frame or NA} one row per (ID, compartment) whose initial condition (the value of that ODE state variable at time 0) should be set explicitly, instead of implicitly deriving it purely from dosing. Required columns: ID, state (a compartment name defined by ODE_text), value (a number, or an expression that may reference a parameter name, mirroring doseData_table's dose column). Set to NA (the default) if every compartment should simply start at 0. \code{doseData_table} and \code{initialConditionData_table} are independent and may be used together, separately, or not at all for a given ID: any dosing events are simulated on top of whatever initial condition is set (or 0, if none is set) for that compartment, exactly as \code{rxode2}'s own \code{inits} argument to \code{solve()} composes with an \code{eventTable}. See \code{\link{make_ShinyCGNM_initialCondition}}.
#' @param runName (default: "") \emph{string} passed through to \code{\link{Cluster_Gauss_Newton_method}}'s \code{runName}; also used to name the generated model function. Punctuation is stripped and spaces are replaced with underscores, matching shinyCGNM's behavior.
#' @param num_minimizersToFind (default: 250) \emph{positive integer} passed through as \code{\link{Cluster_Gauss_Newton_method}}'s \code{num_minimizersToFind}.
#' @param num_iteration (default: 25) \emph{positive integer} passed through as \code{\link{Cluster_Gauss_Newton_method}}'s \code{num_iteration}.
#' @param bootstrap (default: TRUE) \emph{logical} if TRUE the generated script also calls \code{\link{Cluster_Gauss_Newton_Bootstrap_method}} on the fit.
#' @param parallel (default: "none") \emph{"none", "win", or "mac"} if not "none", the generated model function is wrapped for parallel evaluation using doParallel (\code{"win"}) or parallel::mclapply (\code{"mac"}), matching the parallel computation options offered in shinyCGNM.
#' @param includeMiddleOutCode (default: FALSE) \emph{logical} if TRUE, appends a post-hoc middle-out code template (\code{postHoc_likelihood()} plus a \code{plot_profileLikelihood()} call using it) after the main script, matching shinyCGNM's "download post-hoc middle-out code" button. This is a separate mechanism from the built-in \code{MO_weight}/\code{MO_value} columns in \code{parameterInfo_table} (which bake the middle-out constraint into the CGNM search itself): this one applies a user-adjustable constraint to an *already-fit* \code{CGNM_result} instead, only once you edit the generated \code{weight_<ParameterName>} values from their default of 0.
#' @param includeSimulationCode (default: FALSE) \emph{logical} if TRUE, appends a second, simulation-only model function (built the same way as the main one, but skipping \code{Observed_value}) plus a \code{plot_simulationWithCI()} call, matching shinyCGNM's simulation tab. Requires \code{bootstrap = TRUE} (the simulation code needs \code{CGNM_result$bootstrapParameterCombinations} for the confidence band) and \code{simulationTimepoints_table}. Any rows of \code{initialConditionData_table} whose ID matches an ID used in the simulation still apply to it; there is no separate simulation-only initial condition table.
#' @param simulationDoseData_table (default: NA) \emph{data.frame or NA} the dosing regimen to simulate, in the same shape as \code{doseData_table} (and independently validated the same way) &mdash; only used when \code{includeSimulationCode = TRUE}. This can differ from \code{doseData_table} (e.g. a new hypothetical regimen); set to NA if the simulation has no dosing.
#' @param simulationTimepoints_table (default: NA) \emph{data.frame, required when includeSimulationCode = TRUE} the time points/variables to simulate. Required columns: ID, time, Observation_expression (no Observed_value or ResidualError_model, since nothing is being fit here). Unlike \code{observedData_table$Observation_expression}, this one is \strong{not} evaluated as an R expression: the simulation code dumps every compartment/derived (LHS) variable from the ODE and merges by exact name, so each entry here must exactly match one of \code{ODE_text}'s state or LHS variable names (e.g. \code{"C_central"}, not \code{"log10(C_central)"}). If \code{parameterInfo_table} has any individually-varying parameter (\code{VaryByID != 0}), every ID referenced here (and in \code{simulationDoseData_table}) must already be one of the IDs in \code{doseData_table}/\code{observedData_table}, since the bootstrap result only has fitted values for those IDs.
#' @param scriptFileName (default: NA) \emph{NA or string} if not NA, the generated script is additionally written to this file path with \code{writeLines()}.
#' @return \emph{string} the generated R script.
#' @examples
#' \dontrun{
#' ODE_text="
#' d/dt(depot) = -ka*depot
#' d/dt(central) = ka*depot - (CL/V1)*central
#' C_central = central/V1
#' "
#'
#' parameterInfo_table=make_ShinyCGNM_parameterInfo(
#'   ParameterName=c("ka","CL","V1"),
#'   Initial_lower_range=c(0.01,0.01,0.01),
#'   Initial_upper_range=c(100,100,100)
#' )
#'
#' doseData_table=make_ShinyCGNM_doseData(
#'   ID="1", dose=1000, dosing.to="depot", start.time=0, rate=NA,
#'   nbr.doses=1, dosing.interval=NA
#' )
#'
#' observedData_table=make_ShinyCGNM_observationData(
#'   ID="1",
#'   time=c(0.1, 0.2, 0.4, 0.6, 1, 2, 3, 6, 12),
#'   Observation_expression="log10(C_central)",
#'   Observed_value=log10(c(4.91, 8.65, 12.4, 18.7, 24.3, 24.5, 18.4, 4.66, 0.238)),
#'   ResidualError_model=0
#' )
#'
#' script_text=generateCGNM_script(
#'   ODE_text=ODE_text,
#'   parameterInfo_table=parameterInfo_table,
#'   observedData_table=observedData_table,
#'   doseData_table=doseData_table,
#'   runName="oral1cpt"
#' )
#'
#' cat(script_text)
#' eval(parse(text=script_text))
#'
#' ## Same model, but with an explicit initial condition instead of (or in addition
#' ## to) a dose: e.g. central starts at 5 instead of implicitly at 0.
#' initialConditionData_table=make_ShinyCGNM_initialCondition(
#'   ID="1", state="central", value=5
#' )
#' script_text_withInit=generateCGNM_script(
#'   ODE_text=ODE_text,
#'   parameterInfo_table=parameterInfo_table,
#'   observedData_table=observedData_table,
#'   doseData_table=doseData_table,
#'   initialConditionData_table=initialConditionData_table,
#'   runName="oral1cpt_withInit"
#' )
#'
#' ## Same fit, plus a post-hoc middle-out template and a simulation of a denser time
#' ## grid at a higher dose (bootstrap = TRUE is required for the simulation part).
#' ## note: unlike observedData_table, Observation_expression here must be an exact
#' ## compartment/derived-variable name from ODE_text (e.g. "C_central"), not an
#' ## expression like "log10(C_central)" - it is not evaluated, only matched by name.
#' simulationTimepoints_table=make_ShinyCGNM_simulationTimepoints(
#'   ID="1", time=seq(0.1,12,by=0.1), Observation_expression="C_central"
#' )
#' simulationDoseData_table=make_ShinyCGNM_doseData(
#'   ID="1", dose=2000, dosing.to="depot", start.time=0, rate=NA,
#'   nbr.doses=1, dosing.interval=NA
#' )
#' script_text_full=generateCGNM_script(
#'   ODE_text=ODE_text,
#'   parameterInfo_table=parameterInfo_table,
#'   observedData_table=observedData_table,
#'   doseData_table=doseData_table,
#'   runName="oral1cpt_full",
#'   bootstrap=TRUE,
#'   includeMiddleOutCode=TRUE,
#'   includeSimulationCode=TRUE,
#'   simulationDoseData_table=simulationDoseData_table,
#'   simulationTimepoints_table=simulationTimepoints_table
#' )
#' }
#' @export
generateCGNM_script=function(ODE_text,
                              parameterInfo_table,
                              observedData_table,
                              doseData_table=NA,
                              initialConditionData_table=NA,
                              runName="",
                              num_minimizersToFind=250,
                              num_iteration=25,
                              bootstrap=TRUE,
                              parallel="none",
                              includeMiddleOutCode=FALSE,
                              includeSimulationCode=FALSE,
                              simulationDoseData_table=NA,
                              simulationTimepoints_table=NA,
                              scriptFileName=NA){

  if(includeSimulationCode && !bootstrap){
    stop("includeSimulationCode = TRUE requires bootstrap = TRUE: the generated simulation code calls plot_simulationWithCI() using CGNM_result$bootstrapParameterCombinations, which only exists after a bootstrap run.")
  }
  if(includeSimulationCode && (length(simulationTimepoints_table)==1 && is.na(simulationTimepoints_table)[1])){
    stop("includeSimulationCode = TRUE requires simulationTimepoints_table (columns: ID, time, Observation_expression, [Memo]).")
  }

  parallel=match.arg(parallel, c("none","win","mac"))

  odeInfo=.CGNM_compileODE(ODE_text)

  doseData_table=.CGNM_validate_doseTable(doseData_table, odeInfo$stateVariables)
  DOSEparameters=.CGNM_extractDoseParameters(doseData_table)

  initialConditionData_table=.CGNM_validate_initialConditionTable(initialConditionData_table, odeInfo$stateVariables)
  InitialConditionParameters=.CGNM_extractInitialConditionParameters(initialConditionData_table)

  requiredParameterNames=unique(c(odeInfo$ODEparameters, DOSEparameters, InitialConditionParameters))
  parameterInfo_table=.CGNM_validate_parameterInfoTable(parameterInfo_table, requiredParameterNames)

  observedData_table=.CGNM_validate_observedDataTable(observedData_table)

  noObservationIDs=setdiff(unique(c(doseData_table$ID, initialConditionData_table$ID)), unique(observedData_table$ID))
  if(length(noObservationIDs)>0){
    message(paste0("Note: the following ID(s) appear in doseData_table/initialConditionData_table but have no rows in observedData_table, so no simulation output will be generated for them: ",
                    paste(noObservationIDs, collapse=", ")))
  }

  if(includeSimulationCode){
    simulationDoseData_table=.CGNM_validate_doseTable(simulationDoseData_table, odeInfo$stateVariables, label="simulationDoseData_table")
    simulationTimepoints_table=.CGNM_validate_simulationTimepointsTable(simulationTimepoints_table, c(odeInfo$stateVariables, odeInfo$LHSvariables))

    simulationDoseParameters=.CGNM_extractDoseParameters(simulationDoseData_table)
    missingSimParams=setdiff(simulationDoseParameters, parameterInfo_table$ParameterName)
    if(length(missingSimParams)>0){
      stop(paste0("parameterInfo_table is missing row(s) for parameter(s) referenced by simulationDoseData_table: ",
                   paste(missingSimParams, collapse=", ")))
    }

    if(any(parameterInfo_table$VaryByID!=0)){
      fitIDs=unique(c(doseData_table$ID, observedData_table$ID))
      simulationIDs=unique(c(simulationDoseData_table$ID, simulationTimepoints_table$ID))
      newSimulationIDs=setdiff(simulationIDs, fitIDs)
      if(length(newSimulationIDs)>0){
        stop(paste0("parameterInfo_table has individually-varying parameter(s) (VaryByID != 0), but simulationDoseData_table/simulationTimepoints_table reference ID(s) that are not part of the original fit (doseData_table/observedData_table): ",
                     paste(newSimulationIDs, collapse=", "),
                     ". The bootstrap result has no fitted value for a new ID's individually-varying parameter(s), so simulation can only be run for IDs already in the fit. Either restrict the simulation to those IDs, or make the relevant parameter(s) shared (VaryByID = 0)."))
      }
    }
  }

  runName=gsub("[[:punct:]]","", runName)
  runName=gsub(" ","_", runName)

  UseMO=any(parameterInfo_table$MO_weight!=0)

  ll=list(Dose_dat=doseData_table, InitialCondition_dat=initialConditionData_table, parameterInfo_dat=parameterInfo_table, ObservedData_dat=observedData_table, UseMO=UseMO)
  input=list(ODE_text=ODE_text, runNameText=runName)
  rv=list(ODEparameter=odeInfo$ODEparameters, ODEvariable=odeInfo$stateVariables, LHSvariable=odeInfo$LHSvariables, DOSEparameter=DOSEparameters)

  modelCode=puttogether_model_code(input, ll, rv, for_simulation=FALSE)
  runCode=makeCGNM_runCode(parallel=parallel, ll=ll, numIter=num_iteration, numMinimizersTofind=num_minimizersToFind, bootstrap=bootstrap, runName=runName)

  script_text=paste0(modelCode, runCode)

  if(includeMiddleOutCode){
    script_text=paste0(script_text, "\n\n", .makeCodeForPosthoc_middleout(parameterInfo_table))
  }

  if(includeSimulationCode){
    sim_ll=list(Dose_dat=simulationDoseData_table, InitialCondition_dat=initialConditionData_table, parameterInfo_dat=parameterInfo_table, ObservedData_dat=simulationTimepoints_table, UseMO=UseMO)
    simulationModelCode=puttogether_model_code(input, sim_ll, rv, for_simulation=TRUE)
    script_text=paste0(script_text, "\n\n", simulationModelCode, "\n\n", simulation_code_text(runName))
  }

  if(!(length(scriptFileName)==1 && is.na(scriptFileName))){
    writeLines(script_text, scriptFileName)
  }

  return(script_text)
}
