library(shiny)
library(shinyjs)
library(DT)
library(dplyr)
library(shinydashboard)
library(rxode2)
library(mlr3misc)
library(shinybusy)
library(spsComps)
library(CGNM)


js <- '
$(document).on("keyup", function(e) {
  if(e.keyCode == 13){
    Shiny.onInputChange("keyPressed", Math.random());
  }
});
'

callback <- c(
  "var tbl = $(table.table().node());",
  "var id = tbl.closest('.datatables').attr('id');",
  "table.on('autoFill', function(e, datatable, cells){",
  "  var out = Array(cells.length);",
  "  for(var i=0; i<cells.length; ++i){",
  "    var c = cells[i][0];",
  "    var value = c.set === null ? '' : c.set;",
  # null causes problem in R
  "    out[i] = {row: c.index.row+1, col: c.index.column, value: value};",
  # if you want a color for the autofilled cells:
  #"    $(table.cell(c.index.row, c.index.column).node())",
  #"      .css('background-color', 'yellow');",
  "  }",
  "  Shiny.setInputValue(id + '_cells_filled:DT.cellInfo', out);",
  "  table.rows().invalidate();",
  # this updates the column type
  "});"
)


dt_output = function(title, id) {
  fluidRow(column(
    12, h1(paste0('Table ', sub('.*?([0-9]+)$', '\\1', id), ': ', title)),
    hr(), DTOutput(id)
  ))
}
render_dt = function(data, editable = 'cell', server = TRUE, ...) {
  renderDT(data, selection = 'none', server = server, editable = editable, ...)
}


ui <- function(req) {
  dashboardPage(
    dashboardHeader(title = "Shiny CGNM"),
    dashboardSidebar(
      sidebarMenu(
        menuItem(
          "Step1: specify ODE",
          tabName = "ODE_specify",
          icon = icon("glass-whiskey")
        ),
        menuItem(
          "Step2: give dosing information",
          tabName = "DoseInfo",
          icon = icon("glass-whiskey")
        ),
        menuItem(
          "Step3: specify parameter info",
          tabName = "ParameterInfo",
          icon = icon("glass-whiskey")
        ),
        menuItem(
          "Step4: give Observed data",
          tabName = "ObservedData",
          icon = icon("glass-whiskey")
        ),
        menuItem(
          "Step5: generate code",
          tabName = "generateCODE",
          icon = icon("glass-whiskey")
        )
      ),
      wellPanel(
        style = "background: darkgrey",
        h4("Model Variables"),
        textOutput("ODE_lhs"),
        h4("Parameters"),
        textOutput("ODE_para"),
        textOutput("DOSE_para")
      ),
      wellPanel(
        style = "background: grey",
        h4("Load and Save"),
        h6(
          "All GUI inputs can be saved and use later by saving and loading .rds file"
        ),
        fileInput(
          "openSavedStatus",
          label = "Load saved .rds file",
          placeholder = ".rds file"
        ),
        downloadButton(
          "saveCurrentStatus",
          label = "Save model as .rds file",
          icon = shiny::icon("download"),
          style = "color: black; background-color: white; border-color: grey"
        )
      ),
      h6(
        "To see an example model implemented using this GUI you can press the button below. (PBPK model for the DDI between pitavastatin and rifampicin, as published in Yoshikado et al.)"
      ),
      actionButton("loadExample_button", "Load example model")



    ),
    dashboardBody(tabItems(
      tabItem(
        tabName = "ODE_specify",
        add_busy_spinner(spin = "fading-circle"),
        tags$script(js),
        textAreaInput(
          "ODE_text",
          "Specify the system of ODEs (follow rxODE syntax)",
          width = "100%",
          height = "500",
          resize = "vertical"
        ),
        actionButton("compileODE", "Compile ODEs"),
        htmlOutput("ODE_compileMessage"),
      ),

      tabItem(
        tabName = "DoseInfo",
        add_busy_spinner(spin = "fading-circle"),
        numericInput(
          "numDoses",
          "Number of dose events",
          1,
          min = 1,
          step = 1
        ),
      #  dt_output('client-side processing (editable = "cell")', 'dose_dt')
        DTOutput("dose_dt"),
       verbatimTextOutput("table"),
      fileInput("dose_csv", "You can import dose data from a csv file by uploarding it below:"),
      h6("Dose csv file should contain the following columns: ID, dose, dosing.to, start.time and rate."),
      downloadButton('download_dose_csv',"Download dosing information as csv file"),


      ),

      tabItem(
        tabName = "ParameterInfo",
        add_busy_spinner(spin = "fading-circle"),
        DTOutput("parameterInfo_dt"),
        h4(
          "Lower_bound and Upper_bound are optional so can be left empty.  All table entries needs to be numerical values."
        ),
        fileInput("para_csv", "You can import parameter information from a csv file by uploarding it below:"),
        h6("Parameter information csv file should contain the following columns: ParameterName, Initial_lower_range, Initial_upper_range, Lower_bound, Upper_bound, and Unit. Only the rows with relevant ParameterName to the current model will be read in."),
        downloadButton('download_para_csv',"Download parameter information as csv file"),
        #   verbatimTextOutput("parameterInfo_table"),

      ),

      tabItem(
        tabName = "ObservedData",
        add_busy_spinner(spin = "fading-circle"),
        numericInput(
          "numObservations",
          "Number of observations",
          1,
          min = 1,
          step = 1
        ),
        DTOutput("ObservedData_dt"),
        fileInput("obs_csv", "You can import observation data from a csv file by uploarding it below:"),
        h6("Observation csv file should contain the following columns: ID, time, Observation_expression, Observed_value, and Memo."),
        downloadButton('download_obs_csv',"Download observation as csv file"),


        #   verbatimTextOutput("ObservedData_table"),

      ),

      tabItem(
        tabName = "generateCODE",
        add_busy_spinner(spin = "fading-circle"),
        shinyjs::useShinyjs(),
        h3("Conduct test according to the order below."),
        h5("If preceeding test fails most likely the following test will crash the GUI, so if the error is detected fix it before proceeding to the next test."),
        h4("Test 1: check the code can be put together, also check if the code make sense"),
        h6("tip: focusing to see if any unexpected NA values used in the code or not"),
        actionButton("makeCode_button", "Make model code"),
        verbatimTextOutput("test1_output"),
        h4("Test 2: check if your model can be evaluated, if error occurs, fix model code."),
        h6("tip: check online resources on rxode2. If model cannot be evaluated it is a problem independent of CGNM."),
        actionButton("testModelCode_button", "Test model code"),
        verbatimTextOutput("testModelCode_output"),
        h4("Test 3: check if CGNM runs, if error occurs, check all the observation values and initial lower and upper bounds of the parameters."),
        actionButton("testCGNMCode_button", "Test CGNM code"),
        h6("tip: if model evaluation test (Test 2) runs without an error and cannot spot error for CGNM inputs yourself, you may contact the author of CGNM with rds file and generated R-script which can be obtained by pressing the button below."),
        htmlOutput("testCode_output"),
        verbatimTextOutput("testCode_output2"),
        downloadButton("download_CGNMcode_button", label = "Download R script to run CGNM"),
        downloadButton("download_CGNMcode_parallel_mac_button", label = "Download R script to run CGNM in parallel in MAC"),
        downloadButton("download_CGNMcode_parallel_win_button", label = "Download R script to run CGNM in parallel in Windows"),
        downloadButton("download_posthocMiddleout_button", label = "Download R script for posthoc middle out"),

      )

    ))
  )
}

server <- function(input, output, session) {
  ## variable initialization----


  rv<-reactiveValues(ODEparameter=c(),
                    ODEvariable = c(),
                    LHSvariable = c(),
                    DOSEparameter=c(),
                  #  DOSEdata = data.frame(),
                  #  parameterNameVector = c(),
                  #  testCode_text = "",)
  )
  ll<-list(Dose_dat =
                      data.frame(
                        ID = rep(1, 1),
                        dose = NA,
                        dosing.to = NA,
                        start.time = 0,
                        rate = NA
                      ),
                    parameterInfo_dat =
                      data.frame(
                        ParameterName = c(),
                        Initial_lower_range = c(),
                        Initial_upper_range = c(),
                        Lower_bound = c(),
                        Upper_bound = c(),
                        Unit = c()
                      ),
                    ObservedData_dat =
                      data.frame(
                        ID = c(1),
                        time = NA ,
                        Observation_expression = NA,
                        Observed_value = NA,
                        Memo = NA
                      )
                    )
  #ODEparameter = reactiveVal(c())

#  DOSEparameter = reactiveVal(c())



  #odeText <- reactiveVal()

  ## TAB: ODE----

  ### when the compile ODE button being pressed----

  compileODE_func = function(odeText_in) {
    compile_message = "ODE could not be compiled correctly, check the text expression above."
    compile_message = capture.output(withCallingHandlers(
      tryCatch({
        RxODE(odeText_in)#input$ODE_text)
      }, error = function(e) {
        err <<- conditionMessage(e)
      }),
      warning = function(w) {
        warn <<- conditionMessage(w)
      }
    ))


    compileSuccess = tryCatch({
      model = RxODE(odeText_in)
      TRUE
    }, error = function(e) {
      FALSE
      FALSE
    }, warning = function(w) {
      FALSE
    })

    output$ODE_compileMessage <-
      renderUI({
        HTML(paste(unlist(
          strsplit(compile_message, split = "\n", fixed = TRUE)
        ), collapse = '<br/>'))
      })

    if (compileSuccess) {
      rv$ODEvariable <- model$get.modelVars()$state
      rv$LHSvariable <- model$get.modelVars()$lhs


      output$ODE_lhs <-
        renderText(paste(c(rv$ODEvariable, rv$LHSvariable), collapse = ", "))
      rv$ODEparameter=model$get.modelVars()$params
      output$ODE_para <-
        renderText(paste(rv$ODEparameter, collapse = ", "))
    } else{
      rv$ODEvariable <- c()
      rv$LHSvariable <<- c()
      output$ODE_lhs <- renderText("")
      rv$ODEparameter=c()
      output$ODE_para <- renderText("")
    }
  }

  observeEvent(input$compileODE, {
    compileODE_func(input$ODE_text)
  })


  ### OE key pressed for ODE text----
  observeEvent(input[["keyPressed"]], {
    compile_message = tryCatch({
      is.expression(parse(text = input$ODE_text))
      ""
    }, error = function(e) {
      unlist(strsplit(paste(e), split = ": <text>:", fixed = TRUE))[2]
    }, warning = function(w) {
      unlist(strsplit(paste(w), split = ": <text>:", fixed = TRUE))[2]
    })
    output$ODE_compileMessage <- renderText(compile_message)
  })

  ## TAB: Dose info----

  ### OE: number of doses spinbox changed ----

  observeEvent(input$numDoses, {

    if(!is.na(input$numDoses)){
      if(input$numDoses>0){
      if (dim(ll$Dose_dat)[1] > input$numDoses) {
        ll$Dose_dat <<- ll$Dose_dat[seq(1, input$numDoses),]
      }else if (dim(ll$Dose_dat)[1] < input$numDoses) {
        ll$Dose_dat <-
          rbind(ll$Dose_dat,
                data.frame(
                  ID = rep(1, input$numDoses - dim(ll$Dose_dat)[1]),
                  dose = NA,
                  dosing.to = NA,
                  start.time = NA,
                  rate = NA
                ))
      }
      output[["dose_dt"]] <- renderDT({
        datatable(
          ll$Dose_dat,
          editable = list(target = "cell"),
          selection = "none",
          extensions = "AutoFill",
          callback = JS(callback),
          options = list(lengthMenu = list(c(100, -1), c('100','all')), autoFill =  list(horizontal = FALSE))
        )
      }, server = FALSE)
    }}
  })


  ### OE: Dose info data table changed----

  observeEvent(input$dose_dt_cell_edit, {

    if(input$dose_dt_cell_edit$col ==3& !input$dose_dt_cell_edit$value%in% rv$ODEvariable &input$dose_dt_cell_edit$value!=""){
      shinyCatch(stop(paste(
        input$dose_dt_cell_edit$value, "is not a compartment name"
               )))
    }

    ll$Dose_dat <<- editData(ll$Dose_dat, input$dose_dt_cell_edit)

  })

  observeEvent(input$dose_dt_cells_filled, {

    ll$Dose_dat <<- editData(ll$Dose_dat, input$dose_dt_cells_filled)

  })


  observe({
          validEntriesDose = unique(c(ll$Dose_dat$dose, ll$Dose_dat$start.time, ll$Dose_dat$rate))
          shinyCatch(rv$DOSEparameter<<-extract_vars(as.formula(
            paste(
              "testExpression~",
              paste(validEntriesDose, collapse = "+")
            )
          ))$rhs)

          output$DOSE_para = renderText(paste(rv$DOSEparameter, collapse = ", "))
  }
  )
  #
  # Dose_data <- reactive({
  #   info <- rbind(input[["dose_dt_cells_filled"]], input[["dose_dt_cell_edit"]])
  #   if (!is.null(info)) {
  #     info <- unique(info)
  #     info$value[info$value == ""] <- NA
  #     ll$Dose_dat <- editData(ll$Dose_dat, info)
  #
  #     if (info$col == 3) {
  #       if (!info$value %in% c(rv$ODEvariable)) {
  #         shinyCatch(stop(paste(
  #           info$value, "is not a compartment name"
  #         )))
  #       }
  #     }
  #     if (info$col %in% c(2, 4, 5)) {
  #       validEntriesDose = unique(c(ll$Dose_dat$dose, ll$Dose_dat$start.time, ll$Dose_dat$rate))
  #       shinyCatch(rv$DOSEparameter<-extract_vars(as.formula(
  #         paste(
  #           "testExpression~",
  #           paste(validEntriesDose, collapse = "+")
  #         )
  #       ))$rhs)
  #
  #       output$DOSE_para = renderText(paste(rv$DOSEparameter, collapse = ", "))
  #     }
  #   }
  #
  #   output[["table"]] <- renderPrint({
  #     ll$Dose_dat
  #   })
  #   ll$Dose_dat#[(!is.na(Dose_dat$dose)) &
  #   # (Dose_dat$dosing.to %in% ODEvariable), ]
  #
  #
  # })

  #



  ## TAB: parameter info----

  ### RE: parameterInfo_dt cell changed----
  observeEvent(input$parameterInfo_dt_cell_edit, {

    if(input$parameterInfo_dt_cell_edit$col %in% seq(2,4) & is.na(as.numeric(input$parameterInfo_dt_cell_edit$value)) &input$parameterInfo_dt_cell_edit$value!=""){
      shinyCatch(stop(paste(
        input$parameterInfo_dt_cell_edit$value, "is not a number"
      )))
    }

    ll$parameterInfo_dat <<- editData(ll$parameterInfo_dat, input$parameterInfo_dt_cell_edit)

  })

  observeEvent(input$parameterInfo_dt_cells_filled, {

    ll$parameterInfo_dat <<- editData(ll$parameterInfo_dat, input$parameterInfo_dt_cells_filled)

  })

#
#
#   parameterInfo_Data <- reactive({
#     info <-
#       rbind(input[["parameterInfo_dt_cells_filled"]], input[["parameterInfo_dt_cell_edit"]])
#     if (!is.null(info)) {
#       info <- unique(info)
#       info$value[info$value == ""] <- NA
#       ll$parameterInfo_dat <<- editData(ll$parameterInfo_dat, info)
#       if (info$value != "") {
#         if (is.na(as.numeric(info$value))) {
#           shinyCatch(stop(
#             paste(
#               info$value,
#               ":not acceptable input, it needs to be a number."
#             )
#           ))
#         }
#       }
#     }
#
#     ll$parameterInfo_dat
#   })

  ### OBS: when ODEparameter and DOSEparameter changed
  observe({
#    rv$parameterNameVector
    temp_PNV <- (unique(c(rv$ODEparameter, rv$DOSEparameter)))

    if (length(temp_PNV) > 0) {
      if (dim(ll$parameterInfo_dat)[1] == 0) {
        ll$parameterInfo_dat <<-
          data.frame(
            ParameterName = temp_PNV,
            Initial_lower_range = NA,
            Initial_upper_range = NA,
            Lower_bound = 0,
            Upper_bound = NA,
            Unit = NA
          )
      } else{
        ll$parameterInfo_dat <<-
          ll$parameterInfo_dat[ll$parameterInfo_dat$ParameterName %in% temp_PNV,]
        if (dim(ll$parameterInfo_dat)[1] < length(temp_PNV)) {
          ll$parameterInfo_dat <<-
            rbind(
              ll$parameterInfo_dat,
              data.frame(
                ParameterName = temp_PNV[!temp_PNV %in% ll$parameterInfo_dat$ParameterName],
                Initial_lower_range = NA,
                Initial_upper_range = NA,
                Lower_bound = 0,
                Upper_bound = NA,
                Unit = NA
              )
            )
        }
      }

      output[["parameterInfo_dt"]] <- renderDT({
        datatable(
          ll$parameterInfo_dat,
          editable = list(target = "cell", disable = list(columns = c(1))),
          selection = "none",
          extensions = "AutoFill",
          callback = JS(callback),
          options = list(lengthMenu = list(c(100, -1), c('100','all')), autoFill =  list(horizontal = FALSE))
        )
      }, server = FALSE)
    }
  })

  # output[["parameterInfo_table"]] <-
  #   renderPrint({
  #     parameterInfo_Data()
  #   })



  ##TAB: observation----

checkIfEvaluatable=function(textIn){
  shinyCatch({
    eval(parse(text = textIn))

  }, shiny = TRUE)

  a=tryCatch(eval(parse(text = textIn)), error = function(e) NA)

  return(a)

}
  ### REACT: when the ObservedData_dt cell changed----
  observeEvent(input$ObservedData_dt_cell_edit, {
    if(input$ObservedData_dt_cell_edit$col == 4 & input$ObservedData_dt_cell_edit$value!=""){

      print(checkIfEvaluatable(input$ObservedData_dt_cell_edit$value))

    }

    ll$ObservedData_dat <<- editData(ll$ObservedData_dat, input$ObservedData_dt_cell_edit)

  })

  observeEvent(input$ObservedData_dt_cells_filled, {

    ll$ObservedData_dat <<- editData(ll$ObservedData_dat, input$ObservedData_dt_cells_filled)

  })


  ### OE: number of observation spinbox----
  observeEvent(input$numObservations, {
    if (dim(ll$ObservedData_dat)[1] > input$numObservations) {
      ll$ObservedData_dat <<-
        ll$ObservedData_dat[seq(1, input$numObservations),]
    } else if (dim(ll$ObservedData_dat)[1] < input$numObservations) {
      ll$ObservedData_dat <<-
        rbind(
          ll$ObservedData_dat,
          data.frame(
            ID = rep(
              max(ll$ObservedData_dat$ID),
              input$numObservations - dim(ll$ObservedData_dat)[1]
            ),
            time = NA ,
            Observation_expression = NA,
            Observed_value = NA,
            Memo = NA
          )
        )
    }
    output[["ObservedData_dt"]] <- renderDT({
      datatable(
        ll$ObservedData_dat,
        editable = list(target = "cell"),
        selection = "none",
        extensions = "AutoFill",
        callback = JS(callback),
        options = list(lengthMenu = list(c(100, -1), c('100','all')), autoFill =  list(horizontal = FALSE))
      )
    }, server = FALSE)
  })


  ##TAB: code generation----

  makeCodeForPosthoc_middleout=function(){

    paste0("## Use the code below to draw profile likelihood where the likelihood is definied post hoc. Works on the CGNM version 0.7.0 or above.


## Change the values below to be >0 if wish to do middleout.  Set to 0 for the parameters not wishing to do the middle out.
",
      paste(paste0("weight_",ll$parameterInfo_dat$ParameterName,"=0"), collapse = "\n"),"

## Change the values below to be the middle out values. (The values used below initially are the mean of the lower and upper range.)
",
      paste(paste0("middleOutValue_",ll$parameterInfo_dat$ParameterName,"=", (as.numeric(ll$parameterInfo_dat$Initial_lower_range)+as.numeric(ll$parameterInfo_dat$Initial_upper_range))/2), collapse = "\n"),"

postHoc_likelihood=function(CGNM_result, initial=FALSE){

if(initial){
    out=CGNM_result$initialY[,!is.na(CGNM_result$runSetting$targetVector)]-t(matrix(rep(CGNM_result$runSetting$targetVector[!is.na(CGNM_result$runSetting$targetVector)],dim(CGNM_result$initialY)[1]),nrow=dim(CGNM_result$initialY)[2]))

    ",paste(paste0("x",seq(1,length(ll$parameterInfo_dat$ParameterName)),"=CGNM_result$initialX[,",seq(1,length(ll$parameterInfo_dat$ParameterName)),"]"),  collapse="\n    "),"

  }else{
    out=CGNM_result$Y[,!is.na(CGNM_result$runSetting$targetVector)]-t(matrix(rep(CGNM_result$runSetting$targetVector[!is.na(CGNM_result$runSetting$targetVector)],dim(CGNM_result$Y)[1]),nrow=dim(CGNM_result$Y)[2]))

    ",paste(paste0("x",seq(1,length(ll$parameterInfo_dat$ParameterName)),"=CGNM_result$X[,",seq(1,length(ll$parameterInfo_dat$ParameterName)),"]"),  collapse="\n    "),"

}

## If one wishes to do middle out in normal scale make changes to the following code by removing log10
out=cbind(out,
",

    paste0("    weight_",ll$parameterInfo_dat$ParameterName,"*(log10(eval(parse(text=CGNM_result$runSetting$ReparameterizationDef[",seq(1,length(ll$parameterInfo_dat$ParameterName)),"])))-", paste0("log10(middleOutValue_",ll$parameterInfo_dat$ParameterName),"))", collapse = ",\n"),
"
)

return(rowSums((out)^2,na.rm = TRUE))"
    ,
    "
}
middleOutValue_df=data.frame(value=c(",paste(paste0("middleOutValue_", ll$parameterInfo_dat$ParameterName ), collapse=","),"),
                           parameterName=c(\"",paste(ll$parameterInfo_dat$ParameterName , collapse="\",\""),"\"))

middleOutValue_df=middleOutValue_df[c(",paste(paste0("weight_", ll$parameterInfo_dat$ParameterName ), collapse="!=0,"),
"!=0),]

plot_profileLikelihood(CGNM_result,Likelihood_function = postHoc_likelihood)+geom_vline(data=middleOutValue_df,aes(xintercept=value))")

  }

  puttogether_CGNM_code = function() {
    # writing out ODE
    compileODE_func(input$ODE_text)
    #parameterInfo_Data()
    #Dose_data()
    #ObservedData_DATA()

    odeCodeText = paste0("ODE_text=\"\n",
                         input$ODE_text,
                         "\"\n compiledModel=RxODE(ODE_text)")



    #parameterText = paste("x=c(",
    #                      paste(ll$parameterInfo_dat$Initial_lower_range, collapse = ","),
    #                      ")\n")

    # writing nonlinear function

    parameterText = paste0(
      "\nmodel_function=function(x){\n",

    #  "\nmodel_function=function(x, dataSet){\n",
      paste(
        paste0(ll$parameterInfo_dat$ParameterName, "=x[", seq(1, dim(ll$parameterInfo_dat)[1]), "]"),
        collapse = "\n"
      )
    )

    parameterText = paste0(
      "\n\n",
      parameterText,
      "\n\nmodelingParameter=c(",
      paste(paste0(rv$ODEparameter, "=", rv$ODEparameter), collapse = ","),
      ")\n\n"
    )

    unique_IDs = unique(c(ll$ObservedData_dat$ID, ll$Dose_dat$ID))




    dataSet_text = paste0(
      #"\n\ndataSet=data.frame(seq=seq(1,",
      "\n\ndataSet<<-data.frame(seq=seq(1,",

      dim(ll$ObservedData_dat)[1],
      "),
      time=c(",
      paste(ll$ObservedData_dat$time, collapse = ","),
      "),
      ID=c(",
      paste(ll$ObservedData_dat$ID, collapse = ","),
      "),
      Observed_value=c(",
      paste(ll$ObservedData_dat$Observed_value, collapse = ","),
      "),
      Observation_expression=c(\"",
      paste(ll$ObservedData_dat$Observation_expression, collapse = "\",\""),
      "\"))\n")

      dose_obs_Text = paste0("simResult_df=data.frame()\n")

    for (ID_nu in unique_IDs) {
      dose_df = subset(ll$Dose_dat, ID == ID_nu)
      obs_df = subset(ll$ObservedData_dat, ID == ID_nu)

      dose_obs_Text = paste0(dose_obs_Text, "\n\n ## ID: ", ID_nu, "
      ev <- eventTable()")

      if(dim(dose_df)[1]>0){


      for (i in seq(1, dim(dose_df)[1])) {
        if (is.na(dose_df$rate[i])) {
          dose_obs_Text = paste0(
            dose_obs_Text,
            "
      ev$add.dosing(dose = ",
            dose_df$dose[i],
            ", start.time=",
            dose_df$start.time[i],
            ",dosing.to=\"",
            dose_df$dosing.to[i],
            "\")"
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
            ",dosing.to=\"",
            dose_df$dosing.to[i],
            "\")"
          )
        }
      }
      }


      if(dim(obs_df)[1]>0){

        dose_obs_Text = paste0(
          dose_obs_Text,
          "\n      ev$add.sampling(c(",
          paste(sort(as.numeric(unique(obs_df$time))), collapse = ", "),
          "))\n      odeSol=compiledModel$solve(modelingParameter, ev)"
        )

        uniqueObsVariables = unique(obs_df$Observation_expression)

        for (obsVariable_nu in uniqueObsVariables) {
          dose_obs_Text = paste0(
            dose_obs_Text,
            "\n      simResult_df=rbind(simResult_df,data.frame(value=with(data.frame(odeSol), ",
            obsVariable_nu ,
            "), Observation_expression=\"",
            obsVariable_nu,
            "\", time=odeSol[,\"time\"], ID=",
            ID_nu,
            "))"
          )
        }
      }
    }
    dose_obs_Text = paste0(
      dose_obs_Text,
      "

      mergedData=merge(dataSet,simResult_df, all.x = TRUE)
      mergedData=mergedData[order(mergedData$seq),]
      return(mergedData$value)
}"
    )
    #(ParameterName=temp_PNV,Initial_lower_range=1.11, Initial_upper_range=NA, Lower_bound=0, Upper_bound=NA, Unit=NA)


    CGNM_run_text = paste0(
      paste0(
        "
ParaNames=c(\"",
        paste(ll$parameterInfo_dat$ParameterName, collapse = "\",\""),
        "\")
"
      ),

      paste0("UR=c(",
        paste(ll$parameterInfo_dat$Initial_upper_range, collapse = ","),
        ")
"
      ),
      paste0(
        "LR=c(",
        paste(ll$parameterInfo_dat$Initial_lower_range, collapse = ","),
        ")
"
      ),
      paste0(
        "U_bound=c(",
        paste(as.numeric(ll$parameterInfo_dat$Upper_bound), collapse = ","),
        ")
"
      ),
      paste0(
        "L_bound=c(",
        paste(as.numeric(ll$parameterInfo_dat$Lower_bound), collapse = ","),
        ")
"
      ),
      paste0(
        "observation=dataSet$Observed_value
"
      )
    )


    testCode_text <-
      paste0(
        "library(CGNM)\nlibrary(rxode2)\n\n",
        odeCodeText,
        dataSet_text,
        parameterText,
        dose_obs_Text,
        CGNM_run_text,
        "\n\n"
      )

    return(testCode_text)


  }
  # code generation end

##parallelization Code ----

  parallelCode_mac="


library(parallel)

obsLength=length(observation)

## Given CGNM searches through wide range of parameter combination, it can encounter
## parameter combinations that is not feasible to evaluate. This try catch function
## is implemented within CGNM for regular functions but for the matrix functions
## user needs to implement outside of CGNM

modelFunction_tryCatch=function(x){
  out=tryCatch({model_function(x)},
               error=function(cond) {rep(NA, obsLength)}
  )
  return(out)
  }

model_matrix_function=function(x){
  Y_list=mclapply(split(x, rep(seq(1:nrow(x)),ncol(x))), modelFunction_tryCatch,mc.cores = (parallel::detectCores()-1), mc.preschedule = FALSE)

  Y=t(matrix(unlist(Y_list),ncol=length(Y_list)))

  return(Y)
  }

CGNM_result=Cluster_Gauss_Newton_method(model_matrix_function,targetVector = observation, initial_lowerRange = LR, initial_upperRange = UR, lowerBound = L_bound, upperBound = U_bound,ParameterNames = ParaNames)
"


  parallelCode_win="


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
 out=tryCatch({model_function(x)},
              error=function(cond) {rep(NA, obsLength)}
 )
 return(out)
}

model_matrix_function=function(x){
  X=as.matrix(x)

  if(is.matrix(X)){
    Y_list=foreach(i=1:dim(X)[1], .export = c(\"model_function\", \"modelFunction_tryCatch\", \"dataSet\", \"obsLength\", \"compiledModel\"), .packages = c(\"rxode2\"))%dopar%{
      modelFunction_tryCatch(as.numeric(X[i,]))
    }

    Y=t(matrix(unlist(Y_list),ncol=length(Y_list)))

  }else{

   Y= modelFunction_tryCatch(X)
  }

  return(Y)

}

CGNM_result=Cluster_Gauss_Newton_method(model_matrix_function,targetVector = observation, initial_lowerRange = LR, initial_upperRange = UR, lowerBound = L_bound, upperBound = U_bound,ParameterNames = ParaNames)

stopCluster(cl=cluster)
"


#
# parallelCode_mac="
#
#
# library(parallel)
#
# obsLength=length(observation)
#
# ## Given CGNM searches through wide range of parameter combination, it can encounter
# ## parameter combinations that is not feasible to evaluate. This try catch function
# ## is implemented within CGNM for regular functions but for the matrix functions
# ## user needs to implement outside of CGNM
#
# modelFunction_tryCatch=function(x, dataSet){
#   out=tryCatch({model_function(x, dataSet)},
#                error=function(cond) {rep(NA, obsLength)}
#   )
#   return(out)
#   }
#
# model_matrix_function=function(x, dataSet){
#   Y_list=mclapply(split(x, rep(seq(1:nrow(x)),ncol(x))), modelFunction_tryCatch,mc.cores = (parallel::detectCores()-1), mc.preschedule = FALSE)
#
#   Y=t(matrix(unlist(Y_list),ncol=length(Y_list)))
#
#   return(Y)
#   }
#
# CGNM_result=Cluster_Gauss_Newton_method(model_matrix_function,targetVector = observation, initial_lowerRange = LR, initial_upperRange = UR, lowerBound = L_bound, upperBound = U_bound,ParameterNames = ParaNames, dataSet=dataSet)
# "
#
#
# parallelCode_win="
#
#
# library(foreach)
# library(doParallel)
#
# numCoretoUse=detectCores()-1
# registerDoParallel(numCoretoUse)
# cluster=makeCluster(numCoretoUse, type =\"PSOCK\")
# registerDoParallel(cl=cluster)
#
# obsLength=length(observation)
#
# # Given CGNM searches through wide range of parameter combination, it can encounter
# # parameter combinations that is not feasible to evaluate. This try catch function
# # is implemented within CGNM for regular functions but for the matrix functions
# # user needs to implement outside of CGNM
#
# modelFunction_tryCatch=function(x, dataSet){
#  out=tryCatch({model_function(x, dataSet)},
#               error=function(cond) {rep(NA, obsLength)}
#  )
#  return(out)
# }
#
# model_matrix_function=function(x, dataSet){
#   X=as.matrix(x)
#
#   if(is.matrix(X)){
#     Y_list=foreach(i=1:dim(X)[1], .export = c(\"model_function\", \"modelFunction_tryCatch\", \"dataSet\", \"obsLength\", \"compiledModel\", \"dataSet\"), .packages = c(\"rxode2\"))%dopar%{
#       modelFunction_tryCatch(as.numeric(X[i,]), dataSet)
#     }
#
#     Y=t(matrix(unlist(Y_list),ncol=length(Y_list)))
#
#   }else{
#
#    Y= modelFunction_tryCatch(X)
#   }
#
#   return(Y)
#
# }
#
# CGNM_result=Cluster_Gauss_Newton_method(model_matrix_function,targetVector = observation, initial_lowerRange = LR, initial_upperRange = UR, lowerBound = L_bound, upperBound = U_bound,ParameterNames = ParaNames, dataSet=dataSet)
#
# stopCluster(cl=cluster)
# "

### OE: makeCode_button pressed----
observeEvent(input$makeCode_button, {

  output[["test1_output"]] <- renderText({
    paste0(puttogether_CGNM_code(),"\n\n\n",
    makeCodeForPosthoc_middleout())
  })
})

  ### OE: testModelCode_button pressed----
  observeEvent(input$testModelCode_button, {
    #
    # putTogetherModel_message = capture.output(withCallingHandlers(
    #   tryCatch({
    #     testCode_text=puttogether_CGNM_code()
    #   }, error = function(e) {
    #     err <<- conditionMessage(e)
    #   }),
    #   warning = function(w) {
    #     warn <<- conditionMessage(w)
    #   }
    # ))
    #
    #
    # output$testModelCode_output <-
    #   renderUI({
    #     HTML(paste(unlist(
    #       strsplit(putTogetherModel_message, split = "\n", fixed = TRUE)
    #     ), collapse = '<br/>'))
    #   })




    output[["testModelCode_output"]] <-
      renderText({testCode_text=puttogether_CGNM_code()
        eval(parse(text = paste(
          testCode_text, paste("x=c(",
                          paste(ll$parameterInfo_dat$Initial_lower_range, collapse = ","),
                          ")\n\nmodel_function(x)")
        #  ")\n\nmodel_function(x, dataSet)")
        )))
      })
  })


  ### OE: testCGNMCode_button pressed----
  observeEvent(input$testCGNMCode_button, {


    tryCatch({
      withCallingHandlers({
        shinyjs::html("testCode_output", "")
        testCode_text=puttogether_CGNM_code()
        use_testCode = paste0(
          testCode_text,
          "set.seed(0)

    CGNM_result=Cluster_Gauss_Newton_method(model_function,targetVector = observation, num_minimizersToFind = 50, num_iteration=1,initial_lowerRange = LR, initial_upperRange = UR, lowerBound = L_bound, upperBound = U_bound,ParameterNames = ParaNames)"
          #CGNM_result=Cluster_Gauss_Newton_method(model_function,targetVector = observation, num_minimizersToFind = 50, num_iteration=1,initial_lowerRange = LR, initial_upperRange = UR, lowerBound = L_bound, upperBound = U_bound,ParameterNames = ParaNames, dataSet=dataSet)"
        )
        eval(parse(text = use_testCode))
        message("<b>CGNM code verified</b>")

      }
      ,
      message = function(m) {
        shinyjs::html(
          id = "testCode_output",
          add = TRUE,
          html = paste0(m$message, '<br>')
        )
      })
    }, shiny = FALSE)

    #print("here4")

    #
    #
    # out_text = capture.output(capture.output(withCallingHandlers(
    #   tryCatch({
    #     eval(parse(text = use_testCode))
    #     print("CGNM R stcript verified")
    #   }, error = function(e) {
    #     conditionMessage(e)
    #   }),
    #   warning = function(w) {
    #     conditionMessage(w)
    #   }
    # )))
    # output[["testCode_output"]] <-
    #   renderUI({
    #     HTML(paste(unlist(
    #       strsplit(out_text, split = "[1]", fixed = TRUE)
    #     ), collapse = '<br/>'))
    #   })


    # output[["testCode_output2"]] <- renderText({
    #   eval(parse(text = use_testCode))$runSetting$runName
    # })

  })

  ### download CGNM code download handler----
  output$download_CGNMcode_button <- downloadHandler(
    filename = function() {
      paste("CGNM_run_scriptByShinyCGNM-", Sys.Date(), ".R", sep = "")
    },
    content = function(file) {
      testCode_text=puttogether_CGNM_code()
      use_testCode = paste0(
        testCode_text,
        "
    CGNM_result=Cluster_Gauss_Newton_method(model_function,targetVector = observation,initial_lowerRange = LR, initial_upperRange = UR, lowerBound = L_bound, upperBound = U_bound,ParameterNames = ParaNames)

")
      #    CGNM_result=Cluster_Gauss_Newton_method(model_function,targetVector = observation,initial_lowerRange = LR, initial_upperRange = UR, lowerBound = L_bound, upperBound = U_bound,ParameterNames = ParaNames, dataSet=dataSet)

      writeLines(use_testCode, file)
    }
  )

  ### download postHocMiddleOut code download handler----
  output$download_posthocMiddleout_button <- downloadHandler(
    filename = function() {
      paste("CGNM_postHocMiddleOut_run_scriptByShinyCGNM-", Sys.Date(), ".R", sep = "")
    },
    content = function(file) {
      middleOut_code_text=makeCodeForPosthoc_middleout()
      writeLines(middleOut_code_text, file)
    }
  )


  ### download observation file download handler----

  output$download_obs_csv <- downloadHandler(
    filename = function(){paste("CGNM_obsevation_data_file-", Sys.Date(), ".csv", sep = "")},
    content = function(fname){
      write.csv(ll$ObservedData_dat, fname)
    }
  )


  ### download parameter info file download handler----

  output$download_para_csv <- downloadHandler(
    filename = function(){paste("CGNM_parameter_file-", Sys.Date(), ".csv", sep = "")},
    content = function(fname){
      write.csv(ll$parameterInfo_dat, fname)
    }
  )


  ### download parameter info file download handler----

  output$download_dose_csv <- downloadHandler(
    filename = function(){paste("CGNM_dose_file-", Sys.Date(), ".csv", sep = "")},
    content = function(fname){
      write.csv(ll$Dose_dat, fname)
    }
  )


  output$download_CGNMcode_parallel_mac_button <- downloadHandler(
    filename = function() {
      paste("CGNM_run_macParallel_scriptByShinyCGNM-", Sys.Date(), ".R", sep = "")
    },
    content = function(file) {
      testCode_text=puttogether_CGNM_code()
      use_testCode = paste0(
        testCode_text,
        parallelCode_mac)
      writeLines(use_testCode, file)
    }
  )


  output$download_CGNMcode_parallel_win_button <- downloadHandler(
    filename = function() {
      paste("CGNM_run_winParallel_scriptByShinyCGNM-", Sys.Date(), ".R", sep = "")
    },
    content = function(file) {
      testCode_text=puttogether_CGNM_code()
      use_testCode = paste0(
        testCode_text,
        parallelCode_win)
      writeLines(use_testCode, file)
    }
  )



  ## save and load the status ----
  output$saveCurrentStatus <- downloadHandler(
    filename = function() {
      paste("CGNM_shiny_save-", Sys.Date(), ".rds", sep = "")
    },
    content = function(file) {
      saveRDS(
        list(
          ODE_text = input$ODE_text,
          parameterInfo_dat = ll$parameterInfo_dat,
          Dose_dat = ll$Dose_dat,
          ObservedData_dat = ll$ObservedData_dat
        ),
        file = file
      )
    }
  )



  ## Read para csv ----
  input_para_csv_file <- reactive({
    if (is.null(input$para_csv)) {
      return("")
    }

    # actually read the file
    read.csv(file = input$para_csv$datapath, header = TRUE)
  })


  observeEvent(input_para_csv_file(), {
    Dose_dat_read=input_para_csv_file()

    if(is.data.frame(Dose_dat_read)){

      validData=FALSE

      if(sum(names(Dose_dat_read)%in% c("ParameterName","Initial_lower_range","Initial_upper_range","Lower_bound", "Upper_bound","Unit"))==6){
        validData=TRUE
        Dose_dat_read=Dose_dat_read[, c("ParameterName","Initial_lower_range","Initial_upper_range","Lower_bound", "Upper_bound","Unit")]
      }else if (dim(Dose_dat_read)[2]==6){
        validData=TRUE
        names(Dose_dat_read)= c("ParameterName","Initial_lower_range","Initial_upper_range","Lower_bound", "Upper_bound","Unit")
      }else{
        shinyCatch(stop("csv file with dosing information needs to have exactly 6 columns containing ParameterName,Initial_lower_range,Initial_upper_range,Lower_bound, Upper_bound,Unit."))
      }

      if(validData){


        Dose_dat_update=ll$parameterInfo_dat
        overLappingParameters=Dose_dat_read$ParameterName[Dose_dat_read$ParameterName %in% Dose_dat_update$ParameterName]

        for(paraNu in overLappingParameters){
          Dose_dat_update[Dose_dat_update$ParameterName==paraNu,]=Dose_dat_read[Dose_dat_read$ParameterName==paraNu,]
        }



        if (!identical(Dose_dat_update, ll$parameterInfo_dat)) {
          ll$parameterInfo_dat<<- Dose_dat_update
          output[["parameterInfo_dt"]] <- renderDT({
            datatable(
              ll$parameterInfo_dat,
              editable = list(target = "cell"),
              selection = "none",
              extensions = "AutoFill",
              callback = JS(callback),
              options = list(lengthMenu = list(c(100, -1), c('100','all')),autoFill =  list(horizontal = FALSE))
            )
          }, server = FALSE)
        }
      }

    }
  })


  ## Read dose csv ----
  input_dose_csv_file <- reactive({
    if (is.null(input$dose_csv)) {
      return("")
    }

    # actually read the file
    read.csv(file = input$dose_csv$datapath, header = TRUE)
  })


  observeEvent(input_dose_csv_file(), {
    Dose_dat_read=input_dose_csv_file()

    if(is.data.frame(Dose_dat_read)){

      validData=FALSE

      if(sum(names(Dose_dat_read)%in% c("ID","dose","dosing.to","start.time","rate"))==5){
        validData=TRUE
        Dose_dat_read=Dose_dat_read[, c("ID","dose","dosing.to","start.time","rate")]
      }else if (dim(Dose_dat_read)[2]==5){
        validData=TRUE
        names(Dose_dat_read)= c("ID","dose","dosing.to","start.time","rate")
      }else{
        shinyCatch(stop("csv file with dosing information needs to have exactly 5 columns containing ID, dose, dosing.to, start.time, and rate."))
      }

      if(validData){
        updateNumericInput(session, "numDoses", value = dim(Dose_dat_read)[1])

        if (!identical(Dose_dat_read, ll$Dose_dat)) {
          ll$Dose_dat <<-Dose_dat_read

          output[["dose_dt"]] <- renderDT({
            datatable(
              ll$Dose_dat,
              editable = list(target = "cell"),
              selection = "none",
              extensions = "AutoFill",
              callback = JS(callback),
              options = list(lengthMenu = list(c(100, -1), c('100','all')), autoFill =  list(horizontal = FALSE))
            )
          }, server = FALSE)
        }
      }

    }
  })



  ## Read observation csv ----
  input_obs_csv_file <- reactive({
    if (is.null(input$obs_csv)) {
      return("")
    }

    # actually read the file
    read.csv(file = input$obs_csv$datapath, header = TRUE)
  })


  observeEvent(input_obs_csv_file(), {
    Obs_dat_read=input_obs_csv_file()

    if(is.data.frame(Obs_dat_read)){

      validData=FALSE

      if(sum(names(Obs_dat_read)%in% c("ID","time","Observation_expression","Observed_value", "Memo"))==5){
        validData=TRUE
        Obs_dat_read=Obs_dat_read[, c("ID","time","Observation_expression","Observed_value", "Memo")]
      }else if (dim(Obs_dat_read)[2]==5){
        validData=TRUE
        names(Obs_dat_read)= c("ID","time","Observation_expression","Observed_value", "Memo")
      }else{
        shinyCatch(stop("csv file with observation information needs to have exactly 5 columns containing ID, time, Observation_expression, Observed_value, and Memo"))
      }

      if(validData){

        updateNumericInput(session, "numObservations", value = dim(Obs_dat_read)[1])

        if (!identical(Obs_dat_read,ll$ObservedData_dat)) {
          ll$ObservedData_dat <<- Obs_dat_read
          output[["ObservedData_dt"]] <- renderDT({
            datatable(
              ll$ObservedData_dat,
              editable = list(target = "cell"),
              selection = "none",
              extensions = "AutoFill",
              callback = JS(callback),
              options = list(lengthMenu = list(c(100, -1), c('100','all')),autoFill =  list(horizontal = FALSE))
            )
          }, server = FALSE)
        }
      }

    }
  })


  ### Reactive restore file----
  restore_file <- reactive({
    validate(need(input$openSavedStatus, message = FALSE))
    input$openSavedStatus
  })

  ### Reactive to store restored information----
  restored_state <- reactive({
    rs <- readRDS(restore_file()$datapath)
    rs
  })

  ### Restore state button pressed----
  observeEvent(restored_state(), {
    restoreState_function(restored_state())
  })

  ### Restore state function ----
  restoreState_function = function(rs) {
    updateTextAreaInput(session, "ODE_text", value = rs$ODE_text)

    compileODE_func(rs$ODE_text)

    if (!identical(rs$Dose_dat, ll$Dose_dat)) {
      ll$Dose_dat <<- rs$Dose_dat

      output[["dose_dt"]] <- renderDT({
        datatable(
          ll$Dose_dat,
          editable = list(target = "cell"),
          selection = "none",
          extensions = "AutoFill",
          callback = JS(callback),
          options = list(lengthMenu = list(c(100, -1), c('100','all')),autoFill =  list(horizontal = FALSE))
        )
      }, server = FALSE)
    }

    updateNumericInput(session, "numDoses", value = dim(rs$Dose_dat)[1])

    if (!identical(rs$parameterInfo_dat, ll$parameterInfo_dat)) {
      ll$parameterInfo_dat<<- rs$parameterInfo_dat
      output[["parameterInfo_dt"]] <- renderDT({
        datatable(
          ll$parameterInfo_dat,
          editable = list(target = "cell"),
          selection = "none",
          extensions = "AutoFill",
          callback = JS(callback),
          options = list(lengthMenu = list(c(100, -1), c('100','all')),autoFill =  list(horizontal = FALSE))
        )
      }, server = FALSE)
    }

    updateNumericInput(session, "numObservations", value = dim(rs$ObservedData_dat)[1])

    if (!identical(rs$ObservedData_dat,ll$ObservedData_dat)) {
      ll$ObservedData_dat <<- rs$ObservedData_dat
      output[["ObservedData_dt"]] <- renderDT({
        datatable(
          ll$ObservedData_dat,
          editable = list(target = "cell"),
          selection = "none",
          extensions = "AutoFill",
          callback = JS(callback),
          options = list(lengthMenu = list(c(100, -1), c('100','all')),autoFill =  list(horizontal = FALSE))
        )
      }, server = FALSE)
    }
  }



  ### Load example button pressed----
  observeEvent(input$loadExample_button,
               {
                 rs = readRDS("inst/CGNM_shiny_save-example1.rds")
                 restoreState_function(rs)
               })

}

#shinyApp(ui = ui,
#         server = server,
#         enableBookmarking = "server")
