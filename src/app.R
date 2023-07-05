library(shiny)
library(shinydashboard)
library(zip)

# UI
ui <- dashboardPage(
  dashboardHeader(title = "Phenotypic/Genotypic Data Analysis"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data", tabName = "data_tab"),
      menuItem("Model Tuning", tabName = "model_tuning_tab"),
      menuItem("Results", tabName = "results_tab")
    )
  ),
  dashboardBody(
    tabItems(
      # Data tab
      tabItem(
        tabName = "data_tab",
        h2("Data"),
        fluidRow(
          box(
            title = "Upload Files",
            fileInput("phenotypic_data", "Phenotypic Data", multiple = FALSE),
            fileInput("molecular_markers", "Molecular Markers", multiple = FALSE)
          )
        ),
        fluidRow(
          box(
            title = "Data Configuration",
            numericInput("env_col", "Env ID Column", value = 12),
            numericInput("trait_col", "Target Trait Column", value = 3),
            numericInput("gid_col", "Genotype/Variety ID Column", value = 2),
            numericInput("nan_freq", "NaN Threshold Limit %", value = 20),
            numericInput("maf_prop_j", "MAF Proportion", value = 0, min = 0, max = 1)
          )
        ),
        actionButton("run_analysis", "Run Analysis")
      ),
      
      # Model Tuning tab
      tabItem(
        tabName = "model_tuning_tab",
        h2("Model Tuning"),
        fluidRow(
          box(
            title = "Model Options",
            checkboxInput("e_l", "E + L", value = T),
            checkboxInput("e_l_g", "E + L + G", value = T),
            checkboxInput("e_l_g_ge", "E + L + G + GE", value = T)
          ),
          box(
            title = "Other Options",
            checkboxInput("centering", "Centering", value = T),
            checkboxInput("standardization", "Standardization", value = T),
            checkboxInput("weighting", "Weighting", value = FALSE),
            checkboxInput("esc", "ESC", value = FALSE)
          )
        ),
        fluidRow(
          box(
            title = "Model Parameters",
            numericInput("folds", "Folds", value = 5),
            numericInput("n_iter", "Number of Iterations (Optional)", value = 1000),
            numericInput("burn_in", "Burn In (Optional)", value = 100)
          )
        ),
        actionButton("run_analysis", "Run Analysis")
      ),
      
      # Results tab
      tabItem(
        tabName = "results_tab",
        h2("Results"),
        fluidRow(
          box(
            title = "Download Results",
            downloadButton("download_link", "Download Results")
          )
        ),
        fluidRow(
          box(
            title = "View Data",
            uiOutput("data_files")
          )
        )
      )
    )
  )
)
# Server
server <- function(input, output, session) {

  zipReportFiles <- function() {
    # Set the working directory to the "report" folder
    setwd("../output/report")
    
    # Get the file names in the "report" folder
    fileNames <- list.files()
    
    # Create a temporary zip file
    tempZipFile <- tempfile(fileext = ".zip")
    
    # Create the zip file with the files in the "report" folder
    zip(tempZipFile, files = fileNames)
    
    # Return to the original working directory
    setwd("../../src")
    
    # Return the path to the zip file
    return(tempZipFile)
  }
  
  observeEvent(input$run_analysis, {
    # Run analysis code here
    # You can access the user inputs using input$<input_id>
    # For example, input$phenotypic_data will contain the uploaded phenotypic data file
    # input$env_col will contain the value for the env_col input
    # Modify the code below to integrate the user inputs with your existing R code
    
    # Load required packages and modules
    library(zip)
    library(shinyjs)
    library(data.table)
    library(ggplot2)
    library(gridExtra)
    library(BGLR)
    library(matrixStats)
    
    # Source modules
    source("modules/1_data_load.R")
    source("modules/2_matrices.R")
    source("modules/3_cv_prep.R")
    source("modules/4_fit_models.R")
    source("modules/5_output_results.R")
    
    # Read user inputs
    phenos.path <- input$phenotypic_data$datapath
    marker.path <- input$molecular_markers$datapath
    env.col <- input$env_col
    trait.col <- input$trait_col
    gid.col <- input$gid_col
    nan.freq <- input$nan_freq
    prop.maf.j <- input$maf_prop_j / 100
    e.l <- input$e_l
    e.l.g <- input$e_l_g
    e.l.g.ge <- input$e_l_g_ge
    ctr <- input$centering
    std <- input$standardization
    weighting <- input$weighting
    esc <- input$esc
    folds <- input$folds
    nIter <- input$n_iter
    burnIn <- input$burn_in
    
    # Import packages
    library(data.table)
    library(ggplot2)
    library(gridExtra)
    library(BGLR)
    library(matrixStats)
    
    # Import modules
    source("modules/1_data_load.R")
    source("modules/2_matrices.R")
    source("modules/3_cv_prep.R")
    source("modules/4_fit_models.R")
    source("modules/5_output_results.R")
    
    cv1 <- TRUE
    cv2 <- TRUE
    cv0 <- TRUE
    cv00 <- TRUE
    
    ### 1 - Data Load ###
    loaded.data <- loadData(phenos.path, marker.path)
    
    # Create output directory if it doesn't exist
    if (!dir.exists("../output")) { dir.create("../output") }
    
    markers <- loaded.data$markers
    phenos <- loaded.data$phenos
    
    rm(loaded.data)
    
    # Create the NaNs csv files
    createNaNFiles(phenos, markers, nan.freq)  # These files are the Mod1 outputs
    
    ### 2 - G/E Matrices ###
    
    # E matrix
    generateMatrix("../output/E/", phenos, markers = NULL, col.env.id = env.col)
    # G matrix
    generateMatrix("../output/G/", phenos=phenos, markers=markers, 
                   col.env.id = gid.col, prop.maf.j =  prop.maf.j)
    # ZE matrix
    createZMatrix(phenos, env.col, "../output/ZE/")
    # ZL matrix
    createZMatrix(phenos, gid.col, "../output/ZL/")
    
    # Interaction matrix
    g1.file <- '../output/G/G.rda'          # path to matrix file 1
    g2.file <- '../output/E/G.rda'          # path to matrix file 2
    
    generateIntMatrix(g1.file, g2.file, output.path='../output/GE/')
    
    ### 3 - Phenotype data prep ###
    set.seed(1)
    phenos.cv <- phenos
    phenos.cv <- cvPrep(phenos.cv, "../output/cv/", col.id = gid.col, folds = folds,
                        cv1 = cv1, cv2 = cv2)
    phenos.cv <- cvPrep(phenos.cv, "../output/cv/", col.id = trait.col, folds = folds,
                        cv0 = cv0, cv00 = cv00)
    
    ### 4 - Fit models ####
    
    # current structure: cv --> trait --> E+L --> folds --> cols 15:19
    # Output structure: trait (e.g height) --> CV --> fold_n --> predictions.csv 
    
    ab.list <- list()
    ab.list[1] <- '../output/ZE/Z.rda'
    ab.list[2] <- '../output/ZL/Z.rda'
    ab.list[3] <- '../output/G/EVD.rda'  
    ab.list[4] <- '../output/GE/EVD.rda'  
    
    # Find the CV columns
    cv.list <- list(
      cv1 = match("CV1", colnames(phenos.cv)),
      cv2 = match("CV2", colnames(phenos.cv)),
      cv0 = grep("CV0_", colnames(phenos.cv)),
      cv00 = grep("CV00_", colnames(phenos.cv))
    )
    
    set.seed(1)
    for (i in 1:length(cv.list)) {
      cv <- names(cv.list)[i]
      val <- cv.list[i]$cv
      
      # Run E + L
      getPredictions(phenos.cv, cv, trait.col, gid.col, as.numeric(val), ab.list[1:2], 
                     folds = folds, esc = esc, nIter = nIter, burnIn = burnIn)
      # Run E + L + G
      getPredictions(phenos.cv, cv, trait.col, gid.col, as.numeric(val), ab.list[1:3], 
                     folds = folds, esc = esc, nIter = nIter, burnIn = burnIn)
      # Run E + L + G + GE
      getPredictions(phenos.cv, cv, trait.col, gid.col, as.numeric(val), ab.list, 
                     folds = folds, esc = esc, nIter = nIter, burnIn = burnIn)
      
    }
    
    ### 5 - Get Results ###
    getCvResults(phenos.cv, env.col, trait.col)
    
    # Call the zipReportFiles function to create the zip file
    zipPath <- zipReportFiles()
    
    # Display a download link for the zip file
    output$download_link <- downloadHandler(
      filename = "report.zip",
      content = function(file) {
        file.copy(zipPath, file)
      }
    )
    
    # Display notification in UI
    shinyalert::shinyalert(
      title = "Analysis Completed",
      text = "The analysis has finished. Please check the output directory for results.",
      type = "success"
    )

  })
  
  output$data_files <- renderUI({
    # Get the file names in the "output/report" directory
    fileNames <- list.files("../output/report", full.names = TRUE)
    
    # Generate download links for each file
    links <- lapply(fileNames, function(file) {
      fileName <- basename(file)
      downloadLink(file, fileName)
    })
    
    # Return the list of download links
    tagList(links)
  })
  
}

# Run the app
shinyApp(ui, server)
