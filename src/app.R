library(shiny)
library(shinydashboard)
library(zip)
library(shinyalert)
library(data.table)
library(DT)

# UI
ui <- dashboardPage(
  dashboardHeader(title = "Phenotypic Trait Data Analysis"),
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
        )
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
        ),
        fluidRow(
          box(
            title = "Selected Data",
            style = "overflow-y: scroll; max-height: 400px;",
            tableOutput("selected_file_contents")
          )
        )
      )
    )
  )
)

# Server
server <- function(input, output, session) {

  # TODO - Switch all files to temporary directory
  output.path <- file.path(tempdir(), "app_data")
  if (!dir.exists(output.path)) { dir.create(output.path, recursive = TRUE) }

  # TODO: Fix function to only download zip file, not entire directory
  zipReportFiles <- function(output.path) {
    # Set the working directory to the "report" folder
    out.path <- file.path(output.path, "output", "report")

    # Create output directory if it doesn't exist
    if (!dir.exists(out.path)) { dir.create(out.path, recursive = TRUE) }
  
    # Get the file names in the "report" folder
    fileNames <- list.files(out.path)
  
    if (length(fileNames) == 0) {
      # Handle case when no files exist
      print("No files to include in the zip file.")
      return(NULL)
    }
  
    # Create the zip file path
    zipFilePath <- file.path(out.path, "results.zip")
  
    # Create the zip file with the files in the "report" folder
    zip(zipFilePath, files = file.path(out.path, fileNames))
  
    # Return the path to the zip file
    return(zipFilePath)
  }
  
  observeEvent(input$run_analysis, {
    
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
    
    cv1 <- TRUE
    cv2 <- TRUE
    cv0 <- TRUE
    cv00 <- TRUE
    
    ### 1 - Data Load ###
    loaded.data <- loadData(phenos.path, marker.path)
    markers <- loaded.data$markers
    phenos <- loaded.data$phenos
    
    rm(loaded.data)
    
    # Create the NaNs csv files
    createNaNFiles(phenos, markers, nan.freq, output.path)  # These files are the Mod1 outputs
    
    ### 2 - G/E Matrices ###
    
    # E matrix
    generateMatrix(file.path(output.path,"E"), phenos, markers = NULL, col.env.id = env.col)
    # G matrix
    generateMatrix(file.path(output.path,"G"), phenos=phenos, markers=markers, 
                   col.env.id = gid.col, prop.maf.j =  prop.maf.j)
    # ZE matrix
    createZMatrix(phenos, env.col, file.path(output.path,"ZE"))
    # ZL matrix
    createZMatrix(phenos, gid.col, file.path(output.path,"ZL/"))
    
    # Interaction matrix
    g1.file <- file.path(output.path,'G','G.rda')          # path to matrix file 1
    g2.file <- file.path(output.path,'E','G.rda')          # path to matrix file 2
    
    generateIntMatrix(g1.file, g2.file, output.path=file.path(output.path,'GE/'))
    
    ### 3 - Phenotype data prep ###
    set.seed(1)
    phenos.cv <- phenos
    phenos.cv <- cvPrep(phenos.cv, file.path(output.path,'cv'), col.id = gid.col, folds = folds,
                        cv1 = cv1, cv2 = cv2)
    phenos.cv <- cvPrep(phenos.cv, file.path(output.path,'cv'), col.id = trait.col, folds = folds,
                        cv0 = cv0, cv00 = cv00)
    
    ### 4 - Fit models ####
    
    # current structure: cv --> trait --> E+L --> folds --> cols 15:19
    # Output structure: trait (e.g height) --> CV --> fold_n --> predictions.csv 
    
    ab.list <- list()
    ab.list[1] <- file.path(output.path,'ZE/Z.rda')
    ab.list[2] <- file.path(output.path,'ZL/Z.rda')
    ab.list[3] <- file.path(output.path,'G/EVD.rda')
    ab.list[4] <- file.path(output.path,'GE/EVD.rda')
    
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
      getPredictions(phenos.cv, output.path, cv, trait.col, gid.col, as.numeric(val), ab.list[1:2], 
                     folds = folds, esc = esc, nIter = nIter, burnIn = burnIn)
      # Run E + L + G
      getPredictions(phenos.cv, output.path, cv, trait.col, gid.col, as.numeric(val), ab.list[1:3], 
                     folds = folds, esc = esc, nIter = nIter, burnIn = burnIn)
      # Run E + L + G + GE
      getPredictions(phenos.cv, output.path, cv, trait.col, gid.col, as.numeric(val), ab.list, 
                     folds = folds, esc = esc, nIter = nIter, burnIn = burnIn)
      
    }
    
    ### 5 - Get Results ###
    getCvResults(phenos.cv, env.col, trait.col, output.path)

    print("Results collected")
    
    # Call the zipReportFiles function to create the zip file
    zipPath <- zipReportFiles(output.path)

    print("zipPath executed")
    
    # Display a download link for the zip file
    output$download_link <- downloadHandler(
      filename = function() { 
        "results.zip"
      },
      content = function(file) {
        file.copy(zipPath, file)
      }
    )

    print("Download link created")
    
    # Display notification in UI
    shinyalert::shinyalert(
      session = session,
      title = "Analysis Completed",
      text = "The analysis has finished. Please check the output directory for results.",
      type = "success"
    )

    unlink(output.path)

  })
  
  output$data_files <- renderUI({
    # Get the file names in the "output/report" directory
    final.path <- file.path(output.path, "output", "report")
    fileNames <- list.files(final.path, full.names = TRUE)
    
    # Generate file links
    links <- lapply(fileNames, function(file) {
      fileName <- basename(file)
      link <- tagList(
        tags$a(href = "#", onclick = paste0("Shiny.setInputValue('selected_file', '", file, "')"), fileName)
      )
      return(link)
    })
    
    # Return the list of file links
    tagList(links)
  })
  
  observeEvent(input$selected_file, {
    if (!is.null(input$selected_file)) {
      # Read the selected file and render its contents
      data <- fread(input$selected_file)
      output$selected_file_contents <- renderTable({
        data
      })
    }
  })
  
}

# Run the app
shinyApp(ui, server)
