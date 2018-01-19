
###################################Shiny Modelling Tool###########################################
##################################################################################################


##########Imports external R files into the tool#####################################
source("/eas/struktureinheiten/quz/studenten/NishantKumar/shiny_app/FitGLD.R");
source("/eas/struktureinheiten/quz/studenten/NishantKumar/shiny_app/PlotScatter.R");


##########Load R packages required for the tool######################################
require(MASS);
require(Matrix); 
require(gld);
library(shiny)

#### inits

#########returns an object for storing reactive values###############################
CURRENT_DATA <- reactiveValues(
  all=data.frame(),
  plot=data.frame()
)

#################################################################################################
#########UI Script of the tool###################################################################
#################################################################################################
##pageWithSidebar Creates a Shiny UI that contains a header with the application title, a sidebar for input controls, and a main area for output.###
ui <- pageWithSidebar(                     
  ##Creates header with the application title
  headerPanel("CSV Viewer"),       
  ##Creates sidebar for input controls
  sidebarPanel(
    ##Creates a file upload control to upload files
    fileInput('file1', 'Choose CSV File',   
              accept=c('text/csv','text/comma-separated-values,text/plain','.csv')),
    ##Generates a checkbox control
    checkboxInput('header', 'Header', TRUE), 
    ##Creates a set of radio buttons used to select an item from a list
    radioButtons(inputId='sep', label='Column separator',choices=c(Comma=',', Semicolon=';',Tab='\t'), selected=',', inline=TRUE),
    radioButtons('quote', 'Quote',c(None='','Double Quote'='"','Single Quote'="'"),'"', inline=TRUE),
    hr(),
    fluidRow(
      column(6,
             checkboxGroupInput(inputId='Models',label='Select Models',choices=c(Gaussian='Gaussian',Uniform ='Uniform',GLD='GLD'), selected=c(), inline=TRUE))
    ),
    uiOutput(outputId="uiColumnSelection")
  ),
  ##Creates a main area for output Plots and Images
  mainPanel(  
    ##Tabsets used for dividing output into multiple independently viewable sections
    tabsetPanel(
      tabPanel("Plot",plotOutput("plot")),
      tabPanel("Data",tableOutput("contents"))     
    ),
    ## Creates a panel only when a conditional expression is true
    conditionalPanel(
      condition = "input.Models == 'Gaussian'",
      plotOutput('gaussian')
    ),
    conditionalPanel(
      condition = "input.Models == 'GLD'",
      plotOutput('gld')
    )  
  )
)

#### dynamic part of UI###########
uiColumnSelection <- function(CURRENT_DATA) {
  if( any(dim(CURRENT_DATA$all)== 0) ) {
    list(NULL)
  } else {
    list(
         hr(),
         fluidRow(
         column(6,checkboxGroupInput(inputId='Columns',label='Select columns',choices=colnames(CURRENT_DATA$all), selected=colnames(CURRENT_DATA$plot)))
         )
        )
  }
}

################################################################################
################################ server ########################################
################################################################################
server <- function(input,output,session) {
  
  ############## render UI for column selection ################
  output$uiColumnSelection <- renderUI(uiColumnSelection(CURRENT_DATA));
  
  ############## read data from selected file   ################
  observe({
    if( is.null(input$file1) ) {
      ## no input file selected
      CURRENT_DATA$all  <<- data.frame();
      CURRENT_DATA$plot <<- data.frame();
    } else {
      isolate({
        if( file.exists(input$file1$datapath) ) {
          CURRENT_DATA$all <<- read.csv(file=input$file1$datapath, header=input$header, sep=input$sep, quote=input$quote)
          if( ncol(CURRENT_DATA$all) > 4 ) {
            CURRENT_DATA$plot <<- CURRENT_DATA$all[,1:4]
          } else {
            CURRENT_DATA$plot <<- CURRENT_DATA$all
          }
        } else {
          ## file cannot be accessed
          CURRENT_DATA$all <<- data.frame();
          CURRENT_DATA$plot <<- data.frame();
        }
      })
    }
  })
  
  ################ observe checkbox group input #################
  observe({
    if( ! is.null(input$Columns) ) {
 #browser()
 isolate({
        if( length(input$Columns)!=ncol(CURRENT_DATA$plot) ) {
          if( length(input$Columns)==0 ) {
            ## nothing selected
            CURRENT_DATA$plot <<- data.frame()
          } else {
            for( i in 1:length(input$Columns) ) {
              if( i==1 ) {
                TMP_DF <- CURRENT_DATA$all[[input$Columns[i]]];
              } else {
                TMP_DF <- cbind(TMP_DF,CURRENT_DATA$all[[input$Columns[i]]]);
              }
            }
            TMP_DF <- data.frame(TMP_DF); colnames(TMP_DF) <- input$Columns
            CURRENT_DATA$plot <<- TMP_DF
          }
        } else {
          ## selection unchanged --> nothing to do here
        }
      })
    } else {
      ## nothing here since dynamic UI part not yet created
    }
  })
  
  ######### plot ###########
  output$plot <- renderPlot({
    DATA <- list("raw data"=CURRENT_DATA$plot)
    ## extend data depend on model selection
    PlotScatter(DATA=DATA)
  })
  
  
  #########Gaussian########################
  output$gaussian = renderPlot({ 
    DATA <- CURRENT_DATA$plot
    sample_size <- 300;
    MODELS <- list();
    X <- DATA;
    MODELS$covmat <- cov(as.matrix(X));    
    MODELS$Pearson <- cor(X, method="pearson");
    MODELS$Spearman <- cor(X, method="spearman");
    MODELS$means <- apply(DATA,2,mean);
    SAMPLED <- list();
    SAMPLED$Gaussian <- mvrnorm(n=sample_size,mu=MODELS$means,Sigma=MODELS$covmat);
    DATA_G <- list(
      "raw_data" = DATA,  
      "Gaussian model"=SAMPLED$Gaussian
    ) 
    
    if (!is.null(DATA_G[[1]])){
      PlotScatter(DATA=DATA_G)
    }
  })
  
  
  #############GLD#######################
  output$gld = renderPlot({
    DATA <- CURRENT_DATA$plot
    sample_size <- 300;
    MODELS <- list();
    X <- DATA;
    MODELS$covmat <- cov(as.matrix(X));
    MODELS$Pearson <- cor(X, method="pearson");
    MODELS$Spearman <- cor(X, method="spearman");
    MODELS$means <- apply(DATA,2,mean);
    MODELS$gld_coeffs <- apply(DATA,2,FitGLD);
    SAMPLED <- list();
    SAMPLED$Gaussian <- mvrnorm(n=sample_size,mu=MODELS$means,Sigma=MODELS$covmat);
    C <- nearPD(2*sin(pi/6*MODELS$Spearman))$mat;
    Z <- mvrnorm(n=sample_size,mu=rep(0,nrow(C)),Sigma=C);
    U <- pnorm(q=Z, mean=0, sd=1)
    SAMPLED$GLD <- SAMPLED$Gaussian;
    for( i in 1:ncol(SAMPLED$GLD) ) {
      SAMPLED$GLD[,i] <- qgl(p=U[,i], lambda1=MODELS$gld_coeffs[,i])
    }
    DATA_F <- list(
      "raw_data" = DATA,  
      "Gaussian model"=SAMPLED$Gaussian,
      "GLD model"=SAMPLED$GLD
    )
    if (!is.null(DATA_F[[1]])){
      PlotScatter(DATA=DATA_F)
    }
  })
}

####### combines ui and server parts of the code into a functioning app ###########
shinyApp(ui=ui, server=server)
