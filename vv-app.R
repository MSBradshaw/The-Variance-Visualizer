library(shiny)
library(ggplot2)
library(tibble)
library(readr)
library(dplyr)
library(plotrix)

ui <- fluidPage(
  fileInput("file1", "Choose CSV File",
            accept = c(
              "text/csv",
              "text/comma-separated-values,text/plain",
              ".csv",".pileup"), multiple = TRUE
  ),
  selectInput("image_type", "Display Type",
              c('normal','variance','consensus','hybrid'), selected = 2
              ),
  numericInput("front", "Trim from front:", 600, min = 1, max = 10000),
  numericInput("end", "Trim from back:", 600, min = 1, max = 10000),
      htmlOutput("image_area")
  
)
server <- function(input, output) {
  
  source('parsepileup.R')
  
  output$image_area <- renderUI({
    if( is.null(input$file1) ){
      return(NULL)
    }
    names <- input$file1$name %>% str_replace('.pileup','')
    return( HTML( build_image( input$file1$datapath,input$image_type,'',input$front,input$end,names ) ) )
    })
}
shinyApp(ui = ui, server = server)
