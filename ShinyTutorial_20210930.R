#
#
#
# SCript developed from R Ladies East Lansing Shiny App Intro, 30 September 2021
# Slides: https://docs.google.com/presentation/d/1x1XDqkpQ6dt2O6hZNzS_6C0PHysC6pOt-si1qyz_4ZU/edit#slide=id.gf141dd44a3_0_2
# Code: https://github.com/rstudio-education/shiny.rstudio.com-tutorial
#
#
#
##############################################################################
############################### SIMPLE EXAMPLE ###############################
##############################################################################

# library(shiny)
# 
# # add elements to your app as arguments to fluidPage()
# ui <- fluidPage(
#   sliderInput(inputId = "num", 
#               label = "Choose a number", 
#               value = 25, min = 1, max = 100), 
#   
#   plotOutput("hist")
# )
# 
# server <- function(input, output, session) {
#   output$hist <- renderPlot({
#     hist(rnorm(input$num))}) # $hist must match "hist" in UI
#   
# }
# 
# shinyApp(ui=ui, server=server)

################################################################################
############################### SSL DATA ATTEMPT ###############################
################################################################################

library(tidyverse)
library(shiny)

# Read in steller sea lion data 
ssl_orig <- readRDS("../Data_Processed/SSL_UsedAndAvail_WithCovars.rds")  %>% 
  group_by(DeployID, Date) %>% 
  mutate(choice_id = cur_group_id()) %>% 
  ungroup()
ssl <- ssl_orig %>% 
  dplyr::select(DeployID, Used, choice_id, year, month, weekofyear, Bathymetry, DistLand, Dist500m, 
                slope, ssh, eke, wind, fish, ship, sst)
cols <- c("DeployID", "Used", "choice_id", "year", "month")
ssl[cols] <- lapply(ssl[cols], factor)



# add elements to your app as arguments to fluidPage()
ui <- fluidPage(
  # A scroll down list to pick a sea lion
  # A scroll down list to pick a variable
  sidebarPanel(
    selectInput("indl", "Individual", unique(ssl$DeployID)), 
    selectInput('xcol', 'X Variable', c("DeployID", "Used", "year", "month")),
    selectInput('ycol', 'Y Variable', names(ssl)[7:16], 
    selected = names(ssl)[[2]])),
  # A space for a violin plot 
  mainPanel(
    plotOutput('plot1'))
)

server <- function(input, output, session) {
  
  selectedData <- reactive({
    subset(ssl, DeployID ==  input$indl)
  })
    
  output$plot1 <- renderPlot({
    ggplot(selectedData(), aes_string(x=input$xcol, y=input$ycol, fill=input$xcol)) +
      geom_violin() +
      # geom_boxplot(width = 0.1, alpha=0.2, color="grey") +
      geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.2) +
      theme(axis.text.x = element_text(angle = 90))
  }) # $hist must match "hist" in UI
  
}

shinyApp(ui=ui, server=server)
