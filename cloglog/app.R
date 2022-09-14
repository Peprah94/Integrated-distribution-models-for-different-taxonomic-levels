#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Cloglog with Occupancy probability and detection probability"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("det_prob",
                        "Detection probability:",
                        min = 0,
                        max = 1,
                        value = 0.5)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
      cloglog <- function(psi, det_prob){
        ret <- log(-log(1-psi*det_prob))
        return(ret)
      }
      
      cloglog_offset <- function(psi, det_prob){
        ret <- log(-log(1-psi)) - log(det_prob)
        return(ret)
      }
        x    <- seq(0,1,0.1)
        lambda <- cloglog(x, input$det_prob)
        lambda1 <- cloglog_offset(x, input$det_prob)

        # draw the histogram with the specified number of bins
        par(mfrow = c(1,2))
        plot(x, lambda, col = 'darkgray')
        plot(x, lambda1)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
