# App to run CWD Model


library(shiny)

# Packages
if(!require(ggplot2)) install.packages(ggplot2)
if(!require(tidyverse)) install.packages(tidyverse)
if(!require(deSolve)) install.packages(deSolve)


# Imporvements to do:
# - be able to compare models with diferent parametes
# - explain what each parameter does
# - make an reacive value function so plot render is not bulky
# - make the other parameters vary, and set proper ranges
# - make better plot





# Model Specification (user defined function) ------------------------------

SIEmod <- function(times, State, Pars) {
    
    with(as.list(c(State, Pars)), {
        
        dS <-   alpha - S * (beta * I + gamma * E + m)
        dI <-   S * (beta * I + gamma * E) - I * (m + mu)
        dE <-   epsilon * I - tau * E
        
        return(list(c(dS, dI, dE)))
    })
}



# UI ----------------------------------------------------------------------

# Define UI for application that draws a histogram
ui <- fluidPage(
    
    titlePanel("CWD Model"), 
    
    sidebarLayout(
        sidebarPanel(
            position = "right",
            
            sliderInput(inputId = "beta", 
                        label = "beta", 
                        min = 0, 
                        value = 0.002446, 
                        max = 0.005, 
                        step = 0.0001), 
            
            sliderInput(inputId = "mu", 
                        label = "mu", 
                        min = 0, 
                        value = 2.617254, 
                        max = 3,
                        step = 0.05),
            
            sliderInput(inputId = "time", 
                        label = "time", 
                        min = 0, 
                        value = 10, 
                        max = 50,
                        step = 5)
        ),
        
        mainPanel(
            plotOutput("plot1")
        )
    )
)



# Server ------------------------------------------------------------------


# Define server logic required to draw a histogram
server <- function(input, output) {

    # rv <- reactiveValues()
    
    alpha <- 4.48307
    m <- 0.103202
    gamma <- 0.206146
    epsilon <- 0.150344
    tau <- 0.135785
    S0 <- 180
    
    # rv$beta <- 0.002446
    
    # rv$mu <- 2.617254 

    
    
    output$plot1 <- renderPlot({
        parsSIE <- c(alpha = alpha,
                     m = m, 
                     gamma = gamma,
                     epsilon = epsilon,
                     tau = tau,
                     S0 = S0,
                     beta = input$beta,
                     mu =  input$mu)
        
        times_vec <- seq(0, input$time, 1)
        
        out <- as.data.frame(ode(func = SIEmod, y = c(S = S0, I = 0 ,  E =0.0000003),  
                                 parms = parsSIE, times = times_vec))
        
        names(out)[2:4] <- c("Susceptible", "Infectious", "Virus")
        
        # Converting to long format for plots
        outL <- out %>% 
            gather(key = Compartment, value = Number, -time) %>% 
            data.frame()
        
        outL$Compartment <- factor(outL$Compartment, levels = c("Susceptible", "Infectious", "Virus"))
        
        outL <- bind_cols(outL, data.frame(t(matrix(parsSIE))))
        head(outL)
        
        names(outL)[4:ncol(outL)] <- c("alpha", "m", "gamma", "epsilon",
                                       "tau", "S0", "beta", "mu")
        
        # Output <- rbind(Output, outL)
        Output <- outL
        
        SIE_plot <- Output  %>%  filter(Compartment== "Infectious") %>%
            mutate(beta = as.character(round(beta, 6)),
                   mu = as.character(round(mu, 3))) %>%
            ggplot(aes(x = time, y = Number, group = mu, color = mu)) +
            ylab("Number of infectious hosts") + 
            xlab("Time") + geom_line(size = 1) + facet_wrap(~ beta) +
            theme_bw() + theme(text = element_text(size = 20), 
                               strip.background = element_blank())
        SIE_plot
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
