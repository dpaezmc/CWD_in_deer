library(ggplot2)
library(tidyverse)
library(deSolve)

###---Model specification

SIEmod <- function(times, State, Pars) {

     with(as.list(c(State, Pars)), {

      dS <-   alpha - S * (beta * I + gamma * E + m)
      dI <-   S * (beta * I + gamma * E) - I * (m + mu)
      dE <-   epsilon * I - tau * E

	return(list(c(dS, dI, dE)))
	})
 }
#-
###---Definition of parameter values

   alpha <- 4.48307
   beta_vec <- rlnorm(25, log(0.002446), 2) ##Generating 25 beta values from a lognormal with mean equal to estimated value
   mu_vec <-   rlnorm(25, log(2.617254), 0.5) ## Generating 25 mu values from a lognormal with mean equal to estimated value
   m <- 0.103202
   gamma <- 0.206146
   epsilon <- 0.150344
   tau <- 0.135785
   S0 <- 180
 
   ###---Associate every value in beta_vec with every value in mu_vec
   beta_mu_vals <- expand.grid(beta_vec, mu_vec)
   names(beta_mu_vals) <- c("beta", "mu")
#---End of parameter set up
   
   
   ####---Output storage
   Output <- data.frame(matrix(ncol = 11, nrow = 0))
   
   ###--- For loop starts here
   
   for(p in 1:nrow(beta_mu_vals)){
   
   ###---Defining parameter values
   parsSIE    <- c(alpha = alpha,
                   m = m, 
                   gamma = gamma,
                   epsilon = epsilon,
                   tau = tau,
                   S0 = S0,
                   beta = beta_mu_vals[p, 1],
                   mu =   beta_mu_vals[p, 2])
 
   ###---Time scale
   times <- seq(0, 20, 1)

   ###---Solution of diff eqn model
   out <- as.data.frame(ode(func = SIEmod, y = c(S = S0, I = 0 ,  E =0.0000003),  
                            parms = parsSIE, times = times))

   names(out)[2:4] <- c("Susceptible", "Infectious", "Virus")

###---converting to long format for plots
 outL <- out %>% gather(key = Compartment, value = Number, -c(time)) %>% data.frame()

 outL$Compartment <- factor(outL$Compartment, levels = c("Susceptible", "Infectious", "Virus"))

 outL <-  bind_cols(outL, data.frame(t(matrix(parsSIE))))
 head(outL)
 names(outL)[4:ncol(outL)] <- c("alpha", "m", "gamma", "epsilon",
                                "tau", "S0", "beta", "mu")

 Output <- rbind(Output, outL)

 }


###focusing only on infectious group


     SIE_plot <- Output  %>%  filter(Compartment== "Infectious") %>%
        mutate(beta = as.character(round(beta, 6)),
               mu = as.character(round(mu, 3))) %>%
      ggplot(aes(x = time, y = Number, group = mu, color = mu)) +
   ylab("Number of infectious hosts") + 
     xlab("Time") + geom_line(size = 1) + facet_wrap(~ beta) +
     theme_bw() + theme(text = element_text(size = 20), 
                        strip.background = element_blank())
SIE_plot




