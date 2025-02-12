################################################################################
#
# SHINY APP: Generalised Dirichlet Distributions
#
# Implemented by Johann-Friedrich Salzmann (2025)
# 
# Reference: Ascari (2019). https://boa.unimib.it/handle/10281/241237
#
################################################################################

library(pacman)
p_load(shiny)
p_load(shinyjs)
p_load(ggtern)
p_load(gtools)
p_load(dplyr)
p_load(magrittr)
p_load(compositions)
p_load(ggcorrplot2)
p_load(reshape2)
p_load(ggforce)

# Fix for ggcorrplot2 positron connect package installation problem
#remotes::install_github("caijun/ggcorrplot2",force=T)
p_load(rlang)
p_load(dplyr)
source("ggcorrplot.R")
source("ggcorrplot.mixed.R")

# Flexible Dirichlet RNG
fdirichlet_rng = function(n, alpha, p, tau) {
  D = length(alpha)
  ni = rmultinom(1,n,p)[,1]
  sample = NULL
  for(i in 1:D){
    sample = rbind(sample, rdirichlet(ni[i], alpha + tau * diag(D)[i,]))
  }
  sample
}

# Flexible Dirichlet PDF
fdirichlet_pdf = function(n, alpha, p, tau) {
  D = length(alpha)
  mixture = 0
  for(i in 1:D){
    mixture = mixture + p[i] * rdirichlet(n, alpha + tau * diag(D)[i,])
  }
  mixture
}

# Extended Flexible Dirichlet RNG
efdirichlet_rng = function(n, alpha, p, tau) {
  D = length(alpha)
  ni = rmultinom(1,n,p)[,1]
  sample = NULL
  for(i in 1:D){
    sample = rbind(sample, rdirichlet(ni[i], alpha + tau[i] * diag(D)[i,]))
  }
  sample
}

# Extended Flexible Dirichlet PDF
efdirichlet_pdf = function(n, alpha, p, tau) {
  D = length(alpha)
  mixture = 0
  for(i in 1:D){
    mixture = mixture + p[i] * rdirichlet(n, alpha + tau[i] * diag(D)[i,])
  }
  mixture
}

# Double Flexible Dirichlet RNG
dfdirichlet_rng = function(n, alpha, p, tau) {
  D = length(alpha)
  p_flat = c(p)
  ni = rmultinom(1,n,p_flat)[,1]
  nij = matrix(ni,nrow=D)
  sample = NULL
  for(i in 1:D){
    for(j in 1:D){
      sample = rbind(sample, rdirichlet(nij[i,j], alpha + tau * (diag(D)[i,] + diag(D)[j,])))
    }
  }
  sample
}

# Double Flexible Dirichlet PDF
dfdirichlet_pdf = function(n, alpha, p, tau) {
  D = length(alpha)
  mixture = 0
  for(i in 1:D){
    for(j in 1:D){
      mixture = mixture + p[i,j] * rdirichlet(n, alpha + tau * (diag(D)[i,] + diag(D)[j,]))
    }
  }
  mixture
}

# Extended Double Flexible Dirichlet RNG
edfdirichlet_rng = function(n, alpha, p, tau) {
  D = length(alpha)
  p_flat = c(p)
  ni = rmultinom(1,n,p_flat)[,1]
  nij = matrix(ni,nrow=D)
  sample = NULL
  for(i in 1:D){
    for(j in 1:D){
      sample = rbind(sample, rdirichlet(nij[i,j], alpha + tau[i,j] * diag(D)[i,] + tau[j] * diag(D)[j,]))
    }
  }
  sample
}

# Extended Double Flexible Dirichlet PDF
edfdirichlet_pdf = function(n, alpha, p, tau) {
  D = length(alpha)
  mixture = 0
  for(i in 1:D){
    for(j in 1:D){
      mixture = mixture + p[i,j] * rdirichlet(n, alpha + tau[i,j] * diag(D)[i,] + tau[j] * diag(D)[j,])
    }
  }
  mixture
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- Filter(Negate(is.null), c(list(...), plotlist))
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

ui <- fluidPage(
  useShinyjs(),  # Enable shinyjs functionality
  titlePanel("Generalised Dirichlet Distributions"),
  titlePanel(tags$h5(
    "A Shiny App by",
    tags$a("Johann-Friedrich Salzmann",href="https://jfsalzmann.com",target="_blank"),
    "(2025)",
    br(),
    "Implements generalisations of the Dirichlet Class proposed in",
    tags$a("Ascari (2019)", href="https://boa.unimib.it/handle/10281/241237", target="_blank"),
    br(),br()
  )),
  
  sidebarLayout(
    sidebarPanel(
      actionButton("generate", "Generate"),
      selectInput("dist_type", "Choose Distribution:", choices = c("Dirichlet", "Flexible Dirichlet", "Extended Flexible Dirichlet",
                                                                   "Double Flexible Dirichlet", "Extended Double Flexible Dirichlet")),
      sliderInput("n", "Number of Samples:", min = 5, max = 10000, value = 1500),
      textInput("alpha1", "Alpha 1:", "3"),
      textInput("alpha2", "Alpha 2:", "3"),
      textInput("alpha3", "Alpha 3:", "3"),
      
      # tau values (only shown when applicable)
      textInput("tau1", "Tau (Tau 1):", "1"),
      textInput("tau2", "Tau 2:", "1"),
      textInput("tau3", "Tau 3:", "1"),
      
      # p values (only shown when applicable)
      sliderInput("p1", "p1:", min = 0, max = 1, value = 0.33),
      sliderInput("p2", "p2:", min = 0, max = 1, value = 0.33),
      sliderInput("p3", "p3:", min = 0, max = 1, value = 0.33),
      sliderInput("p11", "p11:", min = 0, max = 1, value = 0.11),
      sliderInput("p12", "p12:", min = 0, max = 1, value = 0.11),
      sliderInput("p13", "p13:", min = 0, max = 1, value = 0.11),
      sliderInput("p21", "p21:", min = 0, max = 1, value = 0.11),
      sliderInput("p22", "p22:", min = 0, max = 1, value = 0.11),
      sliderInput("p23", "p23:", min = 0, max = 1, value = 0.11),
      sliderInput("p31", "p31:", min = 0, max = 1, value = 0.11),
      sliderInput("p32", "p32:", min = 0, max = 1, value = 0.11),
      sliderInput("p33", "p33:", min = 0, max = 1, value = 0.11),
      
      textInput("tau11", "Tau 11:", "1"),
      textInput("tau12", "Tau 12:", "1"),
      textInput("tau13", "Tau 13:", "1"),
      textInput("tau21", "Tau 21:", "1"),
      textInput("tau22", "Tau 22:", "1"),
      textInput("tau23", "Tau 23:", "1"),
      textInput("tau31", "Tau 31:", "1"),
      textInput("tau32", "Tau 32:", "1"),
      textInput("tau33", "Tau 33:", "1")
    ),
    
    mainPanel(
      plotOutput("ternaryPlot", height = "auto")
    )
  )
)

server <- function(input, output, session) {
  
  # Dynamic UI: Enable/Disable Inputs Based on Distribution Selection
  observe({
    dist <- input$dist_type
    
    # Disable all by default
    shinyjs::disable("p1"); shinyjs::disable("p2"); shinyjs::disable("p3")
    shinyjs::disable("p11"); shinyjs::disable("p12"); shinyjs::disable("p13")
    shinyjs::disable("p21"); shinyjs::disable("p22"); shinyjs::disable("p23")
    shinyjs::disable("p31"); shinyjs::disable("p32"); shinyjs::disable("p33")
    shinyjs::disable("tau1"); shinyjs::disable("tau2"); shinyjs::disable("tau3")
    shinyjs::disable("tau11"); shinyjs::disable("tau21"); shinyjs::disable("tau31")
    shinyjs::disable("tau12"); shinyjs::disable("tau22"); shinyjs::disable("tau32")
    shinyjs::disable("tau13"); shinyjs::disable("tau23"); shinyjs::disable("tau33")
    
    if (dist == "Flexible Dirichlet" || dist == "Extended Flexible Dirichlet") {
      shinyjs::enable("p1"); shinyjs::enable("p2"); shinyjs::enable("p3")
    }
    
    if (dist == "Double Flexible Dirichlet" || dist == "Extended Double Flexible Dirichlet") {
      shinyjs::enable("p11"); shinyjs::enable("p12"); shinyjs::enable("p13")
      shinyjs::enable("p21"); shinyjs::enable("p22"); shinyjs::enable("p23")
      shinyjs::enable("p31"); shinyjs::enable("p32"); shinyjs::enable("p33")
    }
    
    if (dist == "Flexible Dirichlet" || dist == "Double Flexible Dirichlet") {
      shinyjs::enable("tau1")
    }
    
    if (dist == "Extended Flexible Dirichlet") {
      shinyjs::enable("tau1"); shinyjs::enable("tau2"); shinyjs::enable("tau3")
    }
    
    if (dist == "Extended Double Flexible Dirichlet") {
      shinyjs::enable("tau11"); shinyjs::enable("tau21"); shinyjs::enable("tau31")
      shinyjs::enable("tau12"); shinyjs::enable("tau22"); shinyjs::enable("tau32")
      shinyjs::enable("tau13"); shinyjs::enable("tau23"); shinyjs::enable("tau33")
    }
  })
  
  output$ternaryPlot <- renderPlot({
    input$generate  # Reactivity trigger
    
    # Convert input to numeric
    alpha <- as.numeric(c(input$alpha1, input$alpha2, input$alpha3))
    p <- matrix(NA, nrow = 3, ncol = 3)
    p[1,] <- as.numeric(c(input$p11, input$p12, input$p13))
    p[2,] <- as.numeric(c(input$p21, input$p22, input$p23))
    p[3,] <- as.numeric(c(input$p31, input$p32, input$p33))
    p_vector <- as.numeric(c(input$p1, input$p2, input$p3))
    tau <- as.numeric(c(input$tau1, input$tau2, input$tau3))
    tau_matrix <- matrix(NA, nrow = 3, ncol = 3)
    tau_matrix[1,] <- as.numeric(c(input$tau11, input$tau12, input$tau13))
    tau_matrix[2,] <- as.numeric(c(input$tau21, input$tau22, input$tau23))
    tau_matrix[3,] <- as.numeric(c(input$tau31, input$tau32, input$tau33))
    
    # Choose function based on input
    samples <- switch(input$dist_type,
                      "Dirichlet" = rdirichlet(input$n, alpha),
                      "Flexible Dirichlet" = fdirichlet_rng(input$n, alpha, p_vector, tau[1]),
                      "Extended Flexible Dirichlet" = efdirichlet_rng(input$n, alpha, p_vector, tau),
                      "Double Flexible Dirichlet" = dfdirichlet_rng(input$n, alpha, p, tau[1]),
                      "Extended Double Flexible Dirichlet" = edfdirichlet_rng(input$n, alpha, p, tau_matrix))
    samples_mixed <- switch(input$dist_type,
                     "Flexible Dirichlet" = fdirichlet_pdf(input$n, alpha, p_vector, tau[1]),
                     "Extended Flexible Dirichlet" = efdirichlet_pdf(input$n, alpha, p_vector, tau),
                     "Double Flexible Dirichlet" = dfdirichlet_pdf(input$n, alpha, p, tau[1]),
                     "Extended Double Flexible Dirichlet" = edfdirichlet_pdf(input$n, alpha, p, tau_matrix))
    
    samples <- as.data.frame(samples)
    colnames(samples) <- c("V1", "V2", "V3")
    
    if(input$dist_type != "Dirichlet"){
      samples_mixed <- as.data.frame(samples_mixed)
      colnames(samples_mixed) <- c("V1", "V2", "V3")
      cluster = " Cluster"
    } else {
      cluster = ""
    }
    
    samples_clr <- as_tibble(clr(samples))
    samples_pwlr <- as_tibble(pwlr(samples))
    
    cor_pwlr <- cor(samples_pwlr)
    cor_clr <- cor(samples_clr)
    
    p1 = ggtern(samples, aes(x=V1, y=V3, z=V2)) + 
      geom_point(size=.5, alpha = .8) + 
      geom_density_tern(bdl=.05) +
      stat_density_tern(aes(fill=..level.., alpha=..level..),geom='polygon',bdl=.05) +#now you need to use stat_density_tern
      scale_fill_gradient2(high = "red") +                                    #define the fill color
      guides(color = "none", fill = "none", alpha = "none")  +
      labs(title = paste0("Ternary Plot of ", input$dist_type, cluster, " Realisations")) +
      theme_ggtern()
    
    if(input$dist_type != "Dirichlet"){
      p2 = ggtern(samples_mixed, aes(x=V1, y=V3, z=V2)) + 
        geom_point(size=.5, alpha = .8) + 
        geom_density_tern(bdl=.05) +
        stat_density_tern(aes(fill=..level.., alpha=..level..),geom='polygon',bdl=.05) +#now you need to use stat_density_tern
        scale_fill_gradient2(high = "red") +                                    #define the fill color
        guides(color = "none", fill = "none", alpha = "none")  +
        labs(title = paste0("Ternary Plot of ", input$dist_type, " Mixture Density Realisations")) +
        theme_ggtern()
    } else {
      p2 = NULL
    }
    
    p31 = ggcorrplot.mixed(cor_pwlr,upper="ellipse") + labs(title = "PWLR Corr Plot") + theme(panel.background = element_rect(fill = "grey98"))
    p32 = ggcorrplot.mixed(cor_clr) + labs(title = "CLR Corr Plot") + theme(panel.background = element_rect(fill = "grey98"))
    
    multiplot(p1,p2,p31,p32,cols = 1)
  }, height = 1400,width=700)
}

shinyApp(ui, server)
