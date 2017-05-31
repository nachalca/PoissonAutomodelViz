
#' Shiny visualization app for Poisson auto-models
#'
#' @param sz size of the regular lattice
#' @param M overall mean of the covariate
#' @return shiny app
#' @useDynLib PoissonAutomodelViz
#' @importFrom Rcpp evalCpp
#' @export
#' @examples \dontrun{
#' launchApp()
#' }

launchApp <- function(sz = 15, M = 10) {

  # simulate the counts
  x <- rpois(sz*sz, M)
  X <- scale(x)

  # ---------------------
  # UI file
  # ---------------------
  ui <- shiny::shinyUI( shiny::bootstrapPage(
    #navbarPage("MSQ", id="barnavi",  theme = "bootstrap.min.css",
    shiny::titlePanel("Playing with Winsorized models"),
    shiny::fluidRow(
      shiny::sidebarLayout(
        shiny::sidebarPanel(width=3,
                     #          conditionalPanel(condition="input.conditionedPanels==1",
                     shiny::numericInput(inputId = "R",
                                  label = "Winsorizate:",
                                  min = 15, max = 500, value = 30, step = 5),
                     shiny::numericInput(inputId = "meanpar",label = "Intercept",min=1,max=4,step=1,value=2),
                     shiny::numericInput(inputId = "beta0",
                                  label = "Slope:",
                                  min = 0, max = 20, value = 10, step = .01),
                     shiny::numericInput(inputId = "beta1", label = "Marginal Expectation:",
                                  min = -1, max = 1, value = .1, step = .01),
                     shiny::numericInput(inputId = "etans",
                                  label = "NS Dependence :",
                                  min = -1, max = 1, value = 0, step = .01),
                     shiny::numericInput(inputId = "etaew",
                                  label = "EW Dependence :",
                                  min = -1, max = 1, value = 0, step = .01),
                     shiny::submitButton("Update View")
                     #                         ),
                     #           conditionalPanel(condition="input.conditionedPanels==2",
                     #                            numericInput(inputId = "R2",
                     #                                         label = "Winsorizate:",
                     #                                         min = 15, max = 50, value = 30, step = 5),
                     #                            numericInput(inputId = "beta2",
                     #                                         label = "Marginal Expectation:",
                     #                                         min = 2, max = 20, value = 10, step = 4)
                     #                           )
        ), # close sidebarPanel

        shiny::mainPanel(
          #         tabsetPanel(
          #             tabPanel('Tile plots',tableOutput("tab1"),plotOutput("plot1.ns"),value=1),
          #             tabPanel('Scatter plots',tableOutput("tab2"),plotOutput("plot2"),value=2),
          #               id="conditionedPanels" )

          shiny::h4('Expected Crashes'),
          shiny::tableOutput("tab1"),

          shiny::h4("Tile plot"),
          shiny::plotOutput("plot1.ns"),
          shiny::tableOutput("tab2"),
          shiny::h5('Fit a model ignoring dependence structure'),
          shiny::tableOutput("tab3"),
          shiny::tableOutput("tab4"),

          shiny::h4("Trace plots: mean and moran stat of every iteration"),
          shiny::plotOutput("plot2"),

          shiny::h4("Redction factor statistic"),
          shiny::tableOutput("tab5"),
          shiny::plotOutput("plot3")

        ) # close mainPanel
      ) # close sidebarLayout
    ) # close fluidRow
  ) # close fluidPage
  ) # close shinyUI
  #----------------------------

  # ---------------------
  # Server function
  # ---------------------
  server <- function(input, output,session) {

    Nraw <- cbind(id=1:(sz*sz), sqrgrid.4nbrs(sz))
    NNor  <- rbind(Nraw[seq(1,(sz*sz),2),],Nraw[seq(2,(sz*sz),2),])

    # data for plot 1
    ch1 <- shiny::reactive({ sim_mrfC(meanpar=input$meanpar,bet0=input$beta0, bet1=input$beta1, ens=input$etans, eew=input$etaew, R=input$R,iter=150, a = X, NN=Nraw) })
    ch2 <- shiny::reactive({ sim_mrfC(meanpar=input$meanpar,bet0=input$beta0, bet1=input$beta1, ens=input$etans, eew=input$etaew, R=input$R,iter=150, a = X, NN=Nraw) })
    ch3 <- shiny::reactive({ sim_mrfC(meanpar=input$meanpar,bet0=input$beta0, bet1=input$beta1, ens=input$etans, eew=input$etaew, R=input$R,iter=150, a = X, NN=Nraw) })
    ch4 <- shiny::reactive({ sim_mrfC(meanpar=input$meanpar,bet0=input$beta0, bet1=input$beta1, ens=input$etans, eew=input$etaew, R=input$R,iter=150, a = X, NN=Nraw) })

    d1.trace <- shiny::reactive({  rbind(data.frame(ch=1, ch1()[-1,]),data.frame(ch=2, ch2()[-1,]),data.frame(ch=3, ch3()[-1,]),data.frame(ch=4, ch4()[-1,]) )   })
    d1.sim <- shiny::reactive({  (ch1()[150,]+ch2()[150,]+ch3()[150,]+ch4()[150,])/4  })


    # d1.ns <- reactive({ sim_mrfC(meanpar=input$meanpar,bet0=input$beta0, bet1=input$beta1, ens=input$etans, eew=input$etaew, R=input$R,iter=500, a = X, NN=Nraw) })
    m1_ns <- shiny::reactive({ reshape2::melt( matrix(d1.sim(), ncol=sz, byrow=T) ) })
    mns  <- shiny::reactive({ apply(d1.trace()[,-1], 1, mean) })
    mrn  <- shiny::reactive({ t(apply(d1.trace()[,-1], 1, moran, NN=Nraw)) })
    mrn2  <- shiny::reactive({  moran(d1.sim() , NN=Nraw) })

    #cors <- reactive({ t(apply(d1.trace()[ ,-1], 1, function(z,y) cor( z,y), y=X)) })

    # ignore dependence models
    reg <- shiny::reactive({ glm(y ~ x, family=ifelse(input$meanpar<3,gaussian, poisson), data=data.frame(y=d1.sim(), x=X )) })

    # beta.ci <- reactive({ mns()[500] + c(-1.96, 1.96)*sqrt(mns()[500]/225) })

    d2.1 <- shiny::reactive({ data.frame(iter=rep(2:150,4),chain=rep(1:4,each=149), mean=mns(), moran.ns=mrn()[,1], moran.ew = mrn()[,2]) }) #, cors=cors()
    d2  <- shiny::reactive({ melt(d2.1(), id.vars=c('iter','chain') ) })

    # recomended bound for eta given R and beta
    eta.bnd <- shiny::reactive({  (log(input$R)-log(input$beta0))/(input$R - input$beta0)  })
    #    eta.bnd2 <-reactive({  (log(input$R2)-log(input$beta2))/(input$R2 - input$beta2)  })

    # Reduction factor statistic
    redufun <- function(x, chain) {
      n <- length(x)/4
      W <- mean( tapply(x, chain, var) )
      B <- n*var( tapply(x, chain, mean) )
      sqrt( (1-1/n) + (1/n)*(B/W) )
    }

    redfact <- shiny::reactive({ as.numeric(apply(d1.trace()[ , -1], 2, redufun, chain=d1.trace()[,1])) })

    # Descriptions
    mpdes <- c('Centered Parametrization: Linear effect of covariate, Multiplicative effect of Neighbors.',
               'Identity Link: Linear effect of covariate, Linear effect of Neighbors.',
               'Log-link: Multiplicative effect of covariate, Multiplicative effect of Neighbors.',
               'Multiplicative effect of covariate, Linear effect of Neighbors.')
    #==============================================
    output$plot1.ns <- renderPlot({
      print( ggplot2::qplot(data=m1_ns(), x=Var2, y=Var1, fill=value,geom='tile', main='NS dependence') +
               ggplot2::xlab('EW')+ ggplot2::ylab('NS')+
               ggplot2::scale_fill_gradient2(midpoint=ifelse(input$meanpar<3, input$beta0, exp(input$beta0))) )
    })

    output$tab1 <- shiny::renderTable({ data.frame(MeanPar = input$meanpar, description= mpdes[input$meanpar]) })

    output$tab2 <- shiny::renderTable({ data.frame(bnd = c(eta.bnd(), 3*max( X )), description=c('Recomended bound for eta parameter', '3*max(X): suggested R (Only relevant for meanpar=1)') )  })

    #output$tab3 <- renderTable({ data.frame(x='MLE ignoring dependence',est = summary(reg())$coefficients[,1] , confint(reg()) ) } , digits=4)
    output$tab3 <- shiny::renderTable({ cbind(summary( reg() )$coefficients, confint(reg()) )  })

    output$tab4 <- shiny::renderTable({ data.frame(x='Moran Stat', mrn.ns = mrn2()[1], mrn.ew = mrn2()[2]) }, digits=3)

    output$plot2 <- shiny::renderPlot({
      print( ggplot2::qplot(data=d2(), x=iter,y=value,color=as.factor(chain)) +
               ggplot2::geom_line() + ggplot2::facet_wrap(~variable, scale='free_y', ncol=3)  )
    })

    output$tab5 <- shiny::renderTable({ ( summary( data.frame(redfact())   )  ) })
    output$plot3 <- shiny::renderPlot({
      print( plot(redfact() ) )
    })
  }

  # ---------------------


  shiny::shinyApp(ui, server)
  }
