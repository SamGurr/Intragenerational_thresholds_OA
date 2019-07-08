library(shiny)

server <- shinyServer(function(input, output,session){
  
  vals <- reactiveValues()
  observe({
    vals$x                   <- input$tank_volume*1000 # convert to cells per milliliter
    vals$t                   <- input$target_cell_density # in cells per milliliter
    vals$q                   <- input$how_many_tanks
    vals$f                   <- input$flow_rate*1000 # convert to milliliters
    vals$total_cells         <- vals$x*vals$t
    vals$vol                 <- input$feed_conical
    vals$algae_mix_conc      <- (input$algae_1*(input$perc_algae1/100)) + 
      (input$algae_2*(input$perc_algae2/100)) + 
      (input$algae_3*(input$perc_algae3/100)) # total cell concentration of mixed diet
    vals$assumed_hourly_loss <- (vals$f*60)*vals$t  # in cells per milliliter
    vals$percent_loss_hourly <- (vals$assumed_hourly_loss/vals$total_cells)*100 
    vals$Liters_algae_initial <- ((((vals$total_cells*(input$perc_algae1/100))/input$algae_1)/1000)+
                                    (((vals$total_cells*(input$perc_algae2/100))/input$algae_2)/1000)+
                                    (((vals$total_cells*(input$perc_algae3/100))/input$algae_3)/1000))
    vals$LPH_feed <- ((vals$percent_loss_hourly/100)*vals$Liters_algae_initial)
    
  }) 
  
  output$signature <- renderText({
    paste("App created by Sam J. Gurr - last update 02/04/2019") 
    
  })
  
  output$contactinfo <- renderText({
    paste("Contact: samjgurr@gmail.com") 
    
  })
  
  output$batchtotal <- renderText({
    paste("BATCH FEED: total volume feed per tank (L) =", 
          vals$Liters_algae_initial)
  })
  
  output$batch_algae1 <- renderText({
    paste("BATCH FEED: algae 1 (L per tank) =",
          ((vals$total_cells*(input$perc_algae1/100))/input$algae_1)/1000)
  })
  
  output$batch_algae2 <- renderText({
    paste("BATCH FEED: algae 2 (L per tank) =",
          ((vals$total_cells*(input$perc_algae2/100))/input$algae_2)/1000)
  })
  
  output$batch_algae3 <- renderText({
    paste( "BATCH FEED: algae 3 (L per tank) =",
           ((vals$total_cells*(input$perc_algae3/100))/input$algae_3)/1000)
  })
  
  output$hourly_flow_loss <- renderText({
    paste("Hourly loss from flow through (cells mL) =",
          vals$assumed_hourly_loss)
  })
  
  output$replenish <- renderText({
    paste("Volume to replenish conical/tank-1  hour-1 (L) =",
          (vals$assumed_hourly_loss/vals$algae_mix_conc)/1000) 
  })
  
  output$total_continuous <- renderText({
    paste("CONTINUOUS FEED: total volume mixed algae (L) =",
          ((vals$percent_loss_hourly/100)*vals$Liters_algae_initial)*24*vals$q)
  })
  
  output$cont_algae1 <- renderText({
    paste("CONTINUOUS FEED: algae 1 (L) =",
          (((((vals$percent_loss_hourly/100)*vals$total_cells)*24*vals$q)*(input$perc_algae1/100))/(input$algae_1))/1000)
  })
  
  output$cont_algae2 <- renderText({
    paste("CONTINUOUS FEED: algae 2 (L) =",
          (((((vals$percent_loss_hourly/100)*vals$total_cells)*24*vals$q)*(input$perc_algae2/100))/(input$algae_2))/1000)
  })
  
  output$cont_algae3 <- renderText({
    paste("CONTINUOUS FEED: algae 3 (L) =",
          (((((vals$percent_loss_hourly/100)*vals$total_cells)*24*vals$q)*(input$perc_algae3/100))/(input$algae_3))/1000)
  })
  
  output$dilution <- renderText({
    paste("DILUTION: FILTERED SEAWATER (L) =", 
          ((input$feed_conical/1000)-((vals$assumed_hourly_loss/vals$algae_mix_conc)/1000))*24*vals$q) 
    
  })
  
})



