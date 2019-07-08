library(shiny)

ui <- shinyUI(fluidPage(
  mainPanel(
    numericInput("tank_volume", "tank_volume (L)", 250 ,min = 0, max = 5000),
    numericInput("how_many_tanks", "how_many_tanks of the same volume", 3,min = 1, max = 10),
    numericInput("flow_rate", "flow_rate (LPM)", 0.4,min = 0.00, max = 5.00),
    numericInput("target_cell_density", "target_cell_density (cells per mL)", 25000,min = 0, max = 1000000),
    numericInput("feed_conical", "target mL per hour continuous feed", 1000 ,min = 0, max = 1000),
    numericInput("algae_1", "Algae #1 cells mL-1 (pav, iso, tet, nano, etc.)", 450000,min = 0, max = 2000000),
    numericInput("perc_algae1", "Percent_comp_1", 100,min = 0, max = 100),
    numericInput("algae_2", "Algae #2 cells mL-1  (pav, iso, tet, nano, etc.)", 2000000,min = 0, max = 2000000),
    numericInput("perc_algae2", "Percent_comp_2", 0,min = 0, max = 100),
    numericInput("algae_3", "Algae #3 cells mL-1  (pav, iso, tet, nano, etc.)", 2000000,min = 0, max = 2000000),
    numericInput("perc_algae3", "Percent_comp_3", 0,min = 0, max = 100),
    textOutput("signature"),
    textOutput("contactinfo"),
    textOutput("batchtotal"),
    textOutput("batch_algae1"),
    textOutput("batch_algae2"),
    textOutput("batch_algae3"),
    textOutput("hourly_flow_loss"),
    textOutput("replenish"),
    textOutput("total_continuous"),
    textOutput("cont_algae1"),
    textOutput("cont_algae2"),
    textOutput("cont_algae3"),
    textOutput("dilution"))
))