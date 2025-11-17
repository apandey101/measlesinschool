# ==============================================================================
# MEASLES TRANSMISSION MODEL - SHINY APP 
# ==============================================================================
# 
# Interactive web application for exploring measles transmission dynamics
# in school settings with class-based mixing
#
# Updated on 11/10/2025
# ==============================================================================

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(shinycssloaders)
library(plotly)
library(readxl)
library(shinyBS)
source("sim_cr_v1.R")

# ==============================================================================
# LOAD SCHOOL DATA
# ==============================================================================
# Load your school vaccination data
school_data_raw <- read_xlsx("SC_vaccination_2024_2025.xlsx")
# FILTER OUT SCHOOLS WITH MISSING DATA
school_data <- school_data_raw %>%
  filter(
    !is.na(`Total.Students`),
    !is.na(`Percent.Immunized`),
    !is.na(Region),
    !is.na(County),
    !is.na(`School.Name`),
    `Total.Students` > 0,
    `Percent.Immunized` >= 0,
    `Percent.Immunized` <= 100
  ) %>%
  mutate(Percent.Immunized=Percent.Immunized*100)
school_data$Total.Students <- as.numeric(school_data$Total.Students)

# Check if we have data after filtering
if (nrow(school_data) == 0) {
  stop("No valid school data found after filtering. Please check your CSV file for missing values.")
}

# Fixed class size  
FIXED_CLASS_SIZE <- 25



# ==============================================================================
# UI (User Interface)
# ==============================================================================

ui <- fluidPage(
  
  # Application title
  titlePanel(
    div(
      h2("Measles Transmission Model: School-Based Simulation for South Carolina"),
      h4("Interactive Agent-Based Model", style = "color: gray;")
    )
  ),
  
  # Sidebar with input controls
  sidebarLayout(
    
    # SIDEBAR PANEL - Input Parameters
    sidebarPanel(
      width = 3,
      
      h4("🏫 School Selection"),
      
      selectInput("region",
                  "Select Region:",
                  choices = c("Choose a region..." = "", sort(unique(school_data$Region)))),
      
      selectInput("county",
                  "Select County:",
                  choices = c("First select a region..." = "")),
      
      selectInput("school",
                  "Select School:",
                  choices = c("First select a county..." = "")),
      
      # Display selected school info
      conditionalPanel(
        condition = "input.school != ''",
        wellPanel(
          style = "background-color: #e8f5e9;",
          h5(strong("📋 Selected School"), style = "margin-top: 0; color: #2e7d32;"),
          htmlOutput("school_info_display")
        )
      ),
      
      hr(),
      sliderInput("initial_infected",
                  "Initial Infected Students:",
                  min = 1,
                  max = 10,
                  value = 1,
                  step = 1),
      hr(),
      h4("🚨 Quarantine & Isolation Policy"),
      sliderInput(
        "isolation_delay",
        tagList(
          "Isolation Delay (days):",
          tags$span(
            id = "isolation_delay_help",
            icon("info-circle"),
            class = "text-muted",
            style = "margin-left:6px; cursor: help;"
          )
        ),
        min = 0,
        max = 4,
        value = 1,
        step = 1
      ),
      
      bsTooltip(
        "isolation_delay_help",
        "This delay applies to the first (index) case after rash development. All subsequent cases are isolated with delay upon prodromal symptoms.",
        placement = "right",
        trigger = "hover"
      ),
      
      sliderInput(
        "quarantine_efficacy",
        tagList(
          "Quarantine Efficacy (%)",
          tags$span(
            id = "qe_help",
            icon("info-circle"),
            class = "text-muted",
            style = "margin-left:6px; cursor: help;"
          )
        ),
        min = 0, max = 100, value = 80, step = 10
      ),
      
      bsTooltip(
        "qe_help",
        "Probability (in %) that each unvaccinated classmate of a symptomatic case is successfully quarantined. Vaccinated classmates remain in school.",
        placement = "right",
        trigger = "hover"
      ),
      
      sliderInput("quarantine_duration",
                  "Quarantine Duration (days):",
                  min = 7,
                  max = 28,
                  value = 21,
                  step = 7),
      
      hr(),
      
      
      # Run button
      actionButton("run_sim",
                   "🚀 Run Simulation",
                   class = "btn-primary btn-lg btn-block",
                   style = "margin-top: 10px; margin-bottom: 10px;"),
      
      # Download button
      downloadButton("download_data",
                     "📥 Download Results",
                     class = "btn-success btn-block"),
      
      hr(),
      
      # ==============================================================================
      # ADVANCED OPTIONS - COLLAPSIBLE SECTION
      # ==============================================================================
      tags$details(
        tags$summary(
          style = "cursor: pointer;
                   font-weight: bold;
                   font-size: 16px;
                   color: #1976d2;
                   padding: 12px;
                   background-color: #f5f5f5;
                   border-radius: 4px;
                   border: 2px dashed #1976d2;
                   user-select: none;
                   margin-bottom: 5px;",
          "⚙️ Advanced Options (Click to expand)"
        ),
        
        div(
          style = "padding: 15px;
                   background-color: #fafafa;
                   border: 1px solid #e0e0e0;
                   border-radius: 4px;
                   margin-top: 10px;
                   margin-bottom: 10px;",
          
          # Warning header
          h5(
            strong("⚠️ Override Default Parameters"),
            style = "color: #d32f2f; margin-top: 0; margin-bottom: 10px;"
          ),
          
          # Explanation text
          p(
            style = "font-size: 13px; color: #666; line-height: 1.5; margin-bottom: 15px;",
            "Expand this section only if you want to test different scenarios."
          ),
          
          # Section 1: Transmission Parameters
          div(
            style = "background-color: #fff;
                     padding: 12px;
                     border-radius: 4px;
                     margin-bottom: 15px;
                     border-left: 4px solid #2196F3;",
            h4("🎲 Simulation Settings"),
            
            sliderInput("n_simulations",
                        "Number of Simulations:",
                        min = 10,
                        max = 500,
                        value = 150,
                        step = 10),
            
            sliderInput("n_days",
                        "Maximum Simulation Days:",
                        min = 60,
                        max = 180,
                        value = 100,
                        step = 20),
            
            hr(),
            
            h6(strong("Transmission Parameters"),
               style = "margin-top: 0; color: #1976d2;"),
            
            # Transmission Parameters (replace the old beta sliders)
            sliderInput("c_within",  "Within-class contacts per day:",  min = 0, max = 20, value = 10, step = 1),
            sliderInput("c_between", "Between-class contacts per day:", min = 0, max = 10, value = 5,  step = 1),
            
            sliderInput("p_within",  "Per-contact transmission (within):",  min = 0, max = 1, value = 0.19, step = 0.01),
            sliderInput("p_between", "Per-contact transmission (between):", min = 0, max = 1, value = 0.09, step = 0.01),
            
            tags$small(
              style = "color:#666;",
              "Default parameters leads to R0 of about 15."
            ))
            ,
          
          # Section 2: Vaccination Coverage Override
          div(
            style = "background-color: #fff;
                     padding: 12px;
                     border-radius: 4px;
                     border-left: 4px solid #FF9800;",
            
            h6(strong("Vaccination Coverage Override"),
               style = "margin-top: 0; color: #ef6c00;"),
            
            sliderInput("vaccination_coverage",
                        "Custom Vaccination Coverage (%):",
                        min = 0,
                        max = 100,
                        value = 90,  # Will be updated dynamically when school selected
                        step = 1),
            
            checkboxInput("use_custom_vaccination",
                          "Use custom vaccination coverage (uncheck to use school's rate)",
                          value = FALSE),
            
            tags$small(
              style = "color: #666;",
              "Unchecking this box will use the school's vaccination rate displayed above."
            )
          )
        )
      )
    ),  # END OF SIDEBAR PANEL
    
    # MAIN PANEL - Output
    mainPanel(
      width = 9,
      
      # Tabs for different visualizations
      tabsetPanel(
        type = "tabs",
        
        # =======================
        # Tab 1: Main Results
        # =======================
        tabPanel(
          "📊 Main Results",
          br(),
          
          # Caution banner
          tags$div(
            class = "alert alert-warning",
            style = "background-color: #fff8e1; border-left: 4px solid #ffb300; padding: 12px; border-radius: 4px;",
            h5(
              icon("exclamation-triangle"),
              strong(" Ongoing Work — Interpret with Caution"),
              style = "margin-top: 0; color: #8a6d3b;"
            ),
            p("Results depend on parameter choices and policy assumptions, including isolation timing, quarantine coverage (unvaccinated-only policy), vaccination levels, and school structure. Small changes can lead to large shifts in outcomes.")
          ),
          
          # Summary statistics boxes
          fluidRow(
            column(
              4,
              wellPanel(
                style = "background-color: #e3f2fd;",
                h4(
                  textOutput("attack_rate_text"),
                  style = "color: #1976d2; margin: 0;"
                ),
                tags$p(
                  "Final Outbreak Size",
                  tags$span(
                    id = "attack_rate_help",
                    icon("info-circle"),
                    class = "text-muted",
                    style = "margin-left: 6px; cursor: help;"
                  ),
                  style = "margin: 5px 0 0 0; font-size: 14px;"
                )
              )
            ),
            column(
              4,
              wellPanel(
                style = "background-color: #ffebee;",
                h4(
                  textOutput("total_infected_text"),
                  style = "color: #c62828; margin: 0;"
                ),
                p("Total Infected (Median)", style = "margin: 5px 0 0 0; font-size: 14px;")
              )
            ),
            column(
              4,
              wellPanel(
                style = "background-color: #f3e5f5;",
                h4(
                  textOutput("duration_text"),
                  style = "color: #6a1b9a; margin: 0;"
                ),
                tags$p(
                  "Median Outbreak Duration (days)",
                  tags$span(
                    id = "duration_help",
                    icon("info-circle"),
                    class = "text-muted",
                    style = "margin-left: 6px; cursor: help;"
                  ),
                  style = "margin: 5px 0 0 0; font-size: 14px;"
                )
              )
            )
          ),
          
          # Tooltips for stats
          bsTooltip(
            id = "attack_rate_help",
            title = "Percentage of the total student population that became infected over the course of the outbreak/simulation window (median across all simulations).",
            placement = "right",
            trigger = "hover"
          ),
          bsTooltip(
            id = "duration_help",
            title = "If the median duration equals the simulation horizon, it may be censored. Increase the simulation window to estimate the true median outbreak duration.",
            placement = "right",
            trigger = "hover"
          ),
          
          # Main epidemic curve plot
          withSpinner(
            plotlyOutput("epidemic_curve_plot", height = "600px"),
            type = 6,
            color = "#1976d2"
          ),
          
          br(),
          
          wellPanel(
            style = "background-color: #e8f5e9; border-left: 4px solid #4caf50;",
            h5(
              icon("lightbulb"), strong(" Interpreting These Results"),
              style = "color: #2e7d32; margin-top: 0;"
            ),
            p(
              "These simulations explore ", strong("outbreak uncertainty"),
              ". Identical conditions can produce very different outcomes:",
              style = "color: #2e7d32; margin-bottom: 8px;"
            ),
            tags$ul(
              style = "color: #2e7d32;",
              tags$li(
                strong("The median shows what typically happens"),
                " - the middle outcome when all scenarios are ranked."
              ),
              tags$li(
                strong("95% intervals show likely variation"),
                " - most outbreaks fall within this range."
              ),
              tags$li(
                strong("Outliers occur"),
                " - about 2.5% of outbreaks are worse, 2.5% are better."
              )
            ),
            div(
              style = "background-color: #c8e6c9; padding: 10px; border-radius: 4px; margin-top: 8px;",
              p(
                icon("check-circle"), strong(" Policy Implication:"),
                " Plan for both typical scenarios (median) and worst-case outcomes (upper interval).",
                style = "color: #1b5e20; margin: 0; font-size: 13px;"
              )
            )
          )
        ),
        
        # =======================
        # Tab 2: Infectious States
        # =======================
        tabPanel(
          "🦠 Infectious States",
          br(),
          withSpinner(
            plotOutput("infectious_plot", height = "500px"),
            type = 6,
            color = "#c62828"
          ),
          br(),
          wellPanel(
            h4("About This Plot"),
            p("This plot shows the number of infectious students still attending school over time, combining both
               the Prodromal and Rash phases into a single 'Active Infectious' category. It reflects the individuals
               capable of transmitting infection to others (i.e., not isolated or quarantined)."),
            tags$ul(
              tags$li(strong("Red line:"), "Median number of active infectious students across all simulations."),
              tags$li(strong("Shaded area:"), "95% uncertainty interval around the median."),
              tags$li(strong("Gray dashed lines:"), "Day and magnitude of the peak number of infectious students."),
              tags$li(strong("Purple dotted line:"), "Median outbreak end — the day when 50% of simulated outbreaks have ended.")
            ),
            p(
              strong("Important Note:"),
              " The curve represents only those who remain infectious while present in school settings, not those already isolated or quarantined.
                This is the key measure of potential in-school transmission."
            )
          )
        ),
        
        # =======================
        # Tab 3: Attack Rate Distribution
        # =======================
        tabPanel(
          "📈 Final Outbreak Size Distribution",
          br(),
          
          # Definition banner (tidy + PH friendly)
          tags$div(
            class = "alert alert-info",
            style = "background-color: #e3f2fd; border-left: 4px solid #1976d2; padding: 12px; border-radius: 4px;",
            h5(
              icon("info-circle"), strong(" What is the Final Outbreak Size?"),
              style = "margin-top: 0; color: #0d47a1;"
            ),
            tags$p(
              "Percentage of the total student population that becomes infected by the end of the simulation window.",
              style = "margin-bottom: 0;"
            )
          ),
          
          withSpinner(
            plotOutput("attack_rate_dist_plot", height = "500px"),
            type = 6,
            color = "#2e7d32"
          ),
          
          br(),
          
          wellPanel(
            h4("How to Use This Plot"),
            p("This histogram shows the distribution of final outbreak sizes across all simulations (each bar summarizes many stochastic runs)."),
            tags$ul(
              tags$li(strong("Lower bars to the right:"), " scenarios tend to produce higher overall spread."),
              tags$li(strong("Narrow distribution:"), " outcomes are consistent; policies/parameters dominate over randomness."),
              tags$li(strong("Wide or bimodal:"), " near-threshold behavior — some outbreaks fade, others take off.")
            ),
            tags$p(
              "Lines: ",
              tags$span(strong("Blue solid = median"), " (typical outcome); "),
              tags$span(strong("Orange dashed = mean"), " (sensitive to rare, very large outbreaks).")
            ),
            div(
              style = "background-color: #e8f5e9; padding: 10px; border-radius: 4px; border-left: 4px solid #4caf50; margin-top: 8px;",
              p(
                icon("stethoscope"), strong(" Public Health Takeaway:"),
                " Compare scenarios (e.g., quarantine coverage, isolation timing, vaccination) by how much they shift the median and shrink the spread of outcomes. Plan capacity using the upper tail (worse-but-plausible runs).",
                style = "margin: 0;"
              )
            )
          )
        )
        ,
        
        # =======================
        # Tab 5: About
        # =======================
        tabPanel(
          "ℹ️ About",
          br(),
          h3("Measles Transmission Model"),
          
          h4("Model Description"),
          p("This agent-based model represents measles transmission within a school community. Each student is modeled as an individual with their own vaccination status, infection stage, and classroom assignment. Transmission occurs primarily within classrooms, with a lower rate of between-class mixing. The model focuses on how vaccination coverage, timing of isolation for the index case, and the effectiveness of contact tracing and quarantine of unvaccinated contacts shape outbreak dynamics."),
          
          p("The model runs in daily time steps. When a student becomes infected they progress through a latent (incubation) phase before entering the prodromal (pre-rash infectious) phase and then develop a rash. For the first (index) case, once a rash appears the individual is isolated after a user-specified delay, preventing further school-based transmission. For all subsequent cases, isolation is triggered after the onset of the prodromal phase with user-defined delay."),
          
          p("Quarantine is applied only to unvaccinated classmates identified as exposed to a symptomatic case. Vaccinated classmates are not quarantined in this scenario. Quarantined students are removed from school mixing for the defined quarantine period, then return to their prior immune or susceptible status."),
          
          p(
            strong("Important:"),
            " This model focuses exclusively on school-based transmission and does not model household, community, or social-network spread outside class and school. It is designed to explore the impact of school-based control measures such as isolation timing, quarantine coverage, and vaccination thresholds."
          ),
          
          h4("Key Features"),
          tags$ul(
            tags$li(strong("Two levels of contact:"), "Higher exposure risk within classrooms; lower between-class mixing."),
            tags$li(strong("Vaccination effects:"), "Vaccinated students are less susceptible and, if infected, less infectious to others."),
            tags$li(strong("Targeted interventions:"), "Index case isolation after rash + delay, isolation of subsequent symptomatic students (prodormal onset + delay), and quarantine of unvaccinated exposed classmates only."),
            tags$li(strong("Stochastic simulation:"), "Each simulation run represents a possible outbreak trajectory, enabling exploration of uncertainty in outcomes.")
          ),
          
          h4("Interpretation"),
          p("This model is intended for comparative scenario analysis rather than exact prediction. It allows users to explore how small changes — such as a one-day delay in index case isolation, a drop in vaccine coverage, or incomplete quarantine of contacts — can dramatically alter outbreak size and duration. Running multiple simulations provides a view of the range of possible outcomes and the stochastic variability inherent in outbreak dynamics."),
          
          h4("References & Further Resources"),
          p(
            strong("Measles disease information:"),
            tags$a("US CDC – Measles", href = "https://www.cdc.gov/pinkbook/hcp/table-of-contents/chapter-13-measles.html?CDC_AAref_Val=https://www.cdc.gov/vaccines/pubs/pinkbook/meas.html", target = "_blank"),
            ", ",
            tags$a("ECDC Measles Factsheet", href = "https://www.ecdc.europa.eu/en/measles/facts", target = "_blank")
          ),
          p(
            strong("Other Models:"),
            tags$a("epiENGAGE Measles Outbreak Simulator", href = "https://epiengage-measles.tacc.utexas.edu/", target = "_blank"),
            ", ",
            tags$a("epiworldR Shiny", href = "https://gvegayon.shinyapps.io/epiworldRShiny/", target = "_blank")
          ),
          
          br(),
          hr(),
          p("Model Version 1.0 | Created 2025 | Updated on 11/11/2025", style = "text-align: center; color: gray;")
        )
      ) # end tabsetPanel
    )  # end mainPanel
    
  ) # end sidebarLayout
)   # end fluidPage



# ==============================================================================
# SERVER
# ==============================================================================

server <- function(input, output, session) {
  
  # Reactive values to store results
  results <- reactiveVal(NULL)
  
  # ==============================================================================
  # REACTIVE LOGIC FOR SCHOOL SELECTION
  # ==============================================================================
  
  # Create reactive values for school data
  selected_school_data <- reactiveVal(NULL)
  
  # Update county dropdown when region changes
  observeEvent(input$region, {
    req(input$region != "")
    
    counties <- school_data %>%
      filter(Region == input$region) %>%
      pull(County) %>%
      unique() %>%
      sort()
    
    updateSelectInput(session, "county",
                      choices = c("Select a county..." = "", counties))
    
    # Reset school selection
    updateSelectInput(session, "school",
                      choices = c("First select a county..." = ""))
    selected_school_data(NULL)
  })
  
  # Update school dropdown when county changes
  observeEvent(input$county, {
    req(input$county != "")
    
    schools <- school_data %>%
      filter(Region == input$region, County == input$county) %>%
      pull(`School.Name`) %>%
      sort()
    
    updateSelectInput(session, "school",
                      choices = c("Select a school..." = "", schools))
    
    selected_school_data(NULL)
  })
  
  # Update school data when school is selected
  observeEvent(input$school, {
    req(input$school != "")
    
    school_info <- school_data %>%
      filter(Region == input$region, 
             County == input$county,
             `School.Name` == input$school)
    
    if (nrow(school_info) == 0) {
      showNotification(
        "Error: Selected school not found in database.",
        type = "error",
        duration = 5
      )
      return()
    }
    
    selected_school_data(school_info)
    
    # UPDATE VACCINATION COVERAGE SLIDER TO MATCH SCHOOL'S RATE
    updateSliderInput(session, "vaccination_coverage",
                      value = school_info$`Percent.Immunized`)
    
    # RESET THE CHECKBOX TO USE SCHOOL'S RATE (NOT CUSTOM)
    updateCheckboxInput(session, "use_custom_vaccination", 
                        value = FALSE)
  })
  
  # ==============================================================================
  # REACTIVE VALUE FOR VACCINATION COVERAGE
  # ==============================================================================
  
  # This determines which vaccination rate to use: school's or custom
  vaccination_rate <- reactive({
    req(selected_school_data())
    
    if (isTRUE(input$use_custom_vaccination)) {
      # User checked the box - use custom rate from slider
      return(input$vaccination_coverage)
    } else {
      # User unchecked (or never checked) - use school's rate
      return(selected_school_data()$`Percent.Immunized`)
    }
  })
  
  # Display school information
  output$school_info_display <- renderUI({
    req(selected_school_data())
    
    info <- selected_school_data()
    
    HTML(paste0(
      "<strong>School:</strong> ", info$`School.Name`, "<br/>",
      "<strong>Region:</strong> ", info$Region, "<br/>",
      "<strong>County:</strong> ", info$County, "<br/>",
      "<strong>Students:</strong> ", format(info$`Total.Students`, big.mark = ","), "<br/>",
      "<strong>Vaccination Rate:</strong> ", info$`Percent.Immunized`, "%"
    ))
  })
  
  output$school_size_display <- renderText({
    req(selected_school_data())
    paste(format(selected_school_data()$`Total.Students`, big.mark = ","), "students")
  })
  
  output$vaccination_coverage_display <- renderText({
    req(selected_school_data())
    paste0(selected_school_data()$`Percent.Immunized`, "%")
  })
  
  # Run simulation when button is clicked
  observeEvent(input$run_sim, {
    
    # Validate that a school is selected
    if (is.null(selected_school_data())) {
      showNotification(
        "Please select a school before running simulation.",
        type = "error",
        duration = 5
      )
      return()
    }
    
    # Validate school data
    school_info <- selected_school_data()
    if (is.na(school_info$`Total.Students`) || school_info$`Total.Students` <= 0) {
      showNotification(
        "Error: Selected school has invalid student count.",
        type = "error",
        duration = 5
      )
      return()
    }
    
    if (is.na(school_info$`Percent.Immunized`)) {
      showNotification(
        "Error: Selected school has missing vaccination data.",
        type = "error",
        duration = 5
      )
      return()
    }
    
    # Show notification
    showNotification(
      "Running simulations... This may take 30-60 seconds.",
      duration = NULL,
      id = "sim_notification",
      type = "message"
    )
    
    # Create parameter list from inputs
    sim_params <- list(
      n_simulations = input$n_simulations,
      school_size   = as.integer(school_info$`Total.Students`),
      avg_class_size = FIXED_CLASS_SIZE,
      age_range     = c(5, 17),
      vaccination_coverage = vaccination_rate() / 100,
      
      # Disease parameters
      latent_period = 10,
      prodromal_period = 4,
      rash_period = 4,
      isolation_delay = input$isolation_delay,
      isolation_period = 3,
      vaccine_efficacy = 0.97,
      vaccine_infectiousness_reduction = 0.80,
      prodromal_infectiousness_multiplier = 1.0,
      rash_infectiousness_multiplier     = 0.6,
      
      # ---- NEW: Contact-rate transmission ----
      c_within  = input$c_within,    # e.g., 8–12/day
      c_between = input$c_between,   # e.g., 1–3/day
      p_within  = input$p_within,    # per-contact probability
      p_between = input$p_between,
      
      # Quarantine
      quarantine_contacts = TRUE,
      quarantine_efficacy = input$quarantine_efficacy / 100,
      quarantine_duration = input$quarantine_duration,
      
      # Simulation
      initial_infected = input$initial_infected,
      n_days = input$n_days,
      seed_start = 1234567,
      verbose = FALSE
    )
    
    
    # Run the simulation
    sim_results <- tryCatch({
      do.call(run_multiple_simulations, sim_params)
    }, error = function(e) {
      showNotification(
        paste("Error running simulation:", e$message),
        type = "error",
        duration = 10
      )
      return(NULL)
    })
    
    # Store results
    results(sim_results)
    
    # Remove notification
    removeNotification("sim_notification")
    
    # Show success notification
    showNotification(
      "Simulation complete!",
      type = "message",
      duration = 3
    )
  })
  
  # Summary statistics outputs
  output$attack_rate_text <- renderText({
    if (is.null(results())) return("Run simulation")
    sprintf("%.1f%%", median(results()$summary_stats$attack_rate) * 100)
  })
  
  output$total_infected_text <- renderText({
    if (is.null(results())) return("Run simulation")
    sprintf("%.0f", median(results()$summary_stats$total_infected))
  })
  
  output$duration_text <- renderText({
    if (is.null(results())) return("Run simulation")
    sprintf("%.0f days", median(results()$summary_stats$actual_days))
  })
  
  # Main epidemic curve plot
  output$epidemic_curve_plot <- renderPlotly({
    req(results())
    library(dplyr)
    library(plotly)
    
    df <- results()$all_daily_data
    
    q025 <- function(x) stats::quantile(x, 0.025, na.rm = TRUE, names = FALSE, type = 8)
    q975 <- function(x) stats::quantile(x, 0.975, na.rm = TRUE, names = FALSE, type = 8)
    
    # Summaries by day (medians + 95% intervals)
    summary_data <- df %>%
      group_by(day) %>%
      summarise(
        Infected_median  = stats::median(P + Ra+Iso+QP,   na.rm = TRUE),
        Infected_lower   = q025(P + Ra+Iso+QP),
        Infected_upper   = q975(P + Ra+Iso+QP),
        
        # Recovered_median = stats::median(R,        na.rm = TRUE),
        # Recovered_lower  = q025(R),
        # Recovered_upper  = q975(R),
        .groups = "drop"
      )
    
    # Colors & labels
    state_info <- list(
      
      Infected  = list(color = "#E63946", label = "Infected (Prodromal + Rash)",
                       med = "Infected_median", low = "Infected_lower", up = "Infected_upper")
    )
    
    # Start an empty plotly object
    p <- plot_ly()
    
    # Add ribbons (no legend) then lines (with legend) for each state
    for (st in c("Infected")) {
      info <- state_info[[st]]
      
      # Ribbon (95% CI)
      p <- p %>%
        add_ribbons(
          data = summary_data,
          x = ~day,
          ymin = as.formula(paste0("~", info$low)),
          ymax = as.formula(paste0("~", info$up)),
          fillcolor = info$color,
          line = list(color = "rgba(0,0,0,0)"),
          opacity = 0.25,
          name = info$label,
          legendgroup = st,
          showlegend = FALSE,
          hoverinfo = "skip"
        )
      
      # Median line
      p <- p %>%
        add_lines(
          data = summary_data,
          x = ~day,
          y = as.formula(paste0("~", info$med)),
          name = info$label,
          legendgroup = st,
          showlegend = TRUE,
          line = list(color = info$color, width = 3),
          hovertemplate = paste0(
            "<b>", info$label, "</b>",
            "<br>Day: %{x}",
            "<br>Median: %{y:.1f}",
            "<extra></extra>"
          )
        )
    }
    
    p %>%
      layout(
        title = list(text = "Epidemic Curve"),
        xaxis = list(title = "Days"),
        yaxis = list(title = "Number of Infected Individuals"),
        legend = list(title = list(text = "State")),
        hovermode = "x unified"
      )
  })
  
  
  
  # Infectious states plot
  # Infectious states plot (MEDIAN instead of mean)
  output$infectious_plot <- renderPlot({
    req(results())
    library(dplyr)
    library(ggplot2)
    
    df <- results()$all_daily_data
    
    # helpers to match your main plot's quantile type
    q025 <- function(x) stats::quantile(x, 0.025, na.rm = TRUE, names = FALSE, type = 8)
    q975 <- function(x) stats::quantile(x, 0.975, na.rm = TRUE, names = FALSE, type = 8)
    
    # Summarize by day using MEDIANS and 95% interval
    summary_data <- df %>%
      group_by(day) %>%
      summarise(
        Infected_median = stats::median(P + Ra, na.rm = TRUE),
        Infected_lower  = q025(P + Ra),
        Infected_upper  = q975(P + Ra),
        
        Prodromal_median = stats::median(P,  na.rm = TRUE),
        Rash_median      = stats::median(Ra, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Peak based on the MEDIAN curve
    peak_idx   <- which.max(summary_data$Infected_median)
    peak_day   <- summary_data$day[peak_idx]
    peak_value <- summary_data$Infected_median[peak_idx]
    
    # Percent split at peak (guard against divide-by-zero)
    prod_pct <- if (peak_value > 0) 100 * summary_data$Prodromal_median[peak_idx] / peak_value else 0
    rash_pct <- if (peak_value > 0) 100 * summary_data$Rash_median[peak_idx]      / peak_value else 0
    
    # Median outbreak duration (purple dotted line)
    median_duration <- stats::median(results()$summary_stats$actual_days, na.rm = TRUE)
    
    p <- ggplot(summary_data, aes(x = day, y = Infected_median)) +
      # 95% band
      geom_ribbon(aes(ymin = Infected_lower, ymax = Infected_upper),
                  fill = "#E63946", alpha = 0.3) +
      # median line
      geom_line(color = "#E63946", linewidth = 2) +
      
      # Peak lines
      geom_vline(xintercept = peak_day, linetype = "dashed",
                 color = "gray40", linewidth = 0.8) +
      geom_hline(yintercept = peak_value, linetype = "dashed",
                 color = "gray40", linewidth = 0.8) +
      
      # Median outbreak duration line (purple dotted)
      geom_vline(xintercept = median_duration, linetype = "dotted",
                 color = "#6a1b9a", linewidth = 1.2) +
      
      # Peak annotation
      annotate("text", x = peak_day, y = peak_value * 1.1,
               label = sprintf("Peak: Day %d\n%.0f infectious", peak_day, peak_value),
               size = 4.5, fontface = "bold", color = "#E63946") +
      
      # Duration annotation (kept same style/positioning)
      annotate("text", x = median_duration, y = peak_value * 0.5,
               label = sprintf("Median outbreak\nend: Day %.0f", median_duration),
               size = 4, fontface = "bold", color = "#6a1b9a", hjust = -0.1) +
      
      labs(
        title = "Active Infectious Students in School",
        subtitle = sprintf("Median across simulations"),
        x = "Days",
        y = "Number of Infectious Individuals",
        caption = "Shaded area = 95% interval | Purple dotted line = median outbreak end"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray30"),
        plot.caption = element_text(size = 9, color = "gray50", hjust = 0),
        axis.title = element_text(size = 13, face = "bold"),
        axis.text = element_text(size = 11),
        panel.grid.minor = element_blank()
      )
    
    print(p)
  })
  
  
  
  # Attack rate distribution plot
  output$attack_rate_dist_plot <- renderPlot({
    req(results())
    
    # Create attack rate distribution data
    attack_rates <- results()$summary_stats$attack_rate * 100  # Convert to percentage
    
    # Calculate statistics
    median_ar <- median(attack_rates)
    mean_ar <- mean(attack_rates)
    
    # Create plot
    p <- ggplot(data.frame(attack_rate = attack_rates), aes(x = attack_rate)) +
      geom_histogram(bins = 30, fill = "#4CAF50", color = "white", alpha = 0.7) +
      
      # Add median line
      geom_vline(xintercept = median_ar, linetype = "solid", 
                 color = "#1976d2", linewidth = 1.2) +
      
      # Add mean line
      geom_vline(xintercept = mean_ar, linetype = "dashed", 
                 color = "#FF6F00", linewidth = 1.2) +
      
      # Annotations
      annotate("text", x = median_ar, y = Inf, 
               label = sprintf("Median: %.1f%%", median_ar),
               color = "#1976d2", angle = 90, vjust = 1.5, hjust = 1.1,
               size = 4.5, fontface = "bold") +
      
      annotate("text", x = mean_ar, y = Inf, 
               label = sprintf("Mean: %.1f%%", mean_ar),
               color = "#FF6F00", angle = 90, vjust = -0.5, hjust = 1.1,
               size = 4.5, fontface = "bold") +
      
      labs(
        title = "Distribution of Outbreak Sizes Across Simulations",
        subtitle = sprintf("Based on %d simulations", length(attack_rates)),
        x = "Outbreak Size (%)",
        y = "Number of Simulations",
        caption = "Blue solid = median | Orange dashed = mean"
      ) +
      
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray30"),
        plot.caption = element_text(size = 9, color = "gray50", hjust = 0),
        axis.title = element_text(size = 13, face = "bold"),
        axis.text = element_text(size = 11),
        panel.grid.minor = element_blank()
      )
    
    print(p)
  })
  
  # Detailed statistics tables
  output$attack_rate_table <- renderTable({
    r <- req(results())
    ss <- r$summary_stats
    ar <- ss$attack_rate
    
    # Guard: if attack_rate looks like proportions, scale to percent
    to_pct <- function(x) {
      if (is.null(x) || length(x) == 0) return(numeric(0))
      if (max(x, na.rm = TRUE) <= 1) x * 100 else x
    }
    arp <- to_pct(ar)
    
    fmt1 <- function(x) ifelse(is.finite(x), sprintf("%.1f%%", x), "—")
    
    data.frame(
      Statistic = c("Median", "Mean", "95% CI Lower", "95% CI Upper", "Min", "Max"),
      Value = c(
        fmt1(stats::median(arp, na.rm = TRUE)),
        fmt1(mean(arp, na.rm = TRUE)),
        fmt1(stats::quantile(arp, 0.025, na.rm = TRUE, names = FALSE, type = 8)),
        fmt1(stats::quantile(arp, 0.975, na.rm = TRUE, names = FALSE, type = 8)),
        fmt1(min(arp, na.rm = TRUE)),
        fmt1(max(arp, na.rm = TRUE))
      )
    )
  }, striped = TRUE, hover = TRUE)
  
  output$case_count_table <- renderTable({
    r <- req(results())
    ss <- r$summary_stats
    
    q2.5  <- function(x) stats::quantile(x, 0.025, na.rm = TRUE, names = FALSE, type = 8)
    q97.5 <- function(x) stats::quantile(x, 0.975, na.rm = TRUE, names = FALSE, type = 8)
    med0  <- function(x) sprintf("%.0f", stats::median(x, na.rm = TRUE))
    ci0   <- function(x) sprintf("%.0f - %.0f", q2.5(x), q97.5(x))
    
    data.frame(
      Statistic = c("Total Infected", "Isolated", "Quarantined", "Breakthrough"),
      Median = c(
        med0(ss$total_infected),
        med0(ss$total_isolated),
        med0(ss$total_quarantined),
        med0(ss$breakthrough_infections)
      ),
      `95% CI` = c(
        ci0(ss$total_infected),
        ci0(ss$total_isolated),
        ci0(ss$total_quarantined),
        ci0(ss$breakthrough_infections)
      ),
      check.names = FALSE
    )
  }, striped = TRUE, hover = TRUE)
  
  output$dynamics_table <- renderTable({
    r <- req(results())
    ss <- r$summary_stats
    df <- r$all_daily_data
    
    # Duration summary
    od <- ss$actual_days
    dur_mean <- mean(od, na.rm = TRUE)
    dur_min  <- min(od, na.rm = TRUE)
    dur_max  <- max(od, na.rm = TRUE)
    
    # Helpers
    fmt_day <- function(x) if (is.na(x)) "—" else sprintf("Day %.0f", as.numeric(x))
    fmt_cnt <- function(x) if (is.na(x)) "—" else sprintf("%.0f students", as.numeric(x))
    
    peak_day <- NA_real_
    peak_val <- NA_real_
    
    # Detect "long" vs "wide"
    nm <- names(df)
    has_long_cols <- all(c("P", "Ra") %in% nm)
    has_day <- "day" %in% nm
    
    if (has_long_cols && has_day) {
      # LONG: columns sim(optional), day, P, Ra
      med_by_day <- df |>
        dplyr::mutate(total_inf = P + Ra) |>
        dplyr::group_by(day) |>
        dplyr::summarise(median_total = stats::median(total_inf, na.rm = TRUE), .groups = "drop")
      if (nrow(med_by_day) > 0 && any(is.finite(med_by_day$median_total))) {
        idx <- which.max(med_by_day$median_total)
        peak_day <- med_by_day$day[idx]
        peak_val <- med_by_day$median_total[idx]
      }
    } else if (has_long_cols && !has_day) {
      # LONG but no day: assume row index is day
      vec <- (df$P + df$Ra)
      if (length(vec) > 0 && any(is.finite(vec))) {
        peak_day <- which.max(vec)
        peak_val <- max(vec, na.rm = TRUE)
      }
    } else {
      # WIDE: P_1..P_T and Ra_1..Ra_T (or P.1, Ra.1)
      p_cols  <- grep("^(P[._][0-9]+)$", nm, value = TRUE)
      ra_cols <- grep("^(Ra[._][0-9]+)$", nm, value = TRUE)
      
      if (length(p_cols) > 0 && length(ra_cols) > 0) {
        # Align by position (assume same length/order across P_x and Ra_x)
        tot_mat <- as.matrix(df[, p_cols, drop = FALSE]) +
          as.matrix(df[, ra_cols, drop = FALSE])
        # Median across simulations for each day (column)
        median_curve <- apply(tot_mat, 2, stats::median, na.rm = TRUE)
        if (length(median_curve) > 0 && any(is.finite(median_curve))) {
          peak_day <- which.max(median_curve)
          peak_val <- max(median_curve, na.rm = TRUE)
        }
      }
    }
    
    data.frame(
      Metric = c("Outbreak Duration", "Peak Day (median curve)", "Peak Infectious (median curve)"),
      Value  = c(
        sprintf("%.1f days (%.0f – %.0f)", dur_mean, dur_min, dur_max),
        fmt_day(peak_day),
        fmt_cnt(peak_val)
      ),
      stringsAsFactors = FALSE
    )
  }, striped = TRUE, hover = TRUE)
  
  output$intervention_table <- renderTable({
    req(selected_school_data())
    r <- req(results())
    ss <- r$summary_stats
    
    school_info <- selected_school_data()
    n_school <- as.numeric(school_info$`Total.Students`)
    vax_pct  <- vaccination_rate()  # Use reactive value
    qeff_pct <- as.numeric(input$quarantine_efficacy)
    
    # Basic counts
    n_vaccinated   <- round(n_school * vax_pct / 100)
    n_susceptible  <- n_school - n_vaccinated
    qeff_prop      <- qeff_pct / 100
    
    # Estimate
    med_quarantined <- stats::median(ss$total_quarantined, na.rm = TRUE)
    est_prevented   <- max(0, round(med_quarantined * qeff_prop * 0.30))
    
    data.frame(
      Metric = c("Susceptible Students", "Vaccinated Students", "Quarantine Efficacy",
                 "Cases Prevented by Quarantine (est.)"),
      Value = c(
        sprintf("%s", format(n_susceptible, big.mark = ",")),
        sprintf("%s (%.0f%%)", format(n_vaccinated, big.mark = ","), vax_pct),
        sprintf("%.0f%%", qeff_pct),
        sprintf("%s", format(est_prevented, big.mark = ","))
      ),
      stringsAsFactors = FALSE
    )
  }, striped = TRUE, hover = TRUE)
  
  output$params_table <- renderTable({
    req(selected_school_data())
    school_info <- selected_school_data()
    
    data.frame(
      Category  = c("School", "School", "School",
                    "Vaccination", "Disease", "Disease",
                    "Transmission", "Transmission",
                    "Quarantine", "Quarantine", "Quarantine",
                    "Simulation", "Simulation"),
      Parameter = c("School Size", "Average Class Size", "Number of Classes",
                    "Vaccination Coverage", "Isolation Delay", "Infectious Days at School",
                    "Within-Class Transmission", "Between-Class Transmission",
                    "Quarantine Contacts", "Quarantine Efficacy", "Quarantine Duration",
                    "Initial Infected", "Number of Simulations"),
      Value = c(
        sprintf("%s students", format(as.integer(school_info$`Total.Students`), big.mark = ",")),
        sprintf("%s students", FIXED_CLASS_SIZE),
        sprintf("%s classes", format(ceiling(school_info$`Total.Students` / FIXED_CLASS_SIZE), big.mark = ",")),
        sprintf("%.0f%%", vaccination_rate()),  # Use reactive value
        sprintf("%d days", as.integer(input$isolation_delay)),
        sprintf("%d days (prodromal + early rash)", 3L + as.integer(input$isolation_delay)),
        sprintf("%.3f", as.numeric(input$beta_within_class)),
        sprintf("%.4f", as.numeric(input$beta_between_class)),
        "Yes (classmates)",
        sprintf("%.0f%%", as.numeric(input$quarantine_efficacy)),
        sprintf("%d days", as.integer(input$quarantine_duration)),
        sprintf("%d", as.integer(input$initial_infected)),
        sprintf("%d", as.integer(input$n_simulations))
      ),
      stringsAsFactors = FALSE
    )
  }, striped = TRUE, hover = TRUE)
  
  # Download handler
  output$download_data <- downloadHandler(
    filename = function() {
      paste("measles_simulation_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(results())
      write.csv(results()$summary_stats, file, row.names = FALSE)
    }
  )
}

# ==============================================================================
# RUN APP
# ==============================================================================

shinyApp(ui = ui, server = server)