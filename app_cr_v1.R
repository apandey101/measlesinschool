# ==============================================================================
# MEASLES TRANSMISSION MODEL - SHINY APP (V2.2)
# ==============================================================================
# 
# Interactive web application for exploring measles transmission dynamics
# in school settings with class-based mixing
#
# Update: Comparison of WITH vs WITHOUT intervention scenarios
# Update: UI improvements, plotly graphs, recalibrated parameters
# Update: Separate isolation delays for index case (after rash) and 
#              secondary cases (after prodromal onset)
#
# Updated on 12/11/2025
# ==============================================================================

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(shinycssloaders)
library(plotly)
library(readxl)
library(shinyBS)
source("sim_cr_v1.R")  # V2.3: Fixed quarantine efficacy per unique contact

# ==============================================================================
# LOAD SCHOOL DATA
# ==============================================================================
school_data_raw <- read_xlsx("SC_vaccination_2025_2026.xlsx")
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
  mutate(Percent.Immunized = Percent.Immunized * 100)
school_data$Total.Students <- as.numeric(school_data$Total.Students)

if (nrow(school_data) == 0) {
  stop("No valid school data found after filtering.")
}

FIXED_CLASS_SIZE <- 25

# ==============================================================================
# UI (User Interface)
# ==============================================================================

ui <- fluidPage(
  
  titlePanel(
    div(
      h2("Measles Transmission Model: School-Based Simulation for South Carolina"),
      h4("Interactive Agent-Based Model", style = "color: gray;")
    )
  ),
  
  sidebarLayout(
    
    sidebarPanel(
      width = 3,
      
      h4("🏫 School Selection"),
      
      selectInput("region", "Select Region:",
                  choices = c("Choose a region..." = "", sort(unique(school_data$Region)))),
      
      selectInput("county", "Select County:",
                  choices = c("First select a region..." = "")),
      
      selectInput("school", "Select School:",
                  choices = c("First select a county..." = "")),
      
      conditionalPanel(
        condition = "input.school != ''",
        wellPanel(
          style = "background-color: #e8f5e9;",
          h5(strong("📋 Selected School"), style = "margin-top: 0; color: #2e7d32;"),
          htmlOutput("school_info_display")
        )
      ),
      
      hr(),
      sliderInput("initial_infected", "Initial Infected Students:",
                  min = 1, max = 10, value = 1, step = 1),
      hr(),
      
      h4("🚨 Isolation & Quarantine Policy"),
      
      # Index case isolation delay (after rash)
      sliderInput(
        "isolation_delay_index",
        tagList(
          "Index Case Isolation Delay (days):",
          tags$span(id = "isolation_delay_index_help", icon("info-circle"),
                    class = "text-muted", style = "margin-left:6px; cursor: help;")
        ),
        min = 0, max = 4, value = 3, step = 1
      ),
      bsTooltip("isolation_delay_index_help",
                "Days after rash onset before the first (index) case is isolated. Rash typically appears 4 days after prodromal symptoms begin.",
                placement = "right", trigger = "hover"),
      
      # Secondary case isolation delay (after prodromal onset)
      sliderInput(
        "isolation_delay_secondary",
        tagList(
          "Secondary Case Isolation Delay (days):",
          tags$span(id = "isolation_delay_secondary_help", icon("info-circle"),
                    class = "text-muted", style = "margin-left:6px; cursor: help;")
        ),
        min = 0, max = 8, value = 6, step = 1
      ),
      bsTooltip("isolation_delay_secondary_help",
                "Days after prodromal symptom onset before secondary cases are isolated. Once the index case is identified, subsequent cases can be detected earlier through contact tracing and symptom monitoring.",
                placement = "right", trigger = "hover"),
      
      sliderInput(
        "quarantine_efficacy",
        tagList(
          "Quarantine Efficacy (%)",
          tags$span(id = "qe_help", icon("info-circle"),
                    class = "text-muted", style = "margin-left:6px; cursor: help;")
        ),
        min = 0, max = 100, value = 50, step = 10
      ),
      bsTooltip("qe_help",
                "Probability (in %) that each unvaccinated classmate of a symptomatic case is successfully quarantined. Vaccinated classmates remain in school.",
                placement = "right", trigger = "hover"),
      
      sliderInput("quarantine_duration", "Quarantine Duration (days):",
                  min = 7, max = 28, value = 21, step = 7),
      
      hr(),
      
      actionButton("run_sim", "🚀 Run Simulation",
                   class = "btn-primary btn-lg btn-block",
                   style = "margin-top: 10px; margin-bottom: 10px;"),
      
      downloadButton("download_data", "📥 Download Results",
                     class = "btn-success btn-block"),
      
      hr(),
      
      # Advanced Options
      tags$details(
        tags$summary(
          style = "cursor: pointer; font-weight: bold; font-size: 16px; color: #1976d2;
                   padding: 12px; background-color: #f5f5f5; border-radius: 4px;
                   border: 2px dashed #1976d2; user-select: none; margin-bottom: 5px;",
          "⚙️ Advanced Options (Click to expand)"
        ),
        
        div(
          style = "padding: 15px; background-color: #fafafa; border: 1px solid #e0e0e0;
                   border-radius: 4px; margin-top: 10px; margin-bottom: 10px;",
          
          h5(strong("⚠️ Override Default Parameters"),
             style = "color: #d32f2f; margin-top: 0; margin-bottom: 10px;"),
          
          p(style = "font-size: 13px; color: #666; line-height: 1.5; margin-bottom: 15px;",
            "Expand this section only if you want to test different scenarios."),
          
          div(
            style = "background-color: #fff; padding: 12px; border-radius: 4px;
                     margin-bottom: 15px; border-left: 4px solid #2196F3;",
            h4("🎲 Simulation Settings"),
            
            sliderInput("n_simulations", "Number of Simulations:",
                        min = 10, max = 500, value = 150, step = 10),
            
            sliderInput("n_days", "Maximum Simulation Days:",
                        min = 60, max = 200, value = 150, step = 10),
            
            hr(),
            
            h6(strong("Transmission Parameters"), style = "margin-top: 0; color: #1976d2;"),
            
            sliderInput("c_within", "Within-class contacts per day:",
                        min = 0, max = 30, value = 24, step = 1),
            sliderInput("c_between", "Between-class contacts per day:",
                        min = 0, max = 20, value = 10, step = 1),
            sliderInput("p_within", "Per-contact transmission (within):",
                        min = 0, max = 1, value = 0.0937, step = 0.001),
            sliderInput("p_between", "Per-contact transmission (between):",
                        min = 0, max = 1, value = 0.045, step = 0.001),
            
            tags$small(style = "color:#666;",
                       "Default parameters calibrated for measles R0 ≈ 15.")
          ),
          
          div(
            style = "background-color: #fff; padding: 12px; border-radius: 4px;
                     border-left: 4px solid #FF9800;",
            
            h6(strong("Vaccination Coverage Override"), style = "margin-top: 0; color: #ef6c00;"),
            
            sliderInput("vaccination_coverage", "Custom Vaccination Coverage (%):",
                        min = 0, max = 100, value = 90, step = 1),
            
            checkboxInput("use_custom_vaccination",
                          "Use custom vaccination coverage (uncheck to use school's rate)",
                          value = FALSE),
            
            tags$small(style = "color: #666;",
                       "Unchecking this box will use the school's vaccination rate displayed above.")
          )
        )
      )
    ),
    
    # MAIN PANEL
    mainPanel(
      width = 9,
      
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
        # Tab 2: Intervention Comparison
        # =======================
        tabPanel(
          "📈 Intervention Comparison",
          br(),
          
          tags$div(
            class = "alert alert-warning",
            style = "background-color: #fff8e1; border-left: 4px solid #ffb300; padding: 12px; border-radius: 4px;",
            h5(icon("exclamation-triangle"), strong(" Ongoing Work — Interpret with Caution"),
               style = "margin-top: 0; color: #8a6d3b;"),
            p("Results depend on parameter choices and policy assumptions, including isolation timing, quarantine coverage (unvaccinated-only policy), vaccination levels, and school structure. Small changes can lead to large shifts in outcomes.")
          ),
          
          # Combined summary boxes in a more compact layout
          fluidRow(
            column(
              6,
              wellPanel(
                style = "background-color: #e3f2fd; padding: 15px; margin-bottom: 10px;",
                h5(icon("shield-alt"), strong(" With Intervention"), style = "color: #1976d2; margin-top: 0;"),
                fluidRow(
                  column(4,
                         h4(textOutput("outbreak_size_text_intervention"), style = "color: #1976d2; margin: 5px 0;"),
                         p(
                           tagList(
                             "Final Outbreak Size ",
                             tags$span(id = "fos_int_help", icon("info-circle"),
                                       class = "text-muted", style = "cursor: help;")
                           ),
                           style = "font-size: 12px; margin: 0;"
                         )
                  ),
                  column(4,
                         h4(textOutput("total_infected_text_intervention"), style = "color: #1976d2; margin: 5px 0;"),
                         p("Total Infected", style = "font-size: 12px; margin: 0;")
                  ),
                  column(4,
                         h4(textOutput("duration_text_intervention"), style = "color: #1976d2; margin: 5px 0;"),
                         p(
                           tagList(
                             "Duration (days) ",
                             tags$span(id = "dur_int_help", icon("info-circle"),
                                       class = "text-muted", style = "cursor: help;")
                           ),
                           style = "font-size: 12px; margin: 0;"
                         )
                  )
                )
              )
            ),
            column(
              6,
              wellPanel(
                style = "background-color: #ffebee; padding: 15px; margin-bottom: 10px;",
                h5(icon("virus"), strong(" Without Intervention"), style = "color: #c62828; margin-top: 0;"),
                fluidRow(
                  column(4,
                         h4(textOutput("outbreak_size_text_baseline"), style = "color: #c62828; margin: 5px 0;"),
                         p(
                           tagList(
                             "Final Outbreak Size ",
                             tags$span(id = "fos_base_help", icon("info-circle"),
                                       class = "text-muted", style = "cursor: help;")
                           ),
                           style = "font-size: 12px; margin: 0;"
                         )
                  ),
                  column(4,
                         h4(textOutput("total_infected_text_baseline"), style = "color: #c62828; margin: 5px 0;"),
                         p("Total Infected", style = "font-size: 12px; margin: 0;")
                  ),
                  column(4,
                         h4(textOutput("duration_text_baseline"), style = "color: #c62828; margin: 5px 0;"),
                         p(
                           tagList(
                             "Duration (days) ",
                             tags$span(id = "dur_base_help", icon("info-circle"),
                                       class = "text-muted", style = "cursor: help;")
                           ),
                           style = "font-size: 12px; margin: 0;"
                         )
                  )
                )
              )
            )
          ),
          
          # Tooltips for Final Outbreak Size and Duration
          bsTooltip("fos_int_help",
                    "Percentage of the student population infected by the end of the outbreak. Calculated as (Total Infected / School Size) × 100.",
                    placement = "top", trigger = "hover"),
          bsTooltip("fos_base_help",
                    "Percentage of the student population infected by the end of the outbreak. Calculated as (Total Infected / School Size) × 100.",
                    placement = "top", trigger = "hover"),
          bsTooltip("dur_int_help",
                    "Number of days until outbreak ends (no active infections). If duration equals maximum simulation days, the outbreak may still be ongoing and results are truncated.",
                    placement = "top", trigger = "hover"),
          bsTooltip("dur_base_help",
                    "Number of days until outbreak ends (no active infections). If duration equals maximum simulation days, the outbreak may still be ongoing and results are truncated.",
                    placement = "top", trigger = "hover"),
          
          # Intervention Impact
          fluidRow(
            column(12,
                   wellPanel(
                     style = "background-color: #e8f5e9; border-left: 4px solid #4caf50; padding: 10px;",
                     fluidRow(
                       column(1, icon("chart-line", class = "fa-2x", style = "color: #2e7d32; margin-top: 5px;")),
                       column(11,
                              h4(textOutput("cases_prevented_text"), style = "color: #2e7d32; margin: 0;"),
                              p("Estimated Cases Prevented by Intervention", style = "margin: 0; font-size: 13px;")
                       )
                     )
                   )
            )
          ),
          
          # Comparison plot
          withSpinner(
            plotlyOutput("comparison_curve_plot", height = "500px"),
            type = 6, color = "#1976d2"
          ),
          
          br(),
          
          wellPanel(
            style = "background-color: #f5f5f5; border-left: 4px solid #4caf50;",
            h5(icon("lightbulb"), strong(" Interpreting These Results"), style = "color: #2e7d32; margin-top: 0;"),
            tags$ul(
              style = "color: #333; margin-bottom: 0;",
              tags$li(tags$span("Blue", style = "color: #1976d2; font-weight: bold;"), ": Total infectious students WITH interventions (isolation + contact tracing)."),
              tags$li(tags$span("Red", style = "color: #c62828; font-weight: bold;"), ": Total infectious students WITHOUT interventions (natural epidemic spread)."),
              tags$li("Shaded areas show 95% uncertainty intervals across simulations."),
              tags$li("With interventions, infectious students are isolated and removed from school, reducing the curve.")
            )
          )
        ),
        
        
        
        
        # =======================
        # Tab 3: Infectious States
        # =======================
        
        tabPanel(
          "🦠 Infectious States",
          br(),
          withSpinner(
            plotlyOutput("infectious_plot", height = "500px"),
            type = 6, color = "#c62828"
          ),
          br(),
          wellPanel(
            h4("About This Plot"),
            p("This plot compares the number of infectious students still attending school over time for both scenarios:"),
            tags$ul(
              tags$li(strong("Blue: "), "With Intervention — students are isolated and contacts quarantined."),
              tags$li(strong("Red: "), "Without Intervention — natural epidemic progression."),
              tags$li(strong("Shaded areas: "), "95% uncertainty intervals.")
            ),
            p(strong("Note:"), " The curves represent only those who remain infectious while present in school. With intervention, isolation removes infectious individuals, reducing the peak.")
          )
        ),
        
        
        # =======================
        # Tab 4: About
        # =======================
        tabPanel(
          "ℹ️ About",
          br(),
          h3("Measles Transmission Model"),
          
          h4("Model Description"),
          p("This agent-based model represents measles transmission within a school community. Each student is modeled as an individual with their own vaccination status, infection stage, and classroom assignment. Transmission occurs primarily within classrooms, with a lower rate of between-class mixing. The model focuses on how vaccination coverage, timing of isolation for the index case, and the effectiveness of contact tracing and quarantine of unvaccinated contacts shape outbreak dynamics."),
          
          p("The model runs in daily time steps. When a student becomes exposed they progress through a latent (incubation) phase before entering the prodromal (pre-rash infectious) phase and then develop a rash."),
          
          h5("Isolation Timing"),
          tags$ul(
            tags$li(strong("Index case: "), "The first identified case is isolated after rash onset plus a user-specified delay (0-4 days). Rash typically appears 4 days after prodromal symptoms begin."),
            tags$li(strong("Secondary cases: "), "Once the outbreak is recognized, subsequent cases are isolated based on prodromal symptoms plus a user-specified delay (0-8 days). This allows to reflect enhanced surveillance after the index case is identified.")
          ),
          
          p("Quarantine is applied only to unvaccinated contacts identified as exposed to a symptomatic case. Vaccinated contacts are not quarantined in this scenario. Quarantined students are removed from school mixing for the defined quarantine period, then return to their prior immune or susceptible status."),
          
          p(strong("Important:"), " This model focuses exclusively on school-based transmission and does not model household, community, or social-network spread outside class and school. It is designed to explore the impact of school-based control measures such as isolation timing, quarantine coverage, and vaccination thresholds."),
          
          h4("Key Features"),
          tags$ul(
            tags$li(strong("Two levels of contact:"), " Higher exposure risk within classrooms; lower between-class mixing."),
            tags$li(strong("Vaccination effects:"), " Vaccinated students are less susceptible and, if infected, less infectious to others."),
            tags$li(strong("Separate isolation delays:"), " Different delays for index case (after rash) and secondary cases (after prodromal onset)."),
            tags$li(strong("Targeted interventions:"), " Index case isolation after rash + delay, isolation of subsequent symptomatic students (prodromal onset + delay), and quarantine of unvaccinated exposed contacts only."),
            tags$li(strong("Stochastic simulation:"), " Each simulation run represents a possible outbreak trajectory, enabling exploration of uncertainty in outcomes."),
            tags$li(strong("Scenario comparison:"), " Direct comparison of WITH vs WITHOUT intervention outcomes.")
          ),
          
          h4("Model Assumptions and Limitations"),
          
          h5("Scope"),
          tags$ul(
            tags$li(strong("Single-school model:"), " Simulates transmission within one school only. Household, community, and external transmission are not modeled."),
            tags$li(strong("Closed population:"), " No students enter or leave during the simulation.")
          ),
          
          h5("Contact Structure"),
          tags$ul(
            tags$li(strong("Two-level mixing:"), " Contacts occur within-class (higher intensity) and between-class (lower intensity)."),
            tags$li(strong("Uniform random mixing:"), " Within each contact pool, individuals are selected randomly with equal probability."),
            tags$li(strong("Poisson-distributed contacts:"), " Daily contact counts are drawn from Poisson distributions, introducing stochastic variation.")
          ),
          
          h5("Disease Dynamics"),
          tags$ul(
            tags$li(strong("Fixed stage durations:"), " Latent (10 days), prodromal (4 days), and rash (4 days) periods are fixed values."),
            tags$li(strong("No asymptomatic infections:"), " All infections progress to symptomatic (rash) stage.")
          ),
          
          h5("Vaccination"),
          tags$ul(
            tags$li(strong("Leaky vaccine:"), " Vaccination reduces (but does not eliminate) susceptibility. Breakthrough infections are possible."),
            tags$li(strong("Fixed at initialization:"), " No post-exposure prophylaxis or vaccination campaigns during outbreak."),
            tags$li(strong("Random distribution:"), " Unvaccinated students are randomly distributed, not clustered.")
          ),
          
          h5("Interventions"),
          tags$ul(
            tags$li(strong("Perfect isolation:"), " Isolated cases have zero school contacts; home-based transmission not modeled."),
            tags$li(strong("Retrospective contact tracing:"), " Contacts during the infectious period prior to isolation are traced and quarantined."),
            tags$li(strong("Vaccinated contacts exempt:"), " Only unvaccinated contacts are quarantined."),
            tags$li(strong("Probabilistic quarantine:"), " Quarantine success depends on quarantine efficacy parameter (imperfect compliance).")
          ),
          
          h5("Key Limitations"),
          tags$ul(
            tags$li("Does not capture transmission outside school settings"),
            tags$li("No weekend/holiday effects or behavioral responses to outbreak awareness"),
            tags$li("Uniform mixing may underestimate transmission in friend groups"),
            tags$li("Assumes no clustering of unvaccinated students"),
            tags$li("No post-exposure vaccination modeled")
          ),
          
          h4("Interpretation"),
          p("This model is intended for comparative scenario analysis rather than exact prediction. It allows users to explore how small changes — such as a one-day delay in index case isolation, a drop in vaccine coverage, or incomplete quarantine of contacts — can dramatically alter outbreak size and duration. Running multiple simulations provides a view of the range of possible outcomes and the stochastic variability inherent in outbreak dynamics."),
          
          h4("References & Further Resources"),
          p(
            strong("Measles disease information: "),
            tags$a("US CDC – Measles", href = "https://www.cdc.gov/pinkbook/hcp/table-of-contents/chapter-13-measles.html", target = "_blank"),
            ", ",
            tags$a("ECDC Measles Factsheet", href = "https://www.ecdc.europa.eu/en/measles/facts", target = "_blank")
          ),
          p(
            strong("Other Models: "),
            tags$a("epiENGAGE Measles Outbreak Simulator", href = "https://epiengage-measles.tacc.utexas.edu/", target = "_blank"),
            ", ",
            tags$a("epiworldR Shiny", href = "https://gvegayon.shinyapps.io/epiworldRShiny/", target = "_blank")
          ),
          
          br(),
          hr(),
          p("Model Version 2.2 | Created 2025 | Updated on 12/11/2025", style = "text-align: center; color: gray;"),
          p("For more information about the code and model details, please visit: ",
            tags$a("https://github.com/apandey101/measlesinschool", href = "https://github.com/apandey101/measlesinschool", target = "_blank"),
            style = "text-align: center;")
        )
      )
    )
  )
)


# ==============================================================================
# SERVER
# ==============================================================================

server <- function(input, output, session) {
  
  results_intervention <- reactiveVal(NULL)
  results_baseline <- reactiveVal(NULL)
  selected_school_data <- reactiveVal(NULL)
  
  # Region -> County dropdown
  observeEvent(input$region, {
    req(input$region != "")
    counties <- school_data %>%
      filter(Region == input$region) %>%
      pull(County) %>% unique() %>% sort()
    updateSelectInput(session, "county", choices = c("Select a county..." = "", counties))
    updateSelectInput(session, "school", choices = c("First select a county..." = ""))
    selected_school_data(NULL)
  })
  
  # County -> School dropdown
  observeEvent(input$county, {
    req(input$county != "")
    schools <- school_data %>%
      filter(Region == input$region, County == input$county) %>%
      pull(`School.Name`) %>% sort()
    updateSelectInput(session, "school", choices = c("Select a school..." = "", schools))
    selected_school_data(NULL)
  })
  
  # School selection
  observeEvent(input$school, {
    req(input$school != "")
    school_info <- school_data %>%
      filter(Region == input$region, County == input$county, `School.Name` == input$school)
    if (nrow(school_info) == 0) {
      showNotification("Error: Selected school not found.", type = "error", duration = 5)
      return()
    }
    selected_school_data(school_info)
    updateSliderInput(session, "vaccination_coverage", value = school_info$`Percent.Immunized`)
    updateCheckboxInput(session, "use_custom_vaccination", value = FALSE)
  })
  
  vaccination_rate <- reactive({
    req(selected_school_data())
    if (isTRUE(input$use_custom_vaccination)) input$vaccination_coverage
    else selected_school_data()$`Percent.Immunized`
  })
  
  output$school_info_display <- renderUI({
    req(selected_school_data())
    info <- selected_school_data()
    HTML(paste0(
      "<strong>School:</strong> ", info$`School.Name`, "<br/>",
      "<strong>Region:</strong> ", info$Region, "<br/>",
      "<strong>County:</strong> ", info$County, "<br/>",
      "<strong>Students:</strong> ", format(info$`Total.Students`, big.mark = ","), "<br/>",
      "<strong>Vaccination Rate:</strong> ", round(info$`Percent.Immunized`, 1), "%"
    ))
  })
  
  # Run simulation
  observeEvent(input$run_sim, {
    if (is.null(selected_school_data())) {
      showNotification("Please select a school before running simulation.", type = "error", duration = 5)
      return()
    }
    
    school_info <- selected_school_data()
    if (is.na(school_info$`Total.Students`) || school_info$`Total.Students` <= 0) {
      showNotification("Error: Invalid student count.", type = "error", duration = 5)
      return()
    }
    
    showNotification("Running simulations... This may take a few seconds.",
                     duration = NULL, id = "sim_notification", type = "message")
    
    base_params <- list(
      n_simulations = input$n_simulations,
      school_size = as.integer(school_info$`Total.Students`),
      avg_class_size = FIXED_CLASS_SIZE,
      age_range = c(5, 17),
      vaccination_coverage = vaccination_rate() / 100,
      latent_period = 10,
      prodromal_period = 4,
      rash_period = 4,
      isolation_delay_index = input$isolation_delay_index,
      isolation_delay_secondary = input$isolation_delay_secondary,
      isolation_period = 8,
      vaccine_efficacy = 0.97,
      vaccine_infectiousness_reduction = 0.80,
      prodromal_infectiousness_multiplier = 1.0,
      rash_infectiousness_multiplier = 0.6,
      c_within = input$c_within,
      c_between = input$c_between,
      p_within = input$p_within,
      p_between = input$p_between,
      quarantine_contacts = TRUE,
      quarantine_efficacy = input$quarantine_efficacy / 100,
      quarantine_duration = input$quarantine_duration,
      initial_infected = input$initial_infected,
      n_days = input$n_days,
      seed_start = 1234567,
      verbose = FALSE
    )
    
    # With intervention
    sim_int <- tryCatch({
      params <- base_params
      params$no_intervention <- FALSE
      do.call(run_multiple_simulations, params)
    }, error = function(e) { showNotification(paste("Error:", e$message), type = "error"); NULL })
    
    # Without intervention
    sim_base <- tryCatch({
      params <- base_params
      params$no_intervention <- TRUE
      do.call(run_multiple_simulations, params)
    }, error = function(e) { showNotification(paste("Error:", e$message), type = "error"); NULL })
    
    results_intervention(sim_int)
    results_baseline(sim_base)
    
    removeNotification("sim_notification")
    showNotification("Both simulations complete!", type = "message", duration = 3)
  })
  
  # ===========================================================================
  # TAB 1: Main Results (Intervention Only)
  # ===========================================================================
  
  # Summary stats for first tab (intervention only)
  output$attack_rate_text <- renderText({
    if (is.null(results_intervention())) return("—")
    sprintf("%.1f%%", median(results_intervention()$summary_stats$attack_rate) * 100)
  })
  
  output$total_infected_text <- renderText({
    if (is.null(results_intervention())) return("—")
    sprintf("%.0f", median(results_intervention()$summary_stats$total_infected))
  })
  
  output$duration_text <- renderText({
    if (is.null(results_intervention())) return("—")
    sprintf("%.0f", median(results_intervention()$summary_stats$actual_days))
  })
  
  # Main epidemic curve (intervention only) for Tab 1
  output$epidemic_curve_plot <- renderPlotly({
    req(results_intervention())
    
    df_int <- results_intervention()$all_daily_data
    
    q025 <- function(x) quantile(x, 0.025, na.rm = TRUE)
    q975 <- function(x) quantile(x, 0.975, na.rm = TRUE)
    
    # Total infected per day (P + Ra + Iso + QP)
    summary_int <- df_int %>%
      mutate(total_infected = P + Ra + Iso + QP) %>%
      group_by(day) %>%
      summarise(
        median = median(total_infected, na.rm = TRUE),
        lower  = q025(total_infected),
        upper  = q975(total_infected),
        .groups = "drop"
      )
    
    plot_ly() %>%
      # Intervention ribbon
      add_ribbons(
        data = summary_int, x = ~day, ymin = ~lower, ymax = ~upper,
        fillcolor = "rgba(25, 118, 210, 0.2)", line = list(color = "transparent"),
        name = "95% CI", showlegend = TRUE, hoverinfo = "skip"
      ) %>%
      # Intervention line
      add_lines(
        data = summary_int, x = ~day, y = ~median,
        line = list(color = "#1976d2", width = 3),
        name = "Median (With Intervention)",
        hovertemplate = "<b>With Intervention</b><br>Day: %{x}<br>Median infected: %{y:.1f}<extra></extra>"
      ) %>%
      layout(
        title = list(text = "Epidemic Curve: With Intervention"),
        xaxis = list(title = "Days"),
        yaxis = list(title = "Number of Infected Students"),
        legend = list(orientation = "h", y = -0.15, x = 0.5, xanchor = "center"),
        hovermode = "x unified"
      )
  })
  
  # ===========================================================================
  # TAB 2: Intervention Comparison (With vs Without)
  # ===========================================================================
  
  # Summary stats - Intervention
  output$outbreak_size_text_intervention <- renderText({
    if (is.null(results_intervention())) return("—")
    sprintf("%.1f%%", median(results_intervention()$summary_stats$attack_rate) * 100)
  })
  output$total_infected_text_intervention <- renderText({
    if (is.null(results_intervention())) return("—")
    sprintf("%.0f", median(results_intervention()$summary_stats$total_infected))
  })
  output$duration_text_intervention <- renderText({
    if (is.null(results_intervention())) return("—")
    sprintf("%.0f", median(results_intervention()$summary_stats$actual_days))
  })
  
  # Summary stats - Baseline
  output$outbreak_size_text_baseline <- renderText({
    if (is.null(results_baseline())) return("—")
    sprintf("%.1f%%", median(results_baseline()$summary_stats$attack_rate) * 100)
  })
  output$total_infected_text_baseline <- renderText({
    if (is.null(results_baseline())) return("—")
    sprintf("%.0f", median(results_baseline()$summary_stats$total_infected))
  })
  output$duration_text_baseline <- renderText({
    if (is.null(results_baseline())) return("—")
    sprintf("%.0f", median(results_baseline()$summary_stats$actual_days))
  })
  
  # Cases prevented
  output$cases_prevented_text <- renderText({
    if (is.null(results_intervention()) || is.null(results_baseline())) return("Run simulation to see results")
    med_base <- median(results_baseline()$summary_stats$total_infected)
    med_int <- median(results_intervention()$summary_stats$total_infected)
    prevented <- med_base - med_int
    pct <- ifelse(med_base > 0, (prevented / med_base) * 100, 0)
    sprintf("%.0f cases prevented (%.1f%% reduction)", prevented, pct)
  })
  
  # Comparison epidemic curve (both scenarios) for Tab 2
  output$comparison_curve_plot <- renderPlotly({
    req(results_intervention(), results_baseline())
    
    df_int  <- results_intervention()$all_daily_data
    df_base <- results_baseline()$all_daily_data
    
    q025 <- function(x) quantile(x, 0.025, na.rm = TRUE)
    q975 <- function(x) quantile(x, 0.975, na.rm = TRUE)
    
    # Total infected per day
    summary_int <- df_int %>%
      mutate(total_infected = P + Ra + Iso + QP) %>%
      group_by(day) %>%
      summarise(
        median = median(total_infected, na.rm = TRUE),
        lower  = q025(total_infected),
        upper  = q975(total_infected),
        .groups = "drop"
      )
    
    summary_base <- df_base %>%
      mutate(total_infected = P + Ra + Iso + QP) %>%
      group_by(day) %>%
      summarise(
        median = median(total_infected, na.rm = TRUE),
        lower  = q025(total_infected),
        upper  = q975(total_infected),
        .groups = "drop"
      )
    
    plot_ly() %>%
      # Baseline ribbon
      add_ribbons(
        data = summary_base, x = ~day, ymin = ~lower, ymax = ~upper,
        fillcolor = "rgba(230, 57, 70, 0.2)", line = list(color = "transparent"),
        name = "Without Intervention (95% CI)", legendgroup = "base",
        showlegend = FALSE, hoverinfo = "skip"
      ) %>%
      # Baseline line
      add_lines(
        data = summary_base, x = ~day, y = ~median,
        line = list(color = "#E63946", width = 3),
        name = "Without Intervention", legendgroup = "base",
        hovertemplate = "<b>Without Intervention</b><br>Day: %{x}<br>Median infected: %{y:.1f}<extra></extra>"
      ) %>%
      # Intervention ribbon
      add_ribbons(
        data = summary_int, x = ~day, ymin = ~lower, ymax = ~upper,
        fillcolor = "rgba(25, 118, 210, 0.2)", line = list(color = "transparent"),
        name = "With Intervention (95% CI)", legendgroup = "int",
        showlegend = FALSE, hoverinfo = "skip"
      ) %>%
      # Intervention line
      add_lines(
        data = summary_int, x = ~day, y = ~median,
        line = list(color = "#1976d2", width = 3),
        name = "With Intervention", legendgroup = "int",
        hovertemplate = "<b>With Intervention</b><br>Day: %{x}<br>Median infected: %{y:.1f}<extra></extra>"
      ) %>%
      layout(
        title = list(text = "Epidemic Curve: With vs Without Intervention"),
        xaxis = list(title = "Days"),
        yaxis = list(title = "Number of Infected Students"),
        legend = list(orientation = "h", y = -0.15, x = 0.5, xanchor = "center"),
        hovermode = "x unified"
      )
  })
  
  # ===========================================================================
  # TAB 3: Infectious States (Active at School)
  # ===========================================================================
  
  # Infectious states plot (P + Ra only - active at school)
  output$infectious_plot <- renderPlotly({
    req(results_intervention(), results_baseline())
    
    df_int <- results_intervention()$all_daily_data
    df_base <- results_baseline()$all_daily_data
    
    q025 <- function(x) quantile(x, 0.025, na.rm = TRUE)
    q975 <- function(x) quantile(x, 0.975, na.rm = TRUE)
    
    summary_int <- df_int %>%
      mutate(active = P + Ra) %>%
      group_by(day) %>%
      summarise(median = median(active, na.rm = TRUE),
                lower = q025(active), upper = q975(active), .groups = "drop")
    
    summary_base <- df_base %>%
      mutate(active = P + Ra) %>%
      group_by(day) %>%
      summarise(median = median(active, na.rm = TRUE),
                lower = q025(active), upper = q975(active), .groups = "drop")
    
    plot_ly() %>%
      # Baseline ribbon
      add_ribbons(data = summary_base, x = ~day, ymin = ~lower, ymax = ~upper,
                  fillcolor = "rgba(230, 57, 70, 0.2)", line = list(color = "transparent"),
                  legendgroup = "base", showlegend = FALSE,
                  hoverinfo = "skip") %>%
      # Baseline line
      add_lines(data = summary_base, x = ~day, y = ~median,
                line = list(color = "#E63946", width = 3),
                name = "Without Intervention", legendgroup = "base",
                hovertemplate = "<b>Without Intervention</b><br>Day: %{x}<br>Active Infectious: %{y:.1f}<extra></extra>") %>%
      # Intervention ribbon
      add_ribbons(data = summary_int, x = ~day, ymin = ~lower, ymax = ~upper,
                  fillcolor = "rgba(25, 118, 210, 0.2)", line = list(color = "transparent"),
                  legendgroup = "int", showlegend = FALSE,
                  hoverinfo = "skip") %>%
      # Intervention line
      add_lines(data = summary_int, x = ~day, y = ~median,
                line = list(color = "#1976d2", width = 3),
                name = "With Intervention", legendgroup = "int",
                hovertemplate = "<b>With Intervention</b><br>Day: %{x}<br>Active Infectious: %{y:.1f}<extra></extra>") %>%
      layout(
        title = list(text = "Active Infectious Students at School"),
        xaxis = list(title = "Days"),
        yaxis = list(title = "Number of Active Infectious Students"),
        legend = list(orientation = "h", y = -0.15, x = 0.5, xanchor = "center"),
        hovermode = "x unified"
      )
  })
  
  # ===========================================================================
  # Download Handler
  # ===========================================================================
  
  output$download_data <- downloadHandler(
    filename = function() paste("measles_simulation_", Sys.Date(), ".csv", sep = ""),
    content = function(file) {
      req(results_intervention(), results_baseline())
      df_int <- results_intervention()$summary_stats %>% mutate(scenario = "With_Intervention")
      df_base <- results_baseline()$summary_stats %>% mutate(scenario = "Without_Intervention")
      write.csv(bind_rows(df_int, df_base), file, row.names = FALSE)
    }
  )
}

# ==============================================================================
# RUN APP
# ==============================================================================

shinyApp(ui = ui, server = server)