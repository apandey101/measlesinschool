# Measles-in-School Agent-Based Model
## Overview
This repository contains an agent-based model (ABM) simulating measles transmission within a single school environment. The model allows researchers and decision-makers to explore how school-based control measures—such as isolation timing, quarantine coverage, and vaccination thresholds—shape outbreak dynamics.
The school data on number of students and vaccination rates are taken from South Carolina Department of Public Health and 
represent the percentage of total students who have completed vaccine schedules for all required vaccines, including MMR (measles, mumps, rubella).

Implemented in R with core transmission loops optimized in C++ (via Rcpp), the simulation focuses on the interplay between classroom structure, vaccination gaps, and intervention delays.

## Key Features
Hybrid Architecture: Uses R for data management and UI, with Rcpp for high-speed stochastic transmission and quarantine loops (sim_cr_v1.R).

Granular Mixing: Models two levels of contact risk: high-intensity mixing within fixed classrooms and lower-intensity mixing between classes.

Disease Progression: Tracks agents through Susceptible (S), Exposed (E), Prodromal (P), Rash (Ra) and Recovered (R) stages.

Differentiated Interventions:

Vaccination: Modeled with both susceptibility reduction and reduced infectiousness for breakthrough cases.

Isolation: Distinct logic for the Index case (delayed isolation after rash) vs. subsequent cases (isolation triggered by prodromal onset).

Quarantine: Targeted removal of unvaccinated contacts from the mixing pool.

Stochasticity: Daily contact events are drawn from Poisson distributions, allowing for the exploration of outbreak variability.

## Model Logic & Dynamics
### 1. Agent States & Transitions
The model follows an extended SEIR framework tailored for measles:

S (Susceptible): Unvaccinated individuals.

V (Vaccinated): Individuals with reduced susceptibility.

E (Exposed): Latent phase (non-infectious).

P (Prodromal): Pre-rash infectious phase.

Ra (Rash): Symptomatic infectious phase.

Iso (Isolated): Removed from school mixing.

R (Recovered): Immune.

QS/QE/QP (Quarantined): Unvaccinated contacts removed from mixing.

### 2. Transmission Mechanics
Transmission is contact-based. Each infectious agent (P or Ra) generates contacts daily based on Poisson distributions:

Within-Class: N ~ Poisson(c_within) targeting classmates.

Between-Class: N ~ Poisson(c_between) targeting non-classmates.

Transmission probability depends on the target's status and the infector's stage:

If the target is vaccinated, susceptibility is multiplied by (1 - vaccine_efficacy).

If the infector is vaccinated (breakthrough), infectiousness is multiplied by (1 - vaccine_infectiousness_reduction).

### 3. Intervention Logic
Isolation:

Index Case: Isolated after rash_onset + isolation_delay.

Secondary Cases: Isolated after prodromal_onset + isolation_delay.

Quarantine:

Triggered by symptomatic students.

Eligibility: Unvaccinated contacts (Vaccinated contacts are skipped).

Scope: The current C++ implementation considers all unvaccinated same-school contacts (both classmates and non-classmates) as eligible for quarantine.

Efficacy: Eligible contacts are successfully quarantined with probability quarantine_efficacy.



## --- QUARANTINE LOOP (Contact Tracing) ---
If Unvaccinated Contact is identified:
   S -> QS (Quarantined Susceptible) -> returns to S after duration
   E -> QE (Quarantined Exposed)     -> becomes QP
   P -> QP (Quarantined Prodromal)   -> becomes Iso

## Installation & Requirements
The model requires R and a C++ compiler (standard with RTools on Windows or Xcode on Mac).

## Core simulation requirements
install.packages(c("Rcpp", "dplyr", "ggplot2", "tidyr"))

## Optional: Requirements for the Shiny App UI
install.packages(c("shiny", "plotly", "readxl", "shinyBS", "shinycssloaders"))

## Assumptions & Limitations
Scope: This model represents a single school. Household, community, and social-network transmission outside the school are not modeled.

Mixing: Assumes uniform random mixing within the two contact pools (within-class and between-class).

Quarantine Logic: While the model conceptually emphasizes classmate quarantine, the current code implementation extends quarantine eligibility to all unvaccinated students in the school who are contacts of a case.

Interpretation
This model is designed for comparative scenario analysis (e.g., "How does a 1-day delay in isolation affect total cases?") rather than exact epidemiological prediction. Results should be interpreted as stochastic realizations of potential outbreak trajectories.
