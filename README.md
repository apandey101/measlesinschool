# Measles-in-School Agent-Based Model

## Overview

This repository contains an agent-based model (ABM) simulating measles transmission within a single school environment. The model allows researchers and decision-makers to explore how school-based control measures—such as isolation timing, quarantine coverage, and vaccination thresholds—shape outbreak dynamics.

The school data on number of students and vaccination rates are taken from South Carolina Department of Public Health and represent the percentage of total students who have completed vaccine schedules for all required vaccines, including MMR (measles, mumps, rubella).

Implemented in R with core transmission loops optimized in C++ (via Rcpp), the simulation focuses on the interplay between classroom structure, vaccination gaps, and intervention delays.

## Key Features

- **Hybrid Architecture**: Uses R for data management and UI, with Rcpp for high-speed stochastic transmission and quarantine loops (sim_cr_v1.R).

- **Granular Mixing**: Models two levels of contact risk: high-intensity mixing within fixed classrooms and lower-intensity mixing between classes.

- **Disease Progression**: Tracks agents through Susceptible (S), Exposed (E), Prodromal (P), Rash (Ra) and Recovered (R) stages.

- **Differentiated Interventions**:
  - *Vaccination*: Modeled with both susceptibility reduction and reduced infectiousness for breakthrough cases.
  - *Isolation*: Distinct logic for the Index case (delayed isolation after rash) vs. subsequent cases (isolation triggered by prodromal onset).
  - *Quarantine*: Contact-history-based tracing that retrospectively identifies and quarantines unvaccinated contacts.

- **Stochasticity**: Daily contact events are drawn from Poisson distributions, allowing for the exploration of outbreak variability.

## Model Logic & Dynamics

### 1. Agent States & Transitions

The model follows an extended SEIR framework tailored for measles:

| State | Description |
|-------|-------------|
| S (Susceptible) | Unvaccinated individuals |
| V (Vaccinated) | Individuals with reduced susceptibility |
| E (Exposed) | Latent phase (non-infectious) |
| P (Prodromal) | Pre-rash infectious phase |
| Ra (Rash) | Symptomatic infectious phase |
| Iso (Isolated) | Removed from school mixing |
| R (Recovered) | Immune |
| QS/QE/QP (Quarantined) | Contacts removed from mixing |

### 2. Transmission Mechanics

Transmission is contact-based. Each infectious agent (P or Ra) generates contacts daily based on Poisson distributions:

- **Within-Class**: N ~ Poisson(c_within) targeting classmates.
- **Between-Class**: N ~ Poisson(c_between) targeting non-classmates.

Transmission probability depends on the target's status and the infector's stage:

- If the target is vaccinated, susceptibility is multiplied by (1 - vaccine_efficacy).
- If the infector is vaccinated (breakthrough), infectiousness is multiplied by (1 - vaccine_infectiousness_reduction).

### 3. Intervention Logic

#### Isolation

- **Index Case**: Isolated after `rash_onset + isolation_delay`.
- **Secondary Cases**: Isolated after `prodromal_onset + isolation_delay`.

#### Quarantine (Contact-History-Based Tracing)

The model implements retrospective contact tracing using a rolling contact history.

**Contact Recording**:

- All contacts made by infectious individuals (P or Ra states) are recorded daily, regardless of whether transmission occurred.
- Both within-class and between-class contacts are tracked.
- The contact history maintains a rolling window equal to `isolation_delay` days.

**Quarantine Trigger**:

- Quarantine is initiated when an infectious individual becomes isolated (`newly_isolated` flag).
- The model retrieves all contacts the newly-isolated individual made during the tracking window.

**Quarantine Eligibility**:

| Contact Type | Vaccination Status | Quarantined? |
|--------------|-------------------|--------------|
| Within-class (classmate) | Unvaccinated | Yes |
| Within-class (classmate) | Vaccinated | No |
| Between-class (non-classmate) | Unvaccinated | Yes |
| Between-class (non-classmate) | Vaccinated | No |

**Quarantine State Transitions**:

```
S  → QS (Quarantined Susceptible) → returns to S after quarantine_duration
E  → QE (Quarantined Exposed)     → progresses to QP
P  → QP (Quarantined Prodromal)   → progresses to Iso
```

**Quarantine Efficacy**:

- Eligible contacts are successfully quarantined with probability `quarantine_efficacy`.
- This parameter accounts for imperfect contact tracing (e.g., incomplete recall, contacts not reachable).

## Installation & Requirements

The model requires R and a C++ compiler (standard with RTools on Windows or Xcode on Mac).

```r
# Core simulation requirements
install.packages(c("Rcpp", "dplyr", "ggplot2", "tidyr", "R6"))

# Optional: Requirements for the Shiny App UI
install.packages(c("shiny", "plotly", "readxl", "shinyBS", "shinycssloaders"))
```

## Assumptions & Limitations

- **Scope**: This model represents a single school. Household, community, and social-network transmission outside the school are not modeled.

- **Mixing**: Assumes uniform random mixing within the two contact pools (within-class and between-class).

- **Quarantine Scope**: The contact-history-based quarantine extends to all unvaccinated contacts (both classmates and non-classmates) who had contact with a case during the infectious period prior to isolation. Vaccinated contacts are assumed to be at low risk and are not quarantined.

- **Contact History Window**: The contact tracing lookback period equals `isolation_delay`, representing the window during which the case was infectious but not yet isolated.

## Interpretation

This model is designed for comparative scenario analysis (e.g., "How does a 1-day delay in isolation affect total cases?") rather than exact epidemiological prediction. Results should be interpreted as stochastic realizations of potential outbreak trajectories.
