Measles-in-School Model — Model Logic & Assumptions
=================================================

Overview
--------
This document summarizes the logic, assumptions, and default parameters used by the `measlesinschool` simulation (files: `sim_cr_v1.R`, `app_cr_v1.R`). The model is an agent-based, school-focused measles transmission model with classroom structure, vaccination, targeted quarantine, and isolation interventions. Critical transmission and quarantine routines are implemented in Rcpp/C++ inside `sim_cr_v1.R` for speed.

High-level assumptions
----------------------
- The model represents a single school (students only). Household and community transmission outside of school is not modeled.  
- Students are assigned to classrooms of fixed size (`FIXED_CLASS_SIZE`, default 25).  
- Time advances in daily steps.  
- Transmission occurs via random contact events drawn daily from Poisson distributions (within-class and between-class).  
- Vaccination reduces susceptibility (vaccine efficacy) and, if vaccinated individuals become infected, reduces infectiousness (vaccine infectiousness reduction).  
- Quarantine is applied only to unvaccinated contacts of symptomatic students, and quarantined students are removed from school mixing for the quarantine duration.  
- Isolation removes symptomatic cases from mixing based on policy delays.

Agent states (codes used in the simulation)
-------------------------------------------
- `S` : Susceptible (unvaccinated)  
- `V` : Vaccinated (protected state)  
- `E` : Exposed (latent, not yet infectious)  
- `P` : Prodromal (pre-rash, infectious)  
- `Ra`: Rash (infectious)  
- `Iso`: Isolated (removed from school mixing)  
- `R` : Recovered (immune)  
- `QS`: Quarantined Susceptible (removed from mixing)  
- `QE`: Quarantined Exposed  
- `QP`: Quarantined Prodromal  
- `QV`: Quarantined Vaccinated (kept for completeness; code supports QV release)

State transitions (daily rules)
-------------------------------
- `E` -> `P` after `latent_period` days.  
- `QE` -> `QP` after `latent_period` days (latent progression while quarantined).  
- `P` -> `Ra` after `prodromal_period` days.  
  - For the index case only (`is_index == TRUE`) the timer is reset on progression to `Ra` so isolation delay is measured from rash onset.  
  - For non-index cases, the product is set to `Ra` but the timer is not reset (isolation countdown measured from prodromal onset).  
- Isolation:  
  - Index case: isolate after `rash` + `isolation_delay`.  
  - Non-index cases: isolate when time since `P` onset >= `isolation_delay` (applies to P or Ra stages).  
  - `QP` also transitions to `Iso` after the same isolation delay measured from prodromal onset.  
- `Iso` -> `R` when the isolation duration ends (computed in code).  
- `QS` and `QV` return to `S` or `V` respectively after `quarantine_duration` days.

Transmission model (contact-based)
----------------------------------
- Each infectious (state `P` or `Ra`) who is not isolated or quarantined generates two contact processes daily:  
  - Within-class contacts: Poisson(`c_within`) contacts per day; targets sampled uniformly from other present classmates.  
  - Between-class contacts: Poisson(`c_between`) contacts per day; targets sampled uniformly from other present non-classmates in the school.  
- Per-contact transmission probabilities are `p_within` and `p_between` respectively.  
- Infectivity multipliers are applied by stage: `prodromal_infectiousness_multiplier` and `rash_infectiousness_multiplier`.  
- If the infector is vaccinated, their infectiousness is multiplied by `(1 - vaccine_infectiousness_reduction)`.  
- If the target is vaccinated (`V`), susceptibility is multiplied by `(1 - vaccine_efficacy)` (i.e. reduced susceptibility).  
- Contacts with non-susceptible/non-vulnerable states are "wasted" (no effect).  
- Breakthrough infections (vaccinated but infected) are tracked.

Quarantine & contact tracing logic
----------------------------------
- Triggering events: symptomatic students (states `P`, `Ra`, or `Iso`) act as index cases for contact tracing.  
- Who is eligible: unvaccinated contacts only (vaccinated contacts are skipped).  
- Which contacts are quarantined: (current code behavior)
  - Unvaccinated classmates of symptomatic students (original behavior).  
  - Additionally, unvaccinated same-school non-classmates (between-class contacts) are now also considered eligible and may be quarantined. This expands quarantine beyond classmates to all unvaccinated students in the school.  
  - Note: the UI text previously stated "classmates only" — the code has been updated to include same-school non-classmates. If you want quarantine restricted to classmates only, or separate efficacies for classmates vs between-class contacts, see "Options" below.  
- Application: each eligible (unvaccinated) contact is quarantined with probability `quarantine_efficacy` (applied per contact).  
- State mapping for quarantine:  
  - `S` -> `QS`  
  - `E` -> `QE`  
  - `P` -> `QP`  
- Quarantined individuals are marked `is_quarantined` and removed from mixing (they are not selected as targets for transmission while quarantined).  

Vaccination
-----------
- Vaccination is initialized in `initialize_vaccination()` using `vaccination_coverage` (proportion). Vaccinated students are set to state `V`.  
- Vaccination effects in transmission:  
  - Reduces susceptibility by `vaccine_efficacy` (multiplicative).  
  - Reduces infectiousness of infected vaccinated individuals by `vaccine_infectiousness_reduction`.

Key implementation notes
------------------------
- Core transmission and quarantine loops are implemented in C++ via `Rcpp::sourceCpp(code = '...')` in `sim_cr_v1.R` for performance.  
- The R wrapper functions `school_transmission()` and `apply_quarantine()` call the compiled routines and update the R `population` data.frame accordingly.  
- The simulation is stochastic; multiple runs (default `n_simulations`) are performed and summary statistics (median, 95% intervals) are reported.

Default parameter values (from `app_cr_v1.R` / UI defaults)
------------------------------------------------------------
- School / population:  
  - `FIXED_CLASS_SIZE` = 25 (students per class)  
  - `age_range` = c(5, 17) (used for synthetic population ages)  

- Disease natural history:  
  - `latent_period` = 10 days  
  - `prodromal_period` = 4 days  
  - `rash_period` = 4 days  
  - `isolation_delay` = 1 day (UI slider default)  
  - `isolation_period` = 3 days  

- Intervention / vaccination:  
  - `vaccine_efficacy` = 0.97 (97% reduction in susceptibility)  
  - `vaccine_infectiousness_reduction` = 0.80 (80% reduction in infectiousness)  

- Infectiousness multipliers:  
  - `prodromal_infectiousness_multiplier` = 1.0  
  - `rash_infectiousness_multiplier` = 0.6  

- Contact & transmission parameters (UI defaults):  
  - `c_within` = 10 (mean within-class contacts/day)  
  - `c_between` = 5 (mean between-class contacts/day)  
  - `p_within` = 0.19 (per-contact transmission probability, within class)  
  - `p_between` = 0.09 (per-contact transmission probability, between classes)  

- Quarantine / simulation:  
  - `quarantine_contacts` = TRUE (contacts tracing enabled)  
  - `quarantine_efficacy` = 0.80 (80% — UI slider default 80%)  
  - `quarantine_duration` = 21 days (UI default)  
  - `initial_infected` = 1 (UI default)  
  - `n_simulations` = 150 (advanced option default)  
  - `n_days` = 100 (simulation horizon default)  
  - `seed_start` = 1234567 (deterministic offset for RNG seeds)

Notes on randomness & reproducibility
------------------------------------
- Each simulation call sets the seed using `seed_start` + simulation index; you can set `seed_start` for reproducible batches.  
- Rcpp code uses the R RNG via `R::rpois` and `R::runif` — reproducibility requires consistent R RNG state and seed.

How to run a quick smoke test (R console)
-----------------------------------------
1. Install required packages in R if not already present:  

```r
install.packages(c("Rcpp", "dplyr", "ggplot2", "tidyr"))
# For the Shiny app: install.packages(c("shiny","plotly","readxl","shinyBS","shinycssloaders"))
```

2. In R (or RStudio) source the script and run a single simulation:  

```r
source("sim_cr_v1.R")
params <- list(
  latent_period = 10,
  prodromal_period = 4,
  rash_period = 4,
  isolation_delay = 1,
  isolation_period = 3,
  vaccine_efficacy = 0.97,
  vaccine_infectiousness_reduction = 0.80,
  prodromal_infectiousness_multiplier = 1.0,
  rash_infectiousness_multiplier = 0.6,
  c_within = 10,
  c_between = 5,
  p_within = 0.19,
  p_between = 0.09,
  quarantine_contacts = TRUE,
  quarantine_efficacy = 0.8,
  quarantine_duration = 21
)

pop <- create_school_population(school_size = 100, avg_class_size = 25, age_range = c(5,17))
pop <- initialize_vaccination(pop, vaccination_coverage = 0.9, initial_infected = 1)
res <- run_single_simulation(100, 25, c(5,17), 0.9, params, initial_infected = 1, n_days = 60, seed = 123)
str(res)
```

Known caveats & options
-----------------------
- UI text in `app_cr_v1.R` describes quarantine as "classmates only", but the compiled C++ quarantine implementation currently quarantines all unvaccinated same-school contacts (classmates and non-classmates). The code was intentionally changed to expand quarantine scope; if you prefer the UI wording to match code, update the UI text or revert the quarantine logic.  
- You can add more granular quarantine behaviour (different probabilities for classmates vs between-class contacts) by adding a new parameter (e.g., `quarantine_efficacy_between`) and passing it through the R wrapper into the C++ function.  
- The model assumes uniform mixing within the within-class and between-class pools; it does not model network structure or classrooms with different mixing intensities.  

Extending or modifying the model
--------------------------------
- Change contact structure: implement heterogeneous contact matrices or age-specific mixing by modifying the target-selection logic in `cpp_school_transmission_contacts`.  
- Change quarantine targeting: to restrict quarantine to classmates only, add a condition `if (class_id[j] != symp_class) continue;` inside the C++ `cpp_apply_quarantine` loop.  
- Add household transmission: incorporate an outside-school mixing routine and add a second layer of contacts per agent.

Contact / maintenance
---------------------
- Code author: repository owner. For changes to quarantine logic or parameters, edit `sim_cr_v1.R` (C++ block) and corresponding R wrapper arguments in `apply_quarantine()`.  

Acknowledgements
----------------
- Rcpp used to speed up contact loops.  
- The Shiny app (`app_cr_v1.R`) provides an interactive front-end for parameter exploration.

(End of document)
