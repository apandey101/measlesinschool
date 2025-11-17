# ==============================================================================
# MEASLES TRANSMISSION MODEL WITH RCPP OPTIMIZATION
# ==============================================================================
# Updating the policy...11/07/2025
# This is a C++ optimized version of the measles model
# Key improvements:
# - C++ implementation of transmission function (biggest bottleneck)
# - C++ implementation of quarantine function
# - Optimized data structures
# - Expected speedup: 5-10x faster
#
# ==============================================================================

library(Rcpp)
library(dplyr)
library(ggplot2)
library(tidyr)

# ==============================================================================
# C++ CODE FOR CRITICAL BOTTLENECKS
# ==============================================================================

sourceCpp(code = '
#include <Rcpp.h>
using namespace Rcpp;

// -----------------------------
// Contact-based transmission
// -----------------------------

// [[Rcpp::export]]
List cpp_school_transmission_contacts(
    IntegerVector student_id,
    IntegerVector class_id,
    CharacterVector state,
    LogicalVector is_vaccinated,
    LogicalVector is_isolated,
    LogicalVector is_quarantined,
    double c_within,           // mean contacts/day within class
    double c_between,          // mean contacts/day between classes
    double p_within,           // baseline per-contact transmission prob (within)
    double p_between,          // baseline per-contact transmission prob (between)
    double prodromal_mult,
    double rash_mult,
    double vaccine_reduction,  // reduction in infectiousness if infector vaccinated
    double vaccine_efficacy    // reduction in susceptibility if target vaccinated
) {
  RNGScope scope;
  int n = student_id.size();
  IntegerVector new_exposures;
  IntegerVector breakthrough_cases;
  
  // Find infectious students (P or Ra, not isolated/quarantined)
  std::vector<int> infectious_idx;
  infectious_idx.reserve(n);
  for (int i = 0; i < n; i++) {
    String s = state[i];
    if ((s == "P" || s == "Ra") && !is_isolated[i] && !is_quarantined[i]) {
      infectious_idx.push_back(i);
    }
  }
  
  if (infectious_idx.empty()) {
    return List::create(
      Named("new_exposures") = new_exposures,
      Named("breakthrough_cases") = breakthrough_cases
    );
  }
  
  // For each infectious student
  for (size_t idx = 0; idx < infectious_idx.size(); ++idx) {
    int inf_idx   = infectious_idx[idx];
    int inf_class = class_id[inf_idx];
    bool inf_vacc = is_vaccinated[inf_idx];
    String inf_state = state[inf_idx];
    
    // Determine stage multiplier
    double stage_mult = (inf_state == "Ra") ? rash_mult : prodromal_mult;
    
    // Effective per-contact probabilities for this infector
    double base_p_within  = p_within  * stage_mult * (inf_vacc ? (1.0 - vaccine_reduction) : 1.0);
    double base_p_between = p_between * stage_mult * (inf_vacc ? (1.0 - vaccine_reduction) : 1.0);
    base_p_within  = std::max(0.0, std::min(1.0, base_p_within));
    base_p_between = std::max(0.0, std::min(1.0, base_p_between));
    
    // Build target pools (people present in school; may be S, V, R, etc.)
    std::vector<int> within_targets;
    std::vector<int> between_targets;
    within_targets.reserve(50);
    between_targets.reserve(200);
    
    for (int j = 0; j < n; ++j) {
      if (j == inf_idx) continue;
      if (is_quarantined[j] || is_isolated[j]) continue;  // not present in mixing
      
      if (class_id[j] == inf_class) {
        within_targets.push_back(j);
      } else {
        between_targets.push_back(j);
      }
    }
    
    // Draw number of contacts from Poisson
    int k_within  = (c_within  > 0.0) ? R::rpois(c_within)  : 0;
    int k_between = (c_between > 0.0) ? R::rpois(c_between) : 0;
    
    // WITHIN-CLASS CONTACTS
    for (int c = 0; c < k_within; ++c) {
      if (within_targets.empty()) break;
      int pick = (int) std::floor(R::runif(0.0, (double)within_targets.size()));
      if (pick < 0) pick = 0;
      if (pick >= (int)within_targets.size()) pick = (int)within_targets.size() - 1;
      int j = within_targets[pick];
      
      String t_state = state[j];
      // Only S or V can be infected; others are "wasted" contacts
      if (!(t_state == "S" || t_state == "V")) continue;
      
      double p = base_p_within;
      if (t_state == "V") {
        p *= (1.0 - vaccine_efficacy);  // susceptibility reduction
      }
      p = std::max(0.0, std::min(1.0, p));
      
      if (R::runif(0.0, 1.0) < p) {
        new_exposures.push_back(student_id[j]);
        if (t_state == "V") {
          breakthrough_cases.push_back(student_id[j]);
        }
      }
    }
    
    // BETWEEN-CLASS CONTACTS
    for (int c = 0; c < k_between; ++c) {
      if (between_targets.empty()) break;
      int pick = (int) std::floor(R::runif(0.0, (double)between_targets.size()));
      if (pick < 0) pick = 0;
      if (pick >= (int)between_targets.size()) pick = (int)between_targets.size() - 1;
      int j = between_targets[pick];
      
      String t_state = state[j];
      if (!(t_state == "S" || t_state == "V")) continue;
      
      double p = base_p_between;
      if (t_state == "V") {
        p *= (1.0 - vaccine_efficacy);
      }
      p = std::max(0.0, std::min(1.0, p));
      
      if (R::runif(0.0, 1.0) < p) {
        new_exposures.push_back(student_id[j]);
        if (t_state == "V") {
          breakthrough_cases.push_back(student_id[j]);
        }
      }
    }
  }
  
  return List::create(
    Named("new_exposures") = new_exposures,
    Named("breakthrough_cases") = breakthrough_cases
  );
}


// -----------------------------
// Quarantine logic (unchanged)
// -----------------------------

// [[Rcpp::export]]
List cpp_apply_quarantine(
    IntegerVector student_id,
    IntegerVector class_id,
    CharacterVector state,
    LogicalVector is_quarantined,
    LogicalVector is_vaccinated,        // unchanged
    LogicalVector contacts_traced,      // track which symptomatic sources were already traced
    double quarantine_efficacy
) {
  int n = student_id.size();
  IntegerVector quarantine_ids;
  CharacterVector quarantine_states;

  IntegerVector traced_symptom_ids;

  // Find symptomatic students (trigger tracing) but only those not yet traced
  std::vector<int> symptomatic_idx;
  for (int i = 0; i < n; i++) {
    String s = state[i];
    // Only trigger tracing when the case is isolated (Iso) and not yet traced
    if ((s == "Iso") && !contacts_traced[i]) {
      symptomatic_idx.push_back(i);
    }
  }

  if (symptomatic_idx.size() == 0) {
    return List::create(
      Named("quarantine_ids")   = quarantine_ids,
      Named("quarantine_states")= quarantine_states
    );
  }

  // For each symptomatic student, quarantine UNVACCINATED classmates (existing logic)
  // and also quarantine UNVACCINATED same-school non-classmates (between-class contacts)
  // After processing each symptomatic student, record them as traced (so we do not
  // re-trace the same symptomatic source on subsequent days).
  for (size_t i = 0; i < symptomatic_idx.size(); i++) {
    int symp_idx   = symptomatic_idx[i];
    int symp_class = class_id[symp_idx];

    // record this symptomatic student id to mark as traced by the caller
    traced_symptom_ids.push_back(student_id[symp_idx]);

    for (int j = 0; j < n; j++) {
      if (j == symp_idx) continue;
      if (is_quarantined[j]) continue;

      // Skip vaccinated contacts entirely
      if (is_vaccinated[j]) continue;

      // Determine whether this contact is a classmate or a same-school non-classmate
      bool is_classmate = (class_id[j] == symp_class);
      bool is_same_school_non_classmate = (class_id[j] != symp_class);

      if (!is_classmate && !is_same_school_non_classmate) {
        // In current single-school setup this should not happen, but keep safe check
        continue;
      }

      // Apply quarantine with the same probability for both classmates and
      // same-school non-classmates (unvaccinated only).
      if (R::runif(0, 1) < quarantine_efficacy) {
        String current_state = state[j];
        if (current_state == "S") {
          quarantine_ids.push_back(student_id[j]);
          quarantine_states.push_back("QS");
        } else if (current_state == "E") {
          quarantine_ids.push_back(student_id[j]);
          quarantine_states.push_back("QE");
        } else if (current_state == "P") {
          quarantine_ids.push_back(student_id[j]);
          quarantine_states.push_back("QP");
        }
      }
    }
  }

  return List::create(
    Named("quarantine_ids")   = quarantine_ids,
    Named("quarantine_states")= quarantine_states,
    Named("traced_symptom_ids") = traced_symptom_ids
  );
}

')


# ==============================================================================
# R WRAPPER FUNCTIONS (Call C++ code)
# ==============================================================================

school_transmission <- function(population, params) {
  result <- cpp_school_transmission_contacts(
    student_id      = population$student_id,
    class_id        = population$class_id,
    state           = population$state,
    is_vaccinated   = population$is_vaccinated,
    is_isolated     = population$is_isolated,
    is_quarantined  = population$is_quarantined,
    c_within        = params$c_within,
    c_between       = params$c_between,
    p_within        = params$p_within,
    p_between       = params$p_between,
    prodromal_mult  = params$prodromal_infectiousness_multiplier,
    rash_mult       = params$rash_infectiousness_multiplier,
    vaccine_reduction = params$vaccine_infectiousness_reduction,
    vaccine_efficacy  = params$vaccine_efficacy
  )
  
  if (length(result$new_exposures) > 0) {
    exposure_idx <- match(unique(result$new_exposures), population$student_id)
    population$state[exposure_idx] <- "E"
    population$time_in_state[exposure_idx] <- 0
  }
  if (length(result$breakthrough_cases) > 0) {
    breakthrough_idx <- match(unique(result$breakthrough_cases), population$student_id)
    population$breakthrough_infection[breakthrough_idx] <- TRUE
  }
  population
}



apply_quarantine <- function(population, params) {
  if (params$quarantine_contacts == FALSE || params$quarantine_efficacy == 0) {
    return(population)
  }
  
  result <- cpp_apply_quarantine(
    student_id       = population$student_id,
    class_id         = population$class_id,
    state            = population$state,
    is_quarantined   = population$is_quarantined,
    is_vaccinated    = population$is_vaccinated,   # <-- NEW
    contacts_traced  = population$contacts_traced,
    quarantine_efficacy = params$quarantine_efficacy
  )
  
  if (length(result$quarantine_ids) > 0) {
    for (i in seq_along(result$quarantine_ids)) {
      idx <- which(population$student_id == result$quarantine_ids[i])
      if (length(idx) > 0) {
        population$state[idx] <- result$quarantine_states[i]
        population$is_quarantined[idx] <- TRUE
        population$time_in_state[idx] <- 0
      }
    }
  }

  # Mark symptomatic students whose contacts were traced so we don't re-trace them
  if (!is.null(result$traced_symptom_ids) && length(result$traced_symptom_ids) > 0) {
    traced_idx <- match(unique(result$traced_symptom_ids), population$student_id)
    traced_idx <- traced_idx[!is.na(traced_idx)]
    if (length(traced_idx) > 0) population$contacts_traced[traced_idx] <- TRUE
  }
  
  return(population)
}


# ==============================================================================
# STANDARD R FUNCTIONS (No changes needed)
# ==============================================================================

create_school_population <- function(school_size, avg_class_size, age_range) {
  n_classes <- ceiling(school_size / avg_class_size)
  
  population <- data.frame(
    student_id = 1:school_size,
    class_id = rep(1:n_classes, length.out = school_size),
    age = sample(age_range[1]:age_range[2], school_size, replace = TRUE)
  )
  
  population$state <- "S"
  population$time_in_state <- 0
  population$is_vaccinated <- FALSE
  population$breakthrough_infection <- FALSE
  population$is_isolated <- FALSE
  population$is_quarantined <- FALSE
  population$is_index <- FALSE              # <-- NEW
  population$contacts_traced <- FALSE       # <-- NEW: track whether symptomatic sources were traced
  
  return(population)
}


initialize_vaccination <- function(population, vaccination_coverage, initial_infected) {
  if (vaccination_coverage > 0) {
    n_vaccinated <- round(nrow(population) * vaccination_coverage)
    if (n_vaccinated > 0) {
      vaccinated_ids <- sample(1:nrow(population), n_vaccinated)
      population$is_vaccinated[vaccinated_ids] <- TRUE
      population$state[vaccinated_ids] <- "V"
    }
  }
  
  susceptible_ids <- which(population$state == "S")
  if (length(susceptible_ids) >= initial_infected && initial_infected > 0) {
    infected_ids <- sample(susceptible_ids, initial_infected)
    # DESIGNATE 1 INDEX: start in P, no isolation yet
    index_id <- infected_ids[1]
    population$state[index_id] <- "P"
    population$time_in_state[index_id] <- 0
    population$is_index[index_id] <- TRUE
    
    # Any remaining initial infected (if any) start as E
    if (length(infected_ids) > 1) {
      others <- infected_ids[-1]
      population$state[others] <- "E"
      population$time_in_state[others] <- 0
    }
  } else if (initial_infected > 0) {
    warning("Not enough susceptible individuals")
  }
  
  return(population)
}


update_disease_states <- function(population, params) {
  ## E -> P
  exposed_ready <- which(population$state == "E" &
                           population$time_in_state >= params$latent_period)
  if (length(exposed_ready) > 0) {
    population$state[exposed_ready] <- "P"
    population$time_in_state[exposed_ready] <- 0  # start prodromal clock
  }
  
  ## QE -> QP  (prodromal while quarantined; start prodromal clock)
  qe_ready <- which(population$state == "QE" &
                      population$time_in_state >= params$latent_period)
  if (length(qe_ready) > 0) {
    population$state[qe_ready] <- "QP"
    population$time_in_state[qe_ready] <- 0  # start prodromal clock
  }
  
  ## P -> Ra
  p_ready <- which(population$state == "P" &
                     population$time_in_state >= params$prodromal_period)
  if (length(p_ready) > 0) {
    # Split index vs non-index
    p_idx_index    <- p_ready[ population$is_index[p_ready] ]
    p_idx_nonindex <- p_ready[ !population$is_index[p_ready] ]
    
    # Index: progress to Rash and RESET timer (delay measured from rash)
    if (length(p_idx_index) > 0) {
      population$state[p_idx_index] <- "Ra"
      population$time_in_state[p_idx_index] <- 0
    }
    
    # Non-index: progress to Rash but DO NOT reset timer
    # (keep counting from prodromal onset to enforce "delay from P")
    if (length(p_idx_nonindex) > 0) {
      population$state[p_idx_nonindex] <- "Ra"
      # population$time_in_state[p_idx_nonindex] stays as-is
    }
  }
  
  ## Isolation rules
  # 1) INDEX: isolate AFTER rash + isolation_delay
  rash_index_ready <- which(population$state == "Ra" &
                              population$is_index &
                              population$time_in_state >= params$isolation_delay)
  if (length(rash_index_ready) > 0) {
    population$state[rash_index_ready] <- "Iso"
    population$is_isolated[rash_index_ready] <- TRUE
    population$time_in_state[rash_index_ready] <- 0
  }
  
  # 2) NON-INDEX: isolate when time since P onset >= isolation_delay,
  #    regardless of still in P or already in Ra (we kept the clock running).
  nonindex_delay_ready <- which(!population$is_index &
                                  (population$state == "P" | population$state == "Ra") &
                                  population$time_in_state >= params$isolation_delay)
  if (length(nonindex_delay_ready) > 0) {
    population$state[nonindex_delay_ready] <- "Iso"
    population$is_isolated[nonindex_delay_ready] <- TRUE
    population$time_in_state[nonindex_delay_ready] <- 0
  }
  
  # 3) QP: isolate after the same delay from prodromal onset (no extra wait)
  qp_delay_ready <- which(population$state == "QP" &
                            population$time_in_state >= params$isolation_delay)
  if (length(qp_delay_ready) > 0) {
    population$state[qp_delay_ready] <- "Iso"
    population$is_isolated[qp_delay_ready] <- TRUE
    population$time_in_state[qp_delay_ready] <- 0
  }
  
  ## Iso -> R  (leave as-is; adjust if your definition of isolation duration changes)
  iso_duration <- (params$rash_period - params$isolation_delay) + params$isolation_period
  isolated_ready <- which(population$state == "Iso" &
                            population$time_in_state >= iso_duration)
  if (length(isolated_ready) > 0) {
    population$state[isolated_ready] <- "R"
    population$is_isolated[isolated_ready] <- FALSE
    population$time_in_state[isolated_ready] <- 0
  }
  
  ## QS release
  qs_ready <- which(population$state == "QS" &
                      population$time_in_state >= params$quarantine_duration)
  if (length(qs_ready) > 0) {
    population$state[qs_ready] <- "S"
    population$is_quarantined[qs_ready] <- FALSE
    population$time_in_state[qs_ready] <- 0
  }
  
  ## QV release (kept for completeness)
  qv_ready <- which(population$state == "QV" &
                      population$time_in_state >= params$quarantine_duration)
  if (length(qv_ready) > 0) {
    population$state[qv_ready] <- "V"
    population$is_quarantined[qv_ready] <- FALSE
    population$time_in_state[qv_ready] <- 0
  }
  
  ## Increment time in state
  population$time_in_state <- population$time_in_state + 1
  population
}


run_single_simulation <- function(school_size, avg_class_size, age_range,
                                  vaccination_coverage, params, 
                                  initial_infected, n_days, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  population <- create_school_population(school_size, avg_class_size, age_range)
  population <- initialize_vaccination(population, vaccination_coverage, initial_infected)
  
  state_cols <- c("S", "E", "P", "Ra", "Iso", "R", "V", "QS", "QE", "QP")
  daily_counts <- matrix(0, nrow = n_days, ncol = length(state_cols))
  colnames(daily_counts) <- state_cols
  
  actual_days <- n_days
  outbreak_ended <- FALSE
  
  for (day in 1:n_days) {
    if (!outbreak_ended) {
      # Transmission happens using the current day's present students
      population <- school_transmission(population, params)

      # Progress disease states (E->P, P->Ra, isolate cases, release from quarantine, etc.)
      population <- update_disease_states(population, params)

      # Apply quarantine after states are updated so tracing is triggered when
      # cases have reached the isolated state (`Iso`). This ensures quarantine
      # is triggered at the time of isolation rather than during prodromal phase.
      if (params$quarantine_contacts) {
        population <- apply_quarantine(population, params)
      }
    }
    
    daily_counts[day, ] <- table(factor(population$state, levels = state_cols))
    
    if (!outbreak_ended) {
      active_states <- c("E", "P", "Ra", "QE", "QP", "Iso")
      if (sum(population$state %in% active_states) == 0) {
        actual_days <- day
        outbreak_ended <- TRUE
      }
    }
  }
  
  attack_rate <- sum(population$state %in% c("Iso", "R")) / school_size
  total_infected <- sum(population$state %in% c("P", "Ra", "Iso", "R", "QP"))
  breakthrough_infections <- sum(population$breakthrough_infection)
  total_isolated <- sum(population$state == "Iso") + sum(population$state == "R")
  total_quarantined <- sum(population$state %in% c("QS", "QE", "QP"))
  
  results <- list(
    daily_counts = as.data.frame(daily_counts) %>%
      dplyr::mutate(day = 0:(nrow(daily_counts) - 1)),
    attack_rate = attack_rate,
    total_infected = total_infected,
    total_isolated = total_isolated,
    total_quarantined = total_quarantined,
    breakthrough_infections = breakthrough_infections,
    actual_days = actual_days
  )
  
  results
}

run_multiple_simulations <- function(n_simulations, school_size, avg_class_size,
                                     age_range, vaccination_coverage,
                                     latent_period, prodromal_period, rash_period,
                                     isolation_delay, isolation_period,
                                     vaccine_efficacy, vaccine_infectiousness_reduction,
                                     prodromal_infectiousness_multiplier,
                                     rash_infectiousness_multiplier,
                                     # --- NEW: contact-rate parameters ---
                                     c_within, c_between,
                                     p_within, p_between,
                                     # ------------------------------------
                                     quarantine_contacts, quarantine_efficacy,
                                     quarantine_duration, initial_infected,
                                     n_days, seed_start, verbose) {
  
  params <- list(
    latent_period = latent_period,
    prodromal_period = prodromal_period,
    rash_period = rash_period,
    isolation_delay = isolation_delay,
    isolation_period = isolation_period,
    
    vaccine_efficacy = vaccine_efficacy,
    vaccine_infectiousness_reduction = vaccine_infectiousness_reduction,
    prodromal_infectiousness_multiplier = prodromal_infectiousness_multiplier,
    rash_infectiousness_multiplier = rash_infectiousness_multiplier,
    
    # --- NEW: contact-rate params passed to cpp_school_transmission_contacts ---
    c_within  = c_within,
    c_between = c_between,
    p_within  = p_within,
    p_between = p_between,
    
    # Quarantine
    quarantine_contacts   = quarantine_contacts,
    quarantine_efficacy   = quarantine_efficacy,
    quarantine_duration   = quarantine_duration
  )
  
  if (verbose) {
    cat("=== RCPP-OPTIMIZED MEASLES SIMULATION (contact-rate) ===\n")
    cat(sprintf("Number of simulations: %d\n", n_simulations))
    cat(sprintf("School size: %d students\n", school_size))
    cat("Using c_within, c_between, p_within, p_between\n\n")
  }
  
  all_daily_counts <- list()
  summary_stats <- data.frame(
    sim_id = 1:n_simulations,
    attack_rate = numeric(n_simulations),
    total_infected = numeric(n_simulations),
    total_isolated = numeric(n_simulations),
    total_quarantined = numeric(n_simulations),
    breakthrough_infections = numeric(n_simulations),
    actual_days = numeric(n_simulations)
  )
  
  start_time <- Sys.time()
  
  for (i in 1:n_simulations) {
    sim_seed <- if (!is.null(seed_start)) seed_start + i else NULL
    
    result <- run_single_simulation(
      school_size = school_size,
      avg_class_size = avg_class_size,
      age_range = age_range,
      vaccination_coverage = vaccination_coverage,
      params = params,
      initial_infected = initial_infected,
      n_days = n_days,
      seed = sim_seed
    )
    
    result$daily_counts$sim_id <- i
    all_daily_counts[[i]] <- result$daily_counts
    
    summary_stats$attack_rate[i] <- result$attack_rate
    summary_stats$total_infected[i] <- result$total_infected
    summary_stats$total_isolated[i] <- result$total_isolated
    summary_stats$total_quarantined[i] <- result$total_quarantined
    summary_stats$breakthrough_infections[i] <- result$breakthrough_infections
    summary_stats$actual_days[i] <- result$actual_days
    
    if (verbose && i %% 10 == 0) {
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      est_total <- elapsed / i * n_simulations
      cat(sprintf("Completed %d/%d (%.1f%%) - Est. total time: %.1f sec\n", 
                  i, n_simulations, 100*i/n_simulations, est_total))
    }
  }
  
  all_daily_data <- bind_rows(all_daily_counts)
  
 
  # # ============================================================================
  
  end_time <- Sys.time()
  total_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  if (verbose) {
    cat(sprintf("\n=== COMPLETED IN %.2f SECONDS ===\n", total_time))
    cat(sprintf("Average: %.3f seconds per simulation\n", total_time/n_simulations))
    cat("\n=== SUMMARY STATISTICS ===\n")
    cat(sprintf("Attack Rate: %.1f%% (95%% CI: %.1f%% - %.1f%%)\n",
                100 * median(summary_stats$attack_rate),
                100 * quantile(summary_stats$attack_rate, 0.025),
                100 * quantile(summary_stats$attack_rate, 0.975)))
  }
  
  results <- list(
    all_daily_data = all_daily_data,
    summary_stats = summary_stats,
    params = params,
    school_size = school_size,
    n_simulations = n_simulations,
    n_days = n_days,
    computation_time = total_time
  )
  
  return(results)
}
