# ==============================================================================
# MEASLES TRANSMISSION MODEL WITH IMPROVED QUARANTINE (Contact History)
# ==============================================================================
# Updated:
#
# Added no_intervention option for baseline scenario simulation
# Fixed: Corrected P->Ra time_in_state reset for no_intervention scenario
# Update: Separate isolation delays for index case (after rash) and 
#              secondary cases (after prodromal onset)
# Fixed: Quarantine efficacy now applied per unique contact person (not per
#           contact event), preventing inflated effective quarantine rates
# Fixed: Iso -> R transition now based on remaining infectious period,
#           not a fixed additional period. Added time_since_prodromal tracking.
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
#include <set>
using namespace Rcpp;

// -----------------------------
// Contact-based transmission WITH CONTACT RECORDING
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
  
  // NEW: Track ALL contacts (for quarantine purposes)
  IntegerVector contact_infector_ids;
  IntegerVector contact_target_ids;
  
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
      Named("breakthrough_cases") = breakthrough_cases,
      Named("contact_infector_ids") = contact_infector_ids,
      Named("contact_target_ids") = contact_target_ids
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
    
    // Build target pools (people present in school)
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
      
      // RECORD CONTACT (regardless of transmission outcome)
      contact_infector_ids.push_back(student_id[inf_idx]);
      contact_target_ids.push_back(student_id[j]);
      
      String t_state = state[j];
      // Only S or V can be infected
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
      
      // RECORD CONTACT
      contact_infector_ids.push_back(student_id[inf_idx]);
      contact_target_ids.push_back(student_id[j]);
      
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
    Named("breakthrough_cases") = breakthrough_cases,
    Named("contact_infector_ids") = contact_infector_ids,
    Named("contact_target_ids") = contact_target_ids
  );
}


// -----------------------------
// IMPROVED QUARANTINE: Uses contact history
// FIX V2.3: Apply efficacy per UNIQUE contact, not per contact event
// -----------------------------

// [[Rcpp::export]]
List cpp_apply_quarantine_with_history(
    IntegerVector student_id,
    CharacterVector state,
    LogicalVector is_quarantined,
    LogicalVector is_vaccinated,
    LogicalVector newly_isolated,        // which students just became isolated
    IntegerVector contact_history_infector,  // historical contacts: infector IDs
    IntegerVector contact_history_target,    // historical contacts: target IDs
    double quarantine_efficacy
) {
  int n = student_id.size();
  IntegerVector quarantine_ids;
  CharacterVector quarantine_states;
  
  // Find students who just became isolated
  std::vector<int> isolated_idx;
  for (int i = 0; i < n; i++) {
    if (newly_isolated[i]) {
      isolated_idx.push_back(i);
    }
  }
  
  if (isolated_idx.size() == 0) {
    return List::create(
      Named("quarantine_ids")   = quarantine_ids,
      Named("quarantine_states")= quarantine_states
    );
  }
  
  // For each newly isolated student, find their UNIQUE historical contacts
  for (size_t i = 0; i < isolated_idx.size(); i++) {
    int iso_idx = isolated_idx[i];
    int iso_id  = student_id[iso_idx];
    
    // FIXED: Use a set to collect UNIQUE contact IDs first
    std::set<int> unique_contact_ids;
    
    // Find all unique contacts where this person was the infector
    for (int c = 0; c < contact_history_infector.size(); c++) {
      if (contact_history_infector[c] == iso_id) {
        unique_contact_ids.insert(contact_history_target[c]);
      }
    }
    
    // Now apply quarantine efficacy ONCE per unique contact
    for (std::set<int>::iterator it = unique_contact_ids.begin(); 
         it != unique_contact_ids.end(); ++it) {
      int contact_id = *it;
      
      // Find this contact in the population
      int contact_idx = -1;
      for (int j = 0; j < n; j++) {
        if (student_id[j] == contact_id) {
          contact_idx = j;
          break;
        }
      }
      
      if (contact_idx == -1) continue;
      if (is_quarantined[contact_idx]) continue;
      
      // Skip vaccinated contacts
      if (is_vaccinated[contact_idx]) continue;
      
      // Apply quarantine with efficacy - ONCE per unique contact
      if (R::runif(0, 1) < quarantine_efficacy) {
        String current_state = state[contact_idx];
        if (current_state == "S") {
          quarantine_ids.push_back(contact_id);
          quarantine_states.push_back("QS");
        } else if (current_state == "E") {
          quarantine_ids.push_back(contact_id);
          quarantine_states.push_back("QE");
        } else if (current_state == "P") {
          quarantine_ids.push_back(contact_id);
          quarantine_states.push_back("QP");
        }
      }
    }
  }
  
  return List::create(
    Named("quarantine_ids")   = quarantine_ids,
    Named("quarantine_states")= quarantine_states
  );
}

')


# ==============================================================================
# R WRAPPER FUNCTIONS
# ==============================================================================

school_transmission <- function(population, params, contact_history) {
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
  
  # Process new exposures
  if (length(result$new_exposures) > 0) {
    exposure_idx <- match(unique(result$new_exposures), population$student_id)
    population$state[exposure_idx] <- "E"
    population$time_in_state[exposure_idx] <- 0
    population$time_since_prodromal[exposure_idx] <- NA  # Reset for new exposures
  }
  if (length(result$breakthrough_cases) > 0) {
    breakthrough_idx <- match(unique(result$breakthrough_cases), population$student_id)
    population$breakthrough_infection[breakthrough_idx] <- TRUE
  }
  
  # Store today's contacts in the history
  if (length(result$contact_infector_ids) > 0) {
    contact_history$add_contacts(
      infector_ids = result$contact_infector_ids,
      target_ids = result$contact_target_ids
    )
  }
  
  return(list(population = population, contact_history = contact_history))
}


apply_quarantine <- function(population, params, contact_history) {
  if (params$quarantine_contacts == FALSE || params$quarantine_efficacy == 0) {
    return(population)
  }
  
  # Identify newly isolated individuals (those who JUST transitioned to Iso this step)
  newly_isolated <- population$newly_isolated
  
  if (sum(newly_isolated) == 0) {
    return(population)
  }
  
  # Get contact history for the isolation delay period
  history_contacts <- contact_history$get_all_contacts()
  
  result <- cpp_apply_quarantine_with_history(
    student_id       = population$student_id,
    state            = population$state,
    is_quarantined   = population$is_quarantined,
    is_vaccinated    = population$is_vaccinated,
    newly_isolated   = newly_isolated,
    contact_history_infector = history_contacts$infector_ids,
    contact_history_target   = history_contacts$target_ids,
    quarantine_efficacy = params$quarantine_efficacy
  )
  
  # Apply quarantine
  if (length(result$quarantine_ids) > 0) {
    for (i in seq_along(result$quarantine_ids)) {
      idx <- which(population$student_id == result$quarantine_ids[i])
      if (length(idx) > 0) {
        population$state[idx] <- result$quarantine_states[i]
        population$is_quarantined[idx] <- TRUE
        population$time_in_state[idx] <- 0
        # Note: time_since_prodromal continues to track infectious period in QP
      }
    }
  }
  
  # Reset the newly_isolated flag
  population$newly_isolated <- FALSE
  
  return(population)
}


# ==============================================================================
# CONTACT HISTORY MANAGER (R6 class for cleaner state management)
# ==============================================================================

ContactHistory <- R6::R6Class("ContactHistory",
                              public = list(
                                window_size = NULL,
                                contact_list = NULL,
                                
                                initialize = function(window_size) {
                                  self$window_size <- window_size
                                  self$contact_list <- list()
                                },
                                
                                add_contacts = function(infector_ids, target_ids) {
                                  # Add today's contacts
                                  self$contact_list <- append(self$contact_list, list(
                                    data.frame(
                                      infector_id = infector_ids,
                                      target_id = target_ids,
                                      stringsAsFactors = FALSE
                                    )
                                  ))
                                  
                                  # Keep only the last window_size days
                                  if (length(self$contact_list) > self$window_size) {
                                    self$contact_list <- self$contact_list[(length(self$contact_list) - self$window_size + 1):length(self$contact_list)]
                                  }
                                },
                                
                                get_all_contacts = function() {
                                  if (length(self$contact_list) == 0) {
                                    return(list(infector_ids = integer(0), target_ids = integer(0)))
                                  }
                                  
                                  all_contacts <- bind_rows(self$contact_list)
                                  list(
                                    infector_ids = all_contacts$infector_id,
                                    target_ids = all_contacts$target_id
                                  )
                                },
                                
                                clear = function() {
                                  self$contact_list <- list()
                                }
                              )
)


# ==============================================================================
# STANDARD R FUNCTIONS
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
  population$time_since_prodromal <- NA  # NEW: tracks time since entering P (for Iso -> R)
  population$is_vaccinated <- FALSE
  population$breakthrough_infection <- FALSE
  population$is_isolated <- FALSE
  population$is_quarantined <- FALSE
  population$is_index <- FALSE
  population$newly_isolated <- FALSE
  
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
    index_id <- infected_ids[1]
    population$state[index_id] <- "P"
    population$time_in_state[index_id] <- 0
    population$time_since_prodromal[index_id] <- 0  # Index case starts in P
    population$is_index[index_id] <- TRUE
    
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
  # Reset newly_isolated flag at the start of each update
  population$newly_isolated <- FALSE
  
  # Total infectious period (for Iso -> R calculation)
  total_infectious_period <- params$prodromal_period + params$rash_period
  
  ## E -> P
  exposed_ready <- which(population$state == "E" &
                           population$time_in_state >= params$latent_period)
  if (length(exposed_ready) > 0) {
    population$state[exposed_ready] <- "P"
    population$time_in_state[exposed_ready] <- 0
    population$time_since_prodromal[exposed_ready] <- 0  # Start tracking infectious time
  }
  
  ## QE -> QP
  qe_ready <- which(population$state == "QE" &
                      population$time_in_state >= params$latent_period)
  if (length(qe_ready) > 0) {
    population$state[qe_ready] <- "QP"
    population$time_in_state[qe_ready] <- 0
    population$time_since_prodromal[qe_ready] <- 0  # Start tracking infectious time
  }
  
  ## P -> Ra
  p_ready <- which(population$state == "P" &
                     population$time_in_state >= params$prodromal_period)
  if (length(p_ready) > 0) {
    p_idx_index    <- p_ready[ population$is_index[p_ready] ]
    p_idx_nonindex <- p_ready[ !population$is_index[p_ready] ]
    
    # Index case: reset time_in_state (delay measured from rash onset)
    # But keep time_since_prodromal incrementing
    if (length(p_idx_index) > 0) {
      population$state[p_idx_index] <- "Ra"
      population$time_in_state[p_idx_index] <- 0
      # time_since_prodromal continues (already at 4)
    }
    
    # Non-index cases: behavior depends on intervention scenario
    if (length(p_idx_nonindex) > 0) {
      population$state[p_idx_nonindex] <- "Ra"
      
      # For no_intervention scenario, reset time_in_state so Ra lasts full rash_period
      # For intervention scenario, keep time_in_state for isolation_delay_secondary calculation
      if (params$no_intervention) {
        population$time_in_state[p_idx_nonindex] <- 0
      }
      # time_since_prodromal continues in both cases
    }
  }
  
  # ==========================================================================
  # NO INTERVENTION SCENARIO: Direct Ra -> R after natural rash period
  # ==========================================================================
  if (params$no_intervention) {
    # Natural recovery: Ra -> R after rash_period days (no isolation)
    ra_natural_recovery <- which(population$state == "Ra" &
                                   population$time_in_state >= params$rash_period)
    if (length(ra_natural_recovery) > 0) {
      population$state[ra_natural_recovery] <- "R"
      population$time_in_state[ra_natural_recovery] <- 0
      population$time_since_prodromal[ra_natural_recovery] <- NA
    }
    
    # QP also recovers naturally after full infectious period
    qp_natural_recovery <- which(population$state == "QP" &
                                   !is.na(population$time_since_prodromal) &
                                   population$time_since_prodromal >= total_infectious_period)
    if (length(qp_natural_recovery) > 0) {
      population$state[qp_natural_recovery] <- "R"
      population$is_quarantined[qp_natural_recovery] <- FALSE
      population$time_in_state[qp_natural_recovery] <- 0
      population$time_since_prodromal[qp_natural_recovery] <- NA
    }
    
  } else {
    # ==========================================================================
    # STANDARD INTERVENTION SCENARIO: Isolation rules with separate delays
    # ==========================================================================
    
    ## Isolation rules
    # 1) INDEX CASE: isolate after rash onset + isolation_delay_index
    #    time_in_state was reset to 0 when entering Ra, so this measures time since rash
    rash_index_ready <- which(population$state == "Ra" &
                                population$is_index &
                                population$time_in_state >= params$isolation_delay_index)
    if (length(rash_index_ready) > 0) {
      population$state[rash_index_ready] <- "Iso"
      population$is_isolated[rash_index_ready] <- TRUE
      population$newly_isolated[rash_index_ready] <- TRUE
      population$time_in_state[rash_index_ready] <- 0
      # time_since_prodromal continues for tracking when to recover
    }
    
    # 2) SECONDARY CASES (NON-INDEX): isolate when time since prodromal onset >= isolation_delay_secondary
    #    Note: For non-index, time_in_state was NOT reset at P->Ra transition
    #    So time_in_state represents total time since prodromal onset
    nonindex_delay_ready <- which(!population$is_index &
                                    (population$state == "P" | population$state == "Ra") &
                                    population$time_in_state >= params$isolation_delay_secondary)
    if (length(nonindex_delay_ready) > 0) {
      population$state[nonindex_delay_ready] <- "Iso"
      population$is_isolated[nonindex_delay_ready] <- TRUE
      population$newly_isolated[nonindex_delay_ready] <- TRUE
      population$time_in_state[nonindex_delay_ready] <- 0
      # time_since_prodromal continues for tracking when to recover
    }
    
    # 3) QP: isolate after the secondary case delay (measured from prodromal onset)
    qp_delay_ready <- which(population$state == "QP" &
                              !is.na(population$time_since_prodromal) &
                              population$time_since_prodromal >= params$isolation_delay_secondary)
    if (length(qp_delay_ready) > 0) {
      population$state[qp_delay_ready] <- "Iso"
      population$is_isolated[qp_delay_ready] <- TRUE
      population$is_quarantined[qp_delay_ready] <- FALSE  # No longer quarantined, now isolated
      population$newly_isolated[qp_delay_ready] <- TRUE
      population$time_in_state[qp_delay_ready] <- 0
      # time_since_prodromal continues
    }
  }
  
  ## Iso -> R: Based on TOTAL infectious period, not a fixed isolation duration
  ## Recovery happens when time_since_prodromal >= total infectious period
  if (!params$no_intervention) {
    isolated_ready <- which(population$state == "Iso" &
                              !is.na(population$time_since_prodromal) &
                              population$time_since_prodromal >= total_infectious_period)
    if (length(isolated_ready) > 0) {
      population$state[isolated_ready] <- "R"
      population$is_isolated[isolated_ready] <- FALSE
      population$time_in_state[isolated_ready] <- 0
      population$time_since_prodromal[isolated_ready] <- NA
    }
  }
  
  ## QS release
  qs_ready <- which(population$state == "QS" &
                      population$time_in_state >= params$quarantine_duration)
  if (length(qs_ready) > 0) {
    population$state[qs_ready] <- "S"
    population$is_quarantined[qs_ready] <- FALSE
    population$time_in_state[qs_ready] <- 0
  }
  
  ## QV release
  qv_ready <- which(population$state == "QV" &
                      population$time_in_state >= params$quarantine_duration)
  if (length(qv_ready) > 0) {
    population$state[qv_ready] <- "V"
    population$is_quarantined[qv_ready] <- FALSE
    population$time_in_state[qv_ready] <- 0
  }
  
  ## Increment time in state
  population$time_in_state <- population$time_in_state + 1
  
  ## Increment time_since_prodromal for those in infectious/isolated states
  infectious_states <- c("P", "Ra", "Iso", "QP")
  in_infectious <- which(population$state %in% infectious_states & 
                           !is.na(population$time_since_prodromal))
  if (length(in_infectious) > 0) {
    population$time_since_prodromal[in_infectious] <- 
      population$time_since_prodromal[in_infectious] + 1
  }
  
  population
}


run_single_simulation <- function(school_size, avg_class_size, age_range,
                                  vaccination_coverage, params, 
                                  initial_infected, n_days, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  population <- create_school_population(school_size, avg_class_size, age_range)
  population <- initialize_vaccination(population, vaccination_coverage, initial_infected)
  
  # Initialize contact history tracker
  # Window size should cover the longest delay for contact tracing purposes
  contact_window <- max(params$isolation_delay_index + params$prodromal_period, 
                        params$isolation_delay_secondary)
  contact_history <- ContactHistory$new(window_size = contact_window)
  
  state_cols <- c("S", "E", "P", "Ra", "Iso", "R", "V", "QS", "QE", "QP")
  daily_counts <- matrix(0, nrow = n_days, ncol = length(state_cols))
  colnames(daily_counts) <- state_cols
  
  actual_days <- n_days
  outbreak_ended <- FALSE
  
  for (day in 1:n_days) {
    if (!outbreak_ended) {
      # Transmission (also records contacts)
      trans_result <- school_transmission(population, params, contact_history)
      population <- trans_result$population
      contact_history <- trans_result$contact_history
      
      # Progress disease states (marks newly_isolated flag)
      population <- update_disease_states(population, params)
      
      # Apply quarantine using historical contacts (skipped if no_intervention)
      if (params$quarantine_contacts && !params$no_intervention) {
        population <- apply_quarantine(population, params, contact_history)
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
                                     isolation_delay_index, isolation_delay_secondary,
                                     isolation_period,
                                     vaccine_efficacy, vaccine_infectiousness_reduction,
                                     prodromal_infectiousness_multiplier,
                                     rash_infectiousness_multiplier,
                                     c_within, c_between,
                                     p_within, p_between,
                                     quarantine_contacts, quarantine_efficacy,
                                     quarantine_duration, initial_infected,
                                     n_days, seed_start, verbose,
                                     no_intervention = FALSE) {
  
  params <- list(
    latent_period = latent_period,
    prodromal_period = prodromal_period,
    rash_period = rash_period,
    isolation_delay_index = isolation_delay_index,
    isolation_delay_secondary = isolation_delay_secondary,
    isolation_period = isolation_period,  # Kept for compatibility, but not used for Iso->R
    
    vaccine_efficacy = vaccine_efficacy,
    vaccine_infectiousness_reduction = vaccine_infectiousness_reduction,
    prodromal_infectiousness_multiplier = prodromal_infectiousness_multiplier,
    rash_infectiousness_multiplier = rash_infectiousness_multiplier,
    
    c_within  = c_within,
    c_between = c_between,
    p_within  = p_within,
    p_between = p_between,
    
    quarantine_contacts   = quarantine_contacts,
    quarantine_efficacy   = quarantine_efficacy,
    quarantine_duration   = quarantine_duration,
    
    no_intervention = no_intervention
  )
  
  if (verbose) {
    cat("=== MEASLES SIMULATION ===\n")
    if (no_intervention) {
      cat("*** NO INTERVENTION SCENARIO (Baseline) ***\n")
    } else {
      cat("*** WITH INTERVENTION SCENARIO ***\n")
    }
    cat(sprintf("Number of simulations: %d\n", n_simulations))
    cat(sprintf("School size: %d students\n", school_size))
    if (!no_intervention) {
      cat(sprintf("Index case isolation delay (after rash): %d days\n", isolation_delay_index))
      cat(sprintf("Secondary case isolation delay (after prodromal): %d days\n", isolation_delay_secondary))
      cat(sprintf("Quarantine contacts: %s\n", ifelse(quarantine_contacts, "Yes", "No")))
      cat(sprintf("Quarantine efficacy: %.0f%%\n", quarantine_efficacy * 100))
    }
    cat("\n")
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
    if (!no_intervention) {
      cat(sprintf("Total Quarantined: %.1f (95%% CI: %.1f - %.1f)\n",
                  median(summary_stats$total_quarantined),
                  quantile(summary_stats$total_quarantined, 0.025),
                  quantile(summary_stats$total_quarantined, 0.975)))
    }
  }
  
  results <- list(
    all_daily_data = all_daily_data,
    summary_stats = summary_stats,
    params = params,
    school_size = school_size,
    n_simulations = n_simulations,
    n_days = n_days,
    computation_time = total_time,
    no_intervention = no_intervention
  )
  
  return(results)
}