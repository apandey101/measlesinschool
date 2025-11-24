Measles-in-School Agent-Based ModelOverviewThis repository contains an agent-based model (ABM) simulating measles transmission within a single school environment. The model allows researchers and decision-makers to explore how school-based control measures—such as isolation timing, quarantine coverage, and vaccination thresholds—shape outbreak dynamics.Implemented in R with core transmission loops optimized in C++ (via Rcpp), the simulation focuses on the interplay between classroom structure, vaccination gaps, and intervention delays.Key FeaturesHybrid Architecture: Uses R for data management and UI, with Rcpp for high-speed stochastic transmission and quarantine loops (sim_cr_v1.R).Granular Mixing: Models two levels of contact risk: high-intensity mixing within fixed classrooms and lower-intensity mixing between classes.Realistic Disease Progression: Tracks agents through Susceptible ($S$), Exposed ($E$), Prodromal ($P$), and Rash ($Ra$) stages.Differentiated Interventions:Vaccination: Modeled with both susceptibility reduction and reduced infectiousness for breakthrough cases.Isolation: Distinct logic for the Index case (delayed isolation after rash) vs. subsequent cases (isolation triggered by prodromal onset).Quarantine: Targeted removal of unvaccinated contacts from the mixing pool.Stochasticity: Daily contact events are drawn from Poisson distributions, allowing for the exploration of outbreak variability.Model Logic & Dynamics1. Agent States & TransitionsThe model follows an extended SEIR framework tailored for measles:$S$ (Susceptible): Unvaccinated individuals.$V$ (Vaccinated): Individuals with reduced susceptibility.$E$ (Exposed): Latent phase (non-infectious).$P$ (Prodromal): Pre-rash infectious phase.$Ra$ (Rash): Symptomatic infectious phase.$Iso$ (Isolated): Removed from school mixing.$R$ (Recovered): Immune.$Q$ (Quarantined): Unvaccinated contacts removed from mixing (e.g., $QS$, $QE$, $QP$).2. Transmission MechanicsTransmission is contact-based. Each infectious agent ($P$ or $Ra$) generates contacts daily based on Poisson distributions:Within-Class: $N \sim Poisson(c_{within})$ targeting classmates.Between-Class: $N \sim Poisson(c_{between})$ targeting non-classmates.Transmission probability depends on the target's status and the infector's stage. If the target is vaccinated, susceptibility is multiplied by $(1 - \text{vaccine\_efficacy})$. If the infector is vaccinated (breakthrough), infectiousness is multiplied by $(1 - \text{vaccine\_infectiousness\_reduction})$.3. Intervention LogicIsolation:Index Case: Isolated after $T_{rash} + \text{isolation\_delay}$.Secondary Cases: Isolated after $T_{prodromal\_onset} + \text{isolation\_delay}$.Quarantine:Triggered by symptomatic students ($P, Ra, Iso$).Eligibility: Unvaccinated contacts (Vaccinated contacts are skipped).Scope: The current C++ implementation considers all unvaccinated same-school contacts (both classmates and non-classmates) as eligible for quarantine.Efficacy: Eligible contacts are successfully quarantined with probability $p_{quarantine}$.
graph TD
    %% Define Node Styles
    classDef healthy fill:#e1f5fe,stroke:#01579b,stroke-width:2px;
    classDef exposed fill:#fff9c4,stroke:#fbc02d,stroke-width:2px;
    classDef infectious fill:#ffccbc,stroke:#bf360c,stroke-width:2px;
    classDef isolated fill:#e0e0e0,stroke:#424242,stroke-width:2px,stroke-dasharray: 5 5;
    classDef recovered fill:#c8e6c9,stroke:#2e7d32,stroke-width:2px;

    %% --- Main Population Pools ---
    subgraph General Population
        S(Susceptible S):::healthy
        V(Vaccinated V):::healthy
        R(Recovered R):::recovered
    end

    %% --- Disease Progression Loop ---
    subgraph Disease Progression
        E(Exposed E):::exposed
        P(Prodromal P <br>Pre-rash Infectious):::infectious
        Ra(Rash Ra <br>Symptomatic Infectious):::infectious
    end

    %% --- Intervention States ---
    subgraph Interventions
        Iso(Isolated Iso <br>Removed from mixing):::isolated
        QS(Quarantined Susceptible QS):::isolated
        QE(Quarantined Exposed QE):::exposed
        QP(Quarantined Prodromal QP):::infectious
    end


    %% ==============================
    %% Transitions
    %% ==============================

    %% 1. Infection Events (Transmission)
    S -->|Contact with P or Ra| E
    V -.->|Breakthrough Contact<br>Reduced susceptibility| E

    %% 2. Natural History
    E -->|Latent period| P
    P -->|Prodromal period| Ra
    
    %% 3. Isolation Logic (Removes from P or Ra)
    Ra -- Index Case ONLY:<br>Rash onset + delay --> Iso
    P -- Subsequent Cases:<br>P onset + delay --> Iso
    Ra -- Subsequent Cases:<br>If not yet isolated by P delay --> Iso
    
    Iso -->|Isolation period ends| R

    %% 4. Quarantine Logic (Contact Tracing)
    %% Trigger: Contact with P, Ra, or Iso.
    %% Target: Unvaccinated S, E, or P.
    P -.->|Triggers Tracing| ContactEvent[Contact Event]
    Ra -.->|Triggers Tracing| ContactEvent
    Iso -.->|Triggers Tracing| ContactEvent

    ContactEvent -- Identifies Unvaccinated S --> QS
    ContactEvent -- Identifies E --> QE
    ContactEvent -- Identifies P --> QP

    %% 5. Quarantine Progression
    QS -->|Quarantine duration ends| S
    QE -->|Latent period finishes| QP
    QP -->|P onset + delay triggers isolation| Iso