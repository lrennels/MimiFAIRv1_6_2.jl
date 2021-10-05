# ---------------------------------------------------------------
# Radiative forcing from nitrous oxide
# ---------------------------------------------------------------

# Note: Modified Etminan relationship from Meinshausen et al 2019 https://gmd.copernicus.org/preprints/gmd-2019-222/gmd-2019-222.pdf (Table 3)

@defcomp n2o_forcing begin
    a₂           = Parameter()             # Nitrous oxide forcing coefficient.
    b₂           = Parameter()             # Nitrous oxide forcing coefficient.
    c₂           = Parameter()             # Nitrous oxide forcing coefficient.
    d₂           = Parameter()             # Nitrous oxide forcing coefficient.
    N₂O_pi        = Parameter()             # Initial (pre-industrial) nitrous oxide concentration (ppb).
    #rf_scale_N₂O = Parameter()             # Scaling factor to capture effective radiative forcing uncertainty.
    CO₂          = Parameter(index=[time]) # Atmospheric carbon dioxide concentration (ppm).
    N₂O          = Parameter(index=[time]) # Atmospheric nitrous oxide concentration (ppb).
    CH₄          = Parameter(index=[time]) # Atmospheric methane concentration (ppb).

    rf_n2o       = Variable(index=[time])  # Forcing from atmospheric nitrous oxide concentrations (Wm⁻²).


    function run_timestep(p, v, d, t)

        # Calculate nitrous oxide radiative forcing.
        v.rf_n2o[t] = (p.a₂ * sqrt(p.CO₂[t]) + p.b₂ * sqrt(p.N₂O[t]) + p.c₂ * sqrt(p.CH₄[t]) + p.d₂) * (sqrt(p.N₂O[t]) - sqrt(p.N₂O_pi))
    end
end
