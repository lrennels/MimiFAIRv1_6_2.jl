# ---------------------------------------------------------------
# Radiative forcing from methane
# ---------------------------------------------------------------

# Note: Modified Etminan relationship from Meinshausen et al 2019 https://gmd.copernicus.org/preprints/gmd-2019-222/gmd-2019-222.pdf (Table 3)

@defcomp ch4_forcing begin
    a₃           = Parameter()             # Methane forcing coefficient.
    b₃           = Parameter()             # Methane forcing coefficient.
    d₃           = Parameter()             # Methane forcing coefficient.
    CH₄_pi        = Parameter()             # Initial (pre-industrial) methane concentration (ppb).
    #rf_scale_CH₄ = Parameter()             # Scaling factor to capture effective radiative forcing uncertainty.
    h2o_from_ch4 = Parameter()             # Radiative forcing scaling factor for water capor from methane.
    N₂O          = Parameter(index=[time]) # Atmospheric nitrous oxide concentration (ppb).
    CH₄          = Parameter(index=[time]) # Atmospheric methane concentration (ppb).

    rf_ch4       = Variable(index=[time])  # Forcing from atmospheric methane concentrations (W m⁻²).
    rf_ch4_h2o   = Variable(index=[time])  # Forcing from stratospheric water vapor due to oxidation of CH₄ to HₔO (W m⁻²).


    function run_timestep(p, v, d, t)

        # Calculate methane radiative forcing.
        v.rf_ch4[t] = (p.a₃ * sqrt(p.CH₄[t]) + p.b₃ * sqrt(p.N₂O[t]) + p.d₃) * (sqrt(p.CH₄[t]) - sqrt(p.CH₄_pi))

        # Calculate stratospheric water vapor forcing from oxidation of CH₄ to HₔO.
        v.rf_ch4_h2o[t] = v.rf_ch4[t] * p.h2o_from_ch4
    end
end
