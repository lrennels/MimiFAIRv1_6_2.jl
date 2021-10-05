# ---------------------------------------------------------------
# Radiative forcing from carbon dioxide
# ---------------------------------------------------------------

# Note: Modified Etminan relationship from Meinshausen et al 2019 https://gmd.copernicus.org/preprints/gmd-2019-222/gmd-2019-222.pdf (Table 3)

@defcomp co2_forcing begin
    a₁           = Parameter()             # Carbon dioxide forcing coefficient.
    b₁           = Parameter()             # Carbon dioxide forcing coefficient.
    c₁           = Parameter()             # Carbon dioxide forcing coefficient.
    d₁           = Parameter()             # Carbon dioxide forcing coefficient.
    F2x          = Parameter()              # Radiative forcing from a doubling of CO₂ (W m⁻²).
    N₂O_pi        = Parameter()             # Initial (pre-industrial) nitrous oxide concentration (ppb).
    CO₂_pi        = Parameter()             # Initial (pre-industrial) carbon dioxide concentration (ppm).
    adjust_F2x   = Parameter{Bool}()             # Boolean for whether or not to calculate a scaling scaling term to keep forcing from 2x CO₂ consistent.
    #rf_scale_CO₂ = Parameter()             # Scaling factor to capture effective radiative forcing uncertainty.
    N₂O          = Parameter(index=[time]) # Atmospheric nitrous oxide concentration (ppb).
    CO₂          = Parameter(index=[time]) # Atmospheric carbon dioxide concentration (ppm).

    F2x_scale    = Variable()             # Radiative forcing scaling term (to keep forcing from 2x CO₂ consistent).
    CO₂_α        = Variable(index=[time])  # Radiative forcing scaling coefficient.
    N₂O_α        = Variable(index=[time])  # Radiative forcing scaling coefficient.
    CO₂_α_max        = Variable() # Concentration value where forcing coefficient scaling term reaches its maximum value (ppm).
    rf_co2       = Variable(index=[time])  # Forcing from atmospheric carbon dioxide concentrations (Wm⁻²).


    # Calculate some values required for CO₂ radiative forcing calculations.
    function init(p, v, d)

        # Calculate maximum value for forcing scaling coefficient.
        v.CO₂_α_max = p.CO₂_pi - p.b₁ / (2.0 * p.a₁)

        # "Tune the coefficient of CO₂ forcing to acheive desired F2x, using pre-industrial CO₂ and N₂O. F2x_etminan ~= 3.801."
        if p.adjust_F2x == true

            # Calculate F2x from Etminan radiative forcing equations.
            F2x_etminan = (-2.4e-7 * p.CO₂_pi^2 + 7.2e-4 * p.CO₂_pi - 2.1e-4 * p.N₂O_pi + 5.36) * log(2.0)

            # Calculate scaling term based on user-defined value for F2x.
            v.F2x_scale = p.F2x / F2x_etminan
        else

            # Set F2x_scale to 1.0 if no adjustment.
            v.F2x_scale = 1.0
        end
    end


    function run_timestep(p, v, d, t)

        # Calculate scaling term depending on whether or not CO₂ concentrations are between pre-industrial and maximum values for α.
        if p.CO₂_pi < p.CO₂[t] <= v.CO₂_α_max
            # α when concentrations are between pre-industrial and maximum.
            v.CO₂_α[t] = p.d₁ + p.a₁ * (p.CO₂[t] - p.CO₂_pi) ^ 2 + p.b₁ * (p.CO₂[t] - p.CO₂_pi)

        elseif p.CO₂[t] <= p.CO₂_pi
            # α when concentrations are less than pre-industrial.
            v.CO₂_α[t] = p.d₁
        else
            # α when concentrations are greater than maximum.
            v.CO₂_α[t] = p.d₁ - p.b₁^2 / (4.0 * p.a₁)
        end

        # Calculate α term for N₂O.
        v.N₂O_α[t] = p.c₁ * sqrt(p.N₂O[t])

        # Calculate carbon dioxide radiative forcing (with potential F2x adjustment).
        v.rf_co2[t] = (v.CO₂_α[t] + v.N₂O_α[t]) * log(p.CO₂[t] / p.CO₂_pi) * v.F2x_scale
    end
end
