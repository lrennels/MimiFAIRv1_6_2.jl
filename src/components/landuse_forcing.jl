# --------------------------------------------------
# Land use forcing from surface albedo change.
# --------------------------------------------------

@defcomp landuse_forcing begin

    α_CO₂_land       = Parameter() # Landuse radiative forcing coefficient.
    landuse_emiss    = Parameter(index=[time]) # RCP 'OtherCO2' emissions (GtC yr⁻¹).

    cumulative_emiss = Variable(index=[time])  # Cumulative 'OtherCO2' emissions (GtC).
    forcing_landuse  = Variable(index=[time])  # Forcing from land use change (Wm⁻²).


    function run_timestep(p, v, d, t)

        # Calculate cumulative land use emissions (based on RCP 'OtherCO2' emissions).
        if is_first(t)
            v.cumulative_emiss[t] = p.landuse_emiss[t]
        else
            v.cumulative_emiss[t] = v.cumulative_emiss[t-1] + p.landuse_emiss[t]
        end

        # Caluclate land use forcing (based on regression of non-fossil CO₂ emissions against AR5 land use forcing).
        v.forcing_landuse[t] = p.α_CO₂_land * v.cumulative_emiss[t]
    end
end
