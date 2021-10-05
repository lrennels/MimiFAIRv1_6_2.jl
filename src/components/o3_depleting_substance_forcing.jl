# ---------------------------------------------------------------------
# Radiative forcing from ozone-depleting substances.
# ---------------------------------------------------------------------

@defcomp o3_depleting_substance_forcing begin

    ozone_depleting_substances = Index()                                            # Index for ozone-depleting substances.

    ods_pi                   = Parameter(index=[ozone_depleting_substances])       # Initial (pre-industrial) concentration for ozone-depleting substances (ppt).
    ods_radiative_efficiency = Parameter(index=[ozone_depleting_substances])       # Radiative efficiencies for ozone-depleting substances - from IPCC AR5 WG1 Ch8 SM (Wm⁻²ppb⁻¹).
    conc_ods                 = Parameter(index=[time, ozone_depleting_substances]) # Atmospheric concentrations for ozone-depleting substances (ppt).

    ods_rf_total   = Variable(index=[time])             # Total radiative forcing for all ozone-depleting substances (Wm⁻²).
    ods_rf         = Variable(index=[time, ozone_depleting_substances])  # Individual radiative forcings for each ozone-depleting substances (Wm⁻²).


    function run_timestep(p, v, d, t)

        for g in d.ozone_depleting_substances
            # Caluclate radiative forcing for individual gases.
            # Note: the factor of 0.001 here is because radiative efficiencies are given in Wm⁻²ppb⁻¹ and concentrations of minor gases are in ppt.
            v.ods_rf[t,g] = (p.conc_ods[t,g] - p.ods_pi[g]) * p.ods_radiative_efficiency[g] * 0.001
        end

        # Calculate total radiative forcing as sum of individual gas forcings.
        v.ods_rf_total[t] = sum(v.ods_rf[t,:])
    end
end
