# -----------------------------------------------------
# Concentrations of Ozone-Depleting Substances
# -----------------------------------------------------

@defcomp o3_depleting_substance_cycles begin

    ozone_depleting_substances = Index() # Index for ozone-depleting substances.

    τ_ods                = Parameter(index=[ozone_depleting_substances])                       # Atmospheric (e-folding) lifetime for each gas (years).
    ods_0                = Parameter(index=[ozone_depleting_substances])                       # Initial (pre-industrial) concentration for each gas (ppt).
    emiss2conc_ods       = Parameter(index=[ozone_depleting_substances])                       # Conversion between ppt concentrations and kt emissions.
    emiss_ods            = Parameter(index=[time, ozone_depleting_substances])                 # Emissions for other well-mixed greenhouse gases (kt yr⁻¹).

    conc_ods             = Variable(index=[time, ozone_depleting_substances])                  # Atmospheric concentrations for ozone-depleting substances (ppt).
    #conc_ods                   = Variable(index=[time, ozone_depleting_substances]) # Concentrations for ozone-depleting substances (a subset of 'conc_other_ghg' - used in stratospheric O₃ component).


    function run_timestep(p, v, d, t)

        for g in d.ozone_depleting_substances

            # Set initial concentration values.
            if is_first(t)
                v.conc_ods[t,g] = p.ods_0[g]
            else
                # Calculate concentrations based on simple one-box exponential decay model.
                v.conc_ods[t,g] = v.conc_ods[t-1,g] - v.conc_ods[t-1,g] * (1.0 - exp(-1/p.τ_ods[g])) + 0.5 * (p.emiss_ods[t-1,g] + p.emiss_ods[t,g]) * (1.0/p.emiss2conc_ods[g])
            end
        end

        # For convenience reasons, create a subset variable for ozone-depleting substances (indices 13-28) used in other model components.
        # TODO: Remove this type of indexing (leftover from orifinal FAIR code).
        #v.conc_ods[t,:] = v.conc_other_ghg[t,13:28]
    end
end