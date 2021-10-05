# --------------------------------------------------
# Ozone Radiative Forcing
# --------------------------------------------------

# Note from original Python code: "Calculates total ozone forcing from precursor emissions and
# concentrations based on AerChemMIP and CMIP6 Historical behaviour [Skeie et al. (2020) &
# Thornhill et al. (2021)]. Unlike Stevenson CMIP5 no distinction is made for tropospheric and
# stratospheric. With this formulation, ozone forcing depends on concentrations of CH4, N2O,
# ozone-depleting halogens, and emissions of CO, NVMOC and NOx.

@defcomp o3_forcing begin

    ozone_depleting_substances = Index() # Index for ozone-depleting substances.

    Br                    = Parameter(index=[ozone_depleting_substances])       # Number of bromine atoms in each ozone-depleting species.
    Cl                    = Parameter(index=[ozone_depleting_substances])       # Number of chlorine atoms in each ozone-depleting species.
    FC                    = Parameter(index=[ozone_depleting_substances])       # Fractional stratospheric release for ozone-depleting substances (Reference: Daniel, J. and Velders, G.: A focus on information and options for policymakers, in: Scientific Assessment of Ozone Depletion, WMO, 2011).
    
    ods_pi                = Parameter(index=[ozone_depleting_substances])       # Pre-industrial concentrations for ozone_depleting_substances (ppt).
    conc_ODS              = Parameter(index=[time, ozone_depleting_substances]) # Atmospheric concentrations of ozone_depleting_substances (ppt).
    
    CH₄ 				  = Parameter(index=[time])# Atmospheric methane concentration (ppb).
    N₂O 				  = Parameter(index=[time])# Atmospheric nitrous oxide concentration (ppb).
   # CO₂ 				  = Parameter(index=[time])# Atmospheric carbon dioxide concentration (ppm).
    
    CH₄_pi 				  = Parameter()# Pre-industrial atmospheric methane concentration (ppb).
    N₂O_pi 			      = Parameter()# Pre-industrial atmospheric nitrous oxide concentration (ppb).
   # CO₂_pi 				  = Parameter()# Pre-industrial atmospheric carbon dioxide concentration (ppm).
	total_forcing_O₃_0    = Parameter() # Initial condition for ozone radiative forcing in first model period.

	NOx_emiss_pi      = Parameter() # Pre-industrial nitrogen oxides emissions (MtN yr⁻¹).
    CO_emiss_pi       = Parameter() # Pre-industrial carbon monoxide emissions (MtCO yr⁻¹).
    NMVOC_emiss_pi    = Parameter() # Pre-industrial global non-methane volatile organic compounds emissions (Mt yr⁻¹).
	NOx_emiss         = Parameter(index=[time]) # Nitrogen oxides emissions (MtN yr⁻¹).
    CO_emiss          = Parameter(index=[time]) # Carbon monoxide emissions (MtCO yr⁻¹).
    NMVOC_emiss       = Parameter(index=[time]) # Global non-methane volatile organic compounds emissions (Mt yr⁻¹).
   
    temperature           = Parameter(index=[time]) # Global mean surface temperature anomaly (K).
    feedback 			  = Parameter() # Temperature feedback on ozone radiative forcing (W m⁻²K⁻¹).
   
    Ψ_CH₄                 = Parameter()  # Radiative efficiency coefficient for methane concentrations (W m⁻²ppb⁻¹).
    Ψ_N₂O                 = Parameter()  # Radiative efficiency coefficient for nitrous_oxide concentrations (W m⁻²ppb⁻¹).
    Ψ_ODS                 = Parameter()  # Radiative efficiency coefficient for ozone-depleting substance concentrations in EESC (W m⁻² ppt⁻¹).
    Ψ_CO                  = Parameter()  # Radiative efficiency coefficient for carbon monoxide emissions (W m⁻² (Mt yr⁻¹)⁻¹).
    Ψ_NMVOC               = Parameter()  # Radiative efficiency coefficient for non-methane volatile organic compound emissions (Wm⁻²(Mt yr⁻¹)⁻¹).
    Ψ_NOx                 = Parameter()  # Radiative efficiency coefficient for nitrogen oxides emissions (Wm⁻²(MtN yr⁻¹)⁻¹).


    F_CH₄                 = Variable(index=[time])  # Radiative forcing contribution from methane (Wm⁻²).
    F_N₂O                 = Variable(index=[time])  # Radiative forcing contribution from nitrous_oxide (Wm⁻²).
    F_ODS                 = Variable(index=[time])  # Radiative forcing contribution from ozone-depleting substances (Wm⁻²).
    F_CO                  = Variable(index=[time])  # Radiative forcing contribution from carbon monoxide (Wm⁻²).
    F_NMVOC               = Variable(index=[time])  # Radiative forcing contribution from non-methane volatile organic compounds (Wm⁻²).
    F_NOx                 = Variable(index=[time])  # Radiative forcing contribution from nitrogen oxides (Wm⁻²).
    int_eecs                   = Variable(index=[time, ozone_depleting_substances])  # Intermediate terms to caluclate equivalent effective stratospheric chlorine (EESC) for ozone-depleting substances (just for convenience).
    eecs                       = Variable(index=[time])                              # Equivalent effective stratospheric chlorine (EESC) for ozone-depleting substances.
    total_forcing_O₃           = Variable(index=[time])                              # Radiative forcing from stratospheric ozone (Wm⁻²).


    function run_timestep(p, v, d, t)

        if is_first(t)
            v.total_forcing_O₃[t] = p.total_forcing_O₃_0
        else

            for g in d.ozone_depleting_substances
                # Calculate intermediate variables for EECS (calculations done for each individual ozone-depleting substance).
                v.int_eecs[t,g] = p.Cl[g] * (p.conc_ODS[t,g]-p.ods_pi[g]) * (p.FC[g]/p.FC[1]) + 45 * p.Br[g] * (p.conc_ODS[t,g]-p.ods_pi[g]) * (p.FC[g]/p.FC[1])
            end

            # Calcualte equivalent effective stratospheric chlorine.
            v.eecs[t] = max((p.FC[1]*sum(v.int_eecs[t,:])), 0.0)

            # Calcualte forcing contributions from individual agents.
            v.F_ODS[t]    = p.Ψ_ODS    * v.eecs[t]
            v.F_CH₄[t]    = p.Ψ_CH₄    * (p.CH₄[t] - p.CH₄_pi)
            v.F_N₂O[t]    = p.Ψ_N₂O    * (p.N₂O[t] - p.N₂O_pi)
            #v.F_CO₂[t]    = p.Ψ_CO₂    * (p.CO₂[t] - p.CO₂_pi)
            v.F_CO[t]     = p.Ψ_CO     * (p.CO_emiss[t]    - p.CO_emiss_pi)
            v.F_NMVOC[t]  = p.Ψ_NMVOC  * (p.NMVOC_emiss[t] - p.NMVOC_emiss_pi)
            v.F_NOx[t]    = p.Ψ_NOx    * (p.NOx_emiss[t]   - p.NOx_emiss_pi)

            # Calculate total ozone radiative forcing with (potential) temperature feedback.
            v.total_forcing_O₃[t] = v.F_ODS[t] + v.F_CH₄[t] + v.F_N₂O[t] + v.F_CO[t] + v.F_NMVOC[t] + v.F_NOx[t] + p.feedback * p.temperature[t-1]
        end
    end
end
