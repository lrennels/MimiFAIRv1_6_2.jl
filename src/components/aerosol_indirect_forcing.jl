# --------------------------------------------------
# Direct radiative forcing effect from aerosols.
# --------------------------------------------------

# This follows the `ghan2` indirect forcing calculations from FAIRv1.6.2

@defcomp aerosol_indirect_forcing begin

    ϕ                      = Parameter()             # Scale factor.
    b_SOx                  = Parameter()             # Sensitivity to sulfur oxides emissions.
    b_POM                  = Parameter()             # Sensitivity to black carbon + organic carbon emissions.
    SOx_emiss_pi           = Parameter()             # Pre-industiral sulfur oxides emissions (Mt yr⁻¹).
    BC_emiss_pi            = Parameter()             # Pre-industiral black carbon emissions (Mt yr⁻¹).
    OC_emiss_pi            = Parameter()             # Pre-industrial organic carbon emissions (Mt yr⁻¹).
    rf_scale_aero_indirect = Parameter()             # Scaling factor to capture effective radiative forcing uncertainty.
    SOx_emiss              = Parameter(index=[time]) # Sulfur oxides emissions (MtS yr⁻¹).
    BC_emiss               = Parameter(index=[time]) # Black carbon emissions (Mt yr⁻¹).
    OC_emiss               = Parameter(index=[time]) # Organic carbon emissions (Mt yr⁻¹).

    pd_re       = Variable(index=[time])  # Present-day indirect radiative forcing cntribution from aerosols (Wm⁻²).
    pi_re       = Variable(index=[time])  # Pr-industrial indirect radiative forcing contribution from aerosols (Wm⁻²).
    rf_aero_indirect       = Variable(index=[time])  # Indirect radiative forcing from aerosols (Wm⁻²).


    function run_timestep(p, v, d, t)

        # Calculate present-day forcing contribution.
        v.pd_re[t] = -p.ϕ * log(1.0 + p.SOx_emiss[t] / p.b_SOx + (p.BC_emiss[t] + p.OC_emiss[t]) / p.b_POM)

        # Calculate pre-industrial forcing contribution.
        v.pi_re[t] = -p.ϕ * log(1.0 + p.SOx_emiss_pi / p.b_SOx + (p.BC_emiss_pi + p.OC_emiss_pi) / p.b_POM)

        # Calculate indirect radiative forcing from aerosols (with potential scaling to account for radiative forcing uncertainty).
        v.rf_aero_indirect[t] = (v.pd_re[t] - v.pi_re[t]) * p.rf_scale_aero_indirect
    end
end
