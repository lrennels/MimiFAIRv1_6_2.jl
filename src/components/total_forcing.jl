#--------------------------
# Total raditiave forcing
#--------------------------

@defcomp total_forcing begin

    other_ghg = Index()                                            # Index for other well-mixed greenhouse gases.
    ozone_depleting_substances = Index() # Index for ozone-depleting substances.

    scale_CO₂              = Parameter()                        # Carbon dioxide scaling factor to capture effective radiative forcing uncertainty.
    scale_CH₄              = Parameter()                        # Methane scaling factor to capture effective radiative forcing uncertainty.
    scale_CH₄_H₂O          = Parameter()                        # Stratospheric water vapor from methane oxidation scaling factor to capture effective radiative forcing uncertainty.
    scale_N₂O              = Parameter()                        # Nitrous oxide scaling factor to capture effective radiative forcing uncertainty.
    scale_O₃               = Parameter()                        # Ozone scaling factor to capture effective radiative forcing uncertainty.
    #scale_aerosol_direct   = Parameter()                        # Direct aerosol effect scaling factor to capture effective radiative forcing uncertainty.
    scale_aerosol_indirect = Parameter()                        # Indirect aerosol effect scaling factor to capture effective radiative forcing uncertainty.
    scale_bcsnow           = Parameter()                        # Black carbon on snow scaling factor to capture effective radiative forcing uncertainty.
    scale_landuse          = Parameter()                        # Land-use surface albedo change scaling factor to capture effective radiative forcing uncertainty.
    scale_contrails        = Parameter()                        # Contrails scaling factor to capture effective radiative forcing uncertainty.
    scale_volcanic        = Parameter()                        # Contrails scaling factor to capture effective radiative forcing uncertainty.
    scale_solar        = Parameter()                        # Contrails scaling factor to capture effective radiative forcing uncertainty.

    scale_aerosol_direct_SOx   = Parameter() # Direct aerosol effect (SOx contribution) scaling factor to capture effective radiative forcing uncertainty.
    scale_aerosol_direct_CO_NMVOC    = Parameter() # Direct aerosol effect (CO & NMVOC contribution) scaling factor to capture effective radiative forcing uncertainty.
    #scale_aerosol_direct_NMVOC = Parameter() # Direct aerosol effect (NMVOC contribution) scaling factor to capture effective radiative forcing uncertainty.
    scale_aerosol_direct_NOx_NH3   = Parameter() # Direct aerosol effect (NOx & NH3 contribution) scaling factor to capture effective radiative forcing uncertainty.
    scale_aerosol_direct_BC    = Parameter() # Direct aerosol effect (BC contribution) scaling factor to capture effective radiative forcing uncertainty.
    scale_aerosol_direct_OC    = Parameter() # Direct aerosol effect (OC contribution) scaling factor to capture effective radiative forcing uncertainty.
    #scale_aerosol_direct_NH3   = Parameter() # Direct aerosol effect (NH3 contribution) scaling factor to capture effective radiative forcing uncertainty.
    scale_other_ghg        = Parameter(index=[other_ghg])       # Scaling factor to capture effective radiative forcing uncertainty for other well-mixed greenhouse gases.
    scale_ods              = Parameter(index=[ozone_depleting_substances])       # Scaling factor to capture effective radiative forcing uncertainty for ozone-depeleting substances.



    F_CO₂                     = Parameter(index=[time])            # Radiative forcing from carbon dioxide (Wm⁻²).
    F_CH₄                     = Parameter(index=[time])            # Radiative forcing from methane (Wm⁻²).
    F_CH₄_H₂O                 = Parameter(index=[time])            # Radiative forcing from stratospheric water varpor from methane oxidiation (Wm⁻²).
    F_N₂O                     = Parameter(index=[time])            # Radiative forcing from nitrous oxide (Wm⁻²).
    F_O₃                      = Parameter(index=[time])            # Radiative forcing from ozone (Wm⁻²).
    #F_aerosol_direct          = Parameter(index=[time])            # Radiative forcing from direct aerosol effect (Wm⁻²).
    
        F_aerosol_direct_SOx   = Parameter(index=[time]) # Radiative forcing from direct aerosol effect (SOx contribution) (W m⁻²).
        F_aerosol_direct_CO    = Parameter(index=[time]) # Radiative forcing from direct aerosol effect (CO contribution) (W m⁻²).
        F_aerosol_direct_NMVOC = Parameter(index=[time]) # Radiative forcing from direct aerosol effect (NMVOC contribution) (W m⁻²).
        F_aerosol_direct_NOx   = Parameter(index=[time]) # Radiative forcing from direct aerosol effect (NOx contribution) (W m⁻²).
        F_aerosol_direct_BC    = Parameter(index=[time]) # Radiative forcing from direct aerosol effect (BC contribution) (W m⁻²).
        F_aerosol_direct_OC    = Parameter(index=[time]) # Radiative forcing from direct aerosol effect (OC contribution) (W m⁻²).
        F_aerosol_direct_NH3   = Parameter(index=[time]) # Radiative forcing from direct aerosol effect (NH3 contribution) (W m⁻²).

    F_aerosol_indirect        = Parameter(index=[time])            # Radiative forcing from indirect aerosol effect (Wm⁻²).
    F_bcsnow                  = Parameter(index=[time])            # Radiative forcing from black carbon on snow (Wm⁻²).
    F_landuse                 = Parameter(index=[time])            # Radiative forcing from land-use surface albedo change (Wm⁻²).
    F_contrails               = Parameter(index=[time])            # Radiative forcing from contrails (Wm⁻²).
    F_other_ghg               = Parameter(index=[time, other_ghg]) # Radiative forcing from other well-mixed greenhouse gases (Wm⁻²).
    F_ods                     = Parameter(index=[time, ozone_depleting_substances]) # Radiative forcing from ozone_depleting substances (Wm⁻²).

    F_volcanic                = Parameter(index=[time])            # Radiative forcing from volcanic eruptions (Wm⁻²).
    F_solar                   = Parameter(index=[time])            # Radiative forcing from solar irradience (Wm⁻²).
    F_exogenous               = Parameter(index=[time])            # Radiative forcing from exogenous sources not included in this module (Wm⁻²).

    total_forcing             = Variable(index = [time])           # Total radiative forcing, with individual components scaled by their respective efficacy (Wm⁻²).


    function run_timestep(p, v, d, t)

        # Calculate total radiative forcing as the sum of individual radiative forcings (with potential scaling factor).
        v.total_forcing[t] = p.F_CO₂[t]                 * p.scale_CO₂ +
                             p.F_CH₄[t]                 * p.scale_CH₄ +
                             p.F_CH₄_H₂O[t]             * p.scale_CH₄_H₂O +
                             p.F_N₂O[t]                 * p.scale_N₂O +
                             sum(p.F_other_ghg[t,:]    .* p.scale_other_ghg[:]) +
                             sum(p.F_ods[t,:]          .* p.scale_ods[:]) +
                             p.F_O₃[t]                  * p.scale_O₃ +
                             p.F_aerosol_direct_SOx[t]  * p.scale_aerosol_direct_SOx +
                             (p.F_aerosol_direct_CO[t]  + p.F_aerosol_direct_NMVOC[t]) * p.scale_aerosol_direct_CO_NMVOC +
                             (p.F_aerosol_direct_NOx[t] + p.F_aerosol_direct_NH3[t]) * p.scale_aerosol_direct_NOx_NH3 +
                             p.F_aerosol_direct_BC[t]   * p.scale_aerosol_direct_BC +
                             p.F_aerosol_direct_OC[t]   * p.scale_aerosol_direct_OC +
                             p.F_aerosol_indirect[t]    * p.scale_aerosol_indirect +
                             p.F_bcsnow[t]              * p.scale_bcsnow +
                             p.F_landuse[t]             * p.scale_landuse +
                             p.F_contrails[t]           * p.scale_contrails +
                             p.F_solar[t]               * p.scale_solar +
                             p.F_volcanic[t]            * p.scale_volcanic +
                             p.F_exogenous[t]
    end
end
