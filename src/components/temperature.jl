# -----------------------------------------------------------
# Global Surface Temperature Change
# -----------------------------------------------------------

# This follows the forcing to temperature function described in Geoffroy et al. (2013a, 2013b) and implemented in FAIR v1.6.2.

@defcomp temperature begin

    earth_radius        = Parameter()              # Mean radius of the Earth (m).
    seconds_per_year    = Parameter()              # Number of seconds in a year (s yr⁻¹).
    ocean_heat_exchange = Parameter()              # Heat exchange coefficient between the two ocean layers (W m⁻²K⁻¹).
    deep_ocean_efficacy = Parameter()              # Deep ocean efficacy parameter (unitless).
    dt                  = Parameter()              # Length of timestep (years).
    lambda_global       = Parameter()              # Climate feedback parameter (convention is positive = stable) (W m⁻²K⁻¹).
    T_mix₀              = Parameter(index=[2])     # Initial condition for slow and fast contributions to temperature in mixed layer (K).
    T_deep₀             = Parameter(index=[2])     # Initial condition for slow and fast contributions to temperature in deep ocean (K).
    ocean_heat_capacity = Parameter(index=[2])     # Heat capacities for mixed layer and deep ocean (W m⁻²yr⁻¹).
    forcing             = Parameter(index=[time])  # Total effective radiative forcing (W m⁻²).

    T_mix               = Variable(index=[time,2]) # Slow and fast contributions to temperature in mixed layer (K).
    T_deep              = Variable(index=[time,2]) # Slow and fast contributions to temperature in deep ocean (K).
    c_dtemp             = Variable(index=[time])   # TODO: Fill in description with units.
    heatflux            = Variable(index=[time])   # TODO: Fill in description with units.
    del_ohc             = Variable(index=[time])   # TODO: Fill in description with units.
    factor_lambda_eff   = Variable(index=[time])   # TODO: Fill in description with units.
    lambda              = Variable(index=[time])   # TODO: Fill in description with units.
    lambda_eff          = Variable(index=[time])   # TODO: Fill in description with units.
    ratio               = Variable(index=[time])   # TODO: Fill in description with units.
    int_f               = Variable(index=[time])   # TODO: Fill in description with units.
    int_s               = Variable(index=[time])   # TODO: Fill in description with units.
    T                   = Variable(index=[time])   # Global surface temperature anaomaly (K).

    # Intermediate terms that depend on the sampled ocean heat parameters, but otherwise remain constant in all periods.
    ntoa_joule = Variable() # TODO: Fill in description with units.
    cdeep_p    = Variable() # TODO: Fill in description with units.
    gamma_p    = Variable() # TODO: Fill in description with units.
    g1         = Variable() # TODO: Fill in description with units.
    g2         = Variable() # TODO: Fill in description with units.
    g          = Variable() # TODO: Fill in description with units.
    gstar      = Variable() # TODO: Fill in description with units.
    delsqrt    = Variable() # TODO: Fill in description with units.
    afast      = Variable() # TODO: Fill in description with units.
    aslow      = Variable() # TODO: Fill in description with units.
    cc         = Variable() # TODO: Fill in description with units.
    amix_f     = Variable() # TODO: Fill in description with units.
    amix_s     = Variable() # TODO: Fill in description with units.
    adeep_f    = Variable() # TODO: Fill in description with units.
    adeep_s    = Variable() # TODO: Fill in description with units.
    adf        = Variable() # TODO: Fill in description with units.
    ads        = Variable() # TODO: Fill in description with units.
    exp_f      = Variable() # TODO: Fill in description with units.
    exp_s      = Variable() # TODO: Fill in description with units.


    # Set up some terms that depend on the sampled ocean heat parameters, but otherwise remain constant in all periods.
    # TODO: Add comments for each calculation.
    function init(p, v, d)
        v.ntoa_joule = 4 * π * p.earth_radius^2 * p.seconds_per_year
        v.cdeep_p    = p.ocean_heat_capacity[2] * p.deep_ocean_efficacy
        v.gamma_p    = p.ocean_heat_exchange * p.deep_ocean_efficacy
        v.g1         = (p.lambda_global+v.gamma_p)/p.ocean_heat_capacity[1]
        v.g2         = v.gamma_p/v.cdeep_p
        v.g          = v.g1+v.g2
        v.gstar      = v.g1-v.g2
        v.delsqrt    = sqrt(v.g*v.g - 4*v.g2*p.lambda_global/p.ocean_heat_capacity[1])
        v.afast      = (v.g + v.delsqrt)/2
        v.aslow      = (v.g - v.delsqrt)/2
        v.cc         = 0.5/(p.ocean_heat_capacity[1]*v.delsqrt)
        v.amix_f     = v.cc*(v.gstar+v.delsqrt)
        v.amix_s     = -v.cc*(v.gstar-v.delsqrt)
        v.adeep_f    = -v.gamma_p/(p.ocean_heat_capacity[1]*v.cdeep_p*v.delsqrt)
        v.adeep_s    = -v.adeep_f
        v.adf        = 1/(v.afast*p.dt)
        v.ads        = 1/(v.aslow*p.dt)
        v.exp_f      = exp(-1.0 / v.adf)
        v.exp_s      = exp(-1.0 / v.ads)
    end


    function run_timestep(p, v, d, t)

         if is_first(t)

            # Set initial conditions for slow and fast contributions to temperature in mixed layer and deep ocean.
            v.T_mix[t,:]  = p.T_mix₀
            v.T_deep[t,:] = p.T_deep₀
            v.T[t] = sum(v.T_mix[t,:])

        else

            # TODO: Add comments for each calculation.
            v.int_f[t] = (p.forcing[t-1]*v.adf + p.forcing[t]*(1-v.adf) - v.exp_f*(p.forcing[t-1]*(1+v.adf)-p.forcing[t]*v.adf))/v.afast
            v.int_s[t] = (p.forcing[t-1]*v.ads + p.forcing[t]*(1-v.ads) - v.exp_s*(p.forcing[t-1]*(1+v.ads)-p.forcing[t]*v.ads))/v.aslow

            # Calculate slow and fast contributions to temperature in mixed layer.
            v.T_mix[t,1] = v.exp_f * v.T_mix[t-1,1] + v.amix_f*v.int_f[t]
            v.T_mix[t,2] = v.exp_s*v.T_mix[t-1,2] + v.amix_s*v.int_s[t]

            # Calculate slow and fast contributions to temperature in deep ocean.
            v.T_deep[t,1] = v.exp_f*v.T_deep[t-1,1] + v.adeep_f*v.int_f[t]
            v.T_deep[t,2] = v.exp_s*v.T_deep[t-1,2] + v.adeep_s*v.int_s[t]

            # Calculate global surface temperature anomaly.
            v.T[t] = sum(v.T_mix[t,:])

            # TODO: Add calculation comments.
            v.c_dtemp[t] = p.ocean_heat_capacity[1]*(sum(v.T_mix[t,:])-sum(v.T_mix[t-1,:])) + p.ocean_heat_capacity[2]*(sum(v.T_deep[t,:])-sum(v.T_deep[t-1,:]))

            # TODO: Add calculation comments.
            v.heatflux[t] = v.c_dtemp[t]/p.dt

            # TODO: Add calculation comments.
            v.del_ohc[t]  = v.ntoa_joule * v.c_dtemp[t]

            # TODO: Add calculation comments.
            v.factor_lambda_eff[t] = (p.deep_ocean_efficacy-1.0)*p.ocean_heat_exchange

            # TODO: Add calculation comments.
            if abs(sum(v.T_mix[t,:])) > 1e-6
                v.ratio[t] = (sum(v.T_mix[t,:]) - sum(v.T_deep[t,:])) / sum(v.T_mix[t,:])
                v.lambda_eff[t] = p.lambda_global + v.factor_lambda_eff[t]*v.ratio[t]
            else
                v.lambda_eff[t] = p.lambda_global + v.factor_lambda_eff[t]
            end
        end
    end
end
