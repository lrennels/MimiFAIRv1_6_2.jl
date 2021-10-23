# -----------------------------------------------------------
# Carbon Cycle
# -----------------------------------------------------------

@defcomp co2_cycle begin

    CO₂_0           = Parameter()   # Carbon dioxide concentration in initial model period (ppm).
    CO₂_pi          = Parameter()   # Pre-industrial carbon dioxide concentrations (ppm).
   

       emiss2conc_co2  = Parameter()   # Conversion factor between emissions (GtC) and concentrations (ppm).
    #GU_co2_0        = Parameter()   # Initial model period value for cumulative uptake of agent (unit of E⁻¹).
    r0_co2          = Parameter()   # Pre-industrial 100-year time-integrated airborne fraction.
   # rA_co2          = Parameter()   # Sensitivity of uptake from atmosphere to current atmospheric burden of agent (unit of E⁻¹).
    rT_co2          = Parameter()   # sensitivity of 100-year time-integrated airborne fraction with temperature anomaly (yr K⁻¹).
    rC_co2          = Parameter()   # sensitivity of 100-year time-integrated airborne fraction with atmospheric carbon stock (yr GtC⁻¹).
    τ_co2           = Parameter(index=[4])  # Atmospheric lifetime of gas in iᵗʰ reservior (years).
    a_co2           = Parameter(index=[4])  # Fraction of emissions entering iᵗʰ carbon reservior.
    R0_co2          = Parameter(index=[4])  # Initial model period value for quantity of agent in iᵗʰ atmospheric reservior (unit of E).
    E_co2           = Parameter(index=[time])   # Annual carbon dioxide emissions (GtC yr⁻¹).
    temperature     = Parameter(index=[time]) # Global surface temperature anaomaly (K).
	cumulative_emissions_CO2₀ = Parameter() # Initial cumulative emissions since pre-industrial (GtC).
	airborne_emissions_CO2₀ = Parameter() # Initial total emissions remaining in the atmosphere (GtC).
    iirf_h          = Parameter() # time horizon for time integrated airborne fraction (yr).
    iIRF_max        = Parameter() # maximum allowed value of iirf (yr)
    dt              = Parameter() # timestep length.

    g0_co2          = Variable()   # Constant to set approximation of value for α equal to Millar et al. (2017) numerical solution for iIRF100 carbon cycle parameterization at α=1.
    g1_co2          = Variable()   # Constant to set approximation of gradient for α equal to Millar et al. (2017) numerical solution for iIRF100 carbon cycle parameterization at α=1.
    α_co2           = Variable(index=[time])    # State-dependent multiplicative adjustment coefficient of reservior lifetimes.
    co2             = Variable(index=[time])    # Total atmospheric carbon dioxide concentrations (ppm).
    GA_co2          = Variable(index=[time])    # Atmospheric burden of agent above pre-industrial levels (unit of E).
    GU_co2          = Variable(index=[time])    # Cumulative uptake of agent since model initialization (unit of E⁻¹).
    iIRFT100_co2    = Variable(index=[time])    # 100-year integrated impulse response function (the average airborne fraction over a 100-year period).
    decay_rate_co2  = Variable(index=[time,4])  # Decay rates
    R_co2           = Variable(index=[time,4])  # Quantity of agent in iᵗʰ atmospheric reservior (unit of E).
	cumulative_emissions_CO2 = Variable(index=[time]) # Cumulative emissions since pre-industrial (GtC).
	airborne_emissions_CO2 = Variable(index=[time]) # Total emissions remaining in the atmosphere (GtC).


    function run_timestep(p, v, d, t)

    	if is_first(t)

            # Set initial values.
           # v.GU_co2[t]  = p.GU_co2_0
            v.R_co2[t,:] = p.R0_co2
            v.co2[t]     = p.CO₂_0
			v.cumulative_emissions_CO2[t] = p.cumulative_emissions_CO2₀
			v.airborne_emissions_CO2[t] = p.airborne_emissions_CO2₀

            # Calculate initial burden above pre-industrial values.
            v.GA_co2[t] = sum(v.R_co2[t,:])

            #Initialise simplified carbon cycle parameters
            v.g1_co2 = sum(p.a_co2 .* p.τ_co2 .* (1.0 .- (1.0 .+ p.iirf_h ./ p.τ_co2) .* exp.(-p.iirf_h ./ p.τ_co2)))
            v.g0_co2 = 1.0 / (sinh(sum(p.a_co2 .* p.τ_co2 .* (1.0 .- exp.(-p.iirf_h ./ p.τ_co2))) / v.g1_co2))

        else

            # Calculate iIRF100.
            v.iIRFT100_co2[t] = min(p.r0_co2 + p.rC_co2 * (v.cumulative_emissions_CO2[t-1] - v.airborne_emissions_CO2[t-1]) + p.rT_co2 * p.temperature[t-1], p.iIRF_max)

            # Add a check for negative iIRF100 values (assume they remain fixed at previous value if going negative).
            if v.iIRFT100_co2[t] <= 0.0
                v.iIRFT100_co2[t] =  v.iIRFT100_co2[t-1]
            end

            # Calculate state-dependent lifetime adjustment term based on iIRF100 value.
            v.α_co2[t] = v.g0_co2 * sinh(v.iIRFT100_co2[t] / v.g1_co2)

            # Calculate concentrations in each carbon reservoir.
            for i = 1:4

                # Calculate decay rate for iᵗʰ reservoir.
                v.decay_rate_co2[t,i] = p.dt / (v.α_co2[t] * p.τ_co2[i])

                # Calculate amount of carbon in iᵗʰ reservoir.
#                v.R_co2[t,i] = p.E_co2[t] * p.a_co2[i] / v.decay_rate_co2[t,i] * (1.0 - exp(-v.decay_rate_co2[t,i])) + v.R_co2[t-1,i] * exp(-v.decay_rate_co2[t,i])
                v.R_co2[t,i] = p.E_co2[t-1] / p.emiss2conc_co2 * p.a_co2[i] * v.α_co2[t] * (p.τ_co2[i]/p.dt) * (1.0 - exp(-v.decay_rate_co2[t,i])) + v.R_co2[t-1,i] * exp(-v.decay_rate_co2[t,i])

            end

            # Calcualte atmospheric burden above pre-industiral levels.
            v.GA_co2[t] = sum(v.R_co2[t,:])

  			# Calculate atmospheric CO₂ concentration.
            #v.co2[t] = p.co2_pi + p.emiss2conc_co2 * (v.GA_co2[t-1] + v.GA_co2[t]) / 2.0
            v.co2[t] = p.CO₂_pi + (v.GA_co2[t-1] + v.GA_co2[t]) / 2.0

            # Calculate total emissions remaining in the atmosphere
            v.airborne_emissions_CO2[t] = v.GA_co2[t] * p.emiss2conc_co2

            # Calculate cumulative emissions.
            v.cumulative_emissions_CO2[t] = v.cumulative_emissions_CO2[t-1] + p.E_co2[t]
        end
    end
end
