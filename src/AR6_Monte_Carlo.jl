using JSON

# Create a function that loads the AR6 constrained parameter set and runs a Monte Carlo with FAIRv1.6.2
# Note: Maximum number of parameter samples is 2237
function fair_monte_carlo_ar6(n_samples::Int, ar6_scenario::String = "ssp245")

	# Load FAIR constrained parameter samples from AR6.
	params = JSON.parsefile(joinpath(@__DIR__, "..", "data", "model_data", "fair-1.6.2-wg3-params.json"))

	# Set number of years (just defaulting to 1750-2100 for now).
	n_years = length(1750:2100)

	# Allocate arrays to store results (just temperature and CO₂ concentration for now).
	fair_temperature = zeros(n_years, n_samples)
	fair_co2 = zeros(n_years, n_samples)

	# Calculate index for year 2100 (need to index some values to stop in 2100)
	index_2100 = findfirst(x -> x == 2100, 1750:2100)

	# Create an instnace of MimiFAIRv1.6.2 with IPCC AR6 settings.
	m = get_model(ar6_scenario=ar6_scenario, start_year=1750, end_year=2100)

	for i = 1:n_samples

	   #------------------------------------------------
	   # Update FAIR model parameters with "ith" sample
	   #------------------------------------------------

	   # Carbon cycle
        update_param!(m, :r0_co2, params[i]["r0"])
        update_param!(m, :rT_co2, params[i]["rt"])
        update_param!(m, :rC_co2, params[i]["rc"])

        # Forcing from a doubling of CO₂.
        update_param!(m, :F2x, params[i]["F2x"])

        # Ozone radiative forcing feedback.
        update_param!(m, :feedback, params[i]["ozone_feedback"])

        # "ghan_params" for aerosol indirect forcing effect.
	    update_param!(m, :ϕ,     params[i]["ghan_params"][1])
        update_param!(m, :b_SOx, params[i]["ghan_params"][2])
        update_param!(m, :b_POM, params[i]["ghan_params"][3])

        # Radiative forcing scaling terms (based on ordering of forcing agents in Python code).
	    update_param!(m, :scale_CO₂, params[i]["scale"][1])
        update_param!(m, :scale_CH₄, params[i]["scale"][2])
        update_param!(m, :scale_N₂O, params[i]["scale"][3])
        update_param!(m, :scale_other_ghg, params[i]["scale"][4:15])
        update_param!(m, :scale_ods, params[i]["scale"][16:31])
        update_param!(m, :scale_O₃, params[i]["scale"][32])
        update_param!(m, :scale_CH₄_H₂O, params[i]["scale"][34])
        update_param!(m, :scale_contrails, 0.0) # !!! Default FAIR has contrail forcing switched off. But they sample a scaling term. Not including for now.
        update_param!(m, :scale_aerosol_direct_SOx, params[i]["scale"][36])
        update_param!(m, :scale_aerosol_direct_CO_NMVOC, params[i]["scale"][37])
        update_param!(m, :scale_aerosol_direct_NOx_NH3, params[i]["scale"][38])
        update_param!(m, :scale_aerosol_direct_BC, params[i]["scale"][39])
        update_param!(m, :scale_aerosol_direct_OC, params[i]["scale"][40])
        update_param!(m, :scale_aerosol_indirect, params[i]["scale"][41])
        update_param!(m, :scale_bcsnow, params[i]["scale"][42])
        update_param!(m, :scale_landuse, params[i]["scale"][43])
        update_param!(m, :scale_volcanic, params[i]["scale"][44])
        update_param!(m, :scale_solar, params[i]["scale"][45])

        # Solar radiative forcing.
        update_param!(m, :F_solar, params[i]["F_solar"][1:index_2100])

        # Pre-industrial CO₂ concentration (other concentrations fixed across samples).
        update_param!(m, :CO₂_pi, params[i]["C_pi"][1])

        # Ozone radiative forcing
        update_param!(m, :Ψ_CH₄, params[i]["b_tro3"][1])
        update_param!(m, :Ψ_N₂O, params[i]["b_tro3"][2])
        update_param!(m, :Ψ_ODS, params[i]["b_tro3"][3])
        update_param!(m, :Ψ_CO, params[i]["b_tro3"][4])
        update_param!(m, :Ψ_NMVOC, params[i]["b_tro3"][5])
        update_param!(m, :Ψ_NOx, params[i]["b_tro3"][6])

        # Aerosol direct forcing.
        update_param!(m, :β_SOx, params[i]["b_aero"][1])
        update_param!(m, :β_CO, params[i]["b_aero"][2])
        update_param!(m, :β_NMVOC, params[i]["b_aero"][3])
        update_param!(m, :β_NOx, params[i]["b_aero"][4])
        update_param!(m, :β_BC, params[i]["b_aero"][5])
        update_param!(m, :β_OC, params[i]["b_aero"][6])
        update_param!(m, :β_NH3, params[i]["b_aero"][7])

        # Temperature component
        update_param!(m, :ocean_heat_exchange, params[i]["ocean_heat_exchange"])
        update_param!(m, :deep_ocean_efficacy, params[i]["deep_ocean_efficacy"])
        update_param!(m, :lambda_global, params[i]["lambda_global"])
        update_param!(m, :ocean_heat_capacity, params[i]["ocean_heat_capacity"])

        # Run model
        run(m)

        # Save atmospheric CO₂ concentrations and temperatures resulting from ith parameter sample.
        fair_temperature[:,i] = m[:temperature, :T]
        fair_co2[:,i]         = m[:co2_cycle, :co2]

	end

	# Return temperature and CO₂ projections.
	return fair_temperature, fair_co2
end
