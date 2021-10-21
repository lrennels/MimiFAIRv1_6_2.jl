using Distributions, Dates, Mimi, CSVFiles, DataFrames
import Mimi: SampleStore, add_RV!, add_transform!

"""
    run_mcs(;trials::Int64 = 2237, 
        output_dir::Union{String, Nothing} = nothing, 
        save_trials::Bool = false,
        m::Mimi.Model = get_model())

Return the results of a Monte Carlo Simulation with the defined number of trials
and save data into the `output_dir` folder, optionally also saving trials if 
`save_trials` is set to `true.` If no model is provided, use the default model 
returned by get_model().
"""
function run_mcs(;trials::Int64 = 2237, 
                    output_dir::Union{String, Nothing} = nothing, 
                    save_trials::Bool = false,
                    m::Mimi.Model = get_model())

    trials < 2 && error("Must run `run_mcs` function with a `trials` argument greater than 1 due to a Mimi specification about SampleStores.  TO BE FIXED SOON!")
    
    # Set up output directories
    output_dir = output_dir === nothing ? joinpath(@__DIR__, "../../output/mcs/", "MCS $(Dates.format(now(), "yyyy-mm-dd HH-MM-SS")) MC$trials") : output_dir
    isdir("$output_dir/results") || mkpath("$output_dir/results")

    trials_output_filename = save_trials ?  joinpath("$output_dir/trials.csv") : nothing

    # Get an instance of the mcs
    mcs = get_mcs()

    # run monte carlo trials
    results = run(mcs, m, trials; trials_output_filename = trials_output_filename, results_output_dir = "$output_dir/results")

    return results
end

"""
    get_mcs()

Return a Monte Carlo Simulation definition of type Mimi.SimulationDefinition that
holds all random variables and distributions, as assigned to model component/parameter
pairs, that will be used in a Monte Carlo Simulation.
"""
function get_mcs()
    
    # define the Monte Carlo Simulation
    mcs = @defsim begin
        save(temperature.T, co2_cycle.co2)
    end

    # add the FAIR random variables and transforms - note this could be done within
    # the @defsim macro but we can use the dictionary to make this less verbose

    # Note that if a parameter component is not included in add_transform!, the 
    # parameters are shared model parameters, and each line will create ONE random 
    # variable and assign all component parameters connected to that shared model 
    # parameter to the value taken on by that random variable

    fair_samples_map, fair_samples = get_fair_mcs_params()
    fair_samples_left = deepcopy(fair_samples) # we know we've added everything when this is empty!

    # add and assign all random variables for single dimensional parameters
    for (k,v) in fair_samples_left
        if size(v, 2) == 1 # one column of values
            rv_name = Symbol("rv_$k")
            add_RV!(mcs, rv_name, SampleStore(fair_samples[k][!, 1]))
            add_transform!(mcs, k, :(=), rv_name)
            delete!(fair_samples_left, k)
        end
    end

    # assign one random variable per year with a unique distribution from fair_samples
    # TODO handle dimensions - what years does F_solar include?  Assuming 1750 onwards for 361 years 

    for year in 1750:2100
        rv_name = Symbol("rv_F_solar_$year")
        add_RV!(mcs, rv_name, SampleStore(fair_samples[:F_solar][!,string(year)]))
        add_transform!(mcs, :F_solar, :(=), rv_name, [year])
    end
    delete!(fair_samples_left, :F_solar)

    # Radiative forcing scaling - one distribution per "other" greenhouse gas, and
    # one per ods

    for gas in names(fair_samples[:scale_other_ghg])
        rv_name = Symbol("rv_scale_other_ghg_$(gas)")
        add_RV!(mcs, rv_name, SampleStore(fair_samples[:scale_other_ghg][!, gas]))
        add_transform!(mcs, :scale_other_ghg, :(=), rv_name, [gas])
    end
    delete!(fair_samples_left, :scale_other_ghg)

    for ods in names(fair_samples[:scale_ods])
        rv_name = Symbol("rv_scale_ods_$(ods)")
        add_RV!(mcs, rv_name, SampleStore(fair_samples[:scale_ods][!, ods]))
        add_transform!(mcs, :scale_ods, :(=), rv_name, [ods])
    end
    delete!(fair_samples_left, :scale_ods)

    # ocean_heat_capacity takes an anonymous dim of 2 (deep and mixed, should label 
    # explicilty) - anonymouse dims are named with Ints 1 and 2

    rv_name = Symbol("rv_ocean_heat_capacity_1")
    add_RV!(mcs, rv_name, SampleStore(fair_samples[:ocean_heat_capacity][!, "1"]))
    add_transform!(mcs, :ocean_heat_capacity, :(=), rv_name, [1])

    rv_name = Symbol("rv_ocean_heat_capacity_2")
    add_RV!(mcs, rv_name, SampleStore(fair_samples[:ocean_heat_capacity][!, "2"]))
    add_transform!(mcs, :ocean_heat_capacity, :(=), rv_name, [2])

    delete!(fair_samples_left, :ocean_heat_capacity)

    # check if we've added all FAIR parameters
    isempty(fair_samples_left) ? nothing : error("The following FAIR mcs uncertain parameters has not been added to the simulation: $(keys(fair_samples_left))")
    
    return mcs
end

"""
    get_fair_mcs_params()

Return the FAIR mcs parameters mapped from parameter name to string name, and a dictionary
using the parameter names as keys and a DataFrame holding the values as a value.
"""
function get_fair_mcs_params()

    names_map = get_fair_mcs_params_map()
    params_dict = Dict()
    for (k,v) in names_map
        push!(params_dict, k => load(joinpath(@__DIR__, "..", "..", "data", "mcs_params", "fair_mcs_params_$v.csv")) |> DataFrame)
    end
    return names_map, params_dict

end

"""
    get_fair_mcs_params_map()

Return a dictionary of FAIR elements with the FAIR v1.6.2 parameter name being 
the component, parameter pair and the value being the parameter csv name. 
"""
function get_fair_mcs_params_map()
    return Dict(
        :β_CO                          => "b_aero_CO",
        :scale_CH₄                     => "scale_CH4",
        :F_solar                       => "F_solar",
        :Ψ_CH₄                         => "b_tro3_CH4",
        :scale_N₂O                     => "scale_N2O",
        :CO₂_pi                        => "C_pi",
        :deep_ocean_efficacy           => "deep_ocean_efficacy",
        :scale_bcsnow                  => "scale_bcsnow",
        :scale_aerosol_direct_OC       => "scale_aerosol_direct_OC",
        :b_SOx                         => "ghan_params_SOx",
        :feedback                      => "ozone_feedback",
        :scale_O₃                      => "scale_O3",
        :b_POM                         => "ghan_params_b_POM",
        :r0_co2                        => "r0",
        :β_NH3                         => "b_aero_NH3",
        :lambda_global                 => "lambda_global",
        :scale_landuse                 => "scale_landuse",
        :scale_volcanic                => "scale_volcanic",
        :scale_aerosol_direct_SOx      => "scale_aerosol_direct_SOx",
        :β_NOx                         => "b_aero_NOx",
        :Ψ_N₂O                         => "b_tro3_N2O",
        :ocean_heat_capacity           => "ocean_heat_capacity",
        :β_OC                          => "b_aero_OC",
        :scale_solar                   => "scale_solar",
        :rC_co2                        => "rc",
        :scale_aerosol_direct_BC       => "scale_aerosol_direct_BC",
        :scale_CH₄_H₂O                 => "scale_CH4_H2O",
        :scale_aerosol_indirect        => "scale_aerosol_indirect",
        :scale_ods                     => "scale_ods",
        :Ψ_CO                          => "b_tro3_CO",
        :scale_aerosol_direct_NOx_NH3  => "scale_aerosol_direct_NOx_NH3",
        :scale_other_ghg               => "scale_other_ghg",
        :Ψ_NMVOC                       => "b_tro3_NMVOC",
        :F2x                           => "F2x",
        :β_SOx                         => "b_aero_SOx",
        :β_NMVOC                       => "b_aero_NMVOC",
        :rT_co2                        => "rt",
        :β_BC                          => "b_aero_BC",
        :scale_CO₂                     => "scale_CO2",
        :Ψ_ODS                         => "b_tro3_ODS",
        :scale_aerosol_direct_CO_NMVOC => "scale_aerosol_direct_CO_NMVOC",
        :Ψ_NOx                         => "b_tro3_NOx",
        :ocean_heat_exchange           => "ocean_heat_exchange",
        :ϕ                             => "ghan_params_Pi"
    )
end
