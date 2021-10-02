


# Load required packages.
using CSVFiles, DataFrames, Mimi

# Load helper functions and MimiFAIRv1.6.2 model commponent files.
#include("helper_functions.jl")
include(joinpath("src/components", "temperature_component.jl"))


# Set start and end years.
start_year = 1765
end_year   = 2500

# Load exogenous radiative forcing data from Python run.
python_forcing = DataFrame(load(joinpath(@__DIR__, "data", "rcp85_forcings.csv"), header_exists=false))
python_forcing_total = sum(Matrix(python_forcing), dims=2)[:]

# Load exogenous temperature data from Python run
python_temperature = Matrix(DataFrame(load(joinpath(@__DIR__, "data", "rcp85_temperature.csv"), header_exists=false)))[:]

# Load exogenous CO2 emissions data from Python run.
python_co2_emissions = DataFrame(load(joinpath(@__DIR__, "data", "rcp85_co2Emissions.csv"), header_exists=false))
python_co2_emissions_total = sum(Matrix(python_co2_emissions), dims=2)[:]


    # ---------------------------------------------
    # ---------------------------------------------
    # Initialize Mimi model.
    # ---------------------------------------------
    # ---------------------------------------------
include(joinpath("src/components", "temperature.jl"))
include(joinpath("src/components", "co2_cycle.jl"))

    #TODO = Set up function to load this default data for FAIR (pre-industrial values, etc.)
    co2_pi           = 278.0 # -> 278.05158 is the exact term they have in pre-industiral Python data, but they set CO2 to 278 by default in actual model. 278.05158  # ppm
    m_atmos          = 5.1352e18 # mass of atmosphere, kg
    earth_radius     = 6371000   # m
    seconds_per_year = 60 * 60 * 24 * 365.24219 # Length of tropical year
    molwt_C           = 12.01
    molwt_AIR         = 28.97
    emiss2conc_co2   = m_atmos/1e18*molwt_C/molwt_AIR # Conversion between ppm CO2 and GtC emissions


   # Create a Mimi model.
    m = Model()

    # Set time and gas-grouping indices.
    set_dimension!(m, :time, start_year:end_year)

    # ---------------------------------------------
    # Add components to model
    # ---------------------------------------------
    add_comp!(m, co2_cycle)
    add_comp!(m, temperature)

    # ---------------------------------------------
    # Set component-specific parameters
    # ---------------------------------------------

    # ---- Carbon Cycle ---- #
    set_param!(m, :co2_cycle, :co2_0,  278.001409189439) # From FAIR model run.
    set_param!(m, :co2_cycle, :co2_pi,  co2_pi)
    #set_param!(m, :co2_cycle, :g0_co2,  )
    #set_param!(m, :co2_cycle, :g1_co2,  )
    set_param!(m, :co2_cycle, :iirf_h,  100.0)
    set_param!(m, :co2_cycle, :emiss2conc_co2, emiss2conc_co2)
    set_param!(m, :co2_cycle, :r0_co2, 35.0)
    set_param!(m, :co2_cycle, :rT_co2, 4.165)
    set_param!(m, :co2_cycle, :rC_co2,  0.019)
    set_param!(m, :co2_cycle, :τ_co2, [1000000, 394.4, 36.54, 4.304])
    set_param!(m, :co2_cycle, :a_co2, [0.2173,0.2240,0.2824,0.2763])
    set_param!(m, :co2_cycle, :R0_co2, [0.0003062168651584551, 0.0003156584344017209, 0.0003979550976564552, 0.0003893590420767655]) # From FAIR model run.
    set_param!(m, :co2_cycle, :E_co2, python_co2_emissions_total)
    #set_param!(m, :co2_cycle, :temperature, python_temperature)
    set_param!(m, :co2_cycle, :cumulative_emissions_CO2₀, 0.003)
    set_param!(m, :co2_cycle, :airborne_emissions_CO2₀, 0.0)
    set_param!(m, :co2_cycle, :iIRF_max, 97.0)
    #set_param!(m, :co2_cycle, :dt, 1.0)

    connect_param!(m, :co2_cycle => :temperature, :temperature => :T)

    # ---- Temperature ---- #
    set_param!(m, :temperature, :earth_radius, earth_radius)
    set_param!(m, :temperature, :seconds_per_year, seconds_per_year)
    set_param!(m, :temperature, :ocean_heat_exchange, 0.67)
    set_param!(m, :temperature, :deep_ocean_efficacy, 1.28)
    #set_param!(m, :temperature, :dt, 1.0)
    set_param!(m, :temperature, :lambda_global, 1.18)
    set_param!(m, :temperature, :T_mix₀, [5.30304445e-03, 6.41295867e-05]) # From Python-FAIR model run.
    set_param!(m, :temperature, :T_deep₀, [-0.00013307,  0.00015021]) # From Python-FAIR model run.
    set_param!(m, :temperature, :ocean_heat_capacity, [8.2, 109.0])
    set_param!(m, :temperature, :forcing, python_forcing_total)

    # ---- Parameters Shared Across Multiple Components ---- #
    set_param!(m, :dt, 1.0)


    run(m)

