


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

# Load exogenous concentration data from Python run
python_concentrations = Matrix(DataFrame(load(joinpath(@__DIR__, "data", "rcp85_concentrations.csv"), header_exists=false)))

    # ---------------------------------------------
    # ---------------------------------------------
    # Initialize Mimi model.
    # ---------------------------------------------
    # ---------------------------------------------
include(joinpath("src/components", "ch4_cycle.jl"))
include(joinpath("src/components", "n2o_cycle.jl"))
include(joinpath("src/components", "co2_cycle.jl"))
include(joinpath("src/components", "other_ghg_cycles.jl"))
include(joinpath("src/components", "o3_depleting_cycles.jl"))
include(joinpath("src/components", "o3_forcing.jl"))
include(joinpath("src/components", "aerosol_indirect_forcing.jl"))
include(joinpath("src/components", "temperature.jl"))

    #TODO = Set up function to load this default data for FAIR (pre-industrial values, etc.)
    co2_pi           = 278.0 # -> 278.05158 is the exact term they have in pre-industiral Python data, but they set CO2 to 278 by default in actual model. 278.05158  # ppm
    m_atmos          = 5.1352e18 # mass of atmosphere, kg
    earth_radius     = 6371000   # m
    seconds_per_year = 60 * 60 * 24 * 365.24219 # Length of tropical year
    molwt_C           = 12.01
    molwt_AIR         = 28.97
    emiss2conc_co2   = m_atmos/1e18*molwt_C/molwt_AIR # Conversion between ppm CO2 and GtC emissions


    # Load RCP and other data needed to construct FAIR.
    rcp_emissions, volcano_forcing, solar_forcing, gas_data, gas_fractions, conversions = load_fair_data(start_year, end_year, rcp_scenario)

    # Names of minor greenhouse gases and ozone-depleting substances.
    #This is the original list fro mFAIR1.3, not split up: other_ghg_names = ["CF4", "C2F6", "C6F14", "HFC23", "HFC32", "HFC43_10", "HFC125", "HFC134a", "HFC143a", "HFC227ea", "HFC245fa", "SF6", "CFC_11", "CFC_12", "CFC_113", "CFC_114", "CFC_115", "CARB_TET", "MCF", "HCFC_22", "HCFC_141B", "HCFC_142B", "HALON1211", "HALON1202", "HALON1301", "HALON2402", "CH3BR", "CH3CL"]
    other_ghg_names = ["CF4", "C2F6", "C6F14", "HFC23", "HFC32", "HFC43_10", "HFC125", "HFC134a", "HFC143a", "HFC227ea", "HFC245fa", "SF6"]

    ods_names = ["CFC_11", "CFC_12", "CFC_113", "CFC_114", "CFC_115", "CARB_TET", "MCF", "HCFC_22", "HCFC_141B", "HCFC_142B", "HALON1211", "HALON1202", "HALON1301", "HALON2402", "CH3BR", "CH3CL"]

   # Create a Mimi model.
   m = Model()

    # Set time and gas-grouping indices.
    set_dimension!(m, :time, start_year:end_year)
    set_dimension!(m, :other_ghg, other_ghg_names)
    set_dimension!(m, :ozone_depleting_substances, ods_names)

    # ---------------------------------------------
    # Add components to model
    # ---------------------------------------------
    add_comp!(m, ch4_cycle)
    add_comp!(m, n2o_cycle)
    add_comp!(m, co2_cycle)
    add_comp!(m, other_ghg_cycles)
    add_comp!(m, o3_depleting_cycles)
    add_comp!(m, o3_forcing)
    add_comp!(m, aerosol_indirect_forcing)
    add_comp!(m, temperature)

    # ---------------------------------------------
    # Set component-specific parameters
    # ---------------------------------------------

    # ---- Methane Cycle ---- #
    set_param!(m, :ch4_cycle, :fossil_emiss_CH₄, rcp_emissions.CH4)
    set_param!(m, :ch4_cycle, :natural_emiss_CH₄, rcp_emissions.NaturalCH4)
    set_param!(m, :ch4_cycle, :τ_CH₄, 9.3)
    set_param!(m, :ch4_cycle, :fossil_frac, gas_fractions.ch4_fossil)
    set_param!(m, :ch4_cycle, :oxidation_frac, 0.61)
    set_param!(m, :ch4_cycle, :mol_weight_CH₄, gas_data[gas_data.gas .== "CH4", :mol_weight][1])
    set_param!(m, :ch4_cycle, :mol_weight_C, gas_data[gas_data.gas .== "C", :mol_weight][1])
    set_param!(m, :ch4_cycle, :emiss2conc_ch4, conversions[conversions.gases .== "CH4", :emiss2conc][1])
    set_param!(m, :ch4_cycle, :CH₄_0, gas_data[gas_data.gas .== "CH4", :pi_conc][1])

    # ---- Nitrous Oxide Cycle ---- #
    set_param!(m, :n2o_cycle, :fossil_emiss_N₂O, rcp_emissions.N2O)
    set_param!(m, :n2o_cycle, :natural_emiss_N₂O, rcp_emissions.NaturalN2O)
    set_param!(m, :n2o_cycle, :τ_N₂O, 121.0)
    set_param!(m, :n2o_cycle, :emiss2conc_n2o, conversions[conversions.gases .== "N2O", :emiss2conc][1])
    set_param!(m, :n2o_cycle, :N₂O_0, gas_data[gas_data.gas .== "N2O", :pi_conc][1])


    # ---- Carbon Cycle ---- #
    set_param!(m, :co2_cycle, :CO₂_0,  278.001409189439) # From FAIR model run.
    #set_param!(m, :co2_cycle, :co2_pi,  co2_pi)
    #set_param!(m, :co2_cycle, :g0_co2,  )
    #set_param!(m, :co2_cycle, :g1_co2,  )
    set_param!(m, :co2_cycle, :iirf_h,  100.0)
    #set_param!(m, :co2_cycle, :emiss2conc_co2, emiss2conc_co2)
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

    # ---- Other Well-Mixed Greenhouse Gas Cycles ---- #
    set_param!(m, :other_ghg_cycles, :τ_other_ghg, gas_data[findall((in)(other_ghg_names), gas_data.gas), :lifetimes])
    set_param!(m, :other_ghg_cycles, :emiss_other_ghg, Matrix(rcp_emissions[!,Symbol.(other_ghg_names)]))
    set_param!(m, :other_ghg_cycles, :emiss2conc_other_ghg, conversions[findall((in)(other_ghg_names), conversions.gases), :emiss2conc])

    # ---- Ozone-Depleting Substance Gas Cycles ---- #
    set_param!(m, :o3_depleting_cycles, :τ_ods, gas_data[findall((in)(ods_names), gas_data.gas), :lifetimes])
    set_param!(m, :o3_depleting_cycles, :emiss_ods, Matrix(rcp_emissions[!,Symbol.(ods_names)]))
    set_param!(m, :o3_depleting_cycles, :emiss2conc_ods, conversions[findall((in)(ods_names), conversions.gases), :emiss2conc])

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



    # ---- Ozone Radiative Forcing ---- #
    set_param!(m, :o3_forcing, :Br, gas_data[findall((in)(ods_names), gas_data.gas), :br_atoms])
    set_param!(m, :o3_forcing, :Cl, gas_data[findall((in)(ods_names), gas_data.gas), :cl_atoms])
    set_param!(m, :o3_forcing, :FC, gas_data[findall((in)(ods_names), gas_data.gas), :strat_frac])
    set_param!(m, :o3_forcing, :ODS_pi, gas_data[findall((in)(ods_names), gas_data.gas), :pi_conc])

    set_param!(m, :o3_forcing, :total_forcing_O₃_0, 0.0)

    set_param!(m, :o3_forcing, :NOx_emissions_pi, 0.0)
    set_param!(m, :o3_forcing, :CO_emissions_pi, 0.0)
    set_param!(m, :o3_forcing, :NMVOC_emissions_pi, 0.0)

    set_param!(m, :o3_forcing, :NOx_emissions, rcp_emissions.NOx)
    set_param!(m, :o3_forcing, :CO_emissions, rcp_emissions.CO)
    set_param!(m, :o3_forcing, :NMVOC_emissions, rcp_emissions.NMVOC) 

    set_param!(m, :o3_forcing, :feedback, -0.037)

    set_param!(m, :o3_forcing, :β_CH₄, 2.33379720e-04)
    set_param!(m, :o3_forcing, :β_N₂O, 1.27179106e-03)
    set_param!(m, :o3_forcing, :β_ODS, -6.69347820e-05)
    set_param!(m, :o3_forcing, :β_CO, 1.14647701e-04)
    set_param!(m, :o3_forcing, :β_NMVOC, 5.14366051e-12)
    set_param!(m, :o3_forcing, :β_NOx, 3.78354423e-03)

#connect_param!(m, :o3_forcing => :conc_ODS, :o3_depleting_cycles => :conc_ods)

    #-----------------------------------------------------------
    # THESE NEED TO BE CONNECTIONS BUT READING IN DATA FOR NOW
    #-----------------------------------------------------------
    set_param!(m, :o3_forcing,  :CH₄, python_concentrations[:,2])
    set_param!(m, :o3_forcing,  :N₂O, python_concentrations[:,3])
    #set_param!(m, :o3_forcing,  :CO₂, python_concentrations[:,1])
    set_param!(m, :o3_forcing,  :temperature, python_temperature)
    set_param!(m, :o3_forcing,  :conc_ODS, python_concentrations[:,16:31])



    # ---- Aerosol Indirect Radiative Forcing ---- #
    set_param!(m, :aerosol_indirect_rf, :ϕ, -1.95011431)
    set_param!(m, :aerosol_indirect_rf, :b_SOx, 0.01107147)
    set_param!(m, :aerosol_indirect_rf, :b_POM, 0.01387492)
    #set_param!(m, :aerosol_indirect_rf, :rf_scale_aero_indirect, 0.0)
    #set_param!(m, :aerosol_indirect_rf, :model_years, collect(start_year:end_year))
    set_param!(m, :aerosol_indirect_rf, :SOx_emiss_1765, 1.0)
    set_param!(m, :aerosol_indirect_rf, :BC_OC_emiss_1765, 11.2)
    set_param!(m, :aerosol_indirect_rf, :scale_AR5, true)
    set_param!(m, :aerosol_indirect_rf, :F_1765, -0.3002836449793625)
    set_param!(m, :aerosol_indirect_rf, :F_2011, -1.5236182344467388)


    # ---- Parameters Shared Across Multiple Components ---- #
    set_param!(m, :dt, 1.0)
    set_param!(m, :CH₄_pi, gas_data[gas_data.gas .== "CH4", :pi_conc][1])
    set_param!(m, :emiss2conc_co2, conversions[conversions.gases .== "CO2", :emiss2conc][1])
    set_param!(m, :CO₂_pi, gas_data[gas_data.gas .== "CO2", :pi_conc][1])
    set_param!(m, :N₂O_pi, gas_data[gas_data.gas .== "N2O", :pi_conc][1])
    set_param!(m, :SOx_emiss, rcp_emissions.SOx)
    set_param!(m, :BC_emiss, rcp_emissions.BC)
    set_param!(m, :OC_emiss, rcp_emissions.OC)
    #set_param!(m, :N₂O_0, gas_data[gas_data.gas .== "N2O", :pi_conc][1])
    set_param!(m, :other_ghg_0, gas_data[findall((in)(other_ghg_names), gas_data.gas), :pi_conc])
    set_param!(m, :ods_0, gas_data[findall((in)(ods_names), gas_data.gas), :pi_conc])
    set_param!(m, :fix_pre1850_RCP, false)



#    set_param!(m, :SOx_emiss, rcp_emissions.SOx)
#    set_param!(m, :BC_emiss, rcp_emissions.BC)
#    set_param!(m, :OC_emiss, rcp_emissions.OC)
#    set_param!(m, :NOx_emiss, rcp_emissions.NOx)
#    set_param!(m, :fix_pre1850_RCP, true)
#    set_param!(m, :mol_weight_N, gas_data[gas_data.gas .== "N", :mol_weight][1])
#    set_param!(m, :gtc2ppm, conversions[conversions.gases .== "CO2", :emiss2conc][1])

    run(m)

