

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
include("src/helper_functions.jl")

include(joinpath("src/components", "ch4_cycle.jl"))
include(joinpath("src/components", "n2o_cycle.jl"))
include(joinpath("src/components", "co2_cycle.jl"))
include(joinpath("src/components", "other_ghg_cycles.jl"))
include(joinpath("src/components", "o3_depleting_substance_cycles.jl"))
include(joinpath("src/components", "co2_forcing.jl"))
include(joinpath("src/components", "ch4_forcing.jl"))
include(joinpath("src/components", "n2o_forcing.jl"))
include(joinpath("src/components", "o3_forcing.jl"))
include(joinpath("src/components", "aerosol_direct_forcing.jl"))
include(joinpath("src/components", "aerosol_indirect_forcing.jl"))
include(joinpath("src/components", "other_ghg_forcing.jl"))
include(joinpath("src/components", "o3_depleting_substance_forcing.jl"))
include(joinpath("src/components", "contrails_forcing.jl"))
include(joinpath("src/components", "black_carbon_snow_forcing.jl"))
include(joinpath("src/components", "landuse_forcing.jl"))
include(joinpath("src/components", "total_forcing.jl"))
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
    add_comp!(m, o3_depleting_substance_cycles)
    add_comp!(m, co2_forcing)
    add_comp!(m, ch4_forcing)
    add_comp!(m, n2o_forcing)
    add_comp!(m, o3_forcing)
    add_comp!(m, aerosol_direct_forcing)
    add_comp!(m, aerosol_indirect_forcing)
    add_comp!(m, other_ghg_forcing)
    add_comp!(m, o3_depleting_substance_forcing)
    add_comp!(m, contrails_forcing)
    add_comp!(m, bc_snow_forcing)
    add_comp!(m, landuse_forcing)
    add_comp!(m, total_forcing)
    add_comp!(m, temperature)


### OPTIONAL JUST FOR TROUBLESHOOTIN
#set_param!(m, :temperature, python_temperature)

    # ---------------------------------------------
    # Set component-specific parameters
    # ---------------------------------------------

    # ---- Carbon Cycle ---- #

    set_param!(m, :co2_cycle, :CO₂_0,  278.052989189439) # From FAIR model run.
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
    set_param!(m, :co2_cycle, :E_co2, rcp_emissions.FossilCO2 .+ rcp_emissions.OtherCO2)
    #set_param!(m, :co2_cycle, :temperature, python_temperature)
    set_param!(m, :co2_cycle, :cumulative_emissions_CO2₀, 0.003)
    set_param!(m, :co2_cycle, :airborne_emissions_CO2₀, 0.0)
    set_param!(m, :co2_cycle, :iIRF_max, 97.0)
    #set_param!(m, :co2_cycle, :dt, 1.0)
    
    connect_param!(m, :co2_cycle => :temperature, :temperature => :T)

    #set_param!(m, :co2_cycle, :temperature, python_temperature)

    # ---- Methane Cycle ---- #
    set_param!(m, :ch4_cycle, :fossil_emiss_CH₄, rcp_emissions.CH4)
    set_param!(m, :ch4_cycle, :natural_emiss_CH₄, rcp_emissions.NaturalCH4)
    set_param!(m, :ch4_cycle, :τ_CH₄, 9.3)
    set_param!(m, :ch4_cycle, :fossil_frac, gas_fractions.ch4_fossil)
    set_param!(m, :ch4_cycle, :oxidation_frac, 0.61)
    set_param!(m, :ch4_cycle, :mol_weight_CH₄, gas_data[gas_data.gas .== "CH4", :mol_weight][1])
    set_param!(m, :ch4_cycle, :mol_weight_C, gas_data[gas_data.gas .== "C", :mol_weight][1])
    set_param!(m, :ch4_cycle, :emiss2conc_ch4, conversions[conversions.gases .== "CH4", :emiss2conc][1])
    set_param!(m, :ch4_cycle, :CH₄_0, gas_data[gas_data.gas .== "CH4", :pi_conc_ar6][1])

    # ---- Nitrous Oxide Cycle ---- #
    set_param!(m, :n2o_cycle, :fossil_emiss_N₂O, rcp_emissions.N2O)
    set_param!(m, :n2o_cycle, :natural_emiss_N₂O, rcp_emissions.NaturalN2O)
    set_param!(m, :n2o_cycle, :τ_N₂O, 121.0)
    set_param!(m, :n2o_cycle, :emiss2conc_n2o, conversions[conversions.gases .== "N2O", :emiss2conc][1])
    set_param!(m, :n2o_cycle, :N₂O_0, gas_data[gas_data.gas .== "N2O", :pi_conc_ar6][1])


    # ---- Other Well-Mixed Greenhouse Gas Cycles ---- #
    set_param!(m, :other_ghg_cycles, :τ_other_ghg, gas_data[findall((in)(other_ghg_names), gas_data.gas), :lifetimes])
    set_param!(m, :other_ghg_cycles, :emiss_other_ghg, Matrix(rcp_emissions[!,Symbol.(other_ghg_names)]))
    set_param!(m, :other_ghg_cycles, :emiss2conc_other_ghg, conversions[findall((in)(other_ghg_names), conversions.gases), :emiss2conc])
    set_param!(m, :other_ghg_cycles, :other_ghg_0, gas_data[findall((in)(other_ghg_names), gas_data.gas), :pi_conc_ar6])

    # ---- Ozone-Depleting Substance Gas Cycles ---- #
    set_param!(m, :o3_depleting_substance_cycles, :τ_ods, gas_data[findall((in)(ods_names), gas_data.gas), :lifetimes])
    set_param!(m, :o3_depleting_substance_cycles, :emiss_ods, Matrix(rcp_emissions[!,Symbol.(ods_names)]))
    set_param!(m, :o3_depleting_substance_cycles, :emiss2conc_ods, conversions[findall((in)(ods_names), conversions.gases), :emiss2conc])
    set_param!(m, :o3_depleting_substance_cycles, :ods_0, gas_data[findall((in)(ods_names), gas_data.gas), :pi_conc_ar6])

    # ---- Carbon Dioxide Radiative Forcing ---- #
    set_param!(m, :co2_forcing, :a₁, -2.4785e-07)
    set_param!(m, :co2_forcing, :b₁, 0.00075906)
    set_param!(m, :co2_forcing, :c₁, -0.0021492)
    set_param!(m, :co2_forcing, :d₁, 5.2488)
    set_param!(m, :co2_forcing, :adjust_F2x, true)
    #set_param!(m, :co2_forcing, :N₂O)
    #set_param!(m, :co2_forcing, :CO₂)
    #set_param!(m, :co2_forcing, :rf_scale_CO₂, 1.0)
    #set_param!(m, :co2_forcing, :N₂O, python_concentrations[:,3])
    #set_param!(m, :co2_forcing, :CO₂, python_concentrations[:,1])
    connect_param!(m, :co2_forcing => :CO₂, :co2_cycle => :co2)
    connect_param!(m, :co2_forcing => :N₂O, :n2o_cycle => :N₂O)


    # ---- Methane Radiative Forcing ---- #
    set_param!(m, :ch4_forcing, :a₃, -8.9603e-05)
 set_param!(m, :ch4_forcing, :b₃, -0.00012462)
 set_param!(m, :ch4_forcing, :d₃, 0.045194)
 #set_param!(m, :ch4_forcing, :rf_scale_CH₄, 1.0)
 set_param!(m, :ch4_forcing, :h2o_from_ch4, 0.079047)
  connect_param!(m, :ch4_forcing => :N₂O, :n2o_cycle => :N₂O)
    connect_param!(m, :ch4_forcing => :CH₄, :ch4_cycle => :CH₄)
 #set_param!(m, :ch4_forcing, :CH₄, python_concentrations[:,2])
 #set_param!(m, :ch4_forcing, :N₂O, python_concentrations[:,3])

    # ---- Nitrous Oxide Radiative Forcing ---- #
    set_param!(m, :n2o_forcing, :a₂, -0.00034197)
 set_param!(m, :n2o_forcing, :b₂,  0.00025455)
 set_param!(m, :n2o_forcing, :c₂, -0.00024357)
 set_param!(m, :n2o_forcing, :d₂, 0.12173)
 #set_param!(m, :n2o_forcing, :rf_scale_N₂O, 1.0)
    connect_param!(m, :n2o_forcing => :CO₂, :co2_cycle => :co2)
    connect_param!(m, :n2o_forcing => :N₂O, :n2o_cycle => :N₂O)
    connect_param!(m, :n2o_forcing => :CH₄, :ch4_cycle => :CH₄)

 #set_param!(m, :n2o_forcing, :CO₂, python_concentrations[:,1])
 #set_param!(m, :n2o_forcing, :CH₄, python_concentrations[:,2])
 #set_param!(m, :n2o_forcing, :N₂O, python_concentrations[:,3])


    # ---- Ozone Radiative Forcing ---- #
    set_param!(m, :o3_forcing, :Br, gas_data[findall((in)(ods_names), gas_data.gas), :br_atoms])
    set_param!(m, :o3_forcing, :Cl, gas_data[findall((in)(ods_names), gas_data.gas), :cl_atoms])
    set_param!(m, :o3_forcing, :FC, gas_data[findall((in)(ods_names), gas_data.gas), :strat_frac])
#    set_param!(m, :o3_forcing, :ODS_pi, gas_data[findall((in)(ods_names), gas_data.gas), :pi_conc_ar6])

    set_param!(m, :o3_forcing, :total_forcing_O₃_0, 0.0)

    #set_param!(m, :o3_forcing, :NOx_emissions_pi, 0.0)
    #set_param!(m, :o3_forcing, :CO_emissions_pi, 0.0)
    #set_param!(m, :o3_forcing, :NMVOC_emissions_pi, 0.0)
    #set_param!(m, :o3_forcing, :NOx_emissions, rcp_emissions.NOx)
    #set_param!(m, :o3_forcing, :CO_emissions, rcp_emissions.CO)
    #set_param!(m, :o3_forcing, :NMVOC_emissions, rcp_emissions.NMVOC) 

    set_param!(m, :o3_forcing, :feedback, -0.037)

    set_param!(m, :o3_forcing, :Ψ_CH₄, 2.33379720e-04)
    set_param!(m, :o3_forcing, :Ψ_N₂O, 1.27179106e-03)
    set_param!(m, :o3_forcing, :Ψ_ODS, -6.69347820e-05)
    set_param!(m, :o3_forcing, :Ψ_CO, 1.14647701e-04)
    set_param!(m, :o3_forcing, :Ψ_NMVOC, 5.14366051e-12)
    set_param!(m, :o3_forcing, :Ψ_NOx, 3.78354423e-03)

#connect_param!(m, :o3_forcing => :conc_ODS, :o3_depleting_cycles => :conc_ods)
    
    connect_param!(m, :o3_forcing => :temperature, :temperature => :T)

#    set_param!(m, :o3_forcing,  :temperature, python_temperature)
    #set_param!(m, :o3_forcing,  :conc_ODS, python_concentrations[:,16:31])
  connect_param!(m, :o3_forcing => :conc_ODS, :o3_depleting_substance_cycles => :conc_ods)
  connect_param!(m, :o3_forcing => :N₂O, :n2o_cycle => :N₂O)
    connect_param!(m, :o3_forcing => :CH₄, :ch4_cycle => :CH₄)

    # ---- Aerosol Direct Radiative Forcing ---- #
    set_param!(m, :aerosol_direct_forcing, :β_SOx, -6.2227e-3)
    set_param!(m, :aerosol_direct_forcing, :β_CO, 0.0)
    set_param!(m, :aerosol_direct_forcing, :β_NMVOC, -3.8392e-4)
    set_param!(m, :aerosol_direct_forcing, :β_NOx, -1.16551e-3)
    set_param!(m, :aerosol_direct_forcing, :β_BC, 1.601537e-2)
    set_param!(m, :aerosol_direct_forcing, :β_OC, -1.45339e-3)
    set_param!(m, :aerosol_direct_forcing, :β_NH3, -1.55605e-3)
    #set_param!(m, :aerosol_direct_forcing, :rf_scale_aero_direct, 0.0)


    # ---- Aerosol Indirect Radiative Forcing ---- #
    set_param!(m, :aerosol_indirect_forcing, :ϕ, 0.07334277994353743)#-1.95011431)
    set_param!(m, :aerosol_indirect_forcing, :b_SOx, 3.452849302362568)#0.01107147)
    set_param!(m, :aerosol_indirect_forcing, :b_POM, 33.126485122209154)#0.01387492)
    set_param!(m, :aerosol_indirect_forcing, :rf_scale_aero_indirect, 1.0)
    #set_param!(m, :aerosol_indirect_forcing, :SOx_emiss_pi, 0.0)
    #set_param!(m, :aerosol_indirect_forcing, :BC_emiss_pi, 0.0)
    #set_param!(m, :aerosol_indirect_forcing, :OC_emiss_pi, 0.0)
    #set_param!(m, :aerosol_indirect_forcing, :model_years, collect(start_year:end_year))
    #set_param!(m, :aerosol_indirect_forcing, :SOx_emiss_1765, 1.0)
    #set_param!(m, :aerosol_indirect_forcing, :BC_OC_emiss_1765, 11.2)
    #set_param!(m, :aerosol_indirect_forcing, :scale_AR5, true)
    #set_param!(m, :aerosol_indirect_forcing, :F_1765, -0.3002836449793625)
    #set_param!(m, :aerosol_indirect_forcing, :F_2011, -1.5236182344467388)
    #set_param!(m, :aerosol_indirect_forcing, :fix_pre1850_RCP, false)

    # ---- Other Well-Mixed Greenhouse Gas Radiative Forcings ---- #
    set_param!(m, :other_ghg_forcing, :other_ghg_radiative_efficiency, gas_data[findall((in)(other_ghg_names), gas_data.gas), :rad_eff])
  #  set_param!(m, :other_ghg_forcing, :other_ghg_pi, gas_data[findall((in)(other_ghg_names), gas_data.gas), :pi_conc])
    connect_param!(m, :other_ghg_forcing => :conc_other_ghg, :other_ghg_cycles => :conc_other_ghg)

    # ---- Ozone-Depleting Substance Radiative Forcings ---- #
    set_param!(m, :o3_depleting_substance_forcing, :ods_radiative_efficiency, gas_data[findall((in)(ods_names), gas_data.gas), :rad_eff])
   # set_param!(m, :o3_depleting_substance_forcing, :ods_pi, gas_data[findall((in)(ods_names), gas_data.gas), :pi_conc])
    connect_param!(m, :o3_depleting_substance_forcing => :conc_ods, :o3_depleting_substance_cycles => :conc_ods)

    # ---- Contrails Radiative Forcing ---- #
    set_param!(m, :contrails_forcing, :frac, gas_fractions.nox_aviation)
    set_param!(m, :contrails_forcing, :E_ref_contrails, 2.946)
    set_param!(m, :contrails_forcing, :F_ref_contrails, 0.0448)
    set_param!(m, :contrails_forcing, :ref_is_NO2, true)
    set_param!(m, :contrails_forcing, :mol_weight_NO₂, gas_data[gas_data.gas .== "NO2", :mol_weight][1])
    set_param!(m, :contrails_forcing, :mol_weight_N, gas_data[gas_data.gas .== "N", :mol_weight][1])
    #set_param!(m, :contrails_forcing, :rf_scale_contrails, 0.0) # Default FAIR has contrail forcing switched off.

    # ---- Black Carbon on Snow Radiative Forcing ---- #
    set_param!(m, :bc_snow_forcing, :E_ref_bc, 6.095)
    set_param!(m, :bc_snow_forcing, :F_ref_bc, 0.08)

    # ---- Land Use Change Radiative Forcing ---- #
    set_param!(m, :landuse_forcing, :α_CO₂_land, (-0.2/190))
    set_param!(m, :landuse_forcing, :landuse_emiss, rcp_emissions.OtherCO2)

    set_param!(m, :total_forcing, :scale_CO₂, 1.0)
    set_param!(m, :total_forcing, :scale_CH₄, 1.0)
    set_param!(m, :total_forcing, :scale_CH₄_H₂O, 1.0)
    set_param!(m, :total_forcing, :scale_N₂O, 1.0)
    set_param!(m, :total_forcing, :scale_O₃, 1.0)
    set_param!(m, :total_forcing, :scale_aerosol_indirect, 1.0)
    set_param!(m, :total_forcing, :scale_bcsnow, 1.0)
    set_param!(m, :total_forcing, :scale_landuse, 1.0)
    set_param!(m, :total_forcing, :scale_contrails, 0.0) #!!! Default FAIR has contrail forcing switched off. Set scaling term to 0
    set_param!(m, :total_forcing, :scale_volcanic, 1.0)
    set_param!(m, :total_forcing, :scale_solar, 1.0)
    set_param!(m, :total_forcing, :scale_aerosol_direct_SOx, 1.0)
    set_param!(m, :total_forcing, :scale_aerosol_direct_CO_NMVOC, 1.0)
    set_param!(m, :total_forcing, :scale_aerosol_direct_NOx_NH3, 1.0)
    set_param!(m, :total_forcing, :scale_aerosol_direct_BC, 1.0)
    set_param!(m, :total_forcing, :scale_aerosol_direct_OC, 1.0)
    set_param!(m, :total_forcing, :scale_other_ghg, ones(length(other_ghg_names)))
    set_param!(m, :total_forcing, :scale_ods, ones(length(ods_names)))

    connect_param!(m, :total_forcing => :F_CO₂, :co2_forcing => :rf_co2)
    connect_param!(m, :total_forcing => :F_CH₄, :ch4_forcing => :rf_ch4)
    connect_param!(m, :total_forcing => :F_CH₄_H₂O, :ch4_forcing => :rf_ch4_h2o)
    connect_param!(m, :total_forcing => :F_N₂O, :n2o_forcing => :rf_n2o)
    connect_param!(m, :total_forcing => :F_O₃, :o3_forcing => :total_forcing_O₃)

    connect_param!(m, :total_forcing => :F_aerosol_direct_SOx, :aerosol_direct_forcing => :F_SOx_aero)
    connect_param!(m, :total_forcing => :F_aerosol_direct_CO, :aerosol_direct_forcing => :F_CO_aero)
    connect_param!(m, :total_forcing => :F_aerosol_direct_NMVOC, :aerosol_direct_forcing => :F_NMVOC_aero)
    connect_param!(m, :total_forcing => :F_aerosol_direct_NOx, :aerosol_direct_forcing => :F_NOx_aero)
    connect_param!(m, :total_forcing => :F_aerosol_direct_BC, :aerosol_direct_forcing => :F_BC_aero)
    connect_param!(m, :total_forcing => :F_aerosol_direct_OC, :aerosol_direct_forcing => :F_OC_aero)
    connect_param!(m, :total_forcing => :F_aerosol_direct_NH3, :aerosol_direct_forcing => :F_NH3_aero)

    connect_param!(m, :total_forcing => :F_aerosol_indirect, :aerosol_indirect_forcing => :rf_aero_indirect)

    connect_param!(m, :total_forcing => :F_bcsnow, :bc_snow_forcing => :forcing_BC_snow)
    connect_param!(m, :total_forcing => :F_landuse, :landuse_forcing => :forcing_landuse)
    connect_param!(m, :total_forcing => :F_contrails, :contrails_forcing => :forcing_contrails)
    connect_param!(m, :total_forcing => :F_other_ghg, :other_ghg_forcing => :other_ghg_rf)
    connect_param!(m, :total_forcing => :F_ods, :o3_depleting_substance_forcing => :ods_rf)



    #set_param!(m, :total_forcing, :F_CO₂, python_forcing[:,1])
    #set_param!(m, :total_forcing, :F_CH₄, python_forcing[:,2])
    #set_param!(m, :total_forcing, :F_CH₄_H₂O, , python_forcing[:,34])
    #set_param!(m, :total_forcing, :F_N₂O, , python_forcing[:, 3])
    #set_param!(m, :total_forcing, :F_O₃, , python_forcing[:, 32])
  #  set_param!(m, :total_forcing, :F_aerosol_direct_SOx, python_forcing[:, 36])
  #  set_param!(m, :total_forcing, :F_aerosol_direct_CO, python_forcing[:, 37])
  #  set_param!(m, :total_forcing, :F_aerosol_direct_NMVOC, zeros(length(start_year:end_year))) # need to trick because these should be summed
  #  set_param!(m, :total_forcing, :F_aerosol_direct_NOx, python_forcing[:, 38])
  #  set_param!(m, :total_forcing, :F_aerosol_direct_BC, python_forcing[:, 39])
  #  set_param!(m, :total_forcing, :F_aerosol_direct_OC, python_forcing[:, 40])
  #  set_param!(m, :total_forcing, :F_aerosol_direct_NH3, zeros(length(start_year:end_year)))
  #  set_param!(m, :total_forcing, :F_aerosol_indirect, python_forcing[:, 41])
  #  set_param!(m, :total_forcing, :F_bcsnow, python_forcing[:, 42])
  #  set_param!(m, :total_forcing, :F_landuse, python_forcing[:, 43])
  #  set_param!(m, :total_forcing, :F_contrails, python_forcing[:, 35])
  #  set_param!(m, :total_forcing, :F_other_ghg, python_forcing[:, 4:15])
  #  set_param!(m, :total_forcing, :F_ods, python_forcing[:, 16:31])
    set_param!(m, :total_forcing, :F_volcanic, python_forcing[:, 44])
    set_param!(m, :total_forcing, :F_solar, python_forcing[:, 45])
    set_param!(m, :total_forcing, :F_exogenous, zeros(length(start_year:end_year)))







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
    connect_param!(m, :temperature => :forcing, :total_forcing => :total_forcing)
    #set_param!(m, :temperature, :forcing, sum(Matrix(python_forcing), dims=2)[:])



    # ---- Parameters Shared Across Multiple Components ---- #
    set_param!(m, :dt, 1.0)
    set_param!(m, :emiss2conc_co2, conversions[conversions.gases .== "CO2", :emiss2conc][1])

    set_param!(m, :CH₄_pi, gas_data[gas_data.gas .== "CH4", :pi_conc_ar6][1])
    set_param!(m, :CO₂_pi, gas_data[gas_data.gas .== "CO2", :pi_conc_ar6][1])
    set_param!(m, :N₂O_pi, gas_data[gas_data.gas .== "N2O", :pi_conc_ar6][1])
    set_param!(m, :ods_pi, gas_data[findall((in)(ods_names), gas_data.gas), :pi_conc_ar6])
    set_param!(m, :other_ghg_pi, gas_data[findall((in)(other_ghg_names), gas_data.gas), :pi_conc_ar6])



    set_param!(m, :F2x, 3.71)

    set_param!(m, :SOx_emiss, rcp_emissions.SOx)
    set_param!(m, :BC_emiss, rcp_emissions.BC)
    set_param!(m, :OC_emiss, rcp_emissions.OC)
    set_param!(m, :CO_emiss, rcp_emissions.CO)
    set_param!(m, :NMVOC_emiss, rcp_emissions.NMVOC)
    set_param!(m, :NH3_emiss, rcp_emissions.NH3)
    set_param!(m, :NOx_emiss, rcp_emissions.NOx)

    set_param!(m, :SOx_emiss_pi,   0.0)
    set_param!(m, :CO_emiss_pi,    0.0)
    set_param!(m, :NMVOC_emiss_pi, 0.0)
    set_param!(m, :NOx_emiss_pi,   0.0)
    set_param!(m, :BC_emiss_pi,    0.0)
    set_param!(m, :OC_emiss_pi,    0.0)

#    set_param!(m, :SOx_emiss_pi,   1.22002422)
#    set_param!(m, :CO_emiss_pi,    348.527359)
#    set_param!(m, :NMVOC_emiss_pi, 60.0218262)
#    set_param!(m, :NOx_emiss_pi,   3.87593407)
#    set_param!(m, :BC_emiss_pi,    2.09777075)
#    set_param!(m, :OC_emiss_pi,    15.4476682)


    #set_param!(m, :N₂O_0, gas_data[gas_data.gas .== "N2O", :pi_conc][1])



#    set_param!(m, :SOx_emiss, rcp_emissions.SOx)
#    set_param!(m, :BC_emiss, rcp_emissions.BC)
#    set_param!(m, :OC_emiss, rcp_emissions.OC)
#    set_param!(m, :NOx_emiss, rcp_emissions.NOx)
#    set_param!(m, :fix_pre1850_RCP, true)
#    set_param!(m, :mol_weight_N, gas_data[gas_data.gas .== "N", :mol_weight][1])
#    set_param!(m, :gtc2ppm, conversions[conversions.gases .== "CO2", :emiss2conc][1])



#=
 
    run(m)


julia_rcp85 = deepcopy(m[:temperature, :T])
julia_rcp26 = deepcopy(m[:temperature, :T])

python_rcp85 = deepcopy(python_temperature)
python_rcp26 = Matrix(DataFrame(load(joinpath(@__DIR__, "data", "rcp26_temperature.csv"), header_exists=false)))[:]

R"""
plot(1765:2500, $python_rcp85, type="l", col="red", lwd=4, xlab="Year", ylab="Degrees C", main = "Python vs. Julia FAIR v1.6.2")
lines(1765:2500, $python_rcp26, col="red", lwd=4)
lines(1765:2500, $julia_rcp85, col="blue")
lines(1765:2500, $julia_rcp26, col="blue")
"""
=#
