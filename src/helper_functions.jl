# #----------------------------------------------------------------------------------------------------------------------
# #----------------------------------------------------------------------------------------------------------------------
# This file contains functions and other snippets of code that are used in various calculations for Mimi-FAIR.
# #----------------------------------------------------------------------------------------------------------------------
# #----------------------------------------------------------------------------------------------------------------------



#######################################################################################################################
# LOAD ALL DATA AND PARAMETER VALUES NEEDED FOR Mimi-FAIR
########################################################################################################################
# Description: This function loads and cleans up the Mimi-FAIR data so it can be more easily incorporated into the model.
#
# Function Arguments:
#
#       start_year:      First year to run the model.
#       end_year:        Final year to run the model.
#       rcp_scenario:    A string indicating which RCP scenario to use ("RCP26", "RCP45", "RCP60", & "RCP85").
#----------------------------------------------------------------------------------------------------------------------

function load_fair_data(start_year::Int, end_year::Int, rcp_scenario::String)

    # Calculate indicies to extract RCP data (RCP data spans 1765-2500)
    start_index, end_index = indexin([start_year, end_year], collect(1765:2500))

    # Create vector of names for minor greenhouse gases to loop over.
    other_ghg_names = ["CF4", "C2F6", "C6F14", "HFC23", "HFC32", "HFC43_10", "HFC125", "HFC134a", "HFC143a", "HFC227ea", "HFC245fa", "HFC152a", "HFC236fa", "HFC365mfc", "SF6", "CFC_11", "CFC_12", "CFC_113", "CFC_114", "CFC_115", "CARB_TET", "MCF", "HCFC_22", "HCFC_141B", "HCFC_142B", "HALON1211", "HALON1202", "HALON1301", "HALON2402", "CH3BR", "CH3CL"]
    additional_hfcs = ["HFC152a", "HFC236fa", "HFC365mfc"]

    #---------------------------------------
    # Read in relevant data sets.
    #---------------------------------------

    # RCP scenario emissions.
    rcp_emissions_raw = DataFrame(load(joinpath(@__DIR__, "..", "data", "model_data", rcp_scenario*"_EMISSIONS.csv"), skiplines_begin=36))
    # Natural emissions for methane and nitrous oxide as estimated by FAIR team.
    natural_emissions_raw = DataFrame(load(joinpath(@__DIR__, "..", "data", "model_data", "natural_emissions.csv"), skiplines_begin=3))
    # CMIP6 Solar forcing.
    cmip6_solar_forcing = DataFrame(load(joinpath(@__DIR__, "..", "data", "model_data", "cmip6_solar.csv"), skiplines_begin=6))[start_index:end_index, Symbol("Radiative forcing")]
    # CMIP6 volcanic forcing.
    cmip6_volcano_forcing = DataFrame(load(joinpath(@__DIR__, "..", "data", "model_data", "cmip6_volcanic.csv"), skiplines_begin=8))[start_index:end_index, Symbol("Volcanic ERF")]
    # Fraction of NOx emissions attibutable to aviation (for contrail RF).
    aviation_fraction_raw = DataFrame(load(joinpath(@__DIR__, "..", "data", "model_data", "aviNOx_fraction.csv"), skiplines_begin=4))[:,Symbol(rcp_scenario)]
    # Time-varying shares of anthropogenic methane attribuatable to fossil sources.
    ch4_fossil_frac_raw = DataFrame(load(joinpath(@__DIR__, "..", "data", "model_data", "fossilCH4_fraction.csv"), skiplines_begin=4))[:,Symbol(rcp_scenario)]
    # Information on various gas specieis.
    gas_data = DataFrame(load(joinpath(@__DIR__, "..", "data", "model_data", "fair_ghg_species_data.csv"), skiplines_begin=11))

    #---------------------------------------
    # Emissions
    #---------------------------------------
    emissions = DataFrame()

    emissions.FossilCO2   = rcp_emissions_raw[start_index:end_index, :FossilCO2]
    emissions.OtherCO2    = rcp_emissions_raw[start_index:end_index, :OtherCO2]
    emissions.CH4         = rcp_emissions_raw[start_index:end_index, :CH4]
    emissions.NaturalCH4  = natural_emissions_raw[start_index:end_index, :ch4]
    emissions.N2O         = rcp_emissions_raw[start_index:end_index, :N2O]
    emissions.NaturalN2O  = natural_emissions_raw[start_index:end_index, :n2o]
    emissions.NMVOC       = rcp_emissions_raw[start_index:end_index, :NMVOC]
    emissions.CO          = rcp_emissions_raw[start_index:end_index, :CO]
    emissions.NOx         = rcp_emissions_raw[start_index:end_index, :NOx]
    emissions.SOx         = rcp_emissions_raw[start_index:end_index, :SOx]
    emissions.BC          = rcp_emissions_raw[start_index:end_index, :BC]
    emissions.OC          = rcp_emissions_raw[start_index:end_index, :OC]
    emissions.NH3         = rcp_emissions_raw[start_index:end_index, :NH3]

    # Other greenhouse gases
    for i in other_ghg_names
        if i in additional_hfcs
            # (HFC ADDITION NOTE) we have not pulled these HFCs for the raw RCP 
            # scenarios yet, but we don't need them since we base data on AR6, so
            # just use dummy values to prevent errors -- these missing data will not
            # be used, and if they are Mimi will error when indexing a missing
            emissions[!,Symbol(i)] = fill(missing, length(start_index:end_index))
        else
            emissions[!,Symbol(i)] = rcp_emissions_raw[start_index:end_index, Symbol(i)]
        end
    end

    #---------------------------------------
    # Gas Fractions
    #---------------------------------------
    gas_fractions = DataFrame()

    gas_fractions.nox_aviation = aviation_fraction_raw[start_index:end_index]
    gas_fractions.ch4_fossil = ch4_fossil_frac_raw[start_index:end_index]

    #---------------------------------------
    # Emission to Concentration Conversions
    #---------------------------------------
    emiss_conversions = DataFrame()

    # Set names for concentration conversions.
    emiss_conversions.gases = vcat("CO2", "CH4", "N2O", other_ghg_names)

    # Mass of atmosphere (kg).
    mass_atmos = 5.1352e18

    # Molecular weights of air, CO₂, N₂, CH₄, and other greenhouse gases.
    mol_wt_air    = gas_data[gas_data.gas .== "AIR", :mol_weight][1]
    mol_wt_carbon = gas_data[gas_data.gas .== "C", :mol_weight][1]
    mol_wt_n2     = gas_data[gas_data.gas .== "N2", :mol_weight][1]
    mol_wt_ch4    = gas_data[gas_data.gas .== "CH4", :mol_weight][1]
    mol_wt_others = gas_data[indexin(other_ghg_names, gas_data.gas), :mol_weight]

    # Calculate CO₂ conversion from GtC to ppm.
    emiss2conc_carbon = (mass_atmos / 1.0e18) * (mol_wt_carbon / mol_wt_air)

    # Emission to concentration for CH₄.
    emiss2conc_ch4 = (mass_atmos / 1.0e18) * (mol_wt_ch4 / mol_wt_air)

    # Emission to concentration for N₂O.
    # Note: Use N₂ for N₂O conversion, from FAIR: "Funny units for nitrogen emissions - N₂O is expressed in N₂ equivalent."
    emiss2conc_n2o= (mass_atmos / 1.0e18) * (mol_wt_n2 / mol_wt_air)

    # Conversion factors for other Kyoto Protocol or ozone-depleting gases.
    emiss2conc_others = (mass_atmos / 1.0e18) .* (mol_wt_others ./ mol_wt_air)

    # Combine values conversion values into a single data frame.
    emiss_conversions.emiss2conc = vcat(emiss2conc_carbon, emiss2conc_ch4, emiss2conc_n2o, emiss2conc_others)

    return emissions, cmip6_volcano_forcing, cmip6_solar_forcing, gas_data, gas_fractions, emiss_conversions
end