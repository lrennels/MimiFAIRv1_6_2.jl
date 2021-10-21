# MimiFAIRv1_6_2.jl

This is a work-in-progress repository for a Julia-Mimi implementation of the FAIRv1_6_2 simple climate model.

## Preparing the Software Environment

To add the package to your current environment, run the following command at the julia package REPL:

```julia
pkg> add https://github.com/FrankErrickson/MimiFAIRv1_6_2.jl.git
```

You probably also want to install the Mimi package into your julia environment, so that you can use some of the tools in there:

```julia
pkg> add Mimi
```

## Running the Model

The model uses the Mimi framework and it is highly recommended to read the Mimi  documentation first to understand the code structure. The basic way to access a copy of the default MimiFAIRv2 model and explore the resuts is the following:

```julia
using Mimi 
using MimiFAIRv1_6_2

# Create an instance of MimiFAIRv1_6_2.
m = MimiFAIRv1_6_2.get_model() 

# Run the model.
run(m)

# Access the temperature output.
fair_temps = m[:temperature, :T]

# Explore interactive plots of all the model output.
explore(m)
```

The `get_model()` function currently has the following keyword arguments:  

* `ar6_scenario`: One of the RCMIP scenarios from the original FAIRv2.0 paper. Current options include "ssp119", "ssp126", "ssp245", "ssp370", and"ssp585". The default is "ssp245".  
* `start_year`: The model has an option to be initialized at different time periods, however this is only currently set up to start in 1750.
* `end_year`: The model can be run out to 2500 (the default final year).  

## Running a Monte Carlo Simulation

### Overview

See `mcs/AR6_Monte_Carlo.jl` for the script and details on running a Monte Carlo Simulation, and don't hesitate to contact the developers with any questions, or post an Issue on Github.

This file contains functions to run a Monte Carlo with the AR6 implementation of MimiFAIRv1_6_2 using the constrained parameters from AR6. The primary function, `run_mcs`, loads the constrained parameter samples and then performs the analysis. 

Note that the original parameter samples are stored in `model_data/fair-1.6.2-wg3-params.json` and then are parsed into `mcs_params` using the function `parse_mcs_params` in `mcs/utils.jl`.  This does not need to be repeated, but can be useful for replication and understanding.

### Function Details

The `run_mcs` function is the primary user-facing function provided for the monte carlo simulation and has the signature and function arguments as follows:

```julia
    run_mcs(;trials::Int64 = 2237, 
        output_dir::Union{String, Nothing} = nothing, 
        save_trials::Bool = false,
        m::Mimi.Model = get_model())
```

This function returns the results of a Monte Carlo Simulation with the defined 
number of trials and save data into the `output_dir` folder, optionally also saving 
trials if `save_trials` is set to `true.` If no model is provided, use the default
model returned by get_model(). If an `output_dir` is not provided, data will be
saved to the `output` folder in this repository in a subfolder named based on the 
Date, Time, and number of trials.

Call this function as follows:

```julia
results = MimiFAIRv1_6_2.run_mcs(trials = 100, output_directory = path, save_trials = true)
explore(results)
Mimi.plot(results, :temperature, :T; interactive = true)
```
will run a Monte Carlo simulation with 1000 trials, and return a Mimi.SimulationInstance object that can be `explore`d with a UI (note this is fairly slow at the moment it is under improvement), or display a particular plot for an output variable. Output variable data and trials data will be saved in `path`, or if this isn't provided in the `output` folder in this repository in a subfolder named based on the Date, Time, and number of trials.
\
\
![Illustrative Example of FAIR Temperatures (n = 100)](https://github.com/FrankErrickson/MimiFAIRv1_6_2.jl/blob/main/plot_1.svg)

The output variables, currently temperature and co2, will be saved to the `output_directory` as will all trials values in `trials.csv`.  Adding more variables to output is a matter of augmenting the following section of code.  Feel free to contact the authors with requests on more outputs, or open a PR doing so yourself.
```julia
# define the Monte Carlo Simulation
mcs = @defsim begin
    save(temperature.T, co2_cycle.co2)
end
```
