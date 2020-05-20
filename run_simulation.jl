# run simulation
using PyPlot
using MAT
using Distributions
using StatsBase

include("calc_ext_input_poisson_based.jl")
include("correlation_based_pattern_high_order.jl")
include("calc_ext_input_corr_based.jl")
include("initialize_network.jl")
include("simulation.jl")
include("simulation_new.jl")

simulation_new()
println("Done")
