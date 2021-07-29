module cLNA

using InlineExports
using DifferentialEquations
using LaTeXStrings, Plots, Plots.PlotMeasures
using Statistics
using Distributions
using Parameters: @unpack
using LinearAlgebra: norm

module Sim
	using Distributions
	using StatsBase
	using Random
	using InlineExports
	using Distributed

	include("cLNA_SSA_structs.jl")
	include("cLNA_SSA.jl")
end

module Symb
	using InlineExports
	using Symbolics

	include("cLNA_symb.jl")
end

include("cLNA_structs.jl")
include("cLNA_models.jl")
include("cLNA_tests.jl")

function __init__()
	# Set the backend and canvas size
	pyplot(size=(1000,600))
end

end # module
