module cLNA

using InlineExports
using DifferentialEquations
using LaTeXStrings, Plots, Plots.PlotMeasures
using Statistics
using Distributions
using Parameters: @unpack
using LinearAlgebra: norm, normalize
import Base: *, occursin # For overriding
using BifurcationKit, Setfield
using Interpolations: linear_interpolation

export Sim, Symb, Models, Figures

module Sim
	using Distributions
	# using StatsBase
	using Random
	using LinearAlgebra: Symmetric
	using Combinatorics: combinations
	using InlineExports
	using Distributed
	using Statistics: mean, var
	using UnicodePlots: spy

	include("cLNA_SSA_structs.jl")
	include("cLNA_SSA.jl")
end

module Symb
	using InlineExports
	using Symbolics

	include("cLNA_symb.jl")
end

include("cLNA_structs.jl")
include("cLNA_utils.jl")
include("cLNA_plotter.jl")

module Models
	using ..cLNA
	using InlineExports
	using Distributions
	using Parameters: @unpack

	include("cLNA_models.jl")
	# Additional models to test correlations
	include("cLNA_models_correlations.jl")
end

module Figures
	using ..cLNA # Use exported names from the parent cLNA module
	using InlineExports
	using LaTeXStrings, Measures, Plots
	using Serialization

	function __init__()
		# Set the backend and canvas size
		pyplot(size=(1000,600))
	end

	include("cLNA_figures.jl")
end

end # module
