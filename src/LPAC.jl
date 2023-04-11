module LPAC

using InlineExports
using DifferentialEquations
using LaTeXStrings, Plots, Plots.PlotMeasures
using Statistics
using Distributions
using Parameters: @unpack
using LinearAlgebra: norm, normalize
import Base: *, occursin # For overriding
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

	include("LPAC_SSA_structs.jl")
	include("LPAC_SSA.jl")
end

module Symb
	using InlineExports
	using Symbolics

	include("LPAC_symb.jl")
end

include("LPAC_structs.jl")
include("LPAC_utils.jl")
include("LPAC_solvers.jl")
include("LPAC_plotter.jl")

module Experimental
	using ..LPAC
	using InlineExports
	using Distributions

	include("LPAC_experimental.jl")
end

module Models
	using ..LPAC
	using InlineExports
	using Distributions
	using Parameters: @unpack

	include("LPAC_models.jl")
	# Additional models to test correlations
	include("LPAC_models_correlations.jl")
end

module Figures
	using ..LPAC # Use exported names from the parent LPAC module
	using InlineExports
	using LaTeXStrings, Measures, Plots
	using Serialization

	function __init__()
		# Set the backend and canvas size
		pyplot(size=(1000,600))
	end

	include("LPAC_figures.jl")
end

end # module
