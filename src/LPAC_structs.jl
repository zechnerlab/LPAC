# Structures and types to be used with LPAC

@export struct Model
    parameters::Union{Tuple, Dict{Symbol, T}} where T
    Ω::Number
    Ωc::Number
    # labelsOde::Vector{Symbol}
    # labelsSSA::Vector{Symbol}
    momentsInit::Function
    momentsOde::Function
    ssaSystem::Sim.System
    momentsMapping::Dict{Symbol,Int}
    sigmaMapping::Dict{Symbol, Tuple{Symbol,Symbol}}
end

@export MomentDict = Dict{Symbol, Vector}
@export SigmaDict = Dict{Symbol, Union{Vector,Nothing}}

@export struct Solution
    T::Vector
    M::MomentDict
    σ::SigmaDict
end

@export struct StandardDeviation
    T::Vector
    M::MomentDict
    σ::SigmaDict
end

@export struct SSAsolution
    T
    N
    σN
    M
    σM
end

#eof
