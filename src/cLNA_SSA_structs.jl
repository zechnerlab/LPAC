#=
Code from the sAIC project.
Lorenzo Duso
+ features for the cLNA project by Tommaso Bianucci
=#

# TRANSITION CLASS
# object to define the structure and the rate law of a transition class
@export mutable struct TransitionClass
    rc::Int64                  # Number of reactant compartments
    pc::Int64                  # Number of product compartments
    DeltaN::Int64              # for conveniency = pc - rc
    k::Float64                 # Rate constant
    parameters::Union{Nothing,Float64,Vector{Int64}}
    g::Union{Function,Nothing}         # Reactant compartments kernel (takes a vector of reactant state-vectors)
    pi::Union{Function,Nothing}        # Product compartments distribution
    H::Union{Function,Nothing}         # Total propensity
    fast_sample_reactants!::Union{Function,Nothing} # optional, optimized for specific class
    function TransitionClass(rc::Int64, pc::Int64, k::Float64)
            @assert rc >= 0 "Number of reactant compartments can't be negative!"
            @assert pc >= 0 "Number of product compartments can't be negative!"
            @assert rc <= 2 "The simulator doesn't support more than 2 reactant compartments!"
            @assert pc <= 2 "The simulator doesn't support more than 2 product compartments!"
            @assert k >= 0. "Rate constant can't be negative!"
        return new(rc,pc,pc-rc,k,nothing,nothing,nothing,nothing,nothing)
    end
end

# shortcut definition for one-compartment chemical reaction
## change_vector is the integer update vector of the chemical species
@export function new_chemical_reaction_class(change_vector::Vector{Int64}, rate::Float64)
    class=TransitionClass(1,1,rate)
    class.parameters = change_vector
    class.pi = pi_chemical_reaction!
    return class
end
# shortcut definition for compartment update upon chemical reaction
@export function pi_chemical_reaction!(yc::Vector{Vector{Int64}}, xc::Vector{Vector{Int64}}, change_vector::Vector{Int64})
    for d=1:length(yc[1])
        yc[1][d] = xc[1][d] + change_vector[d]
    end
end


# SYSTEM OBJECT
# object to collect a set of transition classes and speficy model details
# the fields moment_equations and init_moments can be optionally assigned to the system's moment equations (not required for SSA)
@export mutable struct System
    name::String
    n_species::Int64                                  # Number of chemical species
    transition_classes::Vector{TransitionClass}       # Defining model dynamics (as used in SSA simulations)
    MomDict::Dict{Int64,Vector{Int64}}                # link moment index to its gamma exponents
    MomMapping::Dict{Int64,Symbol}
    MomExpMapping::Dict{Symbol,Vector{Int64}}
    moment_equations::Union{Function,Nothing}         # Moment Equations f(dM,M,S,t)
    init_moments::Union{Function,Nothing}             # Based on implementation of Moment Equations
    function System(name::String, n_species::Int64, MomDirectMapping::Vector{Pair{Symbol,Vector{Int64}}})
        for (s, exponents) in MomDirectMapping
           @assert size(exponents,1) == n_species "Moment exponents do not match the dimensionality of the compartment content!"
        end
        # Convert MomDirectMapping into MomDict and MomMapping
        MomDict = Dict([ i=>e for ( i, (s,e) ) in enumerate(MomDirectMapping) ]...)
        MomMapping = Dict([ i=>s for ( i, (s,e) ) in enumerate(MomDirectMapping) ]...)
        MomExpMapping = Dict([ s=>e for (s,e) in MomDirectMapping ]...)
        #
        @assert MomDict[1] == zeros(Int64,n_species) "The moment of index 1 must map to the zero-order moment! (i.e. $(zeros(Int64,n_species))) "
        return new(name,n_species,Vector{TransitionClass}(),MomDict,MomMapping,MomExpMapping,nothing,nothing)
    end
end

"""
    PropensitySet = Vector{Union{Array, Nothing}}

This is the type used for the set of propensity tensors
"""
@export PropensitySet = Vector{Union{Array, Nothing}}

# adds some transition classes to a System object S
@export function add_transition_class(S::System,c...)
    @assert prod(map(class -> typeof(class) == TransitionClass, c))
    S.transition_classes = [S.transition_classes; c...]
end

# computes moments of population state n as given by the index -> exponents dictionary in the system S
@export function compute_moments(S::System, n::Matrix{Int64})
    DD,Ncomp=size(n)
    M=zeros(Int64,length(keys(S.MomDict)))
    for k=1:length(M)
        n_vals=ones(Int64,size(n,2))
        for i=1:length(S.MomDict[k])
            n_vals .*= n[i,:].^S.MomDict[k][i]
        end
        M[k]=sum(n_vals)
    end
    return M
end

@export equal(a, b; tol=1e-6) = abs( a - b ) < tol

@export function compute_propensities(S::System, n::Matrix{Int64}, Mom::Vector{Int64})
    nSpecies, nComp = size(n)
    nTransCl = length(S.transition_classes)
    g = Vector{Union{Array, Nothing}}(nothing, nTransCl)
    G = zeros(nTransCl)
    for (icl, cl) in enumerate(S.transition_classes)
        if isnothing(cl.g) || cl.rc==0
            # If the "g" function is not defined for this class we just leave "nothing"
            continue
        elseif !isnothing(cl.H)
            @debug "Both g() and H() are defined for transition class $icl: using H() for better computational efficiency!"
            continue
        end
        # Initialize a tensor which dimensions are the number of reactant compartments
        dimensions = repeat([nComp], cl.rc)
        gCl = zeros(dimensions...)
        # Now fill it with the propensities for all combinations of reactants
        for (ci, _) in pairs(gCl)
            ci = CartesianIndex(ci)
            I = ci.I
            # Here we want to only populate strictly supra-diagonal elements
            # this is in order to account for each couple of reactants just once
            if !issorted(I) || !allunique(I) 
                            # This is clearly suboptimal, but fine since we run 
                continue    # this function only once at startup
            end
            xc = [ n[:,i] for i in I ]
            val = cl.g(xc, Mom)
            gCl[I...] = val
            G[icl] += val
        end
        g[icl] = gCl
        # # DEBUG: Check for consistency
        # @assert equal( sum(g[icl]), G[icl] ) """
        #     FATAL: G is not consistent with the propensity tensor g!
        #     G[$icl] = $(G[icl])
        #     sum(g[$icl]) = $(sum(g[icl]))
        #     """
    end
    return g, G
end

# check that initial condition matches the system
@export function assert_model(S::System,n0::Matrix{Int64})
    @assert size(n0,1) == S.n_species "Initial condition doesn't match the model dimensionality"
    Mom0=compute_moments(S,n0)
    for (cid, c) in enumerate(S.transition_classes)
        if c.H == nothing && c.g == nothing
            error("Incomplete transition class: specify either the class propensity function 'H' or 'g' ")
        elseif !isnothing(c.H) && !isnothing(c.g)
            @info "Both g() and H() are defined for transition class $cid: using H() for better computational efficiency!"
        end
        if !isnothing(c.H)
            try
                c.H(n0,Mom0)
            catch y
                error("The moment dictionary of the provided model is incomplete!")
            end
        end
        c.rc > 0 && c.g == nothing && c.fast_sample_reactants! == nothing ? error("Incomplete transition class: speficy 'g' or 'fast_sample_reactants!' ") : nothing
        c.pc > 0 && c.pi == nothing ? error("Incomplete transition class: speficy the outcome distribution 'pi' ") : nothing
    end
end

### Default fast sampling functions

"""
fast_sample_generic_unary(r_indices::Vector{Int64}, n::Matrix{Int64}, 
                            propensity::Function, totalPropensity::Number)

Provides a generic fast sampling function.
    propensity is a function that must take the vector of chem content of a single compartment.
"""
function fast_sample_generic_unary!(r_indices::Vector{Int64}, n::Matrix{Int64}, 
                                    Mom::Vector{Int64}, Ncomp::Int64, 
                                    propensity::Function, totalPropensity::Number)
    rv = rand()*totalPropensity
    r_indices[1] = 0
    val = 0.0
    while val < rv
        r_indices[1] += 1
        xc = [ 
                n[:, r_indices[1]], # Layout of n is [species, cell]
             ]
        val += propensity(xc, Mom)
    end
end
"""
fast_sample_generic_unary(r_indices::Vector{Int64}, n::Matrix{Int64}, 
                            propensity::Vector, totalPropensity::Number)

Provides a generic fast sampling function.
    propensity is a vector of pre-computed propensities for each compartment.
"""
function fast_sample_generic_unary!(r_indices::Vector{Int64}, n::Matrix{Int64}, 
                                    Mom::Vector{Int64}, Ncomp::Int64, 
                                    propensity::Vector, totalPropensity::Number)
    rv = rand()*totalPropensity
    r_indices[1] = 1
    val = 1.0*propensity[ r_indices[1] ] # Layout of n is [species, cell]
    while val < rv
        r_indices[1] += 1
        val += propensity[r_indices[1]]
    end
end

"""
fast_sample_generic_binary(r_indices::Vector{Int64}, n::Matrix{Int64}, 
                            propensity::Function, totalPropensity::Number)

Provides a generic fast sampling function.
    propensity is a function that must take the vector of chem content of a single compartment.
"""
function fast_sample_generic_binary!(r_indices::Vector{Int64}, n::Matrix{Int64}, 
                                    Mom::Vector{Int64}, Ncomp::Int64, 
                                    propensity::Function, totalPropensity::Number)
    N = Ncomp
    rv = rand()*totalPropensity
    val = 0.0
    for i=1:N
        for j=i+1:N # NOTE: This requires the propensity to be symmetric!
            r_indices[1] = i
            r_indices[2] = j
            xc = [ 
                    n[:, i], 
                    n[:, j] 
                 ]
            cur = propensity(xc, Mom)
            val += cur # Layout of n is [species, cell]
            if val >= rv
                # Sanity check that we are exiting on correct reactants
                # @assert cur > 0 "BUG: Sampled reactants with 0 propensity!"
                return
                # break
            end
        end
        if val >= rv
            return
            # break
        end
    end
    r_indices[1] = 0
    r_indices[2] = 0
    @show rv totalPropensity propensity #debug
end
"""
fast_sample_generic_binary(r_indices::Vector{Int64}, n::Matrix{Int64}, 
                            propensity::Matrix, totalPropensity::Number)

Provides a generic fast sampling function.
    propensity is a matrix of pre-computed propensities for each compartment couple.
"""
function fast_sample_generic_binary!(r_indices::Vector{Int64}, n::Matrix{Int64}, 
                                    Mom::Vector{Int64}, Ncomp::Int64, 
                                    propensity::Matrix, totalPropensity::Number)
    N = Ncomp
    # @show (N, totalPropensity) #debug
    rv = rand()*totalPropensity
    val = 0.0
    for i=1:N
        for j=i+1:N
            r_indices[1] = i
            r_indices[2] = j
            # r_indices_sampling = r_indices
            # @show r_indices_sampling #debug
            val += propensity[i, j]
            if val >= rv
                # Sanity check that we are exiting on correct reactants
                # @assert propensity[i, j] > 0 "BUG: Sampled reactants with 0 propensity!"
                return
                # break
            end
        end
        if val >= rv
            # Sanity check that we are exiting on correct reactants
            # @assert propensity[r_indices...] > 0 "BUG: Sampled reactants with 0 propensity!"
            return
            # break
        end
    end
    r_indices[1] = 0
    r_indices[2] = 0
    @show rv totalPropensity size(propensity) sum(propensity) val #debug
    spy(propensity) #debug
end

#eof
