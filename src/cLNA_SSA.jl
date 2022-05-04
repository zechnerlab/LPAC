#=
Code from the sAIC project.
Lorenzo Duso
=#

#=
    SAIC_SSA.jl
    This module contains all the SSA machinery to simulate SAIC systems.
=#

"""
    SSA_perturbations(S, n0, durations, changes, Nsamples=100;
                             timestep=0.1, seed=nothing)

Run `Nsamples` SSA trajectories with parameter `changes` (perturbations) at the
given `durations` timepoints.
"""
@export function SSA_perturbations(S::Vector{System}, 
                                    n0::Matrix{Int64},
                                    durations::Vector{Float64}, 
                                    Nsamples::Int64; 
                                    timestep::Float64=durations[1]/100, 
                                    seed::Union{Nothing,Int64}=nothing,
                                    exportRawOutput::Bool=false)
    time_all,mm,t_setpoint,val_setpoint = SSA_perturbations(deepcopy(S),
                                                                n0, durations,
                                                                timestep=timestep,
                                                                asserting=true)
    # Prepare futures array and spawn the processes
    F = []
    for i=1:Nsamples
        fut = @spawnat :any SSA_perturbations(deepcopy(S),
                                                        n0, durations, 
                                                        timestep=timestep,
                                                        asserting=false)
        push!(F, fut)
    end
    # Now prepare to collect the results
    # MM=zeros(length(keys(S.MomDict)),length(time_all))
    # MM2=zeros(size(MM))
    # Take care of the raw output
    n = nothing
    if exportRawOutput
        # The n layout is (Species, Cells)(Trajectories)
        n = [ zeros(Int64, size(n0)) for i=1:Nsamples ]
    end
    #
    function fetchMoments(i::Integer, n; verbose::Bool=false)
        res = fetch(F[i])
        if typeof(res)==RemoteException
            throw(res)
        end
        tt, Mi, ts, vs, n_local = res
        verbose && println("Fetch: SSA simulation # $i")
        if !isnothing(n)
            n[i] = deepcopy(n_local)
        end
        @assert all(Mi .>= 0)
        return Mi
        # MM .+= Mi
        # MM2 .+= Mi.^2
        # [Mi, Mi.^2] # This must be the last line of the loop block for the reduction to work!
    end
    function varm(itr, mean::Array)
        # @assert length(itr) > 1 "Cannot compute variance on a single value!"
        v = zeros(size(mean))
        if length(itr)==1
            return v
        end
        for x in itr
            @. v += (x - mean)^2
        end
        v ./= length(itr)-1
        return v
    end
    momentsGenerator = ( fetchMoments(i, n) for i=1:Nsamples )
    MMav = mean(momentsGenerator)
    # @time begin
    begin
        MMvar = varm(momentsGenerator, MMav)
        momentsGenerator2 = ( fetchMoments(i, nothing; verbose=false).^2 for i=1:Nsamples )
        MM2 = sum(momentsGenerator2)
    end
    # @assert all(MM .>= 0)
    # @assert all(MM2 .>= 0)
    # MMav = MM/Nsamples
    # MMvar = ( MM2 - (MM.^2)/Nsamples ) / (Nsamples-1)
    @assert all(MMav .>= 0) "Moments cannot be negative!"
    @assert all(MMvar .>= 0) "Variances cannot be negative!"
    return time_all, MMav, MMvar, t_setpoint, val_setpoint, n, MM2/(Nsamples)
end

"""
    SSA_perturbations(S, n0, durations, changes;
                             timestep=0.1, seed=nothing, asserting=true)

Run one SSA trajectory with parameter `changes` (perturbations) at the given
`durations` timepoints.
"""
@export function SSA_perturbations(S::Vector{System}, 
                                    n0::Matrix{Int64}, 
                                    durations::Vector{Float64};
                                    timestep::Float64=durations[1]/100, 
                                    seed::Union{Nothing,Int64}=nothing, 
                                    asserting::Bool=true)
    seed!=nothing ? Random.seed!(seed) : nothing
    asserting ? assert_model(S[1],n0) : nothing
    @assert length(durations)==length(S) "Invalid setpoint input"
    LL=length(S)
    n_start = deepcopy(n0)
    n_out   = deepcopy(n0) # We want to be able to export the raw output (to compute PDFs)
    t_start=0.
    time_all=Vector{Float64}()
    moments_all=zeros(Float64,length(compute_moments(S[1],n0)),0) # Avoid Int overflows
    t_setpoint=Vector{Float64}()
    val_setpoint=Vector{Float64}()
    for i=1:LL
        s = S[i]
        timepoints=collect(t_start:timestep:t_start+durations[i])
        push!(t_setpoint,t_start); 
        push!(t_setpoint,timepoints[end]); 
        n_out, mom_out = SSA(s, n_start, timepoints, full_story=false, asserting=false)
        time_all = [time_all;timepoints]
        moments_all = [moments_all mom_out]
        t_start += durations[i]
        n_start = deepcopy(n_out)
    end
    return time_all, moments_all, t_setpoint, val_setpoint, n_out
end

"""
    SSA_perturbations(S, n0, durations, changes, Nsamples=100;
                             timestep=0.1, seed=nothing)

Run `Nsamples` SSA trajectories with parameter `changes` (perturbations) at the
given `durations` timepoints.
"""
@export function SSA_perturbations(S::System, 
                                    n0::Matrix{Int64},
                                    durations::Vector{Float64}, 
                                    changes::Vector{Float64},
                                    Nsamples::Int64; 
                                    timestep::Float64=durations[1]/100, 
                                    seed::Union{Nothing,Int64}=nothing,
                                    exportRawOutput::Bool=false)
    time_all,mm,t_setpoint,val_setpoint = SSA_perturbations(deepcopy(S),
                                                                n0, durations, changes,
                                                                timestep=timestep,
                                                                asserting=true)
    # Prepare futures array and spawn the processes
    F = []
    for i=1:Nsamples
        fut = @spawnat :any SSA_perturbations(deepcopy(S),
                                                        n0, durations, changes,
                                                        timestep=timestep,
                                                        asserting=false)
        push!(F, fut)
    end
    # Now prepare to collect the results
    # MM=zeros(length(keys(S.MomDict)),length(time_all))
    # MM2=zeros(size(MM))
    # Take care of the raw output
    n = nothing
    if exportRawOutput
        # The n layout is (Species, Cells)(Trajectories)
        n = [ zeros(Int64, size(n0)) for i=1:Nsamples ]
    end
    #
    function fetchMoments(i::Integer, n; verbose::Bool=false)
        res = fetch(F[i])
        if typeof(res)==RemoteException
            throw(res)
        end
        tt, Mi, ts, vs, n_local = res
        verbose && println("Fetch: SSA simulation # $i")
        if !isnothing(n)
            n[i] = deepcopy(n_local)
        end
        @assert all(Mi .>= 0)
        return Mi
        # MM .+= Mi
        # MM2 .+= Mi.^2
        # [Mi, Mi.^2] # This must be the last line of the loop block for the reduction to work!
    end
    function varm(itr, mean::Array)
        # @assert length(itr) > 1 "Cannot compute variance on a single value!"
        v = zeros(size(mean))
        if length(itr)==1
            return v
        end
        for x in itr
            @. v += (x - mean)^2
        end
        v ./= length(itr)-1
        return v
    end
    momentsGenerator = ( fetchMoments(i, n) for i=1:Nsamples )
    MMav = mean(momentsGenerator)
    # @time begin
    begin
        MMvar = varm(momentsGenerator, MMav)
        momentsGenerator2 = ( fetchMoments(i, nothing; verbose=false).^2 for i=1:Nsamples )
        MM2 = sum(momentsGenerator2)
    end
    # @assert all(MM .>= 0)
    # @assert all(MM2 .>= 0)
    # MMav = MM/Nsamples
    # MMvar = ( MM2 - (MM.^2)/Nsamples ) / (Nsamples-1)
    @assert all(MMav .>= 0) "Moments cannot be negative!"
    @assert all(MMvar .>= 0) "Variances cannot be negative!"
    return time_all, MMav, MMvar, t_setpoint, val_setpoint, n, MM2/(Nsamples)
end

"""
    SSA_perturbations(S, n0, durations, changes;
                             timestep=0.1, seed=nothing, asserting=true)

Run one SSA trajectory with parameter `changes` (perturbations) at the given
`durations` timepoints.
"""
@export function SSA_perturbations(S::System, 
                                    n0::Matrix{Int64}, 
                                    durations::Vector{Float64}, 
                                    changes::Vector{Float64};
                                    timestep::Float64=durations[1]/100, 
                                    seed::Union{Nothing,Int64}=nothing, 
                                    asserting::Bool=true)
    seed!=nothing ? Random.seed!(seed) : nothing
    asserting ? assert_model(S,n0) : nothing
    @assert length(durations)==length(changes) "Invalid setpoint input"
    LL=length(changes)
    n_start = deepcopy(n0)
    n_out   = deepcopy(n0) # We want to be able to export the raw output (to compute PDFs)
    t_start=0.
    time_all=Vector{Float64}()
    # moments_all=zeros(Int64,length(compute_moments(S,n0)),0)
    moments_all=zeros(Float64,length(compute_moments(S,n0)),0) # Avoid Int overflows
    t_setpoint=Vector{Float64}()
    val_setpoint=Vector{Float64}()
    for i=1:LL
        S.transition_classes[1].k *= changes[i]
        timepoints=collect(t_start:timestep:t_start+durations[i])
        push!(t_setpoint,t_start); 
        # push!(val_setpoint,S.transition_classes[1].k/S.transition_classes[2].k)
        push!(t_setpoint,timepoints[end]); 
        # push!(val_setpoint,S.transition_classes[1].k/S.transition_classes[2].k)
        n_out, mom_out = SSA(S, n_start, timepoints, full_story=false, asserting=false)
        time_all = [time_all;timepoints]
        moments_all = [moments_all mom_out]
        t_start += durations[i]
        n_start = deepcopy(n_out)
    end
    return time_all, moments_all, t_setpoint, val_setpoint, n_out
end


# returns the average and variance of moments of system S with initial condition n0 at times timepoints over Nsamples stochastic realizations
@export function SSA(S::System, n0::Matrix{Int64}, timepoints::Vector{T}, Nsamples::Int64) where T <: Real
    assert_model(S,n0)
    MM=zeros(length(keys(S.MomDict)),length(timepoints))
    MM2=copy(MM)
    MM, MM2 = @distributed (+) for i=1:Nsamples
        println("SSA simulation # $i")
        n,Mi=SSA(S,n0,timepoints,asserting=false)
        # MM .+= Mi
        # MM2 .+= Mi.^2
        @assert all(Mi .>= 0)
        @assert all(Mi.^2 .>= 0)
        [Mi, Mi.^2]
    end
    MMav=MM/Nsamples
    MMvar=MM2/(Nsamples-1)-(MM.^2)/(Nsamples*(Nsamples-1))
    @assert all(MMav .> 0) 
    @assert all(MMvar .> 0) 
    return MMav, MMvar
end



# performs a stochastic simulation of system S with initial condition n0
## if timepoints is a vector, it returns moments at each provided time
## if timepoints is a single value, it returns moments at each event up to final time
## if keyword full_story = true, it returns population state at each timepoint or event,
## else only the final population state is returned
@export function SSA(S::System, n0::Matrix{Int64}, timepoints::Union{Vector{Float64},Float64};
                maxiters::Int64=100000, full_story::Bool=false,
                seed::Union{Nothing,Int64}=nothing,asserting::Bool=true)
    seed!=nothing ? Random.seed!(seed) : nothing
    asserting ? assert_model(S,n0) : nothing
    DD,Ncomp=size(n0)
    CL = length(S.transition_classes)
    n = copy(n0) # n is [species, compartments]
    Mom = compute_moments(S,n)
    g, G = compute_propensities(S, n, Mom) # g is vector of tensors, G is vector of numbers
    rates = [S.transition_classes[i].k for i=1:CL]
    NcompM = [Ncomp]
    if length(timepoints) > 1
        path_flag=false
        TT=length(timepoints)
        MM=zeros(Int64,length(Mom),TT)
        MM[:,1]=Mom
        simtime=timepoints[1]
        full_story ? ( n_story=[zeros(Int64,DD,0) for i=1:TT]; n_story[1]=deepcopy(n0) ) : nothing
    else
        path_flag=true
        reaction_times=zeros(maxiters)
        MMt=zeros(Int64,length(Mom),maxiters)
        MMt[:,1]=Mom
        simtime=0.
        full_story ? ( n_story=[zeros(Int64,DD,0) for i=1:maxiters]; n_story[1]=deepcopy(n0) ) : nothing
    end
    r_indices=zeros(Int64,2)
    xc=[zeros(Int64,DD),zeros(Int64,DD)]
    yc=[zeros(Int64,DD),zeros(Int64,DD)]
    # H_weights=zeros(Int64,CL)
    H_weights=zeros(CL)
    H_classes=zeros(CL)
    time_index=2
    simulation_flag=true
    class_counter=zeros(Int64,length(rates))
    while simulation_flag
        #Evaluate global propensity
        for c=1:CL
            curH = S.transition_classes[c].H
            if !isnothing(curH)
                H_weights[c] = curH(n,Mom)
            else
                H_weights[c] = G[c]
            end
            H_classes[c] = rates[c]*float(H_weights[c])
        end
        Htot = sum(H_classes)
        #Compute time of next reaction
        Htot>0.0 ? simtime -= log(1-rand())/Htot : simtime=Inf  #if total propensity is zero, just end
        if !path_flag
            while timepoints[time_index]<simtime
                MM[:,time_index]=Mom
                full_story ? n_story[time_index]=n[:,1:Mom[1]] : nothing
                time_index+=1
                time_index>TT ? (simulation_flag=false;break) : nothing
            end
        end
        (path_flag && (simtime > timepoints)) ? simulation_flag=false : nothing
        if simulation_flag
            #Compute type and detail of next reaction
            next_class = climbtower(rand()*Htot, H_classes)
            class_counter[next_class]+=1
            if S.transition_classes[next_class].rc > 0
                draw_reactant_indices!(
                    S.transition_classes[next_class].fast_sample_reactants!,
                    r_indices,
                    H_weights[next_class],
                    xc,
                    S,
                    next_class,
                    n,
                    Mom,
                    g[next_class])
            end
            if S.transition_classes[next_class].pc > 0
                draw_product_compartments!(S.transition_classes[next_class].parameters,yc,xc,S,next_class)#,n,Mom)
            end
            # If we filled the whole n matrix...
            if Mom[1] == NcompM[1]
                # ...and if we are about to increase the number of compartments...
                if S.transition_classes[next_class].DeltaN > 0
                    # ...then grow (double) all the data structures that need to be enlarged
                    n = [n zeros(Int64,size(n))]
                    NcompM[1]=size(n,2)
                    grow_propensity_tensors!(g, NcompM[1])
                end
            end
            update_all!(S, next_class, r_indices, xc, yc, n, Mom, g, G)
            if path_flag
                reaction_times[time_index]=simtime
                MMt[:,time_index]=Mom
                full_story ? n_story[time_index]=n[:,1:Mom[1]] : nothing
                time_index+=1
                time_index>maxiters ? ( println("OVERFLOW! Increase maxiters"); simulation_flag=false; break ) : nothing
            end
        end
    end
    #println(class_counter)
    full_story ? nothing : n_story=n[:,1:Mom[1]]
    if path_flag
        time_index-=1
        if full_story
            return  n_story[1:time_index], reaction_times[1:time_index], MMt[:,1:time_index]
        else
            return  n[:,1:Mom[1]], reaction_times[1:time_index], MMt[:,1:time_index]
        end
    else
        if full_story
            return  n_story, MM
        else
            return  n[:,1:Mom[1]], MM
        end
    end
end ## END SSA



#function to compute index of tower sampling from a vector
function climbtower(rr::Float64,vect::Vector{T}) where T <: Real
    i=1
    cumul=vect[1]
    while rr>cumul
        i+=1
        cumul+=vect[i]
    end
    return i
end


# if necessary, expands the state matrix in order to accomodate an increase of compartment number
function check_state_boundary!(n::Matrix{Int64},Ncomp::Int64,NcompM::Vector{Int64})
    if Ncomp == NcompM[1]
        n = [n zeros(Int64,size(n))]
        NcompM[1]=2*Ncomp
    end
end

# if the fast_sample_reactants! function is provided, samples the next reacting compartments efficiently
function draw_reactant_indices!(fast_sample!::Function,
                                r_indices::Vector{Int64},
                                propensity_weight::Number,
                                xc::Vector{Vector{Int64}},
                                S::System,
                                next_class::Int64,
                                n::Matrix{Int64},
                                Mom::Vector{Int64},
                                g::Union{Nothing, Vector, Matrix},
                                )
    r_indices[2] = 0  # necessary to prevent troubles with two compartments ...
    fast_sample!(r_indices, n, Mom)
    @assert r_indices[1]>0 "FATAL: Couldn't sample reactants ($r_indices) for TC=$next_class, total propensity ($propensity_weight) must be wrong!"

    TC = S.transition_classes[next_class]
    for i=1:TC.rc
        for d=1:S.n_species 
            xc[i][d] = n[d, r_indices[i]] 
        end
    end
    # if next_class==5
    #     @assert issorted(r_indices) "r_indices not sorted in class 5"
    #     if ( g[r_indices...] != TC.g(xc, Mom) ) || xc[1][1]<=0 || xc[2][2]<=0
    #         @show next_class r_indices g[r_indices...] TC.g(xc, Mom) xc
    #     end
    #     @assert g[r_indices...] == TC.g(xc, Mom) "Inconsistent propensity in class 5"
    #     @assert xc[1][1]>0
    #     @assert xc[2][2]>0
    # elseif next_class==6
    #     @assert issorted(r_indices) "r_indices not sorted in class 6"
    #     if ( g[r_indices...] != TC.g(xc, Mom) ) || xc[1][2]<=0 || xc[2][1]<=0
    #         @show next_class r_indices g[r_indices...] TC.g(xc, Mom) xc
    #     end
    #     @assert g[r_indices...] == TC.g(xc, Mom) "Inconsistent propensity in class 6"
    #     @assert xc[1][2]>0
    #     @assert xc[2][1]>0
    # end
end

# if the fast_sample_reactants! function is not provided, use the default ones!
function draw_reactant_indices!(fast_sample!::Nothing,
                                r_indices::Vector{Int64},
                                propensity_weight::Number,
                                xc::Vector{Vector{Int64}},
                                S::System,
                                next_class::Int64,
                                n::Matrix{Int64},
                                Mom::Vector{Int64},
                                g::Nothing, # Here we are not using the pre-computed propensities
                                )
    Ncomp = Mom[1]
    DD=size(n,1)
    TC = S.transition_classes[next_class]
    if TC.rc == 1
        fast_sample! = 
            (r_indices::Vector{Int64}, n::Matrix{Int64}, Mom::Vector{Int64}) -> 
            fast_sample_generic_unary!(r_indices, n, Mom, Ncomp, 
                                        TC.g, 
                                        propensity_weight
                                        )
    elseif TC.rc == 2
        fast_sample! = 
            (r_indices::Vector{Int64}, n::Matrix{Int64}, Mom::Vector{Int64}) -> 
            fast_sample_generic_binary!(r_indices, n, Mom, Ncomp, 
                                        TC.g, 
                                        propensity_weight
                                        )
    end
    draw_reactant_indices!(
                            fast_sample!, 
                            r_indices,
                            propensity_weight,
                            xc,
                            S,
                            next_class,
                            n,
                            Mom,
                            g,
                            )
end
function draw_reactant_indices!(fast_sample!::Nothing,
                                r_indices::Vector{Int64},
                                propensity_weight::Number,
                                xc::Vector{Vector{Int64}},
                                S::System,
                                next_class::Int64,
                                n::Matrix{Int64},
                                Mom::Vector{Int64},
                                g::Vector,
                                )
    Ncomp = Mom[1]
    DD=size(n,1)
    TC = S.transition_classes[next_class]
    @assert TC.rc == 1 "ERROR: Number of reactants must be 1"
    fast_sample! = 
        (r_indices::Vector{Int64}, n::Matrix{Int64}, Mom::Vector{Int64}) -> 
        fast_sample_generic_unary!(r_indices, n, Mom, Ncomp, 
                                    g, 
                                    propensity_weight
                                    )
    draw_reactant_indices!(
                            fast_sample!, 
                            r_indices,
                            propensity_weight,
                            xc,
                            S,
                            next_class,
                            n,
                            Mom,
                            g,
                            )
end
function draw_reactant_indices!(fast_sample!::Nothing,
                                r_indices::Vector{Int64},
                                propensity_weight::Number,
                                xc::Vector{Vector{Int64}},
                                S::System,
                                next_class::Int64,
                                n::Matrix{Int64},
                                Mom::Vector{Int64},
                                g::Matrix,
                                )
    Ncomp = Mom[1]
    DD=size(n,1)
    TC = S.transition_classes[next_class]
    @assert TC.rc == 2 "ERROR: Number of reactants must be 2"
    fast_sample! = 
        (r_indices::Vector{Int64}, n::Matrix{Int64}, Mom::Vector{Int64}) -> 
        fast_sample_generic_binary!(r_indices, n, Mom, Ncomp, 
                                    g, 
                                    propensity_weight
                                    )
    draw_reactant_indices!(
                            fast_sample!, 
                            r_indices,
                            propensity_weight,
                            xc,
                            S,
                            next_class,
                            n,
                            Mom,
                            g,
                            )
end

# sample product compartments
function draw_product_compartments!(param::Nothing,yc::Vector{Vector{Int64}},xc::Vector{Vector{Int64}},
                                    S::System,next_class::Int64)
    S.transition_classes[next_class].pi(yc,xc)
end

# sample product compartments, when outcome distribution pi needs some additional parameters
function draw_product_compartments!(param,yc::Vector{Vector{Int64}},xc::Vector{Vector{Int64}},
                                    S::System,next_class::Int64)
    S.transition_classes[next_class].pi(yc,xc,param)
end

# efficient implementation of x^exponent for low indices
function recursive_exponentiation(x::Int64,exponent::Int64)
    if exponent==0
        return 1
    else
        return x*recursive_exponentiation(x,exponent-1)
    end
end

"""
    grow_propensity_tensors!(g::PropensitySet, newN::Int64)

Increase the size of the propensity tensors. To be used when the number of
compartments grows over the current tensor size.
"""
function grow_propensity_tensors!(g::PropensitySet, newN::Int64)
    for (i, gcl) in enumerate(g)
        if isnothing(gcl)
            continue
        end
        valueType = typeof(gcl[1])
        numDimensions = length(size(gcl))
        curNumCompartments = size(gcl, 1)
        @assert newN >= curNumCompartments """
            ERROR: This function can only grow the propensity tensors, not shrink them!
            """
        curGcl = deepcopy(gcl)
        dimensions = repeat([newN], numDimensions)
        g[i] = zeros(valueType, dimensions...)
        curSizeRange = repeat([1:curNumCompartments], numDimensions)
        g[i][curSizeRange...] .= curGcl
    end
end

# updates population state n and moments vector Mom with the drawn class and the reactant and product compartments
function update_all!(S::System, next_class::Int64,
                    r_indices::Vector{Int64},
                    xc::Vector{Vector{Int64}}, yc::Vector{Vector{Int64}},
                    n::Matrix{Int64}, Mom::Vector{Int64}, 
                    g::PropensitySet, G::Vector)
    # @show next_class #debug
    NcompPre=Mom[1]
    MomPre = deepcopy(Mom)
    transitionClass = S.transition_classes[next_class]
    numReactants = transitionClass.rc
    # Here we update the moments, the propensities and the state
    update_moments!(Mom, S, NcompPre, next_class, xc, yc)
    NcompPost=Mom[1]
    newCompartments, changedCompartments, swaps = update_state!(n, S, NcompPre, next_class, r_indices, yc)
    update_propensities!(g, G, S, NcompPre, NcompPost, numReactants, r_indices, n, MomPre, 
        newCompartments, changedCompartments, swaps)
    # _checkMomentsConsistency(Mom, n) #debug
end

function update_moments!(
                            Mom::Vector{Int64}, 
                            S::System, 
                            Ncomp::Int64,
                            transitionClassIndex::Int64, 
                            xc::Vector{Vector{Int64}}, 
                            yc::Vector{Vector{Int64}},
                            )
    transitionClass = S.transition_classes[transitionClassIndex]
    # First we remove the reactants contributions...
    for r=1:transitionClass.rc
        for l=1:length(Mom)
            val=1
            for d=1:S.n_species
                val *= recursive_exponentiation(xc[r][d], S.MomDict[l][d]) # xc[r][d]^S.MomDict[l][d]
            end
            Mom[l] -= val
        end
    end
    # ...then we add the products contributions
    for p=1:transitionClass.pc
        for l=1:length(Mom)
            val=1
            for d=1:S.n_species
                pv = yc[p][d]
                # @assert pv>=0 """
                # FATAL: Species has gone negative in product compartment!
                # Reaction: $(transitionClassIndex)
                # Prod Compartment: $p
                # Species: $d
                # x[$d] = $(pv)
                # """
                val *= recursive_exponentiation(yc[p][d], S.MomDict[l][d])  # yc[p][d]^S.MomDict[l][d]
            end
            Mom[l] += val
        end
    end
end

function update_propensities!(
            g::PropensitySet, G::Vector,
            S::System,
            NcompPre::Int64,
            NcompPost::Int64,
            numReactantsChange::Int64,
            r_indices::Vector{Int64},
            n::Matrix{Int64},
            Mom::Vector{Int64},
            newCompartments::Vector{Int64},
            changedCompartments::Vector{Int64},
            swaps::Vector{Pair{Int64,Int64}},
            )
    # For each transition class:
    # (1) First we want to remove the contributions of all compartments that were
    # either removed or had their content changed.
    # (2) Then we want to update the propensity tensor 'g[transitionClassIndex]'
    # by swapping the compartments that changed their index and updating the 
    # propensities of those that changed their content.
    # (3) Finally we want to add the contributions of all compartments that had their
    # content changed and those that were added.
    removedCompartments = getindex.(swaps, 2)
    remI = vcat(changedCompartments, removedCompartments) # This order ensures that
    addI = vcat(changedCompartments, newCompartments) # both remI,addI are sorted
    if !isempty(removedCompartments) && !(NcompPre in removedCompartments)
        push!(remI, NcompPre)
        push!(addI, removedCompartments...)
    end
    for (transitionClassIndex, transitionClass) in enumerate(S.transition_classes)
        numReactants = transitionClass.rc
        if ( isnothing(g[transitionClassIndex]) 
            || numReactants==0
            )
            continue
        end
        @assert numReactants<=2 "ERROR: Transition classes using 3 or more reactants are not supported!"
        # Step (1): removal of old contributions
        # @show remI #debug
        G[transitionClassIndex] -= _computeCompPropensityContribution(
                                        g[transitionClassIndex],
                                        NcompPre,
                                        remI...,
                                        )
        # Step (2): update of the propensity tensor
        # @show swaps #debug
        # _swapCompartmentIndexInPropensityTensor!.(
        #                                 Ref(g[transitionClassIndex]), 
        #                                 NcompPre,
        #                                 swaps,
        #                                 )
        # @show changedCompartments #debug
        _updatePropensityTensor!.(
                                Ref(g[transitionClassIndex]),
                                Ref(n),
                                Ref(Mom),
                                transitionClass.g,
                                NcompPost,
                                addI,
                                )
        # Step (3): addition of new contributions
        # @show addI #debug
        G[transitionClassIndex] += _computeCompPropensityContribution(
                                        g[transitionClassIndex],
                                        NcompPost,
                                        addI...
                                        )
        # # DEBUG: Check for consistency
        # if NcompPost>1 && (transitionClassIndex in [5,6])
        #     trueG = sum( [ transitionClass.g([n[:,i], n[:,j]], Mom) for i=1:NcompPost for j=i+1:NcompPost ] )
        #     sumg = sum(g[transitionClassIndex][1:NcompPost, 1:NcompPost])
        #     @assert equal(sumg, G[transitionClassIndex]) """
        #         FATAL: G is not consistent with the propensity tensor g!
        #         G[$transitionClassIndex] = $(G[transitionClassIndex])
        #         sum(g[$transitionClassIndex]) = $(sumg)
        #         trueG = $trueG
        #         remI = $remI
        #         addI = $addI
        #         swaps = $swaps
        #         """
        #     @assert equal(sumg, trueG) """
        #         FATAL: the propensity tensor g is not consistent with the true G!
        #         sum(g[$transitionClassIndex]) = $(sumg)
        #         trueG = $trueG
        #         G[$transitionClassIndex] = $(G[transitionClassIndex])
        #         NcompPre = $NcompPre
        #         NcompPost = $NcompPost
        #         remI = $remI
        #         addI = $addI
        #         swaps = $swaps
        #         """
        # end
    end
end

function update_state!(
                        n::Matrix{Int64},
                        S::System, 
                        Ncomp::Int64,
                        transitionClassIndex::Int64,
                        r_indices::Vector{Int64},
                        yc::Vector{Vector{Int64}},
                        )
    transitionClass = S.transition_classes[transitionClassIndex]
    numReactants = transitionClass.rc
    swaps = Vector{Pair{Int64, Int64}}() # This is just cells being moved but unchanged
    changedCompartments = Vector{Int64}() # This is new cells with changed content
    newCompartments = Vector{Int64}()

    # Now manage the {in,de}crease in the number of cells, updating the state
    # matrix accordingly + update the contents according to the reaction
    if transitionClass.DeltaN == -1
        pos_overwrite = max(r_indices[1], r_indices[2])
        n[:, pos_overwrite] .= n[:, Ncomp]
        push!(swaps, Ncomp => pos_overwrite)

        if numReactants == 2
            n[:, r_indices[1]] .= yc[1][:]
            push!(changedCompartments, r_indices[1])
        end

    elseif transitionClass.DeltaN == 0
        for i=1:numReactants
            n[:, r_indices[i]] .= yc[i][:]
            push!(changedCompartments, r_indices[i])
        end

    elseif transitionClass.DeltaN == 1
        n[:,Ncomp+1] .= yc[1][:]
        push!(newCompartments, Ncomp+1)
        if numReactants == 1
            n[:, r_indices[1]] .= yc[2][:]
            push!(changedCompartments, r_indices[1])
        end
    else
        error()
    end
    return newCompartments, changedCompartments, swaps
end

# In case we pass an empty collection, their contribution is 0
function _computeCompPropensityContribution(g::Vector, Ncomp::Int64)
    return 0
end
function _computeCompPropensityContribution(g::Matrix, Ncomp::Int64)
    return 0
end
function _computeCompPropensityContribution(g::Nothing, Ncomp::Int64, I...)
    return 0
end
function _computeCompPropensityContribution(g::Vector, Ncomp::Int64, i::Int64)
    return g[i]
end
function _computeCompPropensityContribution(g::Matrix, Ncomp::Int64, i::Int64)
    return sum( g[1:i-1, i] ) + sum( g[i, i+1:Ncomp] )
end
function _computeCompPropensityContribution(g::Vector, Ncomp::Int64, I...)
    return sum( _computeCompPropensityContribution.(Ref(g), Ncomp, I) )
end
function _computeCompPropensityContribution(g::Matrix, Ncomp::Int64, I...)
    # Here we need to sum all the contributions of each compartment (i.e. upper col-row)
    individualContributions = sum( _computeCompPropensityContribution.(Ref(g), Ncomp, I) )
    # And then to sum all the crossings, where the same values are summed twice, once
    # for each compartment: i.e. g[i,j] is summed both as contribution of 'i' and of 'j'
    doubledContributions = sum( 
                                broadcast( 
                                    x->getindex(g, x...),
                                    combinations(I, 2),
                                    ) 
                                )
    return individualContributions - doubledContributions
end

"""
    _swapCompartmentIndexInPropensityTensor!(g, Ncomp::Int64, i::Int64, j::Int64)

Swap the compartment propensities from position i --> j in the propensity tensor.
"""
function _swapCompartmentIndexInPropensityTensor!(g::Nothing, NcompPre::Int64, i::Int64, j::Int64)
    nothing
end
function _swapCompartmentIndexInPropensityTensor!(g::Vector, NcompPre::Int64, i::Int64, j::Int64)
    g[j] = g[i]
end
function _swapCompartmentIndexInPropensityTensor!(g::Matrix, NcompPre::Int64, i::Int64, j::Int64)
    if i == j
        return
    end
    @assert i > j "ERROR: Only swaps from higher to lower indices are supported! ($i -> $j)"
    g[1:j-1, j] .= g[1:j-1, i]
    g[j, j+1:i-1] .= g[j+1:i-1, i]
    g[j, i:NcompPre-1] .= g[i, i+1:NcompPre]
end
"""
    _swapCompartmentIndexInPropensityTensor!(g, NcompPre::Int64, p::Pair{Int64, Int64})

Swap the compartment propensities from position p.first --> p.second in the propensity tensor.
"""
function _swapCompartmentIndexInPropensityTensor!(
                g::Union{Nothing, Vector, Matrix}, 
                NcompPre::Int64, 
                p::Pair{Int64, Int64},
                )
    _swapCompartmentIndexInPropensityTensor!(g, NcompPre, p.first, p.second)
end

# function _computeSwapLostContribution(g::Union{Nothing, Vector}, P::Pair{Int64, Int64}...)
#     return 0
# end
# function _computeSwapLostContribution(g::Matrix)
#     return 0 # No swap, nothing lost
# end
# function _computeSwapLostContribution(g::Matrix, p::Pair{Int64, Int64})
#     if p.first == p.second
#         return 0
#     end
#     @assert p.second < p.first "FATAL: This swap is ill-posed"
#     return g[p.second, p.first]
# end

function _updatePropensityTensor!(
                g::Nothing, 
                n::Matrix{Int64}, 
                Mom::Vector{Int64},
                propensity::Function, 
                Ncomp::Int64, 
                i::Int64,
                )
    nothing
end
function _updatePropensityTensor!(
                g::Vector, 
                n::Matrix{Int64}, 
                Mom::Vector{Int64},
                propensity::Function, 
                Ncomp::Int64, 
                i::Int64,
                )
    xc = [ n[:,i] ]
    g[i] = propensity(xc, Mom)
end
function _updatePropensityTensor!(
                g::Matrix, 
                n::Matrix{Int64}, 
                Mom::Vector{Int64},
                propensity::Function, 
                Ncomp::Int64, 
                i::Int64,
                )
    # NOTE: propensity MUST be symmetric
    # check = [ propensity([ n[:,i], n[:,j] ], Mom) == propensity([ n[:,j], n[:,i] ], Mom) for j=1:Ncomp ]
    # @assert all(check) "Propensity is not symmetric"

    # propWrapper(j::Int64) = propensity([ n[:,i], n[:,j] ], Mom)
    g[1:i-1, i] .= ( propensity([ n[:,j], n[:,i] ], Mom) for j=1:i-1 )
    g[i, i+1:Ncomp] .= ( propensity([ n[:,i], n[:,j] ], Mom) for j=i+1:Ncomp )
end

function _getTotalPropensity(N::Number, n::Matrix{Int64}, interaction::Function)
    val = 0.0
    for i=1:N
        val += interaction(n[:, i])
    end
    return val
end

function _checkMomentsConsistency(Mom::Vector{Int64}, n::Matrix{Int64})
    # n is [species, compartments]
    N = Mom[1]
    M1 = Mom[2]
    M2 = Mom[3]

    @assert N <= size(n, 2) """
    FATAL: Moment inconsistency detected!
    N is bigger than the state matrix!
    N: $N
    size(n,1): $(size(n, 1))
    """

    M1e = _getTotalPropensity(N, n, x->x[1])
    M2e = _getTotalPropensity(N, n, x->x[1]^2)
    @assert M1 == M1e """
        FATAL: Moment inconsistency detected!
        M1: $M1
        M1e: $M1e
        """
    @assert M2 == M2e """
        FATAL: Moment inconsistency detected!
        M2: $M2
        M2e: $M2e
        """
end

#eof

