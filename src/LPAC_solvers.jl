#=
Solvers: some higher level functions to solve either the differential equations or SSA and produce Solution objects.
=#

function convertSolution(sol::SciMLBase.AbstractODESolution, 
						 MMap::Dict{Symbol,Int}, 
						 σMap::Dict{Symbol,Tuple{Symbol,Symbol}},
						 T::AbstractVector{N} where N<:Number,
						 )::Solution
	u = deepcopy(collect( vcat( [ sol(t)' for t in T ]... ) ))
	M = MomentDict()
	σ = SigmaDict()
	for k in keys(MMap)
		i = MMap[k]
		M[k] = @view u[:,i]
	end
	for k in keys(σMap)
		k2, k1 = σMap[k]
		m2, m = M[k2], M[k1]
		# σ[k] = sqrt.(m2 .- m.^2)
		Var = m2 .- m.^2
		if any(Var .< 0)
			println("::: WARNING: Negative variance estimate for $(k)")
		end
		σ[k] = sqrt.(abs.(Var)) #debug
	end
	σ[:X0] = nothing
	return Solution(T, M, σ)
end
function convertSolution(sol::SciMLBase.AbstractODESolution, 
						 MMap::Dict{Symbol,Int}, 
						 σMap::Dict{Symbol,Tuple{Symbol,Symbol}},
						 )::Solution
	T = sol.t
	# u = deepcopy(collect( hcat(sol.u...)' ))
	return convertSolution(sol, MMap, σMap, T)
end

@export function outOfDomainCheck(u, p=nothing, t=nothing)
	return any(u.<0)
end

@export function outOfDomainCheckVar(u, p=nothing, t=nothing;
						 			MMap::Dict{Symbol,Int}=Dict{Symbol,Int}(), 
									σMap::Dict{Symbol,Tuple{Symbol,Symbol}}=Dict{Symbol,Tuple{Symbol,Symbol}}(),
									)
	for k in keys(σMap)
		k2, k1 = σMap[k]
		Mk1 = u[MMap[k1]]
		Mk2 = u[MMap[k2]]
		m2, m = Mk2, Mk1
		Var = m2 - m^2
		if Var < 0
			return true
		end
	end
	return any(u.<0)
end

@export function LPACsolve(M::Model, u0::Vector{T} where T <: Number;
							T=10.0,
						 	MMap::Dict{Symbol,Int}=Dict{Symbol,Int}(), 
							σMap::Dict{Symbol,Tuple{Symbol,Symbol}}=Dict{Symbol,Tuple{Symbol,Symbol}}(),
							solverFlags...
							)::SciMLBase.AbstractODESolution
	u0Map = [ s=>u0[i] for (s,i) in M.momentsMapping ] #debug
	tspan = (0.0, T)
	p = M.parameters
	prob = ODEProblem(M.momentsOde, u0, tspan, p)
	# solver = Tsit5() # Debug
	solver = AutoTsit5(Rosenbrock23()) # General purpose
	# solver = AutoTsit5(RadauIIA5()) # Implicit
	# @time sol = solve(prob, solver;
	sol = solve(prob, solver;
				isoutofdomain=outOfDomainCheck,
				reltol=1e-6,
				solverFlags...,
				)
	return sol
end
@export function LPACsolve(M::Model;
							T=10.0,
							N0=1., Mpc0=1.,
							# Mpc0bScale=1.,
						 	MMap::Dict{Symbol,Int}=Dict{Symbol,Int}(), 
							σMap::Dict{Symbol,Tuple{Symbol,Symbol}}=Dict{Symbol,Tuple{Symbol,Symbol}}(),
							solverFlags...
							)::SciMLBase.AbstractODESolution
	u0 = M.momentsInit(N0, Mpc0)
	LPACsolve(M, u0;
				T=T, MMap=MMap, σMap=σMap, solverFlags...)
end

@export function SSA_solve(M::Model; T::Number, 
								N0::Number, Mpc0::Number,
								NSSA::Number=100,
								)::Tuple{Solution, Vector}
	S = M.ssaSystem
	nSpecies = S.n_species
	n0 = zeros(Int64, nSpecies, N0) # Rows=species, Cols=cells
	# Get the initial values from the initConditions and set it for the SSA
	Mom0 = M.momentsInit(N0, Mpc0)
	MomMap = M.momentsMapping
	for i=1:nSpecies
		s = :M*("⁰"^(i-1))*"¹"*("⁰"^(nSpecies-i))
		n0[i,:] .= round(Int, Mom0[MomMap[s]] / N0) # Set the init value of each cell to the avg
	end

	# @show n0 #debug

	durations = [float(T)]
	changes = [1.0]
	# tTraj, momTraj, _, _, _ = Sim.SSA_perturbations(deepcopy(S), n0, durations, changes)
	@time t, Moms, Vars, _, _, n, MM2 = Sim.SSA_perturbations(
		deepcopy(S),
		n0,
		durations,
		changes,
		NSSA;
		exportRawOutput=true,
		)

	Moments = MomentDict([ s=>Moms[i,:] for (i,s) in S.MomMapping ]...)
	Sigmas = SigmaDict([ s=>sqrt.(Vars[i,:]) for (i,s) in S.MomMapping ]...)
	Squares = MomentDict([ s=>MM2[i,:] for (i,s) in S.MomMapping ]...)
	Sigmas[:X0] = nothing
	ssaSol = Solution(
				t,
				Moments,
				Sigmas,
				)
	
	return ssaSol, n
end
@export function SSA_solve(M::Vector{Model}; T::Vector{numType} where numType<:Number, 
								N0::Number, Mpc0::Number,
								NSSA::Number=100,
								)::Tuple{Solution, Vector}
	S = [ m.ssaSystem for m in M ]
	nSpecies = S[1].n_species
	n0 = zeros(Int64, nSpecies, N0) # Rows=species, Cols=cells
	# Get the initial values from the initConditions and set it for the SSA
	Mom0 = M[1].momentsInit(N0, Mpc0)
	MomMap = M[1].momentsMapping
	for i=1:nSpecies
		s = :M*("⁰"^(i-1))*"¹"*("⁰"^(nSpecies-i))
		n0[i,:] .= round(Int, Mom0[MomMap[s]] / N0) # Set the init value of each cell to the avg
	end

	# @show n0 #debug

	durations = [float(t) for t in T]
	# tTraj, momTraj, _, _, _ = Sim.SSA_perturbations(deepcopy(S), n0, durations, changes)
	@time t, Moms, Vars, _, _, n, MM2 = Sim.SSA_perturbations(
		S,
		n0,
		durations,
		NSSA;
		exportRawOutput=true,
		)

	Moments = MomentDict([ s=>Moms[i,:] for (i,s) in S[1].MomMapping ]...)
	Sigmas = SigmaDict([ s=>sqrt.(Vars[i,:]) for (i,s) in S[1].MomMapping ]...)
	Squares = MomentDict([ s=>MM2[i,:] for (i,s) in S[1].MomMapping ]...)
	Sigmas[:X0] = nothing
	ssaSol = Solution(
				t,
				Moments,
				Sigmas,
				)
	
	return ssaSol, n
end

#eof
