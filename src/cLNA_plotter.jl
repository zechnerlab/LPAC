#=
Plotter: some plotting-related functions.
=#

const pal = palette(:tab10)

function plotMoment(T, M; σ=nothing, legend=:topleft, ylims=(0.0,Inf), kwargs...)
	p = plot(T, M; 
		ribbon=σ, 
		fillalpha=0.2,
		ylims=ylims, 
		legend=legend,
		# legend=:bottomright,
		# color=pal[1],
		# right_margin=15mm,
		# label=L"N", 
		kwargs...
		)
	xlabel!("Time [a.u.]")
	# ylabel!("Copy number")
end
function plotMoment!(p, T, M; σ=nothing, legend=:topleft, ylims=(0.0,Inf), kwargs...)
	p = plot!(p, T, M; 
		ribbon=σ, 
		fillalpha=0.2,
		ylims=ylims, 
		legend=legend,
		# legend=:bottomright,
		# color=pal[1],
		# right_margin=15mm,
		# label=L"N", 
		kwargs...
		)
	# xlabel!("Time [a.u.]")
	# ylabel!("Copy number")
end
function _plotMoments(sol::Solution, symbols::Symbol...)
	@unpack T, M, σ = sol

	S = collect(symbols)
	s1 = popfirst!(S)
	p = plotMoment(T, M[s1]; σ=safeGet(σ, s1), label=latexstring(s1), 
					color=pal[1], right_margin=15mm, yformatter=:plain)
	i = 2
	for s in S
		p = plotMoment!(p, T, M[s]; σ=safeGet(σ, s), label=latexstring(s), 
					color=pal[1], right_margin=15mm, yformatter=:plain)
		i += 1
	end
	if :X0 in keys(M) && :X0 in keys(σ)
		p = plotMoment!(p, T, M[:X0];
			σ = σ[:X0], 
			label=L"x_0 = \frac{\langle M^1 \rangle}{\langle N \rangle}", 
			color=pal[3],
			linestyle=:dot,
			)
	end
	return p, i
end
function _plotMoments!(p::Union{Plots.Plot, Plots.Subplot}, i::Int, sol::Solution, symbols::Symbol...)
	@unpack T, M, σ = sol

	S = collect(symbols)
	s1 = popfirst!(S)
	p = plotMoment!(p, T, M[s1]; σ=safeGet(σ, s1), label=latexstring(s1), 
						color=pal[i], 
						legend=:topright, 
						ylabel="Total mass", 
						yformatter=:plain)
	i += 1
	for s in S
		p = plotMoment!(p, T, M[s]; σ=safeGet(σ, s), label=latexstring(s), 
						color=pal[i], 
						legend=:topright, 
						ylabel="Total mass", 
						yformatter=:plain)
		i += 1
	end
	return p, i
end

@export function plotMoments(sol::Solution, symbols::Symbol...;
								leftSymbols::Vector{Symbol}=[:N,], 
								# plotX0::Bool=false,
								# plotX0::Bool=true,
								# σX0=nothing
								)
	pMoments, i = _plotMoments(sol, leftSymbols...)
	xlabel!("Time [a.u.]")
	ylabel!("Copy number")
	pMomentsY2 = twinx(pMoments)
	_plotMoments!(pMomentsY2, i, sol, symbols...)
	return pMoments
end

@export function plotMeanComparisons(Tssa, Mssa, Tclna, Mclna; 
					label=L"M", 
					color=pal[1], 
					legend=:bottomright,
					ssaRibbon=nothing,
					)
	p = plotMoment(Tssa, Mssa; 
			color=color, 
			label=L"\langle "*label*L"_{ssa} \rangle",
			legend=legend,
			σ=ssaRibbon,
			)
	p = plotMoment!(p, Tclna, Mclna; 
			color=color, 
			linestyle=:dash,
			linewidth=2.0, 
			label=L"\langle "*label*L"_{approx} \rangle",
			legend=legend,
			)
end
@export function plotVarianceComparisons(Tssa, VarSsa, Tclna, VarClna; 
					label=L"Var", 
					color=pal[1], 
					legend=:bottomright,
					# legend=:topright,
					ssaRibbon=nothing,
					)
	p = plotMoment(Tssa, VarSsa; 
			color=color, 
			label=L"Var("*label*L"_{ssa})",
			legend=legend,
			σ=ssaRibbon,
			)
	p = plotMoment!(p, Tclna, VarClna; 
			color=color, 
			linestyle=:dash,
			linewidth=2.0, 
			label=L"Var("*label*L"_{approx})",
			legend=legend,
			)
end

_var(M2::Number, M1::Number, N::Number) = M2 - ( M1^2 / N)
_var(sol::Solution, M2::Symbol, M1::Symbol, N::Symbol=:N) = @. sol.M[M2] - ( (sol.M[M1])^2 / sol.M[N] )

function correlation(N::Number, M11::Number, M10::Number, M01::Number, M20::Number, M02::Number)
	num = M11 - (M10 * M01 / N)
	denom = sqrt( _var(M20, M10, N) * _var(M02, M01, N) )
	# denom = sqrt( abs( _var(M20, M10, N) * _var(M02, M01, N) ) ) #debug
	return num ./ denom
end
function correlation(sol::Solution, M11::Symbol, M10::Symbol, M01::Symbol, M20::Symbol, M02::Symbol) 
	num = sol.M[M11] .- (sol.M[M10] .* sol.M[M01] ./ sol.M[:N]) 
	# num = sol.M[M11] .- (sol.M[M10] .* sol.M[M01]) 
	denom = sqrt.( _var(sol, M20, M10) .* _var(sol, M02, M01) )
	# denom = sqrt.( abs.( _var(sol, M20, M10) .* _var(sol, M02, M01) ) ) #debug
	# denom = 1 # Un-rescaled
	return num ./ denom
end
@export function plotCorrelation(ssa::Solution, clna::Solution; 
									M11::Symbol=:M¹¹, M10::Symbol=:M¹⁰, M01::Symbol=:M⁰¹,
									M20::Symbol=:M²⁰, M02::Symbol=:M⁰²,
									color=pal[1],
									# legend=:bottomright,
									legend=:topright,
									# legend=false,
									ssaSD::Union{StandardDeviation,Nothing}=nothing,
									)
	l1 = replace(string(M10), "M"=>"x")
	l2 = replace(string(M01), "M"=>"x")
	ssaCorrelation = correlation(ssa, M11, M10, M01, M20, M02)
	ssaCorrelationRibbon = nothing
	if :corr in keys(ssa.M)
		ssaCorrelation = ssa.M[:corr]
		ssaCorrelationRibbon = !isnothing(ssaSD) ? ssaSD.M[:corr] : nothing
	end
	p = plotMoment(ssa.T, ssaCorrelation;
			ribbon=ssaCorrelationRibbon,
			color=color,
			label=latexstring("Corr("*l1*","*l2*")_{ssa}"),
			ylims=(-1.0,1.0),
			# ylims=(-Inf,Inf),
			legend=legend,
			)
	p = plotMoment!(p, clna.T, correlation(clna, M11, M10, M01, M20, M02);
			color=color,
			linestyle=:dash,
			linewidth=2.0, 
			label=latexstring("Corr("*l1*","*l2*")_{approx}"),
			ylims=(-1.0,1.0),
			# ylims=(-Inf,Inf),
			legend=legend,
			)
end

function computeCV(N, M, σN, σM)
	# Coefficient of variation
	cvN = σN ./ N
	cvM = σM ./ M
	CV = hcat(cvN, cvM)
	replace!(CV, NaN=>0.0, Inf=>0.0)
	return CV
end

function _getPlotUpperLimCV(CV)
	# lim = quantile(CV[:], 0.9)
	lim = 2*mean(CV)
	lim *= 1.1
	lim = round(lim, digits=1)
	return lim
end

@export function plotCV(T, N, M, σN, σM)
	CV = computeCV(N, M, σN, σM)
	lim = _getPlotUpperLimCV(CV)
	pCV = plot(T, CV; 
			ylims=(0.0,lim),
			label=[L"CV(N)" L"CV(M^1)"],
			)
	xlabel!("Time [a.u.]")
	return pCV
end

@export function plotSolution(S::Solution, Ω=1.0, Ωc=1.0; Ms::Symbol=:M¹)
	@unpack T, M, σ = S

	pMoments = plotMoments(S, Ms)
	if Ms in keys(σ)
		pCV = plotCV(T, M[:N], M[Ms], σ[:N], σ[Ms])
		return plot(pMoments, pCV)
	else
		return plot(pMoments)
	end
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

function pushSigmaX0!(σ::SigmaDict, M::MomentDict, s2::Symbol, s1::Symbol)
	σ[:X0] = sqrt.( (M[s2] ./ M[:N]) .- (M[s1] ./ M[:N]).^2 ) # Is this the right formula?
	return σ
end
function pushSigmaX0!(S::Solution, s2::Symbol, s1::Symbol)
	pushSigmaX0!(S.σ, S.M, s2, s1)
end


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

@export function cLNAsolve(M::Model;
							T=10.0,
							N0=1., Mpc0=1.,
							# Mpc0bScale=1.,
						 	MMap::Dict{Symbol,Int}=Dict{Symbol,Int}(), 
							σMap::Dict{Symbol,Tuple{Symbol,Symbol}}=Dict{Symbol,Tuple{Symbol,Symbol}}(),
							solverFlags...
							)::SciMLBase.AbstractODESolution
	u0 = M.momentsInit(N0, Mpc0)
	u0Map = [ s=>u0[i] for (s,i) in M.momentsMapping ] #debug
	# @show u0Map #debug
	tspan = (0.0, T)
	p = M.parameters
	prob = ODEProblem(M.momentsOde, u0, tspan, p)
	# solver = Euler() # Debug
	# solver = Midpoint() # Debug
	solver = Tsit5() # Debug
	# solver = AutoTsit5(Rosenbrock23()) # General purpose
	# solver = AutoTsit5(RadauIIA5()) # Implicit
	# @time sol = solve(prob, solver;
	sol = solve(prob, solver;
				isoutofdomain=outOfDomainCheck,
				# isoutofdomain=(u,p,t)->outOfDomainCheckVar(u,p,t; MMap=MMap, σMap=σMap),
				reltol=1e-6,
				# abstol=1e-2,
				solverFlags...,
				)
	return sol
end

# @export function plotSsaSolution(S::SSAsolution)
# 	pMoments = plotMoments(S.T, S.N, S.M, S.σN, S.σM)
# 	pCV = plotCV(S.T, S.N, S.M, S.σN, S.σM)
# 	return plot(pMoments, pCV; layout=(1,2))
# end

@export function SSA_solve(M::Model; T::Number, 
								N0::Number, Mpc0::Number,
								NSSA::Number=100,
								fillalpha::AbstractFloat=0.4,
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

@export function getErrorMeasure(ssaSol::Solution, clnaSol::SciMLBase.AbstractODESolution,
									MMap::Dict{Symbol,Int}, 
						 			σMap::Dict{Symbol,Tuple{Symbol,Symbol}},
									symbols::Symbol...)
	return getErrorMeasure(
				ssaSol,
				convertSolution(clnaSol, MMap, σMap, ssaSol.T),
				symbols...
				)
end
@export function getErrorMeasure(ssa::Solution, clna::Solution,
									symbols::Symbol...)
	T = ssa.T
	errors = [ abs.(clna.M[s] - ssa.M[s]) ./ ssa.M[s] for s in symbols ]
	return mean.(errors)
end

@export function plotMeans(ssa::Solution, clna::Solution, symbols::Symbol...; 
							legend=:bottomright,
							ssaSD::Union{StandardDeviation,Nothing}=nothing,
							)
	sanitizedSymbols = [ s for s in symbols if s in keys(clna.M) && s in keys(ssa.M) ]
	plots = [ 
		plotMeanComparisons(ssa.T, ssa.M[s], clna.T, clna.M[s]; 
							label=latexstring(s), 
							color=pal[i], 
							legend= s==:N ? :bottomright : legend,
							ssaRibbon=isnothing(ssaSD) ? nothing : ssaSD.M[s],
							)
		for (i,s) in enumerate(sanitizedSymbols)
		]
	plot(plots...; layout=(1,length(plots)) )
end
@export function plotVariances(ssa::Solution, clna::Solution, symbols::Symbol...; 
							# legend=:topright,
							legend=:bottomright,
							ssaSD::Union{StandardDeviation,Nothing}=nothing,
							)
	sanitizedSymbols = [ s for s in symbols if s in keys(clna.σ) && s in keys(ssa.σ) ]
	plots = [ 
		plotVarianceComparisons(ssa.T, ssa.σ[s].^2, clna.T, clna.σ[s].^2; 
							label=latexstring(s), color=pal[i],
							legend= s==:N ? :bottomright : legend,
							ssaRibbon=isnothing(ssaSD) ? nothing : ssaSD.σ[s],
							)
		for (i,s) in enumerate(sanitizedSymbols)
		]
	if :M¹¹ in symbols
		push!(plots, 
			plotCorrelation(ssa, clna; M11=:M¹¹, M10=:M¹⁰, M01=:M⁰¹, 
				color=pal[length(sanitizedSymbols)+1],
				ssaSD=ssaSD,
				),
			)
	end
	plot(plots...; layout=(1,length(plots)) )
end

function _plotPopulationHistograms(ssaPool)
	pHistograms = nothing
	if !isnothing(ssaPool)
		histograms=[]
		species = size(ssaPool,1)
		for i=1:species
			h = histogram(ssaPool[i,:]; title="Species $i", legend=false)
			push!(histograms, h)
		end
		pHistograms = nothing
		if species>=2
			h2 = histogram2d(ssaPool[1,:], ssaPool[species,:]; 
						xlabel="Species 1",
						ylabel="Species $species",
						)
			pHistograms = plot(h2, histograms...; layout=(1,species+1))
		else
			pHistograms = plot(histograms...; layout=(1,species))
		end
	end
	return pHistograms
end

function _pow(M::Matrix, e::Vector)
	# @show size(M,2), length(e) #debug
	@assert size(M,2) == length(e)
	newRows = [ (M[i,:] .^ e)' for i=1:size(M,1) ]
	return vcat(newRows...)
end
function _tpow(M::Matrix, e::Vector)
	# @show size(M,1), length(e) #debug
	@assert size(M,1) == length(e)
	newRows = [ (M[:,i] .^ e)' for i=1:size(M,2) ]
	return vcat(newRows...)
end
function _moment(S::Matrix, e::Vector)
	P = _tpow(S,e)
	M = prod(P, dims=2)
	return sum(M)
end

function _plotMomentsHistograms(M::Model, ssaTrajectories, symbols::Symbol...;
									reference=nothing)
	pHistograms = nothing
	if !isnothing(ssaTrajectories)
		histograms=[]
		species = size(ssaTrajectories[1],1)
		Moments = Dict{Symbol,Vector}()
		for (i,s) in enumerate(symbols)
			e = M.ssaSystem.MomExpMapping[s]
			curMom = [ _moment(S, e) for S in ssaTrajectories ]
			Moments[s] = curMom

			h = histogram(Moments[s]; title=string(s), color=pal[i], legend=false)
			if !isnothing(reference)
				ref = reference.M[s][end] # Pick last timepoint
				h = vline!([ref]; color=:black, linewidth=3, legend=false)
			end
			push!(histograms, h)
		end
		pHistograms = nothing
		# if species==2
		# 	M10, M01 = :M¹⁰, :M⁰¹
		# 	h2 = histogram2d(Moments[M10], Moments[M01]; 
		# 				xlabel=string(M10),
		# 				ylabel=string(M01),
		# 				)
		# 	pHistograms = plot(h2, histograms...; layout=(1,1+length(histograms)))
		# else
		# 	pHistograms = plot(histograms...; layout=(1,length(histograms)))
		# end
		pHistograms = plot(histograms...; layout=(1,length(histograms)))
	end
	return pHistograms
end

@export function compareSolutions(
							M::Model, sol::Solution, reference::Solution, symbols::Symbol...;
							Msymb::Symbol = popfirst!(setdiff(symbols, [:N,:N2])),
							rescaleToConcentrations::Bool = false,
							meanlegend=:bottomright,
							varlegend=:bottomright,
							verbose::Bool=true,
							fontscale=1.6,
							ssaSD=nothing,
							ssaPool=nothing,
							ssaTrajectories=nothing,
							)
	if rescaleToConcentrations
		rescaleSolution!(sol, M.Ω, M.Ωc)
		rescaleSolution!(reference, M.Ω, M.Ωc)
	end
	if :M² in keys(M.momentsMapping) && :M¹ in keys(M.momentsMapping)
		pushSigmaX0!(reference, :M², :M¹)
		pushSigmaX0!(sol, :M², :M¹)
	end
	scalefontsizes(fontscale)
	p = nothing
	try
		pCLNA = plotSolution(sol; Ms=Msymb)
		pSSA = plotSolution(reference; Ms=Msymb)
		pMean = plotMeans(reference, sol, symbols...; 
							legend=meanlegend, ssaSD=ssaSD)
		pVar = plotVariances(reference, sol, symbols...; 
							legend=varlegend, ssaSD=ssaSD)

		pHistograms = _plotPopulationHistograms(ssaPool)
		pMHistograms = _plotMomentsHistograms(M, ssaTrajectories, symbols...; reference=sol)

		if !verbose
			p = plot(pMean, pVar; layout=(2,1), size=(800,400))
			if !isnothing(pHistograms)
				# p = plot(pMean, pVar, pHistograms; layout=(3,1), size=(800,600))
				p = plot(pMean, pVar, pHistograms, pMHistograms; layout=(4,1), size=(1000,600))
			end
			# p = plot(pMean, pVar; layout=(2,1))
		elseif sol !=  reference
			p = plot(pSSA, pCLNA, pMean, pVar; layout=(4,1))
		else
			p = plot(pCLNA, pMean, pVar; layout=(3,1))
		end
	finally
		scalefontsizes(1/fontscale)
	end
	return p
end

function _computeSSAStats(
		ssaSols::Vector{Solution})
	T = ssaSols[1].T # TODO: Check that all times agree!
	symbols = keys(ssaSols[1].M)
	_Mvalues = Dict{Symbol, Matrix}()
	_σvalues = Dict{Symbol, Matrix}()
	M = MomentDict()
	σ = SigmaDict()
	Msd = MomentDict()
	σsd = SigmaDict()
	for s in symbols
		# @show s #debug
		_Mvalues[s] = hcat([ sol.M[s] for sol in ssaSols ]...) # Runs are different columns
		_σvalues[s] = hcat([ sol.σ[s] for sol in ssaSols ]...) # Time is on the rows
		M[s] = mean(_Mvalues[s], dims=2)[:,1]
		Msd[s] = sqrt.(var(_Mvalues[s], dims=2))[:,1]
		σ[s] = mean(_σvalues[s], dims=2)[:,1]
		σsd[s] = sqrt.(var(_σvalues[s].^2, dims=2))[:,1] # Here squaring necessary 
		# so that it match the variance (i.e. square) plotting
		# @show maximum(M[s]), maximum(Msd[s]) #debug
		# @show maximum(σ[s])^2, maximum(σsd[s]) #debug
	end
	if :M¹¹ in symbols
		# Then compute the correlation
		N = _Mvalues[:N]
		M11 = _Mvalues[:M¹¹]
		M10 = _Mvalues[:M¹⁰]
		M01 = _Mvalues[:M⁰¹]
		M20 = _Mvalues[:M²⁰]
		M02 = _Mvalues[:M⁰²]
		_correlations = correlation.(N, M11, M10, M01, M20, M02)
		M[:corr] = mean(_correlations, dims=2)[:,1]
		Msd[:corr] = sqrt.(var(_correlations, dims=2))[:,1]
	end

	return Solution(T, M, σ), StandardDeviation(T, Msd, σsd)
end


@export function testAll(M::Model, symbols::Symbol...; 
						T=100.0, NSSA=100, RSSA=10,
						N0=10, Mpc0=10,
						rescaleToConcentrations::Bool=false,
						meanlegend=:bottomright,
						varlegend=:bottomright,
						verbose::Bool=true,
						fontscale=1.6,
						solverFlags...
						)
	# TODO/BEWARE: Rescaling to concentrations does not work well with the ODE solution
	# interpolator: this means that when computing the relative error, completely
	# wrong values are returned by the interpolation!
	N0 = round(Int, M.Ω * N0)
	Mpc0 = round(Int, M.Ωc * Mpc0)
	clnaSol = @time cLNAsolve(M; T=T, N0=N0, Mpc0=Mpc0, MMap=M.momentsMapping, σMap=M.sigmaMapping, solverFlags...)
	# @show typeof(clnaSol)

	Msymb = popfirst!(setdiff(symbols, [:N,:N2]))

	ssaSol = nothing
	
	if NSSA<=0
		clnaSol = convertSolution(clnaSol, M.momentsMapping, M.sigmaMapping)
		# if rescaleToConcentrations
		# 	rescaleSolution!(clnaSol, M.Ω, M.Ωc)
		# end
		# if :M² in keys(M.momentsMapping) && :M¹ in keys(M.momentsMapping)
		# 	pushSigmaX0!(clnaSol, :M², :M¹)
		# end
		# pCLNA = plotSolution(clnaSol; Ms=Msymb)
		# return plot(pCLNA)
		return compareSolutions(M, clnaSol, clnaSol, symbols...;
								Msymb = Msymb,
								rescaleToConcentrations=rescaleToConcentrations,
								meanlegend=meanlegend,
								varlegend=varlegend,
								verbose=verbose,
								fontscale=fontscale,
								)
	else
		ssaSols = Vector{Solution}()
		cells = []
		trajectories = []
		for ir=1:RSSA
			ssaSol, nRaw = SSA_solve(M; T=T, N0=N0, Mpc0=Mpc0, NSSA=NSSA)
			push!(ssaSols, ssaSol)
			Err = getErrorMeasure(ssaSol, 
						clnaSol, M.momentsMapping, M.sigmaMapping,
						symbols...)
			@show Err
			push!(cells, hcat(nRaw...))
			push!(trajectories, nRaw...)
		end
		nPool = hcat(cells...)
		clnaSol = convertSolution(clnaSol, M.momentsMapping, M.sigmaMapping)
		ssaMeanSol, ssaSDSol = _computeSSAStats(ssaSols)
		return compareSolutions(M, clnaSol, ssaMeanSol, symbols...;
								Msymb = Msymb,
								rescaleToConcentrations=rescaleToConcentrations,
								meanlegend=meanlegend,
								varlegend=varlegend,
								verbose=verbose,
								ssaSD=ssaSDSol,
								ssaPool=nPool,
								ssaTrajectories=trajectories,
								)
	end
end

@export function compareModels(M1::Model, M2::Model, symbols::Symbol...; 
						T=100.0,
						N0=10, Mpc0=10,
						rescaleToConcentrations::Bool=false,
						meanlegend=:bottomright,
						varlegend=:bottomright,
						verbose::Bool=true,
						fontscale=1.6,
						solverFlags...
						)
	# TODO/BEWARE: Rescaling to concentrations does not work well with the ODE solution
	# interpolator: this means that when computing the relative error, completely
	# wrong values are returned by the interpolation!
	N0 = round(Int, M1.Ω * N0)
	Mpc0 = round(Int, M1.Ωc * Mpc0)
	clnaSol1 = @time cLNAsolve(M1; T=T, N0=N0, Mpc0=Mpc0, MMap=M1.momentsMapping, σMap=M1.sigmaMapping, solverFlags...)
	clnaSol2 = @time cLNAsolve(M2; T=T, N0=N0, Mpc0=Mpc0, MMap=M2.momentsMapping, σMap=M2.sigmaMapping, solverFlags...)

	Msymb = popfirst!(setdiff(symbols, [:N,:N2]))

	clnaSol1 = convertSolution(clnaSol1, M1.momentsMapping, M1.sigmaMapping)
	clnaSol2 = convertSolution(clnaSol2, M2.momentsMapping, M2.sigmaMapping)
	return compareSolutions(M1, clnaSol1, clnaSol2, symbols...;
							Msymb = Msymb,
							rescaleToConcentrations=rescaleToConcentrations,
							meanlegend=meanlegend,
							varlegend=varlegend,
							verbose=verbose,
							fontscale=fontscale,
							)
end

@export function correlationCheckParamSweep(VM::Vector{Model}, sweepParam::Symbol;
											T=100.0, NSSA=100,
											N0=10, Mpc0=10,
											meanlegend=:bottomright,
											varlegend=:bottomright,
											M11::Symbol=:M¹¹, M10::Symbol=:M¹⁰, M01::Symbol=:M⁰¹,
											M20::Symbol=:M²⁰, M02::Symbol=:M⁰²,
											solverFlags...
											)
	CorrSSA = []
	CorrCLNA = []
	Param = []
	for M in VM
		p = M.parameters[sweepParam]
		push!(Param, p)
		clnaSol = cLNAsolve(M; T=T, N0=N0, Mpc0=Mpc0, MMap=M.momentsMapping, σMap=M.sigmaMapping, solverFlags...)
		clnaSol = convertSolution(clnaSol, M.momentsMapping, M.sigmaMapping)
		ssaSol, nRaw = SSA_solve(M; T=T, N0=N0, Mpc0=Mpc0, NSSA=NSSA)
		nRaw = hcat(nRaw...)
		# Now extract the correlation endpoint
		ssaC = correlation(ssaSol, M11, M10, M01, M20, M02)
		clnaC = correlation(clnaSol, M11, M10, M01, M20, M02)
		push!(CorrSSA, ssaC[end])
		push!(CorrCLNA, clnaC[end])
	end
	plot(Param, CorrSSA; color=pal[1], label="Corr SSA")
	plot!(Param, CorrCLNA; color=pal[1], label="Corr cLNA", linestyle=:dash, linewidth=2)
end
@export function correlationCheckParamSweep(VM::Matrix{Model}, sweepParam1::Symbol, sweepParam2::Symbol;
											T=100.0, NSSA=0,
											N0=10, Mpc0=10,
											meanlegend=:bottomright,
											varlegend=:bottomright,
											M11::Symbol=:M¹¹, M10::Symbol=:M¹⁰, M01::Symbol=:M⁰¹,
											M20::Symbol=:M²⁰, M02::Symbol=:M⁰²,
											solverFlags...
											)
	CorrSSA = similar(VM, Any)
	CorrCLNA = similar(VM, Any)
	Param1 = similar(VM, Any)
	Param2 = similar(VM, Any)
	for (i,M) in enumerate(VM)
		Param1[i] = M.parameters[sweepParam1]
		Param2[i] = M.parameters[sweepParam2]
		clnaSol = cLNAsolve(M; T=T, N0=N0, Mpc0=Mpc0, MMap=M.momentsMapping, σMap=M.sigmaMapping, solverFlags...)
		clnaSol = convertSolution(clnaSol, M.momentsMapping, M.sigmaMapping)
		# Now extract the correlation endpoint
		clnaC = correlation(clnaSol, M11, M10, M01, M20, M02)
		CorrCLNA[i] = clnaC[end]
		if NSSA>0
			ssaSol, nRaw = SSA_solve(M; T=T, N0=N0, Mpc0=Mpc0, NSSA=NSSA)
			# nRaw = hcat(nRaw...)
			ssaC = correlation(ssaSol, M11, M10, M01, M20, M02)
			CorrSSA[i] = ssaC[end]
		end		
	end
	Param1 = Param1[:,1]
	Param2 = Param2[1,:]
	if NSSA>0
		p1 = plot(Param1, Param2, CorrSSA'; st=:surface, title="Corr SSA")
		p2 = plot(Param1, Param2, CorrCLNA'; st=:surface, title="Corr cLNA")
		return plot(p1,p2; layout=(1,2))
	else
		p1 = plot(Param1, Param2, CorrCLNA'; 
						st=:contourf, title="Corr cLNA",
						fillcolor=:curl,
						xlabel=string(sweepParam1), ylabel=string(sweepParam2),
						zlabel="Correlation",
						)
		p1 = plot!(Param1, Param2, CorrCLNA';
						st=:contour,
						levels=[-0.5,0.0,0.5],
						linecolor=[:gray,:black,:gray],
						linewidth=1.5,
						linestyle=:dash, 
						)
		return p1
	end
end

@export function testConvergence(VM::Vector{Model}, symbols::Symbol...; 
						T=100.0,
						N0=10, Mpc0=10,
						)
	# TODO/BEWARE: Rescaling to concentrations does not work well with the ODE solution
	# interpolator: this means that when computing the relative error, completely
	# wrong values are returned by the interpolation!

	RawSolutions = Vector{SciMLBase.AbstractODESolution}()
	Solutions = Vector{Solution}()
	for M in VM
		curN0 = round(Int, M.Ω * N0)
		curMpc0 = round(Int, M.Ωc * Mpc0)
		clnaSol = cLNAsolve(M; T=T, N0=curN0, Mpc0=curMpc0)
		push!(RawSolutions, clnaSol)
	end
	Tout = RawSolutions[end].t
	for (M, RS) in zip(VM, RawSolutions)
		println("Converting solution...") #debug
		clnaSol = convertSolution(RS, M.momentsMapping, M.sigmaMapping, Tout)
		rescaleSolution!(clnaSol, M.Ω, M.Ωc)
		push!(Solutions, clnaSol)
	end

	Msymb = popfirst!(setdiff(symbols, [:N,:N2]))
	_getError(x) = getErrorMeasure(Solutions[end], x, symbols...)
	Err = _getError.(Solutions)
	@show Err
	return compareSolutions(
				VM[1],
				Solutions[1],
				Solutions[end],
				symbols...;
				rescaleToConcentrations=false,
				)
end

@export function bifurcationAnalysis(M::Model; N0=10, Mpc0=10)
	# Ref: https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/tutorialsODE/#Neural-mass-equation-(Hopf-aBS)
	BK = BifurcationKit
	# _norm(x) = norm(x, Inf)
	_norm(x) = norm(x)
	f! = M.momentsOde
	function f(u,p)
		du = similar(u)
		f!(du, u, p, 0)
		return du
	end
	jet = BK.getJet(f; matrixfree=false)
	u0 = M.momentsInit(N0, Mpc0, Mpc0*0.8)
	u0 += randn(length(u0)) # debug/test
	p = NamedTuple(M.parameters)
	Mmap = M.momentsMapping

	opts_br = ContinuationPar(pMin = 2.0, pMax = 2e1,
				# parameters to have a smooth result
				# ds = 0.04, dsmax = 0.05,
				# this is to detect bifurcation points precisely with bisection
				detectBifurcation = 3,
				# Optional: bisection options for locating bifurcations
				nInversion = 8, maxBisectionSteps = 25, nev = 3)
	br, = continuation(jet[1], jet[2], u0, p, (@lens _.kMR), opts_br;
				recordFromSolution = (u, p) -> (M¹⁰ = u[Mmap[:M¹⁰]], M¹¹ = u[Mmap[:M¹¹]], N = u[Mmap[:N]]),
				# tangentAlgo = BorderedPred(),
				plot = true, normC = _norm)
	scene = plot(br, plotfold=false, markersize=3, legend=:topleft)
	return scene,br
end

#eof
