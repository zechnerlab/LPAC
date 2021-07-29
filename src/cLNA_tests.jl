# Some test code to start playing with cLNA

const pal = palette(:tab10)

rescaleCompartments(N, Ω, p=1) = N ./ Ω^p
rescaleContent(M, N, Ωc, p=1) = M ./ (N*Ωc).^p
# rescaleContent(M, N, Ωc, p=1) = M ./ (Ωc).^p
rescaleMixed(MN, N, Ω, Ωc, pN=1, pM=1) = rescaleCompartments(
												rescaleContent(MN, N, Ωc, pM), 
												Ω, pN)

function rescaleClnaSol!(sol, Ω, Ωc)
	N = @view sol[1, :]
	N² = @view sol[2, :]
	M = @view sol[3, :]
	M² = @view sol[4, :]
	NM = @view sol[5, :]

	M .= rescaleContent(M, N, Ωc)
	M² .= rescaleContent(M², N, Ωc, 2)
	NM .= rescaleMixed(NM, N, Ω, Ωc)
	N² .= rescaleCompartments(N², Ω, 2)
	N .= rescaleCompartments(N, Ω)

	return sol
end

function rescaleSsaSol!(sol, Ω, Ωc)
	sol.M .= rescaleContent(sol.M, sol.N, Ωc)
	sol.σM .= rescaleContent(sol.σM, sol.N, Ωc)
	sol.σN .= rescaleCompartments(sol.σN, Ω)
	sol.N .= rescaleCompartments(sol.N, Ω)

	return sol
end

function plotMoment(T, M; σ=nothing, legend=:topleft, kwargs...)
	p = plot(T, M; 
		ribbon=σ, 
		fillalpha=0.2,
		ylims=(0.0,Inf), 
		legend=legend,
		# color=pal[1],
		# right_margin=15mm,
		# label=L"N", 
		kwargs...
		)
end
function plotMoment!(p, T, M; σ=nothing, legend=:topleft, kwargs...)
	p = plot!(p, T, M; 
		ribbon=σ, 
		fillalpha=0.2,
		ylims=(0.0,Inf), 
		legend=legend,
		# color=pal[1],
		# right_margin=15mm,
		# label=L"N", 
		kwargs...
		)
end

@export function plotMoments(T, N, M, σN, σM; plotX0::Bool=false)
	pMoments = plotMoment(T, N; 
		σ=σN,  
		label=L"N", 
		color=pal[1],
		right_margin=15mm)
	if plotX0
		pMoments = plotMoment!(pMoments, T, M./N, 
			label=L"x_0 = \frac{\langle M^1 \rangle}{\langle N \rangle}", 
			color=pal[3])
	end
	xlabel!("Time [a.u.]")
	ylabel!("Copy number")
	pMoments = plotMoment!(twinx(), T, M; 
		σ=σM,  
		yformatter=:plain,
		label=L"M^1", 
		legend=:topright,
		color=pal[2],
		ylabel="Total mass")
	return pMoments
end

@export function plotMeanComparisons(Tssa, Mssa, Tclna, Mclna; 
					label=L"M", 
					color=pal[1], 
					legend=:bottomright,
					)
	p = plotMoment(Tssa, Mssa; 
			color=color, 
			label=L"\langle "*label*L"_{ssa} \rangle",
			legend=legend,
			)
	p = plotMoment!(p, Tclna, Mclna; 
			color=color, 
			linestyle=:dash, 
			label=L"\langle "*label*L"_{clna} \rangle",
			legend=legend,
			)
end
@export function plotVarianceComparisons(Tssa, VarSsa, Tclna, VarClna; 
					label=L"Var", 
					color=pal[1], 
					legend=:bottomright,
					)
	p = plotMoment(Tssa, VarSsa; 
			color=color, 
			label=L"Var("*label*L"_{ssa})",
			legend=legend,
			)
	p = plotMoment!(p, Tclna, VarClna; 
			color=color, 
			linestyle=:dash, 
			label=L"Var("*label*L"_{clna})",
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

@export function plotSolution(sol, Ω=1.0, Ωc=1.0)
	T = sol.t
	u = deepcopy(collect( hcat(sol.u...)' ))
	N = @view u[:,1]
	N² = @view u[:,2]
	M1 = @view u[:,3]
	M1² = @view u[:,4]
	# @show M1²
	NM1 = @view u[:,5]
	# σN = sqrt.( N² .- N.^2 ) #prod version
	# σM = sqrt.( M1² .- M1.^2 ) #prod version
	σN = sqrt.( abs.(N² .- N.^2) ) #debug version
	σM = sqrt.( abs.(M1² .- M1.^2) ) #debug version
	pMoments = plotMoments(T, N, M1, σN, σM)
	pCV = plotCV(T, N, M1, σN, σM)
	plot(pMoments, pCV), M1², T
end

@export function outOfDomainCheck(u, p=nothing, t=nothing)
	N, N², M1, M1² = u[1:4]
	failPositivity = (N<0) || (M1<0)
	failVarianceN = N² - N^2 < 0
	failVarianceM1 = M1² - M1^2 < 0
	return failPositivity || failVarianceN || failVarianceM1
end

@export function cLNAsolve(M::Model;
							T=10.0,
							N0=1., Mpc0=1.)
	u0 = M.momentsInit(N0, Mpc0)
	# @show u0 #debug
	tspan = (0.0, T)
	p = M.parameters
	prob = ODEProblem(M.momentsOde, u0, tspan, p)
	solver = AutoTsit5(Rosenbrock23()) # General purpose
	# solver = AutoTsit5(RadauIIA5()) # Implicit
	@time sol = solve(prob, solver;
				isoutofdomain=outOfDomainCheck,
				reltol=1e-6,
				abstol=1e-2,
				)
	return sol
end

@export function plotSsaSolution(S::SSAsolution)
	pMoments = plotMoments(S.T, S.N, S.M, S.σN, S.σM)
	pCV = plotCV(S.T, S.N, S.M, S.σN, S.σM)
	return plot(pMoments, pCV; layout=(1,2))
end

@export function SSA_solve(M::Model; T::Number, 
								N0::Number, Mpc0::Number,
								NSSA::Number=100,
								fillalpha::AbstractFloat=0.4,
								)
	S = M.ssaSystem
	nSpecies = S.n_species
	n0 = zeros(Int64, nSpecies, N0) # Rows=species, Cols=cells
	n0[1,:] .= Mpc0 # Initialize cell content
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

	M12 = MM2[2,:]
	M1 = Moms[2,:]

	N, M1, σN, σM = ( Moms[1,:], Moms[2,:], sqrt.(Vars[1,:]), sqrt.(M12 .- M1.^2) )

	ssaSol = SSAsolution(t, N, σN, M1, σM)
	
	return ssaSol
end

@export function getErrorMeasure(ssaSol, clnaSol)
	T = ssaSol.T
	Nssa, M1ssa = ssaSol.N, ssaSol.M
	Nclna = [ clnaSol(t)[1] for t in T ]
	M1clna = [ clnaSol(t)[3] for t in T ]
	eN = mean( abs.(Nclna - Nssa) ./ Nssa )
	eM = mean( abs.(M1clna - M1ssa) ./ M1ssa )
	return eN, eM
end

@export function plotMeans(ssaSol, clnaSol)
	Tssa = ssaSol.T
	Nssa, M1ssa = ssaSol.N, ssaSol.M
	Tclna = clnaSol.t
	u = deepcopy(collect( hcat(clnaSol.u...)' ))
	Nclna = @view u[:,1]
	M1clna = @view u[:,3]
	pN = plotMeanComparisons(Tssa, Nssa, Tclna, Nclna; label=L"N", color=pal[1])
	pM = plotMeanComparisons(Tssa, M1ssa, Tclna, M1clna; label=L"M^1", color=pal[2])
	plot(pN, pM; layout=(1,2))
end
@export function plotVariances(ssaSol, clnaSol)
	Tssa = ssaSol.T
	σNssa, σMssa = ssaSol.σN.^2, ssaSol.σM.^2
	Tclna = clnaSol.t
	u = deepcopy(collect( hcat(clnaSol.u...)' ))
	Nclna = @view u[:,1]
	N2clna = @view u[:,2]
	M1clna = @view u[:,3]
	M12clna = @view u[:,4]
	σNclna = N2clna .- Nclna.^2
	σMclna = M12clna .- M1clna.^2

	pN = plotVarianceComparisons(Tssa, σNssa, Tclna, σNclna; label=L"N", color=pal[1])
	pM = plotVarianceComparisons(Tssa, σMssa, Tclna, σMclna; label=L"M^1", color=pal[2])
	plot(pN, pM; layout=(1,2))
end

@export function testAll(M::Model; 
						T=100.0, NSSA=100,
						N0=10, Mpc0=10,
						rescaleToConcentrations::Bool=false
						)
	N0 = round(Int, M.Ω * N0)
	Mpc0 = round(Int, M.Ωc * Mpc0)
	clnaSol = cLNAsolve(M; T=T, N0=N0, Mpc0=Mpc0)
	if rescaleToConcentrations
		rescaleClnaSol!(clnaSol, M.Ω, M.Ωc)
	end
	pCLNA, M1², Tclna = plotSolution(clnaSol)
	if NSSA<=0
		return plot(pCLNA)
	else
		ssaSol = SSA_solve(M; T=T, N0=N0, Mpc0=Mpc0, NSSA=NSSA)
		if rescaleToConcentrations
			rescaleSsaSol!(ssaSol, M.Ω, M.Ωc)
		end
		pSSA = plotSsaSolution(ssaSol)
		pMean = plotMeans(ssaSol, clnaSol)
		pVar = plotVariances(ssaSol, clnaSol)
		eN, eM = getErrorMeasure(ssaSol, clnaSol)
		@show eN
		@show eM
		return plot(pSSA, pCLNA, pMean, pVar; layout=(4,1))
	end
end

#eof
