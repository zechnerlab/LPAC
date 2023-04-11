#=
Experimental: experimental utility functions not used for final figures or results.
=#

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
	lpacSol = @time LPACsolve(M; T=T, N0=N0, Mpc0=Mpc0, MMap=M.momentsMapping, σMap=M.sigmaMapping, solverFlags...)
	# @show typeof(lpacSol)

	Msymb = popfirst!(setdiff(symbols, [:N,:N2]))

	ssaSol = nothing
	
	if NSSA<=0
		lpacSol = convertSolution(lpacSol, M.momentsMapping, M.sigmaMapping)
		return compareSolutions(M, lpacSol, lpacSol, symbols...;
								Msymb = Msymb,
								rescaleToConcentrations=rescaleToConcentrations,
								meanlegend=meanlegend,
								varlegend=varlegend,
								verbose=verbose,
								fontscale=fontscale,
								)
	else
		nPool, trajectories, ssaMeanSol, ssaSDSol = runSSASimulations(
				M, symbols...;
				T=T, NSSA=NSSA, RSSA=RSSA, N0=N0, Mpc0=Mpc0, lpacSol=lpacSol)
		lpacSol = convertSolution(lpacSol, M.momentsMapping, M.sigmaMapping)
		return compareSolutions(M, lpacSol, ssaMeanSol, symbols...;
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
	lpacSol1 = @time LPACsolve(M1; T=T, N0=N0, Mpc0=Mpc0, MMap=M1.momentsMapping, σMap=M1.sigmaMapping, solverFlags...)
	lpacSol2 = @time LPACsolve(M2; T=T, N0=N0, Mpc0=Mpc0, MMap=M2.momentsMapping, σMap=M2.sigmaMapping, solverFlags...)

	Msymb = popfirst!(setdiff(symbols, [:N,:N2]))

	lpacSol1 = convertSolution(lpacSol1, M1.momentsMapping, M1.sigmaMapping)
	lpacSol2 = convertSolution(lpacSol2, M2.momentsMapping, M2.sigmaMapping)
	return compareSolutions(M1, lpacSol1, lpacSol2, symbols...;
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
	CorrLPAC = []
	Param = []
	for M in VM
		p = M.parameters[sweepParam]
		push!(Param, p)
		lpacSol = LPACsolve(M; T=T, N0=N0, Mpc0=Mpc0, MMap=M.momentsMapping, σMap=M.sigmaMapping, solverFlags...)
		lpacSol = convertSolution(lpacSol, M.momentsMapping, M.sigmaMapping)
		ssaSol, nRaw = SSA_solve(M; T=T, N0=N0, Mpc0=Mpc0, NSSA=NSSA)
		nRaw = hcat(nRaw...)
		# Now extract the correlation endpoint
		ssaC = correlation(ssaSol, M11, M10, M01, M20, M02)
		lpacC = correlation(lpacSol, M11, M10, M01, M20, M02)
		push!(CorrSSA, ssaC[end])
		push!(CorrLPAC, lpacC[end])
	end
	plot(Param, CorrSSA; color=pal[1], label="Corr SSA")
	plot!(Param, CorrLPAC; color=pal[1], label="Corr LPAC", linestyle=:dash, linewidth=2)
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
	CorrLPAC = similar(VM, Any)
	Param1 = similar(VM, Any)
	Param2 = similar(VM, Any)
	for (i,M) in enumerate(VM)
		Param1[i] = M.parameters[sweepParam1]
		Param2[i] = M.parameters[sweepParam2]
		lpacSol = LPACsolve(M; T=T, N0=N0, Mpc0=Mpc0, MMap=M.momentsMapping, σMap=M.sigmaMapping, solverFlags...)
		lpacSol = convertSolution(lpacSol, M.momentsMapping, M.sigmaMapping)
		# Now extract the correlation endpoint
		lpacC = correlation(lpacSol, M11, M10, M01, M20, M02)
		CorrLPAC[i] = lpacC[end]
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
		p2 = plot(Param1, Param2, CorrLPAC'; st=:surface, title="Corr LPAC")
		return plot(p1,p2; layout=(1,2))
	else
		p1 = plot(Param1, Param2, CorrLPAC'; 
						st=:contourf, title="Corr LPAC",
						fillcolor=:curl,
						xlabel=string(sweepParam1), ylabel=string(sweepParam2),
						zlabel="Correlation",
						)
		p1 = plot!(Param1, Param2, CorrLPAC';
						st=:contour,
						levels=[-0.5,0.0,0.5],
						linecolor=[:gray,:black,:gray],
						linewidth=1.5,
						linestyle=:dash, 
						)
		return p1
	end
end

#eof
