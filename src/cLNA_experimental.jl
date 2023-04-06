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
	clnaSol = @time cLNAsolve(M; T=T, N0=N0, Mpc0=Mpc0, MMap=M.momentsMapping, σMap=M.sigmaMapping, solverFlags...)
	# @show typeof(clnaSol)

	Msymb = popfirst!(setdiff(symbols, [:N,:N2]))

	ssaSol = nothing
	
	if NSSA<=0
		clnaSol = convertSolution(clnaSol, M.momentsMapping, M.sigmaMapping)
		return compareSolutions(M, clnaSol, clnaSol, symbols...;
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
				T=T, NSSA=NSSA, RSSA=RSSA, N0=N0, Mpc0=Mpc0, clnaSol=clnaSol)
		clnaSol = convertSolution(clnaSol, M.momentsMapping, M.sigmaMapping)
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

#eof
