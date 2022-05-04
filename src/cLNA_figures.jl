#=
cLNA Figures: code for defining the figures for the paper.
=#

@export function NestedBirthDeathQuadratic(;
						T=50.0, NSSA=100, RSSA=1,
						N0=10, Mpc0=1,
						fontscale=1.6,
						size=(503,526).*0.75,
						savepath="Figures",
						)
	M = Models.IEqCFBDq_new(; kC=0,kF=0)
	solverFlags = []
	symbols = [:N, :M¹]
	# testAll(M, :N, :M¹; 
	# 			NSSA=NSSA, RSSA=RSSA, 
	# 			T=T, N0=N0, Mpc0=Mpc0, 
	# 			rescaleToConcentrations=false, 
	# 			meanlegend=:bottomright, varlegend=:bottomright, 
	# 			reltol=1e-3, 
	# 			verbose=false)
	
	### Compute the solutions
	N0 = round(Int, M.Ω * N0)
	Mpc0 = round(Int, M.Ωc * Mpc0)
	clnaSolRaw = @time cLNAsolve(M; T=T, N0=N0, Mpc0=Mpc0, MMap=M.momentsMapping, σMap=M.sigmaMapping, solverFlags...)
	Msymb = popfirst!(setdiff(symbols, [:N,:N2]))
	nPool, trajectories, ssaMeanSol, ssaSDSol = runSSASimulations(
				M, symbols...;
				T=T, NSSA=NSSA, RSSA=RSSA, N0=N0, Mpc0=Mpc0, clnaSol=clnaSolRaw)
	clnaSol = cLNA.convertSolution(clnaSolRaw, M.momentsMapping, M.sigmaMapping)
	
	### Plot the figure
	scalefontsizes(fontscale)
	p = nothing
	try
        local lmargin, tmargin, rmargin, bmargin = -8mm, -2mm, 2mm, -4mm
        local tx, ty = -10mm, -2mm

		pmN = plotMeanComparisons(ssaMeanSol.T, ssaMeanSol.M[:N],
								clnaSol.T, clnaSol.M[:N];
								# ssaRibbon=ssaSDSol.M[:N],
								ssaRibbon=ssaMeanSol.σ[:N],
								clnaRibbon=clnaSol.σ[:N],
								color=cLNA.pal[1],
								label=L"N",
								legend=:bottomright,
								title="Number of compartments ("*latexstring("N")*")",
								ylabel="Abundance",
								xlabel=nothing,
								ylims=(0,39.9),
								)
		pmN = plot!(pmN; # Floating label for the panel
            title=L"(a)",
            titlepos=:left,
            grid=false, showaxis=false, ticks=false,
            inset=(1, bbox(tx,ty, 20mm,20mm, :top, :left)),
            subplot=2,
            bg_inside=nothing,
            )
		pmM = plotMeanComparisons(ssaMeanSol.T, ssaMeanSol.M[:M¹],
								clnaSol.T, clnaSol.M[:M¹];
								# ssaRibbon=ssaSDSol.M[:M¹],
								ssaRibbon=ssaMeanSol.σ[:M¹],
								clnaRibbon=clnaSol.σ[:M¹],
								color=cLNA.pal[2],
								label=L"M^1",
								legend=:bottomright,
								title="Number of molecules ("*latexstring("M^1")*")",
								ylabel="Abundance",
								xlabel="Time",
								ylims=(0,399),
								)
		pmM = plot!(pmM; # Floating label for the panel
            title=L"(b)",
            titlepos=:left,
            grid=false, showaxis=false, ticks=false,
            inset=(1, bbox(tx,ty, 20mm,20mm, :top, :left)),
            subplot=2,
            bg_inside=nothing,
            )
		# pMean = plot(pmN, pmM; layout=(1,2))
		# pMean = plotMeans(ssaMeanSol, clnaSol, symbols...; 
		# 					legend=false, ssaSD=ssaSDSol)

		# pvN = plotVarianceComparisons(ssaMeanSol.T, ssaMeanSol.σ[:N].^2,
		# pN = plotVarianceComparisons!(pmN, ssaMeanSol.T, ssaMeanSol.σ[:N].^2,
		# 						clnaSol.T, clnaSol.σ[:N].^2;
		# 						ssaRibbon=ssaSDSol.σ[:N],
		# 						color=cLNA.pal[1],
		# 						legend=false,
		# 						title=nothing,
		# 						ylabel="Variance",
		# 						xlabel="Time [a.u.]",
		# 						)
		# # pvM = plotVarianceComparisons(ssaMeanSol.T, ssaMeanSol.σ[:M¹].^2,
		# pM = plotVarianceComparisons!(pmM, ssaMeanSol.T, ssaMeanSol.σ[:M¹].^2,
		# 						clnaSol.T, clnaSol.σ[:M¹].^2;
		# 						ssaRibbon=ssaSDSol.σ[:M¹],
		# 						color=cLNA.pal[2],
		# 						legend=false,
		# 						title=nothing,
		# 						ylabel=nothing,
		# 						xlabel="Time [a.u.]",
		# 						)
		# pVar = plot(pvN, pvM; layout=(1,2))
		# pVar = plotVariances(ssaMeanSol, clnaSol, symbols...; 
		# 					legend=false, ssaSD=ssaSDSol)

		# p = plot(pMean, pVar; layout=(2,1), size=size)
		# p = plot(pN, pM; layout=(2,1))
		p = plot(pmN, pmM; layout=(2,1), size=size, 
						   left_margin=lmargin, 
						   top_margin=tmargin, 
						   right_margin=rmargin,
						   bottom_margin=bmargin,
						   )
		savefig(savepath*"/nbdq_NM1.pdf")
		savefig(savepath*"/nbdq_NM1.png")
	finally
		scalefontsizes(1/fontscale)
	end
	return p
end

@export function SAIC(;
						T=200.0, NSSA=100, RSSA=1,
						N0=25, Mpc0=0,
						fontscale=1.6,
						size=(503,526).*0.75,
						savepath="Figures",
						)
	Ms = [ Models.SAIC(; kref=1e1), Models.SAIC(; kref=4e1), Models.SAIC(; kref=2e1) ]
	solverFlags = []
	symbols = [:N, :M¹⁰⁰, :M⁰⁰¹, :M⁰¹⁰]
	
	### Compute the solutions
	N0 = round(Int, Ms[1].Ω * N0)
	Mpc0 = round(Int, Ms[1].Ωc * Mpc0)
	u0 = Ms[1].momentsInit(N0, Mpc0)
	solutions = Vector{Solution}()
	@time for M in Ms
		clnaSolRaw = @time cLNAsolve(M, u0; T=T, MMap=M.momentsMapping, σMap=M.sigmaMapping, solverFlags...)
		clnaSol = cLNA.convertSolution(clnaSolRaw, M.momentsMapping, M.sigmaMapping)
		push!(solutions, clnaSol)
		u0 = deepcopy(clnaSolRaw.u[end])
	end
	clnaSol = solcat(solutions...)
	nPool, trajectories, ssaMeanSol, ssaSDSol = runSSASimulations(
					Ms, symbols...;
					T=T, NSSA=NSSA, RSSA=RSSA, N0=N0, Mpc0=Mpc0)

	#TODO: Save the solutions data and allow for reloading them.
	
	### Plot the figure
	scalefontsizes(fontscale)
	p = nothing
	try
        local lmargin, tmargin, rmargin, bmargin = -8mm, -2mm, 2mm, -4mm
        local tx, ty = -10mm, -2mm

		pmN = plotMeanComparisons(ssaMeanSol.T, ssaMeanSol.M[:N],
								clnaSol.T, clnaSol.M[:N];
								ssaRibbon=ssaMeanSol.σ[:N],
								clnaRibbon=clnaSol.σ[:N],
								clnaLabel=false,
								color=cLNA.pal[1],
								label=L"N",
								legend=:bottomright,
								title="Number of compartments ("*latexstring("N")*")",
								ylabel="Abundance",
								xlabel=nothing,
								# ylims=(0,39.9),
								)
		pmN = plot!(pmN; # Floating label for the panel
            title=L"(a)",
            titlepos=:left,
            grid=false, showaxis=false, ticks=false,
            inset=(1, bbox(tx,ty, 20mm,20mm, :top, :left)),
            subplot=2,
            bg_inside=nothing,
            )
		pmM = plotMeanComparisons(ssaMeanSol.T, ssaMeanSol.M[:M¹⁰⁰],
								clnaSol.T, clnaSol.M[:M¹⁰⁰];
								ssaRibbon=ssaMeanSol.σ[:M¹⁰⁰],
								clnaRibbon=clnaSol.σ[:M¹⁰⁰],
								clnaLabel=false,
								color=cLNA.pal[2],
								# label=L"M^{1,0,0}",
								label=L"Z_1",
								legend=:topright,
								title="Chemical content",
								ylabel="Abundance",
								xlabel="Time",
								# ylims=(0,399),
								)
		pmM = plotMeanComparisons!(pmM, ssaMeanSol.T, ssaMeanSol.M[:M⁰¹⁰],
								clnaSol.T, clnaSol.M[:M⁰¹⁰];
								ssaRibbon=ssaMeanSol.σ[:M⁰¹⁰],
								clnaRibbon=clnaSol.σ[:M⁰¹⁰],
								clnaLabel=false,
								color=cLNA.pal[3],
								# label=L"M^{0,1,0}",
								label=L"Z_2",
								legend=:topright,
								# ylims=(0,399),
								)
		pmM = plotMeanComparisons!(pmM, ssaMeanSol.T, ssaMeanSol.M[:M⁰⁰¹],
								clnaSol.T, clnaSol.M[:M⁰⁰¹];
								ssaRibbon=ssaMeanSol.σ[:M⁰⁰¹],
								clnaRibbon=clnaSol.σ[:M⁰⁰¹],
								clnaLabel=false,
								color=cLNA.pal[4],
								# label=L"M^{0,0,1}",
								label=L"Q_2",
								legend=:topright,
								# ylims=(0,399),
								)
		pmM = plot!(pmM; # Floating label for the panel
            title=L"(b)",
            titlepos=:left,
            grid=false, showaxis=false, ticks=false,
            inset=(1, bbox(tx,ty, 20mm,20mm, :top, :left)),
            subplot=2,
            bg_inside=nothing,
            )
		p = plot(pmN, pmM; layout=(2,1), size=size, 
						   left_margin=lmargin, 
						   top_margin=tmargin, 
						   right_margin=rmargin,
						   bottom_margin=bmargin,
						   )
		savefig(savepath*"/SAIC.pdf")
		savefig(savepath*"/SAIC.png")
	finally
		scalefontsizes(1/fontscale)
	end
	return p
end

"""
	MutualRepression()

TODO: 
	- mult=2; kb=10*1e1; kF=100*5*1e-3; N=50*mult; M = Models.MutualRepression(kb=kb,kd=kb*1e-2*0.1, kMR0=30, kMR1=30,kF=kF,kE=(2/(N-1))*kF); f = M.momentsOde; m0 = M.momentsInit(10,10); p = M.parameters; m = similar(m0); @time f(m,m0,p,0.0); testAll(M, :N, :M¹⁰, :M⁰¹, :M¹¹, :M²⁰, :M⁰²; NSSA=8*4, RSSA=4, T=50, N0=25*mult, Mpc0=0, rescaleToConcentrations=false, meanlegend=false, varlegend=false, reltol=1e-8, verbose=false, fontscale=0.5)
	- mult=2; kb=10*1e1; kF=10*5*1e-3; N=50*mult; M = Models.MutualRepression(kb=kb,kd=kb*1e-2*0.1, kMR0=30, kMR1=30,kF=kF,kE=(2/(N-1))*kF); f = M.momentsOde; m0 = M.momentsInit(10,10); p = M.parameters; m = similar(m0); @time f(m,m0,p,0.0); testAll(M, :N, :M¹⁰, :M⁰¹, :M¹¹, :M²⁰, :M⁰²; NSSA=8*4, RSSA=4, T=50, N0=25*mult, Mpc0=0, rescaleToConcentrations=false, meanlegend=false, varlegend=false, reltol=1e-8, verbose=false, fontscale=0.5)
	- mult=2; kb=10*1e1; kF=1*5*1e-3; N=50*mult; M = Models.MutualRepression(kb=kb,kd=kb*1e-2*0.1, kMR0=30, kMR1=30,kF=kF,kE=(2/(N-1))*kF); f = M.momentsOde; m0 = M.momentsInit(10,10); p = M.parameters; m = similar(m0); @time f(m,m0,p,0.0); testAll(M, :N, :M¹⁰, :M⁰¹, :M¹¹, :M²⁰, :M⁰²; NSSA=8*4, RSSA=4, T=50, N0=25*mult, Mpc0=0, rescaleToConcentrations=false, meanlegend=false, varlegend=false, reltol=1e-8, verbose=false, fontscale=0.5)
	- mult=2; kb=10*1e1; kF=0.1*5*1e-3; N=50*mult; M = Models.MutualRepression(kb=kb,kd=kb*1e-2*0.1, kMR0=30, kMR1=30,kF=kF,kE=(2/(N-1))*kF); f = M.momentsOde; m0 = M.momentsInit(10,10); p = M.parameters; m = similar(m0); @time f(m,m0,p,0.0); testAll(M, :N, :M¹⁰, :M⁰¹, :M¹¹, :M²⁰, :M⁰²; NSSA=8*4, RSSA=4, T=50, N0=25*mult, Mpc0=0, rescaleToConcentrations=false, meanlegend=false, varlegend=false, reltol=1e-8, verbose=false, fontscale=0.5)
"""
@export function MutualRepression(;
						T=nothing, NSSA=128, RSSA=1,
						N0=25, Mpc0=0,
						fontscale=1.6,
						size=(3.75*503,1.25*526).*0.75,
						savepath="Figures",
						aspect_ratio=:none,
						speed=:mid, #∈[:slow,:mid,:fast,:faster]
						)
	kb = 100
	S = Dict(:slow=>0.1, :mid=>1, :fast=>10, :faster=>100)
	# kF = 0.1*5e-3
	# kF = 1.0*5e-3
	# kF = 10*5e-3
	# kF = 100*5e-3
	kF = S[speed]*5e-3
	T = isnothing(T) ? 500/S[speed] : T
	N = 50 # Asymptotic level
	M = Models.MutualRepression(; kb=kb, kd=1e-3*kb, 
								  kMR0=30, kMR1=30,
								  kF=kF, kE=(2/(N-1))*kF,
								  )
	solverFlags = []
	symbols = [:N, :M¹⁰, :M²⁰, :M¹¹] #TODO: amend
	
	### Compute the solutions
	N0 = round(Int, M.Ω * N0)
	Mpc0 = round(Int, M.Ωc * Mpc0)
	clnaSolRaw = @time cLNAsolve(M; T=T, N0=N0, Mpc0=Mpc0, MMap=M.momentsMapping, σMap=M.sigmaMapping, solverFlags...)
	Msymb = popfirst!(setdiff(symbols, [:N,:N2]))
	nPool, trajectories, ssaMeanSol, ssaSDSol = runSSASimulations(
				M, symbols...;
				T=T, NSSA=NSSA, RSSA=RSSA, N0=N0, Mpc0=Mpc0, clnaSol=clnaSolRaw)
	clnaSol = cLNA.convertSolution(clnaSolRaw, M.momentsMapping, M.sigmaMapping)
	
	### Plot the figure
	scalefontsizes(fontscale)
	p = nothing
	try
        local lmargin, tmargin, rmargin, bmargin = -4mm, -2mm, -2mm, -4mm
        # local lmargin, tmargin, rmargin, bmargin = 0mm, 0mm, 0mm, 0mm
        local tx, ty = -10mm, -2mm

		pmN = plotMeanComparisons(ssaMeanSol.T, ssaMeanSol.M[:N],
								clnaSol.T, clnaSol.M[:N];
								# ssaRibbon=ssaSDSol.M[:N],
								ssaRibbon=ssaMeanSol.σ[:N],
								clnaRibbon=clnaSol.σ[:N],
								color=cLNA.pal[1],
								label=L"N",
								legend=:bottomright,
								title="Number of compartments ("*latexstring("N")*")",
								ylabel="Abundance",
								xlabel=nothing,
								# ylims=(0,39.9),
								aspect_ratio=aspect_ratio,
								)
		pmN = plot!(pmN; # Floating label for the panel
            title=L"(a)",
            titlepos=:left,
            grid=false, showaxis=false, ticks=false,
            inset=(1, bbox(tx,ty, 20mm,20mm, :top, :left)),
            subplot=2,
            bg_inside=nothing,
            )
		pmM10 = plotMeanComparisons(ssaMeanSol.T, ssaMeanSol.M[:M¹⁰],
								clnaSol.T, clnaSol.M[:M¹⁰];
								# ssaRibbon=ssaSDSol.M[:M¹⁰],
								ssaRibbon=ssaMeanSol.σ[:M¹⁰],
								clnaRibbon=clnaSol.σ[:M¹⁰],
								color=cLNA.pal[2],
								label=L"M^{1,0}",
								legend=:bottomright,
								title="Number of molecules ("*latexstring("M^{1,0}")*")",
								ylabel="Abundance",
								xlabel="Time",
								# ylims=(0,399),
								aspect_ratio=aspect_ratio,
								)
		pmM10 = plot!(pmM10; # Floating label for the panel
            title=L"(b)",
            titlepos=:left,
            grid=false, showaxis=false, ticks=false,
            inset=(1, bbox(tx,ty, 20mm,20mm, :top, :left)),
            subplot=2,
            bg_inside=nothing,
            )
		pmM11 = plotMeanComparisons(ssaMeanSol.T, ssaMeanSol.M[:M¹¹],
								clnaSol.T, clnaSol.M[:M¹¹];
								# ssaRibbon=ssaSDSol.M[:M¹¹],
								ssaRibbon=ssaMeanSol.σ[:M¹¹],
								# clnaRibbon=clnaSol.σ[:M¹¹],
								color=cLNA.pal[3],
								label=L"M^{1,1}",
								legend=:bottomright,
								title=latexstring("M^{1,1}")*" moment",
								ylabel="Abundance",
								xlabel="Time",
								# ylims=(0,399),
								aspect_ratio=aspect_ratio,
								)
		pmM11 = plot!(pmM11; # Floating label for the panel
            title=L"(c)",
            titlepos=:left,
            grid=false, showaxis=false, ticks=false,
            inset=(1, bbox(tx,ty, 20mm,20mm, :top, :left)),
            subplot=2,
            bg_inside=nothing,
            )
		pmMcorr = plotCorrelation(ssaMeanSol, clnaSol,
								ssaSD=ssaSDSol,
								color=cLNA.pal[4],
								# label=L"Corr(x^{1,0}, x^{0,1})",
								legend=:topright,
								title="Correlation coefficient",
								# ylabel="Abundance",
								xlabel="Time",
								)
		pmMcorr = plot!(pmMcorr; # Floating label for the panel
            title=L"(d)",
            titlepos=:left,
            grid=false, showaxis=false, ticks=false,
            inset=(1, bbox(tx,ty, 20mm,20mm, :top, :left)),
            subplot=2,
            bg_inside=nothing,
            )
		pmHistN = plotMomentsHistograms(M, trajectories, :N; 
								color=cLNA.pal[1], 
								reference=clnaSol,
								aspect_ratio=aspect_ratio,
								)
		pmHistN = plot!(pmHistN; # Floating label for the panel
            title=L"(e)",
            titlepos=:left,
            grid=false, showaxis=false, ticks=false,
            inset=(1, bbox(tx,ty, 20mm,20mm, :top, :left)),
            subplot=2,
            bg_inside=nothing,
            )
		pmHistM10 = plotMomentsHistograms(M, trajectories, :M¹⁰; 
								color=cLNA.pal[2], 
								reference=clnaSol,
								aspect_ratio=aspect_ratio,
								)
		pmHistM10 = plot!(pmHistM10; # Floating label for the panel
            title=L"(f)",
            titlepos=:left,
            grid=false, showaxis=false, ticks=false,
            inset=(1, bbox(tx,ty, 20mm,20mm, :top, :left)),
            subplot=2,
            bg_inside=nothing,
            )
		pmHistM11 = plotMomentsHistograms(M, trajectories, :M¹¹; 
								color=cLNA.pal[3], 
								reference=clnaSol,
								aspect_ratio=aspect_ratio,
								)
		pmHistM11 = plot!(pmHistM11; # Floating label for the panel
            title=L"(g)",
            titlepos=:left,
            grid=false, showaxis=false, ticks=false,
            inset=(1, bbox(tx,ty, 20mm,20mm, :top, :left)),
            subplot=2,
            bg_inside=nothing,
            )
		pmHist2d = plot2dPopulationHistogram(nPool;
								# aspect_ratio=:none,
								# aspect_ratio=0.8,
								aspect_ratio=aspect_ratio,
								# ssaSD=ssaSDSol,
								# color=cLNA.pal[4],
								# # label=L"Corr(x^{1,0}, x^{0,1})",
								# legend=:topright,
								# title="Correlation coefficient",
								# # ylabel="Abundance",
								# xlabel="Time",
								)
		pmHist2d = plot!(pmHist2d; # Floating label for the panel
            title=L"(h)",
            titlepos=:left,
            grid=false, showaxis=false, ticks=false,
            inset=(1, bbox(tx,ty, 20mm,20mm, :top, :left)),
            subplot=2,
            bg_inside=nothing,
            )
		p = plot(pmN, pmM10, pmM11, pmMcorr, 
				 pmHistN, pmHistM10, pmHistM11, pmHist2d; layout=(2,4), size=size, 
						   left_margin=lmargin, 
						   top_margin=tmargin, 
						   right_margin=rmargin,
						   bottom_margin=bmargin,
						   )
		savefig(savepath*"/MutualRepression_$(speed).pdf")
		savefig(savepath*"/MutualRepression_$(speed).png")
	finally
		scalefontsizes(1/fontscale)
	end
	return p
end

#eof
