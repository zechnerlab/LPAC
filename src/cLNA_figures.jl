#=
cLNA Figures: code for defining the figures for the paper.
=#

@export function BinaryBirthDeathCoagulation(;
						T=50.0, NSSA=100, RSSA=1,
						N0=10, Mpc0=1,
						fontscale=1.6,
						size=(1.9*503,2*526).*0.75,
						savepath="Figures",
						name="bBDC_NM1",
						readFromDump::Bool=true,
						tightLayout::Bool=true,
						palette::ColorPalette=cLNA.pal_custom,
						)
	M = Models.IEqCFBDq_new(; kC=0.01,kF=0,kE=0)
	solverFlags = []
	symbols = [:N, :M¹]
	
	# Save the solutions data and allow for reloading them.
	dumpFName="$(savepath)/$(name).jser"
	solDump = nothing
	# Make sure not to attempt to read an unexisting dump
	readFromDump = readFromDump && isfile(dumpFName)
	if !readFromDump
		println("> $(name)::Solving moment equations...")
		### Compute the solutions
		N0 = round(Int, M.Ω * N0)
		Mpc0 = round(Int, M.Ωc * Mpc0)
		clnaSolRaw = @time cLNAsolve(M; T=T, N0=N0, Mpc0=Mpc0, MMap=M.momentsMapping, σMap=M.sigmaMapping, solverFlags...)
		Msymb = popfirst!(setdiff(symbols, [:N,:N2]))
		println("> $(name)::Running SSA simulations...")
		nPool, trajectories, ssaMeanSol, ssaSDSol = runSSASimulations(
					M, symbols...;
					T=T, NSSA=NSSA, RSSA=RSSA, N0=N0, Mpc0=Mpc0, clnaSol=clnaSolRaw)
		clnaSol = cLNA.convertSolution(clnaSolRaw, M.momentsMapping, M.sigmaMapping)
		println("> $(name)::Serializing...")
		solDump = [clnaSol, nPool, trajectories, ssaMeanSol, ssaSDSol]
        @time serialize(dumpFName, solDump)
    else
		println("> $(name)::De-serializing...")
    	@time solDump = deserialize(dumpFName)
    end
    # Unpack dump
    clnaSol, nPool, trajectories, ssaMeanSol, ssaSDSol = solDump
	
	### Plot the figure
	println("> $(name)::Plotting...")
	scalefontsizes(fontscale)
	p = nothing
	try
        local lmargin, tmargin, rmargin, bmargin = 0mm, 5mm, 0mm, 0mm
        local tx, ty = -10mm, -2mm

		pmN = plotMeanComparisons(ssaMeanSol.T, ssaMeanSol.M[:N],
								clnaSol.T, clnaSol.M[:N];
								# ssaRibbon=ssaSDSol.M[:N],
								ssaRibbon=ssaMeanSol.σ[:N],
								clnaRibbon=clnaSol.σ[:N],
								color=palette[2],
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
								color=palette[3],
								label=L"M^1",
								legend=:bottomright,
								title="Number of molecules ("*latexstring("M^1")*")",
								ylabel="Abundance",
								xlabel="Time",
								ylims=(0,415),
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
		# First save as non-tight anyway, avoid strange issues
		savefig(p, savepath*"/$(name).pdf")
		savefig(p, savepath*"/$(name).png")
		if tightLayout
			_savefigTight(p, savepath*"/$(name).pdf")
			_savefigTight(p, savepath*"/$(name).png")
		end
	finally
		scalefontsizes(1/fontscale)
	end
	return p
end

@export function BinaryBirthDeathCoagulationHistograms(;
						T=50.0, NSSA=100, RSSA=1,
						N0=10, Mpc0=1,
						fontscale=1.4,
						size=(2*503,526*0.75).*0.75,
						savepath="Figures",
						name="bBDC_NM1",
						readFromDump::Bool=true,
						tightLayout::Bool=true,
						palette::ColorPalette=cLNA.pal_custom,
						fillalpha=0.3,
						)
	M = Models.IEqCFBDq_new(; kC=0.01,kF=0,kE=0)
	solverFlags = []
	symbols = [:N, :M¹]
	
	# Save the solutions data and allow for reloading them.
	dumpFName="$(savepath)/$(name).jser"
	solDump = nothing
	# Make sure not to attempt to read an unexisting dump
	readFromDump = readFromDump && isfile(dumpFName)
	if !readFromDump
		println("> $(name)::Solving moment equations...")
		### Compute the solutions
		N0 = round(Int, M.Ω * N0)
		Mpc0 = round(Int, M.Ωc * Mpc0)
		clnaSolRaw = @time cLNAsolve(M; T=T, N0=N0, Mpc0=Mpc0, MMap=M.momentsMapping, σMap=M.sigmaMapping, solverFlags...)
		Msymb = popfirst!(setdiff(symbols, [:N,:N2]))
		println("> $(name)::Running SSA simulations...")
		nPool, trajectories, ssaMeanSol, ssaSDSol = runSSASimulations(
					M, symbols...;
					T=T, NSSA=NSSA, RSSA=RSSA, N0=N0, Mpc0=Mpc0, clnaSol=clnaSolRaw)
		clnaSol = cLNA.convertSolution(clnaSolRaw, M.momentsMapping, M.sigmaMapping)
		println("> $(name)::Serializing...")
		solDump = [clnaSol, nPool, trajectories, ssaMeanSol, ssaSDSol]
        @time serialize(dumpFName, solDump)
    else
		println("> $(name)::De-serializing...")
    	@time solDump = deserialize(dumpFName)
    end
    # Unpack dump
    clnaSol, nPool, trajectories, ssaMeanSol, ssaSDSol = solDump
	
	### Plot the figure
	println("> $(name)::Plotting...")
	scalefontsizes(fontscale)
	p = nothing
	try
        local lmargin, tmargin, rmargin, bmargin = 2mm, 10mm, 0mm, 0mm
        local tx, ty = -10mm, -2mm

        pPH = plot1dPopulationHistogram(nPool; 
        				clnaSol=clnaSol,
        				color=palette[3], 
        				linecolor=palette[3], 
						fillalpha=fillalpha,
        				s1Label=L"x_1",
        				)
		pPH = plot!(pPH; # Floating label for the panel
            title=L"(a)",
            titlepos=:left,
            grid=false, showaxis=false, ticks=false,
            inset=(1, bbox(tx,ty, 20mm,20mm, :top, :left)),
            subplot=2,
            bg_inside=nothing,
            )

		pHN = plotMomentsHistogram(M, trajectories, :N;
							color=palette[2],
							linecolor=palette[2],
							fillalpha=fillalpha,
							reference=clnaSol,
							title="",
							# showaxis=:x, 
							xlabel=L"N",
							)
		pHN = plot!(pHN; # Floating label for the panel
            title=L"(b)",
            titlepos=:left,
            grid=false, showaxis=false, ticks=false,
            inset=(1, bbox(tx,ty, 20mm,20mm, :top, :left)),
            subplot=2,
            bg_inside=nothing,
            )

		pHM = plotMomentsHistogram(M, trajectories, :M¹;
							color=palette[3],
							linecolor=palette[3],
							fillalpha=fillalpha,
							reference=clnaSol,
							title="",
							# showaxis=:x,
							xlabel=L"M^1",
							)
		pHM = plot!(pHM; # Floating label for the panel
            title=L"(c)",
            titlepos=:left,
            grid=false, showaxis=false, ticks=false,
            inset=(1, bbox(tx,ty, 20mm,20mm, :top, :left)),
            subplot=2,
            bg_inside=nothing,
            )

		p = plot(pPH, pHN, pHM; layout=(1,3), size=size, 
						   left_margin=lmargin, 
						   top_margin=tmargin, 
						   right_margin=rmargin,
						   bottom_margin=bmargin,
						   )
		savefig(p, savepath*"/$(name)_histograms.pdf")
		savefig(p, savepath*"/$(name)_histograms.png")
		if tightLayout
			_savefigTight(p, savepath*"/$(name)_histograms.pdf")
			_savefigTight(p, savepath*"/$(name)_histograms.png")
		end
	finally
		scalefontsizes(1/fontscale)
	end
	return p
end

@export function SAIC(;
						T=[200.0,200.0,200.0], NSSA=8, RSSA=1,
						N0=25, Mpc0=0,
						fontscale=1.6,
						size=(1.9*503,2*526).*0.75,
						savepath="Figures",
						name="SAIC",
						readFromDump::Bool=true,
						tightLayout::Bool=true,
						palette::ColorPalette=cLNA.pal_custom,
						)
	Ms = [ Models.SAIC(; kref=1e1), Models.SAIC(; kref=4e1), Models.SAIC(; kref=2e1) ]
	setpoints = [1000, 4000, 2000]
	solverFlags = []
	symbols = [:N, :M¹⁰⁰, :M⁰⁰¹, :M⁰¹⁰]

	if typeof(T) <: Number
		T = [ T for M in Ms ]
	end
	
	# Save the solutions data and allow for reloading them.
	dumpFName="$(savepath)/$(name).jser"
	solDump = nothing
	# Make sure not to attempt to read an unexisting dump
	readFromDump = readFromDump && isfile(dumpFName)
	if !readFromDump
		println("> $(name)::Solving moment equations...")
		### Compute the solutions
		N0 = round(Int, Ms[1].Ω * N0)
		Mpc0 = round(Int, Ms[1].Ωc * Mpc0)
		u0 = Ms[1].momentsInit(N0, Mpc0)
		solutions = Vector{Solution}()
		@time for (i,M) in enumerate(Ms)
			clnaSolRaw = @time cLNAsolve(M, u0; T=T[i], MMap=M.momentsMapping, σMap=M.sigmaMapping, solverFlags...)
			clnaSol = cLNA.convertSolution(clnaSolRaw, M.momentsMapping, M.sigmaMapping)
			push!(solutions, clnaSol)
			u0 = deepcopy(clnaSolRaw.u[end])
		end
		clnaSol = solcat(solutions...)
		println("> $(name)::Running SSA simulations...")
		nPool, trajectories, ssaMeanSol, ssaSDSol = runSSASimulations(
						Ms, symbols...;
						T=T, NSSA=NSSA, RSSA=RSSA, N0=N0, Mpc0=Mpc0)
		println("> $(name)::Serializing...")
		solDump = [clnaSol, nPool, trajectories, ssaMeanSol, ssaSDSol]
        @time serialize(dumpFName, solDump)
    else
		println("> $(name)::De-serializing...")
    	@time solDump = deserialize(dumpFName)
    end
    # Unpack dump
    clnaSol, nPool, trajectories, ssaMeanSol, ssaSDSol = solDump
	
	### Plot the figure
	println("> $(name)::Plotting...")
	scalefontsizes(fontscale)
	p = nothing
	try
        local lmargin, tmargin, rmargin, bmargin = 0mm, 5mm, 0mm, 0mm
        local tx, ty = -10mm, -2mm

		pmN = plotMeanComparisons(ssaMeanSol.T, ssaMeanSol.M[:N],
								clnaSol.T, clnaSol.M[:N];
								ssaRibbon=ssaMeanSol.σ[:N],
								clnaRibbon=clnaSol.σ[:N],
								clnaLabel=false,
								color=palette[2],
								label=L"N",
								legend=:bottomright,
								title="Number of compartments ("*latexstring("N")*")",
								ylabel="Abundance",
								xlabel=nothing,
								# ylims=(0,39.9),
								)
		## Shaded gray background
        pmN = vspan!(pmN, cumsum(T)[1:2], color="#585858", alpha=0.1, label=false)

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
								color=palette[3],
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
								color=palette[1],
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
								color=palette[4],
								# label=L"M^{0,0,1}",
								label=L"Q",
								legend=:topright,
								# ylims=(0,399),
								)
		## Setpoint
		SP = _getSetpointValues(clnaSol.T, T, setpoints)
		pmM = plot!(pmM, clnaSol.T, SP, color="black", alpha=0.65, linestyle=:dot, label=L"Q^*")
		## Shaded gray background
        pmM = vspan!(pmM, cumsum(T)[1:2], color="#585858", alpha=0.1, label=false)

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
		# First save as non-tight anyway, avoid strange issues
		savefig(p, savepath*"/$(name).pdf")
		savefig(p, savepath*"/$(name).png")
		if tightLayout
			_savefigTight(p, savepath*"/$(name).pdf")
			_savefigTight(p, savepath*"/$(name).png")
		end
	finally
		scalefontsizes(1/fontscale)
	end
	return p
end
@export function SAICHistograms(;
						T=[200.0,200.0,200.0], NSSA=8, RSSA=1,
						N0=25, Mpc0=0,
						fontscale=1.4,
						size=(4*503,526*0.75).*0.75,
						savepath="Figures",
						name="SAIC",
						readFromDump::Bool=true,
						tightLayout::Bool=true,
						palette::ColorPalette=cLNA.pal_custom,
						fillalpha=0.3,
						)
	Ms = [ Models.SAIC(; kref=1e1), Models.SAIC(; kref=4e1), Models.SAIC(; kref=2e1) ]
	setpoints = [1000, 4000, 2000]
	solverFlags = []
	symbols = [:N, :M¹⁰⁰, :M⁰⁰¹, :M⁰¹⁰]

	if typeof(T) <: Number
		T = [ T for M in Ms ]
	end
	
	# Save the solutions data and allow for reloading them.
	dumpFName="$(savepath)/$(name).jser"
	solDump = nothing
	# Make sure not to attempt to read an unexisting dump
	readFromDump = readFromDump && isfile(dumpFName)
	if !readFromDump
		println("> $(name)::Solving moment equations...")
		### Compute the solutions
		N0 = round(Int, Ms[1].Ω * N0)
		Mpc0 = round(Int, Ms[1].Ωc * Mpc0)
		u0 = Ms[1].momentsInit(N0, Mpc0)
		solutions = Vector{Solution}()
		@time for (i,M) in enumerate(Ms)
			clnaSolRaw = @time cLNAsolve(M, u0; T=T[i], MMap=M.momentsMapping, σMap=M.sigmaMapping, solverFlags...)
			clnaSol = cLNA.convertSolution(clnaSolRaw, M.momentsMapping, M.sigmaMapping)
			push!(solutions, clnaSol)
			u0 = deepcopy(clnaSolRaw.u[end])
		end
		clnaSol = solcat(solutions...)
		println("> $(name)::Running SSA simulations...")
		nPool, trajectories, ssaMeanSol, ssaSDSol = runSSASimulations(
						Ms, symbols...;
						T=T, NSSA=NSSA, RSSA=RSSA, N0=N0, Mpc0=Mpc0)
		println("> $(name)::Serializing...")
		solDump = [clnaSol, nPool, trajectories, ssaMeanSol, ssaSDSol]
        @time serialize(dumpFName, solDump)
    else
		println("> $(name)::De-serializing...")
    	@time solDump = deserialize(dumpFName)
    end
    # Unpack dump
    clnaSol, nPool, trajectories, ssaMeanSol, ssaSDSol = solDump
	
	### Plot the figure
	println("> $(name)::Plotting...")
	scalefontsizes(fontscale)
	p = nothing
	try
        local lmargin, tmargin, rmargin, bmargin = 2mm, 10mm, 0mm, 0mm
        local tx, ty = -10mm, -2mm
        local M = Ms[end]

        # Species = [Z1,Z2,Q]

        pHZ1 = plotMomentsHistogram(M, trajectories, :M¹⁰⁰;
							color=palette[3],
							linecolor=palette[3],
							fillalpha=fillalpha,
							reference=clnaSol,
							title="",
							# showaxis=:x,
							xlabel=L"Z_1",
							)
		pHZ1 = plot!(pHZ1; # Floating label for the panel
            title=L"(a)",
            titlepos=:left,
            grid=false, showaxis=false, ticks=false,
            inset=(1, bbox(tx,ty, 20mm,20mm, :top, :left)),
            subplot=2,
            bg_inside=nothing,
            )
		pHZ2 = plotMomentsHistogram(M, trajectories, :M⁰¹⁰;
							color=palette[1],
							linecolor=palette[1],
							fillalpha=fillalpha,
							reference=clnaSol,
							title="",
							# showaxis=:x,
							xlabel=L"Z_2",
							)
		pHZ2 = plot!(pHZ2; # Floating label for the panel
            title=L"(b)",
            titlepos=:left,
            grid=false, showaxis=false, ticks=false,
            inset=(1, bbox(tx,ty, 20mm,20mm, :top, :left)),
            subplot=2,
            bg_inside=nothing,
            )
		pHN = plotMomentsHistogram(M, trajectories, :N;
							color=palette[2],
							linecolor=palette[2],
							fillalpha=fillalpha,
							reference=clnaSol,
							title="",
							# showaxis=:x,
							xlabel=L"N",
							)
		pHN = plot!(pHN; # Floating label for the panel
            title=L"(c)",
            titlepos=:left,
            grid=false, showaxis=false, ticks=false,
            inset=(1, bbox(tx,ty, 20mm,20mm, :top, :left)),
            subplot=2,
            bg_inside=nothing,
            )
		pPHq = plot1dPopulationHistogram(nPool; 
        				clnaSol=clnaSol,
        				color=palette[4], #TBC
        				linecolor=palette[4],
						fillalpha=fillalpha,
        				s1=3,
        				s1Label=L"q",
        				M=:M⁰⁰¹,
        				)
		pPHq = plot!(pPHq; # Floating label for the panel
            title=L"(d)",
            titlepos=:left,
            grid=false, showaxis=false, ticks=false,
            inset=(1, bbox(tx,ty, 20mm,20mm, :top, :left)),
            subplot=2,
            bg_inside=nothing,
            )
		pHQ = plotMomentsHistogram(M, trajectories, :M⁰⁰¹;
							color=palette[4],
							linecolor=palette[4],
							fillalpha=fillalpha,
							reference=clnaSol,
							title="",
							# showaxis=:x,
							xlabel=L"Q",
							)
		pHQ = plot!(pHQ; # Floating label for the panel
            title=L"(e)",
            titlepos=:left,
            grid=false, showaxis=false, ticks=false,
            inset=(1, bbox(tx,ty, 20mm,20mm, :top, :left)),
            subplot=2,
            bg_inside=nothing,
            )

		p = plot(pHZ1, pHZ2, pHN, pPHq, pHQ; layout=(1,5), size=size, 
						   left_margin=lmargin, 
						   top_margin=tmargin, 
						   right_margin=rmargin,
						   bottom_margin=bmargin,
						   )
		savefig(p, savepath*"/$(name)_histograms.pdf")
		savefig(p, savepath*"/$(name)_histograms.png")
		if tightLayout
			_savefigTight(p, savepath*"/$(name)_histograms.pdf")
			_savefigTight(p, savepath*"/$(name)_histograms.png")
		end
	finally
		scalefontsizes(1/fontscale)
	end
	return p
end

@export function MutualRepression(;
						T=nothing, NSSA=128, RSSA=1,
						N0=25, Mpc0=0,
						N=2*N0, #Asymptotic level
						kb=100, kd=1e-3*kb,
						kMR=30,
						zeta=1e-1/kb,
						fontscale=1.4,
						size=(4*503,526).*0.75,
						savepath="Figures",
						name="MutualRepression",
						readFromDump::Bool=true,
						aspect_ratio=:none,
						speed=:fast, #∈[:stop,:slow,:mid,:fast,:faster]
						tightLayout::Bool=true,
						tightCorrLegend::Bool=false,
						savePlots::Bool=true,
						startPanel::Char='a',
						isFirstRow::Bool=startPanel=='a',
						isLastRow::Bool=true,
						palette::ColorPalette=cLNA.pal_custom3,
						histofillalpha=0.3,
						solverFlags...
						)
	# S = Dict(:stop=>0, :slow=>0.1, :mid=>1, :fast=>10, :faster=>100)
	S = Dict(:stop=>0, :slow=>1e0, :mid=>1e1, :fast=>1e2, :faster=>1e3, :fastest=>1e4)
	kF = S[speed]*5e-3
	@show speed kb/kF #debug
	rate = kb/kF
	rate = isinf(rate) ? string(rate) : round(Int,rate)
	T = isnothing(T) ? 500/S[speed] : T
	M = Models.MutualRepression(; 
								  kb=kb, kd=kd, 
								  kMR0=kMR, kMR1=kMR,
								  kF=kF, kE=(2/(N-1))*kF,
								  )

	symbols = [:N, :M¹⁰, :M²⁰, :M¹¹] #TODO: amend
	
	# Save the solutions data and allow for reloading them.
	dumpFName="$(savepath)/$(name)_$(speed).jser"
	solDump = nothing
	# Make sure not to attempt to read an unexisting dump
	readFromDump = readFromDump && isfile(dumpFName)
	if !readFromDump
		println("> $(name)::Solving moment equations...")
		### Compute the solutions
		N0 = round(Int, M.Ω * N0)
		Mpc0 = round(Int, M.Ωc * Mpc0)
		clnaSolRaw = @time cLNAsolve(M; T=T, N0=N0, Mpc0=Mpc0, MMap=M.momentsMapping, σMap=M.sigmaMapping, solverFlags...)
		Msymb = popfirst!(setdiff(symbols, [:N,:N2]))
		println("> $(name)::Running SSA simulations...")
		nPool, trajectories, ssaMeanSol, ssaSDSol = runSSASimulations(
					M, symbols...;
					T=T, NSSA=NSSA, RSSA=RSSA, N0=N0, Mpc0=Mpc0, clnaSol=clnaSolRaw)
		clnaSol = cLNA.convertSolution(clnaSolRaw, M.momentsMapping, M.sigmaMapping)
		println("> $(name)::Serializing...")
		solDump = [clnaSol, nPool, trajectories, ssaMeanSol, ssaSDSol]
        @time serialize(dumpFName, solDump)
    else
		println("> $(name)::De-serializing...")
    	@time solDump = deserialize(dumpFName)
    end
    # Unpack dump
    clnaSol, nPool, trajectories, ssaMeanSol, ssaSDSol = solDump
	
	### Plot the figure
	println("> $(name)::Plotting...")
	scalefontsizes(fontscale)
	panel::Char=startPanel
	p = nothing
	pmMcorr = nothing
	try
        # local lmargin, tmargin, rmargin, bmargin = -8mm, -2mm, -2mm, -4mm
        local lmargin, tmargin, rmargin, bmargin = 4mm, 0mm, 0mm, 0mm
        local tx, ty = -0.15, -0.05
        local ix, iy, iw, ih = 0.05, 0.05, 0.4, 0.4
        local jx, jy, jw, jh = 0.1, 0.025, 0.5, 0.46
	    if tightCorrLegend
		    jx, jy, jw, jh = 0.09, 0.025, 0.45, 0.33
		end

		pmN = plotMeanComparisons(ssaMeanSol.T, ssaMeanSol.M[:N],
								clnaSol.T, clnaSol.M[:N];
								# ssaRibbon=ssaSDSol.M[:N],
								ssaRibbon=ssaMeanSol.σ[:N],
								clnaRibbon=clnaSol.σ[:N],
								color=palette[1],
								label=L"N",
								legend=:bottomright,
								title=isFirstRow ? L"N" : nothing,
								ylabel="Abundance",
								xlabel=isLastRow ? "Time" : nothing,
								# ylims=(0,39.9),
								aspect_ratio=aspect_ratio,
								)
		pmN = plot!(pmN; # Floating label for the panel
            title=latexstring("($(panel))"),
            titlepos=:left,
            grid=false, showaxis=false, ticks=false,
            inset=(1, bbox(tx,ty, 20mm,20mm, :top, :left)),
            subplot=2,
            bg_inside=nothing,
            )
		if speed != :stop
			pmN = plot!(pmN; # Floating box for moment histogram
	            inset_subplots=[ (1, bbox(ix,iy,iw,ih, :bottom, :left)) ], # bbox(x, y, width, height, origin...)
	            subplot=3,
	            bg_inside=nothing,
	            ticks=nothing,
	            )
	        pmN = plotMomentsHistogram!(pmN[3], M, trajectories, :N; #...and fill the box with the histogram
							color=palette[1], 
							linecolor=palette[1], 
							fillalpha=histofillalpha,
							reference=clnaSol,
							aspect_ratio=aspect_ratio,
							# ylabel="Count",
							# xlabel="Value",
							# ticks=false,
							title="",
							showaxis=:x,
							)
	    end
	    panel+=1

		pmM10 = plotMeanComparisons(ssaMeanSol.T, ssaMeanSol.M[:M¹⁰],
								clnaSol.T, clnaSol.M[:M¹⁰];
								# ssaRibbon=ssaSDSol.M[:M¹⁰],
								ssaRibbon=ssaMeanSol.σ[:M¹⁰],
								clnaRibbon=clnaSol.σ[:M¹⁰],
								color=palette[2],
								label=L"M^{1,0}",
								legend=:bottomright,
								title=isFirstRow ? L"M^{1,0}" : nothing,
								# ylabel="Abundance",
								xlabel=isLastRow ? "Time" : nothing,
								# ylims=(0,399),
								aspect_ratio=aspect_ratio,
								)
		pmM10 = plot!(pmM10; # Floating label for the panel
            title=latexstring("($(panel))"),
            titlepos=:left,
            grid=false, showaxis=false, ticks=false,
            inset=(1, bbox(tx,ty, 20mm,20mm, :top, :left)),
            subplot=2,
            bg_inside=nothing,
            )
		pmM10 = plot!(pmM10; # Floating box for moment histogram
            inset_subplots=[ (1, bbox(ix,iy,iw,ih, :bottom, :left)) ], # bbox(x, y, width, height, origin...)
            subplot=3,
            bg_inside=nothing,
            ticks=nothing,
            )
        plotMomentsHistogram!(pmM10[3], M, trajectories, :M¹⁰; #...and fill the box with the histogram
						color=palette[2], 
						linecolor=palette[2], 
						fillalpha=histofillalpha,
						reference=clnaSol,
						aspect_ratio=aspect_ratio,
						# ylabel="Count",
						# xlabel="Value",
						# ticks=false,
						title="",
						showaxis=:x,
						)
        panel+=1

		pmM11 = plotMeanComparisons(ssaMeanSol.T, ssaMeanSol.M[:M¹¹],
								clnaSol.T, clnaSol.M[:M¹¹];
								# ssaRibbon=ssaSDSol.M[:M¹¹],
								ssaRibbon=ssaMeanSol.σ[:M¹¹],
								# clnaRibbon=clnaSol.σ[:M¹¹],
								color=palette[3],
								label=L"M^{1,1}",
								legend=:bottomright,
								title=isFirstRow ? L"M^{1,1}" : nothing,
								# ylabel="Abundance",
								xlabel=isLastRow ? "Time" : nothing,
								# ylims=(0,399),
								aspect_ratio=aspect_ratio,
								)
		pmM11 = plot!(pmM11; # Floating label for the panel
            title=latexstring("($(panel))"),
            titlepos=:left,
            grid=false, showaxis=false, ticks=false,
            inset=(1, bbox(tx,ty, 20mm,20mm, :top, :left)),
            subplot=2,
            bg_inside=nothing,
            )
		pmM11 = plot!(pmM11; # Floating box for moment histogram
            inset_subplots=[ (1, bbox(ix,iy,iw,ih, :bottom, :left)) ], # bbox(x, y, width, height, origin...)
            subplot=3,
            bg_inside=nothing,
            ticks=nothing,
            )
        plotMomentsHistogram!(pmM11[3], M, trajectories, :M¹¹; #...and fill the box with the histogram
						color=palette[3], 
						linecolor=palette[3], 
						fillalpha=histofillalpha,
						reference=clnaSol,
						aspect_ratio=aspect_ratio,
						# ylabel="Count",
						# xlabel="Value",
						# ticks=false,
						title="",
						showaxis=:x,
						)
        panel+=1

		pmMcorr = plotCorrelation(ssaMeanSol, clnaSol,
								ssaSD=ssaSDSol,
								color=palette[4],
								legend=tightCorrLegend ? :topleft : :bottomright,
								tightLegend=tightCorrLegend,
								title=isFirstRow ? "Correlation" : nothing,
								ylabel="Correlation",
								xlabel=isLastRow ? "Time" : nothing,
								)
		pmMcorr = plot!(pmMcorr; # Floating label for the panel
            title=latexstring("($(panel))"),
            titlepos=:left,
            grid=false, showaxis=false, ticks=false,
            inset=(1, bbox(tx,ty, 20mm,20mm, :top, :left)),
            subplot=2,
            bg_inside=nothing,
            )
		pmMcorr = plot!(pmMcorr; # Floating box for moment histogram
            inset_subplots=[ (1, bbox(jx,jy,jw,jh, :top, :right)) ], # bbox(x, y, width, height, origin...)
            subplot=3,
            bg_inside=nothing,
            ticks=nothing,
            )
		if !tightCorrLegend
			plot2dPopulationHistogram!(pmMcorr[3], nPool;
				aspect_ratio=aspect_ratio,
				clnaSol=clnaSol,
				)
		else
			plot2dPopulationHistogram!(pmMcorr[3], nPool;
				aspect_ratio=aspect_ratio,
				clnaSol=clnaSol,
				colorbar_ticks=nothing, # remove ticks and numbers from the colorbar #test
				)
		end

		p = plot(pmN, pmM10, pmM11, pmMcorr, 
				 ; layout=(1,4), 
				 # ; layout=(2,2), 
				 		   size=size, 
						   left_margin=lmargin, 
						   top_margin=tmargin, 
						   right_margin=rmargin,
						   bottom_margin=bmargin,
						   )
		if savePlots
			# First save as non-tight anyway, avoid strange issues
			savefig(p, savepath*"/$(name)_$(speed)_$(rate).pdf")
			savefig(p, savepath*"/$(name)_$(speed)_$(rate).png")
			if tightLayout
				_savefigTight(p, savepath*"/$(name)_$(speed)_$(rate).pdf")
				_savefigTight(p, savepath*"/$(name)_$(speed)_$(rate).png")
			end
		end
	finally
		scalefontsizes(1/fontscale)
	end
	return p
end
@export function MutualRepressionCorrelations(;
						T=nothing, NSSA=128, RSSA=1,
						N0=25, Mpc0=0,
						N=2*N0, #Asymptotic level
						kb=100, kd=1e-3*kb,
						kMR=30,
						zeta=1e-1/kb,
						fontscale=1.4,
						size=(4*503,526).*0.75,
						savepath="Figures",
						name="MutualRepression",
						readFromDump::Bool=true,
						aspect_ratio=:none,
						speeds=[:stop,:slow,:mid,:fast,:faster],
						tightLayout::Bool=true,
						tightCorrLegend::Bool=false,
						savePlots::Bool=true,
						startPanel::Char='e',
						palette::ColorPalette=cLNA.pal_custom3,
						solverFlags...
						)
	S = Dict(:stop=>0, :slow=>1e0, :mid=>1e1, :fast=>1e2, :faster=>1e3, :fastest=>1e4)
	P = Dict()
	panel::Char=startPanel
	for speed in speeds
		kF = S[speed]*5e-3
		@show speed kb/kF #debug
		rate = kb/kF
		rate = isinf(rate) ? string(rate) : round(Int,rate)
		T = isnothing(T) ? 500/S[speed] : T
		M = Models.MutualRepression(; 
									  kb=kb, kd=kd, 
									  kMR0=kMR, kMR1=kMR,
									  kF=kF, kE=(2/(N-1))*kF,
									  )

		symbols = [:N, :M¹⁰, :M²⁰, :M¹¹] #TODO: amend
		
		# Save the solutions data and allow for reloading them.
		dumpFName="$(savepath)/$(name)_$(speed).jser"
		solDump = nothing
		# Make sure not to attempt to read an unexisting dump
		readFromDump = readFromDump && isfile(dumpFName)
		if !readFromDump
			println("> $(name)::Solving moment equations...")
			### Compute the solutions
			N0 = round(Int, M.Ω * N0)
			Mpc0 = round(Int, M.Ωc * Mpc0)
			clnaSolRaw = @time cLNAsolve(M; T=T, N0=N0, Mpc0=Mpc0, MMap=M.momentsMapping, σMap=M.sigmaMapping, solverFlags...)
			Msymb = popfirst!(setdiff(symbols, [:N,:N2]))
			println("> $(name)::Running SSA simulations...")
			nPool, trajectories, ssaMeanSol, ssaSDSol = runSSASimulations(
						M, symbols...;
						T=T, NSSA=NSSA, RSSA=RSSA, N0=N0, Mpc0=Mpc0, clnaSol=clnaSolRaw)
			clnaSol = cLNA.convertSolution(clnaSolRaw, M.momentsMapping, M.sigmaMapping)
			println("> $(name)::Serializing...")
			solDump = [clnaSol, nPool, trajectories, ssaMeanSol, ssaSDSol]
	        @time serialize(dumpFName, solDump)
	    else
			println("> $(name)::De-serializing...")
	    	@time solDump = deserialize(dumpFName)
	    end
	    # Unpack dump
	    clnaSol, nPool, trajectories, ssaMeanSol, ssaSDSol = solDump
		
		### Plot the figure
		println("> $(name)::Plotting...")
		scalefontsizes(fontscale)
		pmMcorr = nothing
		try
	        # local lmargin, tmargin, rmargin, bmargin = -8mm, -2mm, -2mm, -4mm
	        local lmargin, tmargin, rmargin, bmargin = 0mm, 0mm, 0mm, 0mm
	        local tx, ty = -0.15, -0.05
	        local ix, iy, iw, ih = 0.05, 0.05, 0.4, 0.4
	        # local jx, jy, jw, jh = 0.1, 0.025, 0.5, 0.46
	        local jx, jy, jw, jh = 0.1, 0.025, 0.45, 0.33

			pmMcorr = plotCorrelation(ssaMeanSol, clnaSol,
									ssaSD=ssaSDSol,
									color=palette[4],
									# legend=tightCorrLegend ? :topleft : :bottomright,
									legend=false,
									tightLegend=tightCorrLegend,
									# title="Correlation",
									title=latexstring("k_b/k_D = $rate"),
									ylabel=panel==startPanel ? "Correlation" : nothing,
									xlabel="Time",
									)
			pmMcorr = plot!(pmMcorr; # Floating label for the panel
	            title=latexstring("($(panel))"),
	            titlepos=:left,
	            grid=false, showaxis=false, ticks=false,
	            inset=(1, bbox(tx,ty, 20mm,20mm, :top, :left)),
	            subplot=2,
	            bg_inside=nothing,
	            )
			if speed ∈ [:fastest]
				pmMcorr = plot!(pmMcorr; # Floating box for moment histogram
		            inset_subplots=[ (1, bbox(jx,jy,jw,jh, :bottom, :right)) ], # bbox(x, y, width, height, origin...)
		            subplot=3,
		            bg_inside=nothing,
		            ticks=nothing,
		            )
			else
				pmMcorr = plot!(pmMcorr; # Floating box for moment histogram
		            inset_subplots=[ (1, bbox(jx,jy,jw,jh, :top, :right)) ], # bbox(x, y, width, height, origin...)
		            subplot=3,
		            bg_inside=nothing,
		            ticks=nothing,
		            )
			end
			# @show length(nPool) #debug
			# @show extrema(nPool[1,:]) #debug
			# @show extrema(nPool[2,:]) #debug
			plot2dPopulationHistogram!(pmMcorr[3], nPool;
				aspect_ratio=aspect_ratio,
				clnaSol=clnaSol,
				# s1Label=nothing, s2Label=nothing,
				# colorbar = false,
				colorbar_ticks=nothing, # see https://github.com/JuliaPlots/Plots.jl/issues/3174#issuecomment-806313344
				)

			P[speed] = pmMcorr
			panel+=1
		finally
			scalefontsizes(1/fontscale)
		end
	end # for speed in speeds

	println("> $(name)::Grouping plots...")
	scalefontsizes(fontscale)
	Pvec = [ P[s] for s in speeds ]
	p = nothing
	try
		# local lmargin, tmargin, rmargin, bmargin = -8mm, -2mm, -2mm, -4mm
        local lmargin, tmargin, rmargin, bmargin = 0mm, 0mm, 0mm, 0mm
        local tx, ty = -0.15, -0.05
        local ix, iy, iw, ih = 0.05, 0.05, 0.4, 0.4
        local jx, jy, jw, jh = 0.1, 0.025, 0.5, 0.46
		p = plot(Pvec..., 
			 ; layout=(1,length(Pvec)), 
			 		   size=size, 
					   left_margin=lmargin, 
					   top_margin=tmargin, 
					   right_margin=rmargin,
					   bottom_margin=bmargin,
					   )
		if savePlots
			# First save as non-tight anyway, avoid strange issues
			savefig(p, savepath*"/$(name)_Correlations.pdf")
			savefig(p, savepath*"/$(name)_Correlations.png")
			if tightLayout
				_savefigTight(p, savepath*"/$(name)_Correlations.pdf")
				_savefigTight(p, savepath*"/$(name)_Correlations.png")
			end
		end
	finally
		scalefontsizes(1/fontscale)
	end
	return p
end
@export function MutualRepressionCorrelationsCombine(;
						size=(4*503,526).*0.75,
						savepath="Figures",
						name="MutualRepression",
						fontscale=1.4,
						tightLayout::Bool=true,
						savePlots::Bool=true,
						titlefontsize=0.9*14*fontscale, # Plots.jl default value is 14*fontscale
						row2scale=0.80,
						speed=:fast,
						opts...
						)
	p1 = MutualRepression(; 
						N0=25, N=50, kMR=10, T=50, NSSA=128, 
						size=size, tightLayout=tightLayout, readFromDump=true,
						fontscale=fontscale, savepath=savepath, name=name,
						savePlots=savePlots,
						speed=speed,
						opts...
						)
	p2 = MutualRepressionCorrelations(; 
						N0=25, N=50, kMR=10, T=50, NSSA=128, 
						size=(1,row2scale).*size, tightLayout=tightLayout, readFromDump=true,
						fontscale=fontscale, savepath=savepath, name=name,
						savePlots=savePlots,
						opts...
						)
	p = nothing
	scalefontsizes(fontscale)
	try
		# local lmargin, tmargin, rmargin, bmargin = -8mm, -2mm, -2mm, -4mm
        local lmargin, tmargin, rmargin, bmargin = 0mm, 0mm, 0mm, 6mm
        local newsize = (1,1+row2scale) .* size
        @show newsize
		p = plot(p1, p2
			 ; layout=(2,1), 
			 		   size=newsize, 
					   left_margin=lmargin, 
					   top_margin=tmargin, 
					   right_margin=rmargin,
					   bottom_margin=bmargin,
					   titlefontsize=titlefontsize, #rescale all titles in subplots
					   )
		if savePlots
			# First save as non-tight anyway, avoid strange issues
			savefig(p, savepath*"/$(name)_$(speed)_CorrelationsCombine.pdf")
			savefig(p, savepath*"/$(name)_$(speed)_CorrelationsCombine.png")
			if tightLayout
				_savefigTight(p, savepath*"/$(name)_$(speed)_CorrelationsCombine.pdf")
				_savefigTight(p, savepath*"/$(name)_$(speed)_CorrelationsCombine.png")
			end
		end
	finally
		scalefontsizes(1/fontscale)
	end
	return p
end

function _savefigTight(p, fname)
	p.o["canvas"][:print_figure](fname,bbox_inches="tight")
end

@export function MutualRepressionCorrelationComparison(;
						T=nothing, NSSA=128, RSSA=1,
						N0=25, Mpc0=0,
						N=2*N0, #Asymptotic level
						kb=100, kd=1e-3*kb,
						kMR=30,
						zeta=1e-1/kb,
						fontscale=1.4,
						size=(4*503,526).*0.75,
						# size=(2*503,526).*0.75,
						savepath="Figures",
						name="MutualRepression",
						readFromDump::Bool=true,
						speeds=[:stop,:slow,:mid,:fast,:faster],
						aspect_ratio=:none,
						tightLayout::Bool=true,
						titlefontsize=1.0*14*fontscale, # Plots.jl default value is 14*fontscale
						savePlots::Bool=true,
						solverFlags...
						)
	ps = []
	panel::Char='a'
	for speed in speeds
		local p= MutualRepression(;
			T=T, NSSA=NSSA, RSSA=RSSA, N0=N0, Mpc0=Mpc0, N=N,
			kb=kb, kd=kd, kMR=kMR, zeta=zeta, fontscale=fontscale, size=size,
			savepath=savepath, name=name, readFromDump=readFromDump, aspect_ratio=aspect_ratio,
			speed=speed, tightLayout=tightLayout,
			tightCorrLegend=true, savePlots=false,
			startPanel=panel,
			isLastRow=speed==speeds[end],
			solverFlags...
			)
		push!(ps, p)
		panel+=4
	end
	l = length(ps)
    local lmargin, tmargin, rmargin, bmargin = 0mm, 0mm, 0mm, 15mm
	p = plot(ps...; layout=(l, 1), 
			size=size.*(1,l),
			left_margin=lmargin, 
			top_margin=tmargin, 
			right_margin=rmargin,
			bottom_margin=bmargin,
			titlefontsize=titlefontsize, #rescale all titles in subplots
					   )
	if savePlots
		# First save as non-tight anyway, avoid strange issues
		savefig(p, savepath*"/$(name)_Comparison.pdf")
		savefig(p, savepath*"/$(name)_Comparison.png")
		if tightLayout
			_savefigTight(p, savepath*"/$(name)_Comparison.pdf")
			_savefigTight(p, savepath*"/$(name)_Comparison.png")
		end
	end
	return p
end

@export function MutualRepressionBimodal(;
						T=nothing, NSSA=128, RSSA=1,
						N0=25, Mpc0=5,
						N=2*N0, #Asymptotic level
						kb=100,
						low=5, high=50,
						zeta=1e-1/kb,
						fontscale=1.4,
						size=(4*503,526).*0.75,
						savepath="Figures",
						name="MutualRepressionBimodal",
						readFromDump::Bool=true,
						aspect_ratio=:none,
						speed=:midfast, #∈[:stop,:slow,:mid,:midfast,:fast,:faster]
						νspeed=:midfast, #∈[:stop,:slow,:mid,:midfast,:fast,:faster]
						tightLayout::Bool=true,
						tightCorrLegend::Bool=true,
						savePlots::Bool=true,
						solverFlags...
						)
	# S = Dict(:stop=>0, :slow=>0.1, :mid=>1, :fast=>10, :faster=>100)
	S = Dict(:stop=>0, :slow=>1e0, :mid=>1e1, :midfast=>2.5e1, :fast=>5e1, :faster=>1e2, :fastest=>1e3)
	νS = Dict(:stop=>0, :slow=>1e-4, :mid=>5e-4, :midfast=>1e-3, :fast=>5e-3, :faster=>1e-2, :fastest=>1e-1)
	kF = S[speed]*1e-2
	kν = νS[νspeed]
	@show speed kb/kF #debug
	@show νspeed #debug
	rate = kb/kF
	rate = isinf(rate) ? string(rate) : round(Int,rate)
	T = isnothing(T) ? 500/S[speed] : T
	M = Models.GeneSwitchAnnihilating(;
								N=N,
								n=4.5,
								a=low, b=high,
								β=kb,
								k=0.3,
								ν=kb*kν,
								kD=kF,
								)

	symbols = [:N, :M¹⁰, :M²⁰, :M¹¹] #TODO: amend
	
	# Save the solutions data and allow for reloading them.
	dumpFName="$(savepath)/$(name)_$(speed)_nu$(νspeed).jser"
	@show dumpFName #debug
	solDump = nothing
	# Make sure not to attempt to read an unexisting dump
	readFromDump = readFromDump && isfile(dumpFName)
	@show readFromDump #debug
	if !readFromDump
		println("> $(name)::Solving moment equations...")
		### Compute the solutions
		N0 = round(Int, M.Ω * N0)
		Mpc0 = round(Int, M.Ωc * Mpc0)
		clnaSolRaw = @time cLNAsolve(M; T=T, N0=N0, Mpc0=Mpc0, MMap=M.momentsMapping, σMap=M.sigmaMapping, solverFlags...)
		Msymb = popfirst!(setdiff(symbols, [:N,:N2]))
		println("> $(name)::Running SSA simulations...")
		nPool, trajectories, ssaMeanSol, ssaSDSol = runSSASimulations(
					M, symbols...;
					T=T, NSSA=NSSA, RSSA=RSSA, N0=N0, Mpc0=Mpc0, clnaSol=clnaSolRaw)
		clnaSol = cLNA.convertSolution(clnaSolRaw, M.momentsMapping, M.sigmaMapping)
		println("> $(name)::Serializing...")
		solDump = [clnaSol, nPool, trajectories, ssaMeanSol, ssaSDSol]
        @time serialize(dumpFName, solDump)
    else
		println("> $(name)::De-serializing...")
    	@time solDump = deserialize(dumpFName)
    end
    # Unpack dump
    clnaSol, nPool, trajectories, ssaMeanSol, ssaSDSol = solDump
	
	### Plot the figure
	println("> $(name)::Plotting...")
	scalefontsizes(fontscale)
	p = nothing
	try
        # local lmargin, tmargin, rmargin, bmargin = -8mm, -2mm, -2mm, -4mm
        local lmargin, tmargin, rmargin, bmargin = 4mm, 0mm, 0mm, 0mm
        local tx, ty = -0.15, -0.05
        local ix, iy, iw, ih = 0.05, 0.05, 0.4, 0.4
        # local jx, jy, jw, jh = 0.1, 0.025, 0.5, 0.46
	    local jx, jy, jw, jh = 0.05, 0.025, 0.45, 0.33

		pmN = plotMeanComparisons(ssaMeanSol.T, ssaMeanSol.M[:N],
								clnaSol.T, clnaSol.M[:N];
								# ssaRibbon=ssaSDSol.M[:N],
								ssaRibbon=ssaMeanSol.σ[:N],
								clnaRibbon=clnaSol.σ[:N],
								color=cLNA.pal_custom3[1],
								label=L"N",
								legend=:bottomright,
								title=L"N",
								ylabel="Abundance",
								xlabel="Time",
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
		if speed != :stop
			pmN = plot!(pmN; # Floating box for moment histogram
	            inset_subplots=[ (1, bbox(ix,iy,iw,ih, :bottom, :left)) ], # bbox(x, y, width, height, origin...)
	            subplot=3,
	            bg_inside=nothing,
	            ticks=nothing,
	            )
	        pmN = plotMomentsHistogram!(pmN[3], M, trajectories, :N; #...and fill the box with the histogram
							color=cLNA.pal_custom3[1], 
							reference=clnaSol,
							aspect_ratio=aspect_ratio,
							# ylabel="Count",
							# xlabel="Value",
							# ticks=false,
							title="",
							showaxis=:x,
							)
	    end

		pmM10 = plotMeanComparisons(ssaMeanSol.T, ssaMeanSol.M[:M¹⁰],
								clnaSol.T, clnaSol.M[:M¹⁰];
								# ssaRibbon=ssaSDSol.M[:M¹⁰],
								ssaRibbon=ssaMeanSol.σ[:M¹⁰],
								clnaRibbon=clnaSol.σ[:M¹⁰],
								color=cLNA.pal_custom3[2],
								label=L"M^{1,0}",
								legend=:topleft,
								title=L"M^{1,0}",
								# ylabel="Abundance",
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
		pmM10 = plot!(pmM10; # Floating box for moment histogram
            inset_subplots=[ (1, bbox(ix,iy,iw,ih, :bottom, :left)) ], # bbox(x, y, width, height, origin...)
            subplot=3,
            bg_inside=nothing,
            ticks=nothing,
            )
        plotMomentsHistogram!(pmM10[3], M, trajectories, :M¹⁰; #...and fill the box with the histogram
						color=cLNA.pal_custom3[2], 
						reference=clnaSol,
						aspect_ratio=aspect_ratio,
						# ylabel="Count",
						# xlabel="Value",
						# ticks=false,
						title="",
						showaxis=:x,
						)

		pmM11 = plotMeanComparisons(ssaMeanSol.T, ssaMeanSol.M[:M¹¹],
								clnaSol.T, clnaSol.M[:M¹¹];
								# ssaRibbon=ssaSDSol.M[:M¹¹],
								ssaRibbon=ssaMeanSol.σ[:M¹¹],
								# clnaRibbon=clnaSol.σ[:M¹¹],
								color=cLNA.pal_custom3[3],
								label=L"M^{1,1}",
								legend=:topleft,
								title=L"M^{1,1}",
								# ylabel="Abundance",
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
		pmM11 = plot!(pmM11; # Floating box for moment histogram
            inset_subplots=[ (1, bbox(ix,iy,iw,ih, :bottom, :left)) ], # bbox(x, y, width, height, origin...)
            subplot=3,
            bg_inside=nothing,
            ticks=nothing,
            )
        plotMomentsHistogram!(pmM11[3], M, trajectories, :M¹¹; #...and fill the box with the histogram
						color=cLNA.pal_custom3[3], 
						reference=clnaSol,
						aspect_ratio=aspect_ratio,
						# ylabel="Count",
						# xlabel="Value",
						# ticks=false,
						title="",
						showaxis=:x,
						)

		pmMcorr = plotCorrelation(ssaMeanSol, clnaSol,
								ssaSD=ssaSDSol,
								color=cLNA.pal_custom3[4],
								legend=tightCorrLegend ? :topleft : :bottomright,
								tightLegend=tightCorrLegend,
								title="Correlation",
								ylabel="Correlation",
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
		pmMcorr = plot!(pmMcorr; # Floating box for moment histogram
            inset_subplots=[ (1, bbox(jx,jy,jw,jh, :top, :right)) ], # bbox(x, y, width, height, origin...)
            subplot=3,
            bg_inside=nothing,
            ticks=nothing, # Uncomment to remove ticks from 2d histogram
            )
		plot2dPopulationHistogram!(pmMcorr[3], nPool;
			aspect_ratio=aspect_ratio,
			colorbar_ticks=nothing, # see https://github.com/JuliaPlots/Plots.jl/issues/3174#issuecomment-806313344
			clnaSol=clnaSol,
			)

		p = plot(pmN, pmM10, pmM11, pmMcorr, 
				 ; layout=(1,4), 
				 # ; layout=(2,2), 
				 		   size=size, 
						   left_margin=lmargin, 
						   top_margin=tmargin, 
						   right_margin=rmargin,
						   bottom_margin=bmargin,
						   )
		if savePlots
			# First save as non-tight anyway, avoid strange issues
			savefig(p, savepath*"/$(name)_$(speed)_$(rate)_nu$(νspeed)_$(kν).pdf")
			savefig(p, savepath*"/$(name)_$(speed)_$(rate)_nu$(νspeed)_$(kν).png")
			if tightLayout
				_savefigTight(p, savepath*"/$(name)_$(speed)_$(rate)_nu$(νspeed)_$(kν).pdf")
				_savefigTight(p, savepath*"/$(name)_$(speed)_$(rate)_nu$(νspeed)_$(kν).png")
			end
		end
	finally
		scalefontsizes(1/fontscale)
	end
	return p
end

@export function generateAllFigures(; savepath=".", readFromDump::Bool=true)
	Figures.BinaryBirthDeathCoagulation(; NSSA=100, RSSA=1, savepath=savepath, readFromDump=readFromDump)
	Figures.BinaryBirthDeathCoagulationHistograms(; NSSA=100, RSSA=1, savepath=savepath, readFromDump=readFromDump)
	Figures.SAIC(; T=200, NSSA=8, RSSA=1, savepath=savepath, readFromDump=readFromDump)
	Figures.SAICHistograms(; T=200, NSSA=8, RSSA=1, savepath=savepath, readFromDump=readFromDump)
	# return nothing
	for s in reverse([:stop,:slow,:mid,:fast,:faster,:fastest])
		Figures.MutualRepression(; speed=s, N0=25, N=50, kMR=10, T=50, NSSA=128, 
			savepath="$(savepath)/kMR10_Lite2", readFromDump=readFromDump, tightLayout=true, fontscale=1.4)
	end
	Figures.MutualRepressionCorrelationsCombine(; N0=25, N=50, kMR=10, T=50, NSSA=128, 
		savepath="$(savepath)/kMR10_Lite2", readFromDump=true, tightLayout=true, fontscale=1.4)
	Figures.MutualRepressionCorrelationComparison(; N0=25, N=50, kMR=10, T=50, NSSA=128, 
		savepath="$(savepath)/kMR10_Lite2", readFromDump=true, tightLayout=true, fontscale=1.4)
	return nothing
end

#eof
