#=
Plotter: some plotting-related functions.
=#

const pal = palette(:tab10) # Non color-blind/CVD friendly :(
const pal_oi = palette(:okabe_ito) # Friendly :)
const pal_tl = palette(:tol_light) # Friendly :)
const pal_tb = palette(:tol_bright) # Friendly :)
const pal_custom = palette( [ pal_tb[1], pal_tb[2], pal_tb[3], pal_oi[6], pal_oi[1] ] )
const pal_custom2 = palette( [ pal_tb[1], pal_oi[3], pal_tb[3], pal_oi[6], pal_oi[1] ] )
const pal_custom3 = palette( [ pal_tb[2], pal_oi[3], pal_tb[3], pal_oi[6], pal_oi[1] ] )

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

@export function plotMeanComparisons!(p, Tssa, Mssa, Tclna, Mclna; 
					label=L"M", 
					color=pal[1], 
					legend=:bottomright,
					ssaRibbon=nothing,
					clnaRibbon=nothing,
					clnaLabel::Bool=true,
					plotFlags...
					)
	p = plotMoment!(p, Tssa, Mssa; 
			color=color, 
			label=clnaLabel ? L"\langle "*label*L"_{ssa} \rangle" : label,
			legend=legend,
			σ=ssaRibbon,
			plotFlags...
			)
	p = plotMoment!(p, Tclna, Mclna; 
			color=color, 
			linestyle=:dash,
			linewidth=2.0, 
			label=clnaLabel ? L"\langle "*label*L"_{pred} \rangle" : false,
			legend=legend,
			plotFlags...
			)
	@show label
	μerr = mean( norm1error(Tclna, Mclna, Tssa, Mssa; relative=true) )
	@show μerr
	if !isnothing(clnaRibbon)
		#plot ribbon extrema as dotted lines
		p = plot!(p, Tclna, Mclna .- clnaRibbon;
					color=color,
					linestyle=:dot,
					linewidth=2.0, 
					label=clnaLabel ? L"\sigma("*label*L")_{pred}" : false,
					plotFlags...)
		p = plot!(p, Tclna, Mclna .+ clnaRibbon;
					color=color,
					linestyle=:dot,
					linewidth=2.0, 
					label=false,
					plotFlags...)
		σerr = mean( norm1error(Tclna, clnaRibbon, Tssa, ssaRibbon; relative=true) )
		@show σerr
	end
	return p
end
@export function plotMeanComparisons(Tssa, Mssa, Tclna, Mclna; 
					label=L"M", 
					color=pal[1], 
					legend=:bottomright,
					ssaRibbon=nothing,
					clnaRibbon=nothing,
					clnaLabel::Bool=true,
					plotFlags...
					)
	p = plot()
	plotMeanComparisons!(p, Tssa, Mssa, Tclna, Mclna;
			label=label, color=color, legend=legend, 
			ssaRibbon=ssaRibbon, clnaRibbon=clnaRibbon, clnaLabel=clnaLabel, plotFlags...)
end

@export function plotVarianceComparisons(Tssa, VarSsa, Tclna, VarClna; 
					label=L"Var", 
					color=pal[1], 
					legend=:bottomright,
					# legend=:topright,
					ssaRibbon=nothing,
					plotFlags...
					)
	p = plotMoment(Tssa, VarSsa; 
			color=color, 
			label=L"Var("*label*L"_{ssa})",
			legend=legend,
			σ=ssaRibbon,
			plotFlags...
			)
	p = plotMoment!(p, Tclna, VarClna; 
			color=color, 
			linestyle=:dash,
			linewidth=2.0, 
			label=L"Var("*label*L"_{pred})",
			legend=legend,
			plotFlags...
			)
end
@export function plotVarianceComparisons!(p, Tssa, VarSsa, Tclna, VarClna; 
					label=L"Var", 
					color=pal[1], 
					legend=:bottomright,
					# legend=:topright,
					ssaRibbon=nothing,
					plotFlags...
					)
	p2 = twinx(p)
	plotMoment!(p2, Tssa, VarSsa; 
			color=color, 
			label=L"Var("*label*L"_{ssa})",
			legend=legend,
			σ=ssaRibbon,
			linestyle=:dot,
			fillstyle=:|,
			plotFlags...
			)
	plotMoment!(p2, Tclna, VarClna; 
			color=color, 
			linestyle=:dash,
			linewidth=2.0, 
			label=L"Var("*label*L"_{pred})",
			legend=legend,
			plotFlags...
			)
end

_cov(M11::Number, M10::Number, M01::Number, N::Number) = M11 - (M10*M01 / N)
_cov(sol::Solution, M11::Symbol, M10::Symbol, M01::Symbol, N::Symbol=:N) = @. sol.M[M11] - (sol.M[M10]*sol.M[M01] / sol.M[N])
_var(M2::Number, M1::Number, N::Number) = _cov(M2, M1, M1, N)
_var(sol::Solution, M2::Symbol, M1::Symbol, N::Symbol=:N) = _cov(sol, M2, M1, M1, N)
_std(M2::Number, M1::Number, N::Number) = sqrt( _cov(M2, M1, M1, N) )
# _std(sol::Solution, M2::Symbol, M1::Symbol, N::Symbol=:N) = sqrt.( _cov(sol, M2, M1, M1, N) )
_std(sol::Solution, M2::Symbol, M1::Symbol, N::Symbol=:N) = sqrt.( abs.( _cov(sol, M2, M1, M1, N) ) ) #debug

function correlation(N::Number, M11::Number, M10::Number, M01::Number, M20::Number, M02::Number)
	# num = M11 - (M10 * M01 / N)
	# denom = sqrt( _var(M20, M10, N) * _var(M02, M01, N) )
	# # denom = sqrt( abs( _var(M20, M10, N) * _var(M02, M01, N) ) ) #debug
	num = _cov(M11, M10, M01, N)
	denom = _std(M20, M10, N) * _std(M02, M01, N)
	return num / denom
end
function correlation(sol::Solution, M11::Symbol, M10::Symbol, M01::Symbol, M20::Symbol, M02::Symbol, N::Symbol=:N) 
	# num = sol.M[M11] .- (sol.M[M10] .* sol.M[M01] ./ sol.M[:N]) 
	# # num = sol.M[M11] .- (sol.M[M10] .* sol.M[M01]) 
	# denom = sqrt.( _var(sol, M20, M10) .* _var(sol, M02, M01) )
	# # denom = sqrt.( abs.( _var(sol, M20, M10) .* _var(sol, M02, M01) ) ) #debug
	num = _cov(sol, M11, M10, M01, N)
	denom = _std(sol, M20, M10, N) .* _std(sol, M02, M01, N)
	return num ./ denom
end
@export function plotCorrelation(ssa::Solution, clna::Solution; 
									M11::Symbol=:M¹¹, M10::Symbol=:M¹⁰, M01::Symbol=:M⁰¹,
									M20::Symbol=:M²⁰, M02::Symbol=:M⁰²,
									color=pal[1],
									title=nothing,
									# legend=:bottomright,
									legend=:topright,
									# legend=false,
									tightLegend::Bool=false,
									ssaSD::Union{StandardDeviation,Nothing}=nothing,
									# xlabel="Time [a.u.]",
									# ylabel=nothing,
									plotFlags...
									)
	l1 = replace(string(M10), "M"=>"x")
	l2 = replace(string(M01), "M"=>"x")
	ssaCorrelation = correlation(ssa, M11, M10, M01, M20, M02)
	ssaCorrelationRibbon = nothing
	if :corr in keys(ssa.M)
		ssaCorrelation = ssa.M[:corr]
		# @show ssaSD.M[:corr]
		# @show all(isnan.(ssaSD.M[:corr]))
		ssaCorrelationRibbon = !( isnothing(ssaSD) || all(isnan.(ssaSD.M[:corr])) ) ? ssaSD.M[:corr] : nothing
	end
	clnaCorrelation = correlation(clna, M11, M10, M01, M20, M02)

	# Sanitize correlations to replace any NaN with a 0.0
	ssaCorrelation[ isnan.(ssaCorrelation) ] .= 0.0
	clnaCorrelation[ isnan.(clnaCorrelation) ] .= 0.0

	@show (M10,M01)
	corrErr = mean( norm1error(clna.T, clnaCorrelation, ssa.T, ssaCorrelation; relative=false) )
	@show corrErr

	p = plotMoment(ssa.T, ssaCorrelation;
			ribbon=ssaCorrelationRibbon,
			color=color,
			label=tightLegend ? latexstring("\\rho_{ssa}") : latexstring("\\rho("*l1*","*l2*")_{ssa}"),
			ylims=(-1.0,1.0),
			# ylims=(-Inf,Inf),
			legend=legend,
			title=title,
			plotFlags...
			)
	p = plotMoment!(p, clna.T, clnaCorrelation;
			color=color,
			linestyle=:dash,
			linewidth=2.0, 
			label=tightLegend ? latexstring("\\rho_{pred}") : latexstring("\\rho("*l1*","*l2*")_{pred}"),
			ylims=(-1.0,1.0),
			# ylims=(-Inf,Inf),
			legend=legend,
			plotFlags...
			)
	return p
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

function pushSigmaX0!(σ::SigmaDict, M::MomentDict, s2::Symbol, s1::Symbol)
	σ[:X0] = sqrt.( (M[s2] ./ M[:N]) .- (M[s1] ./ M[:N]).^2 ) # Is this the right formula?
	return σ
end
function pushSigmaX0!(S::Solution, s2::Symbol, s1::Symbol)
	pushSigmaX0!(S.σ, S.M, s2, s1)
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
							# legend= s==:N ? :bottomright : legend,
							legend= s==:N ? false : legend,
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
							ssaRibbon=isnothing(ssaSD) ? nothing : ssaSD.σ[s], # This is ok like this (no squaring necessary)
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

_ep(clnaSol::Solution, M::Symbol) = clnaSol.M[M][end] / clnaSol.M[:N][end]

@export function plot1dPopulationHistogram(ssaPool; 
											s1=1, 
											s1Label="Species $s1",											
											# aspect_ratio=0.8,
											clnaSol::Union{Solution,Nothing}=nothing, M=:M¹, 
											# todo: replace the moments with a function that returns the right symbol based on species index
											plotFlags...
											)
	p = plot()
	plot1dPopulationHistogram!(p, ssaPool;
		s1=s1, s1Label=s1Label, 
		# aspect_ratio=aspect_ratio,
		clnaSol=clnaSol, M=M, plotFlags...)
end
@export function plot1dPopulationHistogram!(p, ssaPool; 
											s1=1, 
											s1Label=latexstring("x_$s1"),
											# aspect_ratio=0.8,
											clnaSol::Union{Solution,Nothing}=nothing, M=:M¹,
											# todo: replace the moments with a function that returns the right symbol based on species index
											plotFlags...
											)
	@assert !isnothing(ssaPool)
	species = size(ssaPool, 1)
	@assert species>=s1 "Species with index $s1 does not exist!"
	h2 = histogram!(p, ssaPool[s1,:]; 
						xlabel=s1Label,
						# aspect_ratio=aspect_ratio,
						# xlims=(0,Inf), ylims=(0,Inf),
						plotFlags...
						)
	if !isnothing(clnaSol)
		# Place a vertical line where the Taylor expansion point is!
		ep = _ep(clnaSol, M)
		h2 = vline!(p, [ep,];
				color=:black, 
				linewidth=3, 
				legend=false,
				)
				# plotFlags...)
	end
	return h2
end

@export function plot2dPopulationHistogram(ssaPool; 
											s1=1, s2=2, 
											s1Label="Species $s1", s2Label="Species $s2",
											aspect_ratio=0.8,
											clnaSol::Union{Solution,Nothing}=nothing, M10=:M¹⁰, M01=:M⁰¹, 
											# todo: replace the moments with a function that returns the right symbol based on species index
											plotFlags...
											)
	p = plot()
	plot2dPopulationHistogram!(p, ssaPool;
		s1=s1, s2=s2, s1Label=s1Label, s2Label=s2Label, aspect_ratio=aspect_ratio,
		clnaSol=clnaSol, M10=M10, M01=M01, plotFlags...)
end
@export function plot2dPopulationHistogram!(p, ssaPool; 
											s1=1, s2=2, 
											s1Label=latexstring("x_$s1"), s2Label=latexstring("x_$s2"),
											aspect_ratio=0.8,
											clnaSol::Union{Solution,Nothing}=nothing, M10=:M¹⁰, M01=:M⁰¹, 
											# todo: replace the moments with a function that returns the right symbol based on species index
											plotFlags...
											)
	@assert !isnothing(ssaPool)
	species = size(ssaPool, 1)
	@assert species>=2
	histogram2d!(p, ssaPool[s1,:], ssaPool[s2,:]; 
						xlabel=s1Label,
						ylabel=s2Label,
						# xscale=:log10, yscale=:log10, # Uncomment this for log scale
						aspect_ratio=aspect_ratio,
						xlims=(0,Inf), ylims=(0,Inf),
						plotFlags...
						)
	if !isnothing(clnaSol)
		# Place a star where the Taylor expansion point is!
		ep = ( _ep(clnaSol, M10), _ep(clnaSol, M01) )
		h2 = scatter!(p, [ep,];
				markersize=20, 
				markerstrokewidth=1.5,
				markerstrokecolor=:black, # Border color
				markercolor=:white, 	  # Fill color
				markershape=:star5, 
				label=false,
				)
	end
	return h2
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
			h2 = plot2dPopulationHistogram(ssaPool; s1=1, s2=species)
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

@export function plotMomentsHistogram(M::Model, ssaTrajectories, symbol::Symbol;
									reference=nothing,
									color=nothing,
									paletteIndex=1,
									plotFlags...)
	p = plot()
	return plotMomentsHistogram!(p, M, ssaTrajectories, symbol;
									reference=reference,
									color=color,
									paletteIndex=paletteIndex,
									plotFlags...)
end

@export function plotMomentsHistogram!(p, M::Model, ssaTrajectories, symbol::Symbol;
									reference=nothing,
									color=nothing,
									paletteIndex=1,
									plotFlags...)
	pHistogram = nothing
	if !isnothing(ssaTrajectories)
		species = size(ssaTrajectories[1],1)
		e = M.ssaSystem.MomExpMapping[symbol]
		curMom = [ _moment(S, e) for S in ssaTrajectories ]

		h = histogram!(p, curMom; 
				title=string(symbol), 
				color=isnothing(color) ? pal[paletteIndex] : color, 
				legend=false, 
				plotFlags...)
		if !isnothing(reference)
			ref = reference.M[symbol][end] # Pick last timepoint
			h = vline!(p, [ref]; 
					color=:black, 
					linewidth=3, 
					legend=false,
					)
					# plotFlags...)
		end
		pHistogram = h
	end
	return pHistogram
end

@export function plotMomentsHistograms(M::Model, ssaTrajectories, symbols::Symbol...;
									reference=nothing,
									color=nothing,
									plotFlags...)
	pHistograms = nothing
	if !isnothing(ssaTrajectories)
		histograms=[]
		species = size(ssaTrajectories[1],1)
		Moments = Dict{Symbol,Vector}()
		for (i,s) in enumerate(symbols)
			e = M.ssaSystem.MomExpMapping[s]
			curMom = [ _moment(S, e) for S in ssaTrajectories ]
			Moments[s] = curMom

			h = histogram(Moments[s]; title=string(s), color=isnothing(color) ? pal[i] : color, legend=false)
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
		pHistograms = plot(histograms...; layout=(1,length(histograms)), plotFlags...)
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
	pHistograms = nothing
	try
		pCLNA = plotSolution(sol; Ms=Msymb)
		pSSA = plotSolution(reference; Ms=Msymb)
		pMean = plotMeans(reference, sol, symbols...; 
							legend=meanlegend, ssaSD=ssaSD)
		pVar = plotVariances(reference, sol, symbols...; 
							legend=varlegend, ssaSD=ssaSD)

		pHistograms = _plotPopulationHistograms(ssaPool)
		pMHistograms = plotMomentsHistograms(M, ssaTrajectories, symbols...; reference=sol)

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

@export function runSSASimulations(M::Union{Model, Vector{Model}}, symbols::Symbol...; 
						T=100.0, NSSA=100, RSSA=10,
						N0=10, Mpc0=10,
						clnaSol=nothing,
						)
	ssaSols = Vector{Solution}()
	cells = []
	trajectories = []
	for ir=1:RSSA
		ssaSol, nRaw = SSA_solve(M; T=T, N0=N0, Mpc0=Mpc0, NSSA=NSSA)
		push!(ssaSols, ssaSol)
		if !isnothing(clnaSol)
			Err = getErrorMeasure(ssaSol, 
						clnaSol, M.momentsMapping, M.sigmaMapping,
						symbols...)
			@show Err
		end
		push!(cells, hcat(nRaw...))
		push!(trajectories, nRaw...)
	end
	nPool = hcat(cells...)
	ssaMeanSol, ssaSDSol = _computeSSAStats(ssaSols)
	return nPool, trajectories, ssaMeanSol, ssaSDSol
end

#eof
