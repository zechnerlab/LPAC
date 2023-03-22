#=
Utils: functions to manipulate data structures and solutions
=#

@export *(s::Symbol, a::String) = Symbol(string(s)*a)
@export occursin(needle::Union{Symbol,String,Number}, haystack::Symbol) = Base.occursin(string(needle), string(haystack))

@export function safeGet(d::Dict{K,V}, k::K)::Union{Nothing,V} where {K,V}
	try
		return d[k]
	catch KeyError
		return nothing
	end
end

rescaleCompartments(N, Ω, p=1) = N ./ Ω^p
rescaleContent(M, N, Ωc, p=1) = M ./ (N*Ωc).^p
# rescaleContent(M, N, Ωc, p=1) = M ./ (Ωc).^p
rescaleMixed(MN, N, Ω, Ωc, p=1) = rescaleCompartments(
												rescaleContent(MN, N, Ωc, p), 
												Ω, p)

# function rescaleClnaSol!(sol, Ω, Ωc)
# 	N = @view sol[1, :]
# 	N² = @view sol[2, :]
# 	M = @view sol[3, :]
# 	M² = @view sol[4, :]
# 	NM = @view sol[5, :]

# 	M .= rescaleContent(M, N, Ωc)
# 	M² .= rescaleContent(M², N, Ωc, 2)
# 	NM .= rescaleMixed(NM, N, Ω, Ωc)
# 	N² .= rescaleCompartments(N², Ω, 2)
# 	N .= rescaleCompartments(N, Ω)

# 	return sol
# end

# function rescaleSsaSol!(sol, Ω, Ωc)
# 	sol.M .= rescaleContent(sol.M, sol.N, Ωc)
# 	sol.σM .= rescaleContent(sol.σM, sol.N, Ωc)
# 	sol.σN .= rescaleCompartments(sol.σN, Ω)
# 	sol.N .= rescaleCompartments(sol.N, Ω)

# 	return sol
# end

function rescaleSolution!(sol::Solution, Ω, Ωc)
	N = sol.M[:N]
	SymbolsM = keys(sol.M)
	Symbolsσ = keys(sol.σ)
	# First rescale the other moments
	for s in setdiff(SymbolsM, [:N,:N2])
		if isnothing(sol.M[s])
			continue
		end
		p = occursin(2, s)||occursin("¹¹", s) ? 2 : 1
		if !occursin(:N, s)
			sol.M[s] .= rescaleContent(sol.M[s], N, Ωc, p)
		else
			sol.M[s] .= rescaleMixed(sol.M[s], N, Ω, Ωc, p)
		end
	end
	for s in setdiff(Symbolsσ, [:N,:N2])
		if isnothing(sol.σ[s])
			continue
		end
		p = occursin(2, s)||occursin("¹¹", s) ? 2 : 1
		if !occursin(:N, s)
			sol.σ[s] .= rescaleContent(sol.σ[s], N, Ωc, p)
		else
			sol.σ[s] .= rescaleMixed(sol.σ[s], N, Ω, Ωc, p)
		end
	end
	# Then rescale the compartments
	sol.M[:N] .= rescaleCompartments(sol.M[:N], Ω)
	sol.σ[:N] .= rescaleCompartments(sol.σ[:N], Ω)
	if :N2 in SymbolsM
		sol.M[:N2] .= rescaleCompartments(sol.M[:N], Ω, 2)
	end
	if :N2 in Symbolsσ
		sol.σ[:N2] .= rescaleCompartments(sol.σ[:N], Ω, 2)
	end
	return sol
end

@export function solcat(S::Solution...)
	Tcat = deepcopy(S[1].T)
	Mcat = deepcopy(S[1].M)
	σcat = deepcopy(S[1].σ)
	
	for s in S[2:end]
		T, M, σ = s.T, s.M, s.σ
		skip = 0
		if T[1] < Tcat[end]	# If time of current solution started from zero
			T .+= Tcat[end] # then offset it by the end time of the previous one
		end
		if T[1] == Tcat[end] # If the time of cur solution starts exactly at the end time of prev
			skip = 1
		end
		Tcat = vcat(Tcat, T[skip:end])
		for k in keys(Mcat)
			@assert k in keys(M) "Moment $k not available in all solutions that are being concatenated!"
			Mcat[k] = vcat(Mcat[k], M[k][skip:end])
		end
		for k in keys(σcat)
			@assert k in keys(σ) "σ of moment $k not available in all solutions that are being concatenated!"
			if !isnothing(σcat[k])
				σcat[k] = vcat(σcat[k], σ[k][skip:end])
			end
		end
	end
	return Solution(Tcat, Mcat, σcat)
end

"""
	norm1error(X,Y, Xref,Yref; 
				relative::Bool=false)

Compute the norm-1 error of the X,Y data series against the reference Xref,Yref.
The data series is linearly interpolated at the reference points (Xref).
"""
@export function norm1error(X,Y, Xref,Yref; relative::Bool=false)
	try
		I = linear_interpolation(X,Y)
		diff = I.(Xref) .- Yref
		err = diff
		if relative
			err = diff ./ Yref
			err[isnan.(err)] .= 0.0
		end
		return abs.( err )
	catch e
		println("Warning: Cannot compute the norm1error due to exception (probably ODE solution failed), returning NaN")
		# @show e
		return NaN
	end
end

#eof
