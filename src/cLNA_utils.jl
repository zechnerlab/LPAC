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

#eof
