#=
	Symbolic manipulations (quick&dirty) for deriving the cLNA
=#

@export struct Reaction
	ΔN::Num
	ΔM::Num
	h::Num
end

@export function linearize(f::Num, x::Num, x₀::Num)
	D = Differential(x)
	fx₀ = substitute(f, Dict(x=>x₀))
	dfx₀ = substitute(D(f), Dict(x=>x₀))
	return fx₀ + dfx₀*(x - x₀)
end

@export function derive_dN(reactions, x, moments)
	@variables _N, _X
	substitutions = [ _N*_X^(i-1)=>moments[i] for i=length(moments):-1:1 ]
	append!(substitutions, [ moments[1]*_X^(i-1)=>moments[i] for i=length(moments):-1:1 ])
	dN = Num(0) # Initialize
	for R in reactions
		cur = R.ΔN * R.h * _N
		cur = substitute(cur, Dict(x=>_X))
		@show cur
		dN += simplify(cur; rewriter=substitute(cur, Dict(substitutions)))
	end
	return dN
end

#eof
