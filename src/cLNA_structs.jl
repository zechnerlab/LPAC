# Structures and types to be used with cLNA

# struct Parameters
#     kI::Number
#     kE::Number
#     kb::Number
#     ke::Number
#     λ::Number
# end

struct Model
    parameters::Union{Tuple, Dict{Symbol, T}} where T
    Ω::Number
    Ωc::Number
    momentsInit::Function
    momentsOde::Function
    ssaSystem::Sim.System
end

struct SSAsolution
    T
    N
    σN
    M
    σM
end

#eof
