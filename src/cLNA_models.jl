#=
Models to be used with the cLNA code, comparing cLNA vs. SSA
=#

##################### IntakeExitCoagulation²FragmentationBirthDeath

@export function IECpqFBD(; 
                            kI=5.0, kE=0.1, 
                            # kC=5e-6, kF=5e-2, 
                            kC=5e-5, kF=5e-3, 
                            kb=10.0, kd=0.1, 
                            λ=10.0, 
                            Ω=1.0, Ωc=1.0,
                            )
    kI, kE, kC, kF, kb, kd, λ, Ω, Ωc = float.( (kI, kE, kC, kF, kb, kd, λ, Ω, Ωc) )
    kI *= Ω; kC /= Ω*Ωc^2
    kb *= Ωc; kF /= Ωc
    p = Dict(:kI=>kI, :kE=>kE, :kC=>kC, :kF=>kF, :kb=>kb, :kd=>kd, :λ=>λ, 
                :Ω=>Ω, :Ωc=>Ωc,
                )
    ssaSys = getIECpqFBDSSASystem(; kI=kI, kE=kE, kC=kC, kF=kF, kb=kb, kd=kd, λ=λ)
    return Model(p, Ω, Ωc, IECpFBD_initial, IECpqFBD_ODEs, ssaSys)
end

function getIECpqFBDSSASystem(;
                            kI = 0.0, kE = 0.0,  # Intake/Exit rates for cells
                            kC = 0.0, kF = 0.0,  # Coag/Frag rates for cells
                            kb = 0.0, kd = 0.0,  # Birth/Death rates for X1
                            λ = 0,  # Parameter of the poisson distrib on intake
                            )
    S = Sim.System("IECpqFBD", 1, # X1
            Dict(1 =>[0], # Number of compartments
                 2 =>[1], # M¹
                 3 =>[2], # M²
                )
            )

    # cells' reaction network
    prodX1 = Sim.new_chemical_reaction_class([1], kb)
    prodX1.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[1]
    prodX1.fast_sample_reactants! = fast_sample_uniform_cell

    deathX1 = Sim.new_chemical_reaction_class([-1], kd)
    deathX1.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[2]
    deathX1.fast_sample_reactants! = fast_sample_mass_x1

    # Compartment reactions
    cellIntake = Sim.TransitionClass(0, 1, kI)
    cellIntake.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> 1
    # This is how the output compartment (yc) of the transition looks like
    cellIntake.pi = function(yc::Vector{Vector{Int64}}, xc::Vector{Vector{Int64}})
                        yc[1][1] = rand(Poisson(λ))
                    end

    cellExit = Sim.TransitionClass(1, 0, kE) # This is a cell death
    cellExit.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[1] # Depends on num cells
    cellExit.fast_sample_reactants! = fast_sample_uniform_cell

    coagulation = Sim.TransitionClass(2, 1, kC)
    # coagulation.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[3]*(Mom[1]-1) # NM2-M2
    coagulation.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[3]*(Mom[1]-2) + Mom[2]^2 # M²(N-2) + (M¹)²
    coagulation.fast_sample_reactants! = 
        (r_indices::Vector{Int64}, n::Matrix{Int64}, Mom::Vector{Int64}) -> 
            fast_sample_generic_binary(r_indices, n, Mom, 
                                        (x,y)->(x[1]+y[1])^2, 
                                        ()->coagulation.H(n, Mom)
                                        )
    coagulation.pi = function(yc::Vector{Vector{Int64}}, xc::Vector{Vector{Int64}})
                            yc[1][:] .= xc[1][:] .+ xc[2][:]
                        end

    fragmentation = Sim.TransitionClass(1, 2, kF)
    fragmentation.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[2]
    fragmentation.fast_sample_reactants! = fast_sample_mass_x1
    fragmentation.pi = function(yc::Vector{Vector{Int64}}, xc::Vector{Vector{Int64}})
                        yc[1][1] = rand(0:xc[1][1])
                        yc[2][1] = xc[1][1] - yc[1][1]
                    end

    Sim.add_transition_class(S, prodX1, deathX1, 
        cellIntake, cellExit, coagulation, fragmentation)

    return S
end

function IECpqFBD_ODEs(dM, M, parameters, t)
    c6 = parameters[:λ] # λ
    c2 = parameters[:kC] # kC
    c1 = parameters[:kE] # kE
    c3 = parameters[:kF] # kF
    c0 = parameters[:kI] # kI
    c4 = parameters[:kb] # kb
    c5 = parameters[:kd] # kd
    # Number of Compartments (N)
    dM[1] = c0+c3*M[3]-c1*M[1]-2*c2*M[3]^2+2*c2*M[3]^2/M[1]
    # N^2
    dM[2] = (M[1]^2*(c0+c1*M[1]+c3*M[3]-2*c1*M[2]+2*c2*M[3]^2+2*c3*M[5]+2*c0*M[1]-8*c2*M[5]*M[3]+4*c2*M[3]^2*M[1])-4*c2*M[3]^2*M[2]+2*c2*(-M[3]+4*M[5])*M[1]*M[3])/M[1]^2
    # Total Mass
    dM[3] = c6*c0+c4*M[1]-c1*M[3]-c5*M[3]
    # M^2
    dM[4] = c6*c0+c0*c6^2+c4*M[1]+c5*M[3]-2*c1*M[4]-2*c5*M[4]+2*c4*M[5]+c1*M[3]^2/M[1]+2*c6*c0*M[3]
    # N*M
    dM[5] = c6*c0+c1*M[3]+c3*M[4]+c0*M[3]+c4*M[2]-c5*M[5]-2*c1*M[5]+2*c2*M[3]^3+c6*c0*M[1]-4*c2*M[4]*M[3]-2*c2*M[3]^2*M[5]/M[1]^2+4*c2*M[4]*M[3]/M[1]
    return
end

function IECpqFBD_ODEs_partiallyClosed(dM, M, parameters, t)
    c4 = parameters[:λ] # λ
    c1 = parameters[:kC] # kC
    c6 = parameters[:kE] # kE
    c2 = parameters[:kF] # kF
    c5 = parameters[:kI] # kI
    c0 = parameters[:kb] # kb
    c3 = parameters[:kd] # kd
    # Number of Compartments (N)
    dM[1] = (M[1]^2*(c5+c2*M[3]-c6*M[1])+2*c1*M[3]^2*M[2]+2*c1*(-2*M[5]+M[3])*M[1]*M[3])/M[1]^2
    # N^2
    dM[2] = (M[1]^3*(c5+c6*M[1]+c2*M[3]-2*c6*M[2]+2*c2*M[5]+2*c5*M[1])-4*c1*M[3]^2*(M[1]^2-2*M[2])*M[2]+2*c1*M[1]^2*(-M[3]+6*M[5])*M[3]+2*c1*(-8*M[5]-3*M[3]+4*M[1]*M[3])*M[2]*M[1]*M[3])/M[1]^3
    # Total Mass
    dM[3] = c4*c5+c0*M[1]-c6*M[3]-c3*M[3]
    # M^2
    dM[4] = c4*c5+c5*c4^2+c0*M[1]+c3*M[3]-2*c6*M[4]-2*c3*M[4]+2*c0*M[5]+c6*M[3]^2/M[1]+2*c4*c5*M[3]
    # N*M
    dM[5] = (M[1]^3*(c4*c5+c6*M[3]+c2*M[4]+c5*M[3]+c0*M[2]-c3*M[5]-2*c6*M[5]+c4*c5*M[1])-2*c1*M[3]^2*M[5]*M[1]+2*c1*M[3]^2*(2*M[5]-M[1]*M[3])*M[2]+4*c1*M[1]^2*(-2*M[5]+M[1]*M[3]+M[3])*M[4])/M[1]^3
    return
end

##################### IntakeExitCoagulation²FragmentationBirthDeath

@export function IECpFBD(; kI=5.0, kE=0.1, kC=0.005/5, kF=0.005, kb=10.0, kd=5*0.1, λ=10.0, 
                            Ω=1.0, Ωc=1.0,
                            )
    kI, kE, kC, kF, kb, kd, λ, Ω, Ωc = float.( (kI, kE, kC, kF, kb, kd, λ, Ω, Ωc) )
    kI *= Ω; kC /= Ω*Ωc
    kb *= Ωc; kF /= Ωc
    p = Dict(:kI=>kI, :kE=>kE, :kC=>kC, :kF=>kF, :kb=>kb, :kd=>kd, :λ=>λ, 
                :Ω=>Ω, :Ωc=>Ωc,
                )
    ssaSys = getIECpFBDSSASystem(; kI=kI, kE=kE, kC=kC, kF=kF, kb=kb, kd=kd, λ=λ)
    return Model(p, Ω, Ωc, IECpFBD_initial, IECpFBD_ODEs, ssaSys)
end

function getIECpFBDSSASystem(;
                            kI = 0.0, kE = 0.0,  # Intake/Exit rates for cells
                            kC = 0.0, kF = 0.0,  # Coag/Frag rates for cells
                            kb = 0.0, kd = 0.0,  # Birth/Death rates for X1
                            λ = 0,  # Parameter of the poisson distrib on intake
                            )
    S = Sim.System("IECpFBD", 1, # X1
            Dict(1 =>[0], # Number of compartments
                 2 =>[1], # M¹
                 3 =>[2], # M²
                )
            )

    # cells' reaction network
    prodX1 = Sim.new_chemical_reaction_class([1], kb)
    prodX1.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[1]
    prodX1.fast_sample_reactants! = fast_sample_uniform_cell

    deathX1 = Sim.new_chemical_reaction_class([-1], kd)
    deathX1.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[2]
    deathX1.fast_sample_reactants! = fast_sample_mass_x1

    # Compartment reactions
    cellIntake = Sim.TransitionClass(0, 1, kI)
    cellIntake.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> 1
    # This is how the output compartment (yc) of the transition looks like
    cellIntake.pi = function(yc::Vector{Vector{Int64}}, xc::Vector{Vector{Int64}})
                        yc[1][1] = rand(Poisson(λ))
                    end
    # cell_division.fast_sample_reactants! = fast_sample_cell_div

    cellExit = Sim.TransitionClass(1, 0, kE) # This is a cell death
    cellExit.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[1] # Depends on num cells
    cellExit.fast_sample_reactants! = fast_sample_uniform_cell

    coagulation = Sim.TransitionClass(2, 1, kC)
    coagulation.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[2]*(Mom[1]-1) # NM1-M1
    coagulation.fast_sample_reactants! = 
        (r_indices::Vector{Int64}, n::Matrix{Int64}, Mom::Vector{Int64}) -> 
            fast_sample_generic_binary(r_indices, n, Mom, (x,y)->x[1]+y[1], ()->coagulation.H(n, Mom))
    coagulation.pi = function(yc::Vector{Vector{Int64}}, xc::Vector{Vector{Int64}})
                            yc[1][:] .= xc[1][:] .+ xc[2][:]
                        end

    fragmentation = Sim.TransitionClass(1, 2, kF)
    fragmentation.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[2]
    fragmentation.fast_sample_reactants! = fast_sample_mass_x1
    fragmentation.pi = function(yc::Vector{Vector{Int64}}, xc::Vector{Vector{Int64}})
                        yc[1][1] = rand(0:xc[1][1])
                        yc[2][1] = xc[1][1] - yc[1][1]
                    end

    Sim.add_transition_class(S, prodX1, deathX1, 
        cellIntake, cellExit, coagulation, fragmentation)

    return S
end

function IECpFBD_ODEs(dM, M, parameters, t)
    c5 = parameters[:λ] # λ
    c2 = parameters[:kC] # kC
    c4 = parameters[:kE] # kE
    c3 = parameters[:kF] # kF
    c0 = parameters[:kI] # kI
    c6 = parameters[:kb] # kb
    c1 = parameters[:kd] # kd
    # Number of Compartments (N)
    dM[1] = c0+c2*M[3]+c3*M[3]-c4*M[1]-c2*M[1]*M[3]
    # N^2
    dM[2] = c0+c4*M[1]+c3*M[3]-c2*M[3]-2*c4*M[2]+2*c2*M[5]+2*c3*M[5]+2*c0*M[1]+c2*M[1]*M[3]-2*c2*M[2]*M[3]-2*c2*M[5]*M[1]+2*c2*M[1]^2*M[3]
    # Total Mass
    dM[3] = c5*c0+c6*M[1]-c4*M[3]-c1*M[3]
    # M^2
    dM[4] = c5*c0+c0*c5^2+c6*M[1]+c1*M[3]-2*c4*M[4]-2*c1*M[4]+2*c6*M[5]+c4*M[3]^2/M[1]+2*c5*c0*M[3]
    # N*M
    dM[5] = c5*c0+c2*M[4]+c4*M[3]+c3*M[4]+c0*M[3]+c6*M[2]-c1*M[5]-2*c4*M[5]+c5*c0*M[1]+c2*M[3]^2*M[1]-c2*M[4]*M[1]-c2*M[5]*M[3]
    return
end

# initialize expected moments vector
function IECpFBD_initial(N0, Mpc0)
    M0 = N0*Mpc0
    M = zeros(5)
    # Number of Compartments (N)
    M[1] = N0 # initial value for Moment(0,), please specify!
    # N^2
    M[2] = M[1]^2
    # Total Mass
    M[3] = M0 # initial value for Moment(1,), please specify!
    # M^2
    M[4] = M[3]^2
    # N*M
    M[5] = M[1]*M[3]
    return M
end

##################### IntakeExitCoagulationFragmentation²BirthDeath

@export function IECFqBD(; kI=5.0, kE=0.1, kC=0.005/5, kF=0.005, kb=10.0, kd=5*0.1, λ=10.0, 
                            Ω=1.0, Ωc=1.0,
                            )
    kI, kE, kC, kF, kb, kd, λ, Ω, Ωc = float.( (kI, kE, kC, kF, kb, kd, λ, Ω, Ωc) )
    kI *= Ω; kC /= Ω
    kb *= Ωc; kF /= Ωc^2
    p = Dict(:kI=>kI, :kE=>kE, :kC=>kC, :kF=>kF, :kb=>kb, :kd=>kd, :λ=>λ, 
                :Ω=>Ω, :Ωc=>Ωc,
                )
    ssaSys = getIECFqBDSSASystem(; kI=kI, kE=kE, kC=kC, kF=kF, kb=kb, kd=kd, λ=λ)
    return Model(p, Ω, Ωc, IECFqBD_initial, IECFqBD_ODEs, ssaSys)
end

function getIECFqBDSSASystem(;
                            kI = 0.0, kE = 0.0,  # Intake/Exit rates for cells
                            kC = 0.0, kF = 0.0,  # Coag/Frag rates for cells
                            kb = 0.0, kd = 0.0,  # Birth/Death rates for X1
                            λ = 0,  # Parameter of the poisson distrib on intake
                            )
    S = Sim.System("IECFqBD", 1, # X1
            Dict(1 =>[0], # Number of compartments
                 2 =>[1], # M¹
                 3 =>[2], # M²
                )
            )

    # cells' reaction network
    prodX1 = Sim.new_chemical_reaction_class([1], kb)
    prodX1.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[1]
    prodX1.fast_sample_reactants! = fast_sample_uniform_cell

    deathX1 = Sim.new_chemical_reaction_class([-1], kd)
    deathX1.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[2]
    deathX1.fast_sample_reactants! = fast_sample_mass_x1

    # Compartment reactions
    cellIntake = Sim.TransitionClass(0, 1, kI)
    cellIntake.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> 1
    # This is how the output compartment (yc) of the transition looks like
    cellIntake.pi = function(yc::Vector{Vector{Int64}}, xc::Vector{Vector{Int64}})
                        yc[1][1] = rand(Poisson(λ))
                    end
    # cell_division.fast_sample_reactants! = fast_sample_cell_div

    cellExit = Sim.TransitionClass(1, 0, kE) # This is a cell death
    cellExit.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[1] # Depends on num cells
    cellExit.fast_sample_reactants! = fast_sample_uniform_cell

    coagulation = Sim.TransitionClass(2, 1, kC)
    coagulation.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[1]*(Mom[1]-1)/2
    coagulation.fast_sample_reactants! = fast_sample_uniform_2cells
    coagulation.pi = function(yc::Vector{Vector{Int64}}, xc::Vector{Vector{Int64}})
                            yc[1][:] .= xc[1][:] .+ xc[2][:]
                        end

    fragmentation = Sim.TransitionClass(1, 2, kF)
    fragmentation.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> div(Mom[3]-Mom[2], 2)
    fragmentation.fast_sample_reactants! = fast_sample_x1_square
    fragmentation.pi = function(yc::Vector{Vector{Int64}}, xc::Vector{Vector{Int64}})
                        yc[1][1] = rand(0:xc[1][1])
                        yc[2][1] = xc[1][1] - yc[1][1]
                    end

    Sim.add_transition_class(S, prodX1, deathX1, 
        cellIntake, cellExit, coagulation, fragmentation)

    return S
end

function IECFqBD_ODEs(dM, M, parameters, t)
    c0 = parameters[:λ] # λ
    c5 = parameters[:kC] # kC
    c6 = parameters[:kE] # kE
    c2 = parameters[:kF] # kF
    c1 = parameters[:kI] # kI
    c4 = parameters[:kb] # kb
    c3 = parameters[:kd] # kd
    # Number of Compartments (N)
    dM[1] = c1+1/2*c5*M[1]-c6*M[1]-1/2*c5*M[1]^2-1/2*c2*M[3]+1/2*c2*M[3]^2/M[1]
    # N^2
    dM[2] = c1+c5*M[1]^3+c5*M[2]+c6*M[1]+1/2*c5*M[1]^2-c2*M[5]-2*c6*M[2]+2*c1*M[1]-1/2*c5*M[1]-1/2*c2*M[3]+1/2*c2*M[3]^2/M[1]-2*c5*M[2]*M[1]-c2*M[3]^2*M[2]/M[1]^2+2*c2*M[5]*M[3]/M[1]
    # Total Mass
    dM[3] = c0*c1+c4*M[1]-c6*M[3]-c3*M[3]
    # M^2
    dM[4] = c0*c1+c1*c0^2+c4*M[1]+c3*M[3]-2*c6*M[4]-2*c3*M[4]+2*c4*M[5]+c6*M[3]^2/M[1]+2*c0*c1*M[3]
    # N*M
    dM[5] = c0*c1+c6*M[3]+c1*M[3]+c4*M[2]+1/2*c5*M[5]-c3*M[5]-2*c6*M[5]-1/2*c2*M[4]+c0*c1*M[1]+1/2*c5*M[1]^2*M[3]-c5*M[5]*M[1]+c2*M[4]*M[3]/M[1]-1/2*c2*M[3]^2*M[5]/M[1]^2
    return
end

# initialize expected moments vector
function IECFqBD_initial(N0, Mpc0)
    M0 = N0*Mpc0
    M = zeros(5)
    # Number of Compartments (N)
    M[1] = N0 # initial value for Moment(0,), please specify!
    # N^2
    M[2] = M[1]^2
    # Total Mass
    M[3] = M0 # initial value for Moment(1,), please specify!
    # M^2
    M[4] = M[3]^2
    # N*M
    M[5] = M[1]*M[3]
    return M
end

##################### IntakeExitCoagulationFragmentationBirthDeath³

@export function IECFBDc(; kI=5.0, kE=0.1, kC=0.005, kF=0.005, kb=10.0, kd=0.1, λ=10.0, 
                            Ω=1.0, Ωc=1.0,
                            )
    kI, kE, kC, kF, kb, kd, λ, Ω, Ωc = float.( (kI, kE, kC, kF, kb, kd, λ, Ω, Ωc) )
    kI *= Ω; kC /= Ω
    kb *= Ωc; kd /= Ωc^2; kF /= Ωc
    p = Dict(:kI=>kI, :kE=>kE, :kC=>kC, :kF=>kF, :kb=>kb, :kd=>kd, :λ=>λ, 
                :Ω=>Ω, :Ωc=>Ωc,
                )
    ssaSys = getIECFBDcSSASystem(; kI=kI, kE=kE, kC=kC, kF=kF, kb=kb, kd=kd, λ=λ)
    return Model(p, Ω, Ωc, IECFBDq_initial, IECFBDc_ODEs, ssaSys)
end

function getIECFBDcSSASystem(;
                            kI = 0.0, kE = 0.0,  # Intake/Exit rates for cells
                            kC = 0.0, kF = 0.0,  # Coag/Frag rates for cells
                            kb = 0.0, kd = 0.0,  # Birth/Death rates for X1
                            λ = 0,  # Parameter of the poisson distrib on intake
                            )
    S = Sim.System("IECFBDc", 1, # X1
            Dict(1 =>[0], # Number of compartments
                 2 =>[1], # M¹
                 3 =>[2], # M²
                 4 =>[3], # M³
                )
            )

    # cells' reaction network
    prodX1 = Sim.new_chemical_reaction_class([1], kb)
    prodX1.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[1]
    prodX1.fast_sample_reactants! = fast_sample_uniform_cell

    deathX1 = Sim.new_chemical_reaction_class([-3], kd)
    deathX1.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> div(Mom[4] - 3*Mom[3] + 2*Mom[2], 6)
    deathX1.fast_sample_reactants! = fast_sample_x1_cube

    # Compartment reactions
    cellIntake = Sim.TransitionClass(0, 1, kI)
    cellIntake.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> 1
    # This is how the output compartment (yc) of the transition looks like
    cellIntake.pi = function(yc::Vector{Vector{Int64}}, xc::Vector{Vector{Int64}})
                        yc[1][1] = rand(Poisson(λ))
                    end
    # cell_division.fast_sample_reactants! = fast_sample_cell_div

    cellExit = Sim.TransitionClass(1, 0, kE) # This is a cell death
    cellExit.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[1] # Depends on num cells
    cellExit.fast_sample_reactants! = fast_sample_uniform_cell

    coagulation = Sim.TransitionClass(2, 1, kC)
    coagulation.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[1]*(Mom[1]-1)/2
    coagulation.fast_sample_reactants! = fast_sample_uniform_2cells
    coagulation.pi = function(yc::Vector{Vector{Int64}}, xc::Vector{Vector{Int64}})
                            yc[1][:] .= xc[1][:] .+ xc[2][:]
                        end

    fragmentation = Sim.TransitionClass(1, 2, kF)
    fragmentation.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[2]
    fragmentation.fast_sample_reactants! = fast_sample_mass_x1
    fragmentation.pi = function(yc::Vector{Vector{Int64}}, xc::Vector{Vector{Int64}})
                        yc[1][1] = rand(0:xc[1][1])
                        yc[2][1] = xc[1][1] - yc[1][1]
                    end

    Sim.add_transition_class(S, prodX1, deathX1, 
        cellIntake, cellExit, coagulation, fragmentation)

    return S
end

function IECFBDc_ODEs(dM, M, parameters, t)
    c0 = parameters[:kC] # k_C, please specify!
    c5 = parameters[:kC] # kC
    c6 = parameters[:kE] # kE
    c2 = parameters[:kF] # kF
    c1 = parameters[:kI] # kI
    c4 = parameters[:kb] # kb
    c3 = parameters[:kd] # kd
    # Number of Compartments (N)
    dM[1] = c1+c2*M[3]+1/2*c5*M[1]-c6*M[1]-1/2*c5*M[1]^2
    # N^2
    dM[2] = c1+c5*M[1]^3+c5*M[2]+c6*M[1]+c2*M[3]+1/2*c5*M[1]^2-2*c6*M[2]+2*c2*M[5]+2*c1*M[1]-1/2*c5*M[1]-2*c5*M[2]*M[1]
    # Total Mass
    dM[3] = c0*c1+c4*M[1]-c6*M[3]-c3*M[3]-1/2*c3*M[3]^3/M[1]^2+3/2*c3*M[3]^2/M[1]
    # M^2
    dM[4] = 1/2*(2*M[1]^3*(c0*c1+c1*c0^2+c4*M[1]-2*c6*M[4]-2*c3*M[4]+2*c4*M[5]+3*c3*M[3]+2*c0*c1*M[3])+M[1]^2*(-9*c3*M[3]+2*c6*M[3]+12*c3*M[4])*M[3]+4*c3*M[3]^3*M[5]+3*c3*M[3]^2*(-2*M[4]-2*M[5]+M[3])*M[1])/M[1]^3
    # N*M
    dM[5] = 1/2*(M[1]^3*(c5*M[5]-4*c6*M[5]-2*c3*M[5]+2*c0*c1+2*c6*M[3]+2*c2*M[4]+2*c1*M[3]+2*c4*M[2]+c5*M[1]^2*M[3]-2*c5*M[5]*M[1]+2*c0*c1*M[1])+2*c3*M[3]^3*M[2]-3*c3*M[3]^2*(M[2]+M[5])*M[1]+6*c3*M[1]^2*M[5]*M[3])/M[1]^3
    return
end

##################### IntakeExitCoagulationFragmentationBirthDeath²

@export function IECFBDq(; kI=5.0, kE=0.1, kC=0.005, kF=0.005, kb=10.0, kd=0.1, λ=10.0, 
                            Ω=1.0, Ωc=1.0,
                            )
    kI, kE, kC, kF, kb, kd, λ, Ω, Ωc = float.( (kI, kE, kC, kF, kb, kd, λ, Ω, Ωc) )
    kI *= Ω; kC /= Ω
    kb *= Ωc; kd /= Ωc; kF /= Ωc
    p = Dict(:kI=>kI, :kE=>kE, :kC=>kC, :kF=>kF, :kb=>kb, :kd=>kd, :λ=>λ, 
                :Ω=>Ω, :Ωc=>Ωc,
                )
    ssaSys = getIECFBDqSSASystem(; kI=kI, kE=kE, kC=kC, kF=kF, kb=kb, kd=kd, λ=λ)
    return Model(p, Ω, Ωc, IECFBDq_initial, IECFBDq_ODEs, ssaSys)
end

function getIECFBDqSSASystem(;
                            kI = 0.0, kE = 0.0,  # Intake/Exit rates for cells
                            kC = 0.0, kF = 0.0,  # Coag/Frag rates for cells
                            kb = 0.0, kd = 0.0,  # Birth/Death rates for X1
                            λ = 0,  # Parameter of the poisson distrib on intake
                            )
    S = Sim.System("IECFBDq", 1, # X1
            Dict(1 =>[0], # Number of compartments
                 2 =>[1], # M¹
                 3 =>[2], # M²
                )
            )

    # cells' reaction network
    prodX1 = Sim.new_chemical_reaction_class([1], kb)
    prodX1.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[1]
    prodX1.fast_sample_reactants! = fast_sample_uniform_cell

    deathX1 = Sim.new_chemical_reaction_class([-2], kd)
    deathX1.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> div(Mom[3]-Mom[2], 2)
    deathX1.fast_sample_reactants! = fast_sample_x1_square

    # Compartment reactions
    cellIntake = Sim.TransitionClass(0, 1, kI)
    cellIntake.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> 1
    # This is how the output compartment (yc) of the transition looks like
    cellIntake.pi = function(yc::Vector{Vector{Int64}}, xc::Vector{Vector{Int64}})
                        yc[1][1] = rand(Poisson(λ))
                    end
    # cell_division.fast_sample_reactants! = fast_sample_cell_div

    cellExit = Sim.TransitionClass(1, 0, kE) # This is a cell death
    cellExit.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[1] # Depends on num cells
    cellExit.fast_sample_reactants! = fast_sample_uniform_cell

    coagulation = Sim.TransitionClass(2, 1, kC)
    coagulation.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[1]*(Mom[1]-1)/2
    coagulation.fast_sample_reactants! = fast_sample_uniform_2cells
    coagulation.pi = function(yc::Vector{Vector{Int64}}, xc::Vector{Vector{Int64}})
                            yc[1][:] .= xc[1][:] .+ xc[2][:]
                        end

    fragmentation = Sim.TransitionClass(1, 2, kF)
    fragmentation.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[2]
    fragmentation.fast_sample_reactants! = fast_sample_mass_x1
    fragmentation.pi = function(yc::Vector{Vector{Int64}}, xc::Vector{Vector{Int64}})
                        yc[1][1] = rand(0:xc[1][1])
                        yc[2][1] = xc[1][1] - yc[1][1]
                    end

    Sim.add_transition_class(S, prodX1, deathX1, 
        cellIntake, cellExit, coagulation, fragmentation)

    return S
end

function IECFBDq_ODEs(dM, M, parameters, t)
    c0 = parameters[:kC] # k_C, please specify!
    c1 = parameters[:kd] # k_d, please specify!
    c2 = parameters[:λ] # \lambda, please specify!
    c3 = parameters[:kI] # k_I, please specify!
    c4 = parameters[:kb] # k_b, please specify!
    c5 = parameters[:kE] # k_E, please specify!
    c6 = parameters[:kF] # k_F, please specify!
    # Number of Compartments (N)
    dM[1] = c3+c6*M[3]+1/2*c0*M[1]-c5*M[1]-1/2*c0*M[1]^2
    # N^2
    dM[2] = c3+c0*M[1]^3+c0*M[2]+c5*M[1]+c6*M[3]+1/2*c0*M[1]^2-2*c5*M[2]+2*c6*M[5]+2*c3*M[1]-1/2*c0*M[1]-2*c0*M[2]*M[1]
    # Total Mass
    dM[3] = c2*c3+c4*M[1]+c1*M[3]-c5*M[3]-c1*M[3]^2/M[1]
    # M^2
    dM[4] = (M[1]^2*(c2*c3+c3*c2^2+c4*M[1]-2*c5*M[4]-2*c1*M[3]+2*c4*M[5]+2*c1*M[4]+2*c2*c3*M[3])+(c5*M[3]-4*c1*M[4]+2*c1*M[3])*M[1]*M[3]+2*c1*M[3]^2*M[5])/M[1]^2
    # N*M
    dM[5] = c2*c3+c5*M[3]+c6*M[4]+c3*M[3]+c4*M[2]+c1*M[5]+1/2*c0*M[5]-2*c5*M[5]+c2*c3*M[1]+1/2*c0*M[1]^2*M[3]-c0*M[5]*M[1]+c1*M[3]^2*M[2]/M[1]^2-2*c1*M[5]*M[3]/M[1]
    return
end

# initialize expected moments vector
function IECFBDq_initial(N0, Mpc0)
    M0 = N0*Mpc0
    M = zeros(5)
    # Number of Compartments (N)
    M[1] = N0 # initial value for Moment(0,), please specify!
    # N^2
    M[2] = M[1]^2
    # Total Mass
    M[3] = M0 # initial value for Moment(1,), please specify!
    # M^2
    M[4] = M[3]^2
    # N*M
    M[5] = M[1]*M[3]
    return M
end

##################### IntakeExitCoagulationFragmentationBirthDeath

@export function IECFBD(; kI=5.0, kE=0.1, kC=0.005, kF=0.005, kb=10.0, kd=5*0.1, λ=10.0, 
                            Ω=1.0, Ωc=1.0,
                            )
    kI, kE, kC, kF, kb, kd, λ, Ω, Ωc = float.( (kI, kE, kC, kF, kb, kd, λ, Ω, Ωc) )
    kI *= Ω; kC /= Ω
    kb *= Ωc; kF /= Ωc
    p = Dict(:kI=>kI, :kE=>kE, :kC=>kC, :kF=>kF, :kb=>kb, :kd=>kd, :λ=>λ, 
                :Ω=>Ω, :Ωc=>Ωc,
                )
    ssaSys = getIECFBDSSASystem(; kI=kI, kE=kE, kC=kC, kF=kF, kb=kb, kd=kd, λ=λ)
    return Model(p, Ω, Ωc, IECFBD_initial, IECFBD_ODEs, ssaSys)
end

function getIECFBDSSASystem(;
                            kI = 0.0, kE = 0.0,  # Intake/Exit rates for cells
                            kC = 0.0, kF = 0.0,  # Coag/Frag rates for cells
                            kb = 0.0, kd = 0.0,  # Birth/Death rates for X1
                            λ = 0,  # Parameter of the poisson distrib on intake
                            )
    S = Sim.System("IECFBD", 1, # X1
            Dict(1 =>[0], # Number of compartments
                 2 =>[1], # M¹
                 3 =>[2], # M²
                )
            )

    # cells' reaction network
    prodX1 = Sim.new_chemical_reaction_class([1], kb)
    prodX1.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[1]
    prodX1.fast_sample_reactants! = fast_sample_uniform_cell

    deathX1 = Sim.new_chemical_reaction_class([-1], kd)
    deathX1.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[2]
    deathX1.fast_sample_reactants! = fast_sample_mass_x1

    # Compartment reactions
    cellIntake = Sim.TransitionClass(0, 1, kI)
    cellIntake.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> 1
    # This is how the output compartment (yc) of the transition looks like
    cellIntake.pi = function(yc::Vector{Vector{Int64}}, xc::Vector{Vector{Int64}})
                        yc[1][1] = rand(Poisson(λ))
                    end
    # cell_division.fast_sample_reactants! = fast_sample_cell_div

    cellExit = Sim.TransitionClass(1, 0, kE) # This is a cell death
    cellExit.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[1] # Depends on num cells
    cellExit.fast_sample_reactants! = fast_sample_uniform_cell

    coagulation = Sim.TransitionClass(2, 1, kC)
    coagulation.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[1]*(Mom[1]-1)/2
    coagulation.fast_sample_reactants! = fast_sample_uniform_2cells
    coagulation.pi = function(yc::Vector{Vector{Int64}}, xc::Vector{Vector{Int64}})
                            yc[1][:] .= xc[1][:] .+ xc[2][:]
                        end

    fragmentation = Sim.TransitionClass(1, 2, kF)
    fragmentation.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[2]
    fragmentation.fast_sample_reactants! = fast_sample_mass_x1
    fragmentation.pi = function(yc::Vector{Vector{Int64}}, xc::Vector{Vector{Int64}})
                        yc[1][1] = rand(0:xc[1][1])
                        yc[2][1] = xc[1][1] - yc[1][1]
                    end

    Sim.add_transition_class(S, prodX1, deathX1, 
        cellIntake, cellExit, coagulation, fragmentation)

    return S
end

function IECFBD_ODEs(dM, M, parameters, t)
    c6 = parameters[:λ] # λ
    c5 = parameters[:kC] # kC
    c0 = parameters[:kE] # kE
    c1 = parameters[:kF] # kF
    c2 = parameters[:kI] # kI
    c4 = parameters[:kb] # kb
    c3 = parameters[:kd] # kd
    # Number of Compartments (N)
    dM[1] = c2+c1*M[3]+1/2*c5*M[1]-c0*M[1]-1/2*c5*M[1]^2
    # N^2
    dM[2] = c2+c5*M[1]^3+c5*M[2]+c0*M[1]+c1*M[3]+1/2*c5*M[1]^2-2*c0*M[2]+2*c1*M[5]+2*c2*M[1]-1/2*c5*M[1]-2*c5*M[2]*M[1]
    # Total Mass
    dM[3] = c6*c2+c4*M[1]-c0*M[3]-c3*M[3]
    # M^2
    dM[4] = c6*c2+c2*c6^2+c4*M[1]+c3*M[3]-2*c0*M[4]-2*c3*M[4]+2*c4*M[5]+c0*M[3]^2/M[1]+2*c6*c2*M[3]
    # N*M
    dM[5] = c6*c2+c0*M[3]+c1*M[4]+c2*M[3]+c4*M[2]+1/2*c5*M[5]-c3*M[5]-2*c0*M[5]+c6*c2*M[1]+1/2*c5*M[1]^2*M[3]-c5*M[5]*M[1]
    return
end

# initialize expected moments vector
function IECFBD_initial(N0, Mpc0)
    M0 = N0*Mpc0
    M = zeros(5)
    # Number of Compartments (N)
    M[1] = N0 # initial value for Moment(0,), please specify!
    # N^2
    M[2] = M[1]^2
    # Total Mass
    M[3] = M0 # initial value for Moment(1,), please specify!
    # M^2
    M[4] = M[3]^2
    # N*M
    M[5] = M[1]*M[3]
    return M
end

#####################

@export function nBDq(; kI=5.0, kE=0.1, kb=10.0, kd=0.01, λ=10.0, Ω=1.0, Ωc=1.0)
    kI, kE, kb, kd, λ, Ω, Ωc = float.( (kI, kE, kb, kd, λ, Ω, Ωc) )
    kI *= Ω
    kb *= Ωc; kd /= Ωc
    p = (kI, kE, kb, kd, λ, λ, 1.0, 1.0)
    ssaSys = getNBDqSSASystem(; kI=kI, kE=kE, kb=kb, kd=kd, λ=λ)
    return Model(p, Ω, Ωc, NDBqInit, NBDqODE, ssaSys)
end

NDBqInit(N, Mpc) = [N, N^2, N*Mpc, (N*Mpc)^2, N*N*Mpc, N*Mpc*Mpc]


function getNBDqSSASystem(;
                            kI = 0.0, kE = 0.0,  # Intake/Exit rates for cells
                            kb = 0.0, kd = 0.0,  # Birth/Death rates for X1
                            λ = 0,  # Parameter of the poisson distrib on intake
                            )
    S = Sim.System("NBDq", 1, # X1
            Dict(1 =>[0], # Number of compartments
                 2 =>[1], # M¹
                 3 =>[2], # M²
                )
            )

    # cells' reaction network
    prodX1 = Sim.new_chemical_reaction_class([1], kb)
    prodX1.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[1]
    prodX1.fast_sample_reactants! = fast_sample_uniform_cell

    deathX1 = Sim.new_chemical_reaction_class([-2], kd)
    deathX1.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> div(Mom[3]-Mom[2], 2)
    deathX1.fast_sample_reactants! = fast_sample_x1_square

    # Compartment reactions
    cellIntake = Sim.TransitionClass(0, 1, kI)
    cellIntake.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> 1
    # This is how the output compartment (yc) of the transition looks like
    cellIntake.pi = function(yc::Vector{Vector{Int64}}, xc::Vector{Vector{Int64}})
                        yc[1][1] = rand(Poisson(λ))
                    end
    # cell_division.fast_sample_reactants! = fast_sample_cell_div

    cellExit = Sim.TransitionClass(1, 0, kE) # This is a cell death
    cellExit.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[1] # Depends on num cells
    cellExit.fast_sample_reactants! = fast_sample_uniform_cell

    Sim.add_transition_class(S, prodX1, deathX1, cellIntake, cellExit)

    return S
end

@export function NBDqODE(du, u, p, t)
    # All the following are implied to be in average
    N, N², M1, M1², NM1, M2 = u
    kI, kE, kb, kd, μI, σ²I = p

    # @show (M1², NM1) #debug

    # This is by considering x(x-1)/2 as propensity
    dN   = kI - kE*N
    dN²  = kI*(1 + 2*N) + kE*(N - 2*N²)
    dM1  = kI*μI - kE*M1 + kb*N - kd*( (M1^2 / N) - M1 )
    dM1² = ( 
           kI*(σ²I + μI^2 + 2*μI*M1)
         # - kE*( 2*M1² - (M1^2 / N)) # Approximation of M2
         - kE*( 2*M1² - M2) # Actually using M2
         + kb*(N + 2*NM1)
         # + kd*( (1 - 2*M1/N)*M1² + (2 + NM1/N)*(M1^2)/N - 2*M1 )
         # + kd*(
         #      -( (2 * M1 * M1² / N) - (M1^2 * NM1 / N^2) )
         #      + M1²
         #      + 2*M2 # Actually using M2
         #      # + 2*(M1^2 / N) # Approx M2
         #      -2*M1
         #  )
         + kd*(
                -4*M1²*M1/N
                +2*NM1*(M1^2)/(N^2)
                +2*M1²
                +2*(M1^2)/N
                -2*M1
            )
        )
    dNM1 = ( kI*( (N+1)*μI + M1 )
         - kE*(2*NM1 - M1)
         + kb*N²
         # + kd*( 2* M1^2 * N²/ N^2 - 4*M1*NM1/N + 2*NM1 )
         +kd*( -2*M1*NM1/N + N²*(M1^2)/(N^2) + NM1 )
         )
    dM2  = ( 
           kI*(σ²I + μI)
         - kE*M2
         + kb*(N + 2*M1)
         + kd*( (2 - 4*M1/N)*M2 + (2*(M1^2)/(N^2) + 4*M1/N - 2)*M1 - 2*(M1^2)/N ) 
         # + kd*(
         #      -2*(M1^3 / N^2)
         #      +4*M2
         #      -2*M1
         #  )
        )

    # u .= [ N, N², M1, M1², NM1, M2 ]
    du .= [ dN, dN², dM1, dM1², dNM1, dM2 ]
end

### SUPPORTING FUNCTIONS FOR EFFICIENT SIMULATION ###
function fast_sample_first(r_indices::Vector{Int64}, n::Matrix{Int64}, Mom::Vector{Int64})
    r_indices[1] = 1
end

function fast_sample_uniform_cell(r_indices::Vector{Int64}, n::Matrix{Int64}, Mom::Vector{Int64})
    r_indices[1] = rand(1:Mom[1])
end

function fast_sample_uniform_2cells(r_indices::Vector{Int64}, n::Matrix{Int64}, Mom::Vector{Int64})
    r_indices[1] = rand(1:Mom[1])
    r_indices[2] = rand(1:Mom[1])
    while r_indices[1] == r_indices[2]
        r_indices[2] = rand(1:Mom[1])
    end
    sort!(r_indices) # This is necessary!
end

function fast_sample_mass_x1(r_indices::Vector{Int64}, n::Matrix{Int64}, Mom::Vector{Int64})
    N = Mom[1]
    H = Mom[2]
    rv = rand()*H
    r_indices[1] = 1
    val = 1.0*n[1,1]
    while val < rv
        r_indices[1] += 1
        val += n[1, r_indices[1]]
    end
    @assert r_indices[1]<=N """
        FATAL: Reactant ($(r_indices[1])) out of range ($N)!
        H = $H
        rv = $rv
        val = $val
        fast_sample_mass_x1
        """
end

# interaction_square(x) = x<=1 ? 0 : div(x*(x-1),2)
interaction_square(x) = x*(x-1)/2
interaction_cube(x) = x*(x-1)*(x-2)/6

function _getTotalPropensity(N, n, interaction)
    val = 0.0
    for i=1:N
        val += interaction(n[1, i])
    end
    return val
end

function fast_sample_x1_square(r_indices::Vector{Int64}, n::Matrix{Int64}, Mom::Vector{Int64})
    N = Mom[1]
    H = div(Mom[3]-Mom[2], 2)
    rv = rand()*H
    r_indices[1] = 1
    val = 1.0*interaction_square(n[1,1]) # Layout of n is [species, cell]
    while val < rv
        r_indices[1] += 1
        x = n[1, r_indices[1]]
        @assert x >= 0 "FATAL: The state ($x) of cell ($(r_indices[1])/$(N)) must be NON-NEGATIVE!"
        val += interaction_square(x)
    end
    @assert r_indices[1]<=N """
        FATAL: Reactant ($(r_indices[1])) out of range ($N)!
        H = $H
        rv = $rv
        val = $val
        fast_sample_x1_square
        """
end

function fast_sample_x1_cube(r_indices::Vector{Int64}, n::Matrix{Int64}, Mom::Vector{Int64})
    N = Mom[1]
    H = div(Mom[4] - 3*Mom[3] + 2*Mom[2], 6)
    rv = rand()*H
    r_indices[1] = 1
    val = 1.0*interaction_cube(n[1,1]) # Layout of n is [species, cell]
    while val < rv
        r_indices[1] += 1
        x = n[1, r_indices[1]]
        @assert x >= 0 "FATAL: The state ($x) of cell ($(r_indices[1])/$(N)) must be NON-NEGATIVE!"
        val += interaction_cube(x)
    end
    @assert r_indices[1]<=N """
        FATAL: Reactant ($(r_indices[1])) out of range ($N)!
        H = $H
        rv = $rv
        val = $val
        fast_sample_x1_cube
        """
end

"""
fast_sample_generic_1cell(r_indices::Vector{Int64}, n::Matrix{Int64}, 
                            propensity::Function, totalPropensity::Number)

Provides a generic fast sampling function.
    propensity is a function that must take the vector of chem content of a single compartment.
"""
function fast_sample_generic_unary(r_indices::Vector{Int64}, n::Matrix{Int64}, Mom::Vector{Int64}, 
                                    propensity::Function, totalPropensity::Function)
    rv = rand()*totalPropensity()
    r_indices[1] = 1
    val = 1.0*propensity(n[:, r_indices[1]]) # Layout of n is [species, cell]
    while val < rv
        r_indices[1] += 1
        val += propensity(n[:, r_indices[1]])
    end
end

function fast_sample_generic_binary(r_indices::Vector{Int64}, n::Matrix{Int64}, Mom::Vector{Int64}, 
                                    propensity::Function, totalPropensity::Function)
    N = Mom[1]
    rv = rand()*totalPropensity()
    val = 0.0
    for i=1:N
        for j=i+1:N
            r_indices[1] = i
            r_indices[2] = j
            val += propensity(n[:, r_indices[1]], n[:, r_indices[2]]) # Layout of n is [species, cell]
            if val >= rv
                break
            end
        end
        if val >= rv
            break
        end
    end
end

function fast_sample_generic_binary_check(r_indices::Vector{Int64}, n::Matrix{Int64}, Mom::Vector{Int64}, 
                                    propensity::Function, totalPropensity::Function)
    N = Mom[1]
    rv = rand()*totalPropensity()
    val = 0.0
    He = 0.0
    for i=1:N
        for j=i+1:N
            r_indices[1] = i
            r_indices[2] = j
            lv = propensity(n[:, r_indices[1]], n[:, r_indices[2]]) # Layout of n is [species, cell]
            He += lv
            if val >= rv
                # break
                nothing
            else
                val += lv
            end
        end
        if val >= rv
            # break
            nothing
        end
    end
    @assert He == totalPropensity() """
    ERROR: Total propensity is not consistent!
    H: $(totalPropensity())
    He: $He
    fast_sample_generic_binary
    """
end

#eof
