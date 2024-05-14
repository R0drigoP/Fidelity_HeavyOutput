using LinearAlgebra
using Statistics
using DelimitedFiles
using Random
using StatsBase

include(joinpath(@__DIR__,"libs.jl"))

#------------------

function main()

    L::UInt32 = parse(UInt32, ARGS[1]) # nb of qubits
    T::UInt32 = parse(UInt32, ARGS[2]) # max nb of layers
    It::UInt32 = parse(UInt32, ARGS[3]) # nb of iterations to average

    alpha = parse(Float64, ARGS[4])
    prob = parse(Float64, ARGS[5])

    #--------------------------

    FtAv = Array{Float64}(undef, T)
    errors = Array{Float64}(undef, T)
    Ft = Array{Float64}(undef, T, It)

    #initial state
    psi_aux::Vector = zeros(2^L)
    psi_aux[1] = 1


    for iter in 1:It
        
        psi0::Vector{ComplexF64} = copy(psi_aux)
        psi0_t::Vector{ComplexF64} = copy(psi_aux)
        
        for t in 1:T
            
            psi0, psi0_t = random_faulty_perm(psi0, psi0_t, L, prob)
            
            for q1 in 1:2:L-1
                
                U = rand_CUE(4)
                H = rand_GUE(4)
                
                psi0, psi0_t = interaction_2qubit(q1, q1+1, psi0, psi0_t, L, U, H, alpha)
                
            end # End of q1 loop

            #F = heavy_output(psi0, psi0_t, L)

            F = abs( inner_product(conj(psi0), psi0_t) )^2 #fidelity
            Ft[t, iter] = F
            
        end # End of time loop
        
    end # End of iteration loop

    for i in 1:T
        FtAv[i] = mean(Ft[i,:])
        errors[i] = std(Ft[i,:])/sqrt(It)
    end

    #dirF = string("<file directory>")
    #dirError = string("<file directory")

    #writedlm(dirF, FtAv, ',')
    #writedlm(dirError, errors, ',')

    display(FtAv)
end

main()
