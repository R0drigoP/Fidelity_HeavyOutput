using LinearAlgebra
using Statistics
using Random
using Permutations

include(joinpath(@__DIR__,"free_fermions.jl"))
include(joinpath(@__DIR__,"auxiliary.jl"))

# Auxiliary functions

function BinToDec(bin::Array{Int})::UInt
    N::UInt = length(bin)
    sum::UInt = 0
    for i in 1:N
        sum += bin[i]*2^(N-i)
    end

    return sum
end


function inner_product(x::Vector, y::Vector)
    if length(x) != length(y)
        ArgumentError("Vectors must have the same length to compute inner product")
    end

    sum = 0

    for i in 1:length(x)
        sum += x[i]*y[i]
    end

    return sum
end

# Permutation algorithm for 1D connectivity
function brick_sort(m::Int, perm0::Array{Int,1}, start_shift::Int = 0)
    perm = copy(perm0)
    layer_counter = 0
    was_sorted1 = false
    was_sorted2 = false
    swap_list = []

    while true
        if start_shift == 0
            was_sorted1 = true
            for i in 1:div(m,2)
                if perm[2*i-1] > perm[2*i]
                    perm[2*i-1], perm[2*i] = perm[2*i], perm[2*i-1]
                    was_sorted1 = false
                    push!(swap_list, [2*i-1, 2*i])
                end
            end

            if was_sorted1 && was_sorted2
                break
            elseif !was_sorted1
                layer_counter += 1
            end
        end

        was_sorted2 = true
        for i in 1:div(m,2) + m % 2 - 1
            if perm[2*i] > perm[2*i+1]
                perm[2*i], perm[2*i+1] = perm[2*i+1], perm[2*i]
                was_sorted2 = false
                push!(swap_list, [2*i, 2*i+1])
            end
        end

        if was_sorted1 && was_sorted2
            break
        elseif !was_sorted2
            layer_counter += 1
        end
        start_shift = 0
    end

    return layer_counter, reverse(swap_list)
end

function brick_sort1D(m::Int, perm::Array{Int,1})
    layer_n1, swap_list1 = brick_sort(m, perm, 0)
    layer_n2, swap_list2 = brick_sort(m, perm, 1)

    if layer_n1 < layer_n2
        return layer_n1, swap_list1
    else
        return layer_n2, swap_list2
    end
end


# Permutates qubits q1 with q2.
# psi0 is the ideal state, psi0_t is faulty state, L is number of qubits
# flag = 0 permute psi_t, flag = 1 don't permute psi_t (faulty)

function perm_2qubits(q1, q2, psi0, psi0_t, L, flag = 0)

    if q1 == q2
        return psi0, psi0_t
    end
    
    nums::Vector{Int} = zeros(L)
    pos::Vector{Int} = zeros(4)
    
    psi = copy(psi0)
    psi_t = copy(psi0_t)

    
    for i in 0:2^(L-2) - 1
        Bi = digits(i, base=2, pad=L-2)|>reverse #-> this gives correct binary but not needed
        
        # j loop generates all 4 numbers and puts them in nums vector in the position of the qubits being altered
        for j in 0:3
            Bj = digits(j, base=2, pad=2)|>reverse
            
            nums[q1] = Bj[1]
            nums[q2] = Bj[2]
            
            # k loop fills nums vector
            p::UInt = 1
            
            for k in 1:L
                if k != q1 && k != q2
                    nums[k] = Bi[p]
                    p += 1
                end
            end

            pos[j+1] = BinToDec(nums) + 1
            
        end # End of j loop

        vec = [ psi[pos[1]], psi[pos[3]], psi[pos[2]], psi[pos[4]] ]
        
        if flag == 0
            vect = [ psi_t[pos[1]], psi_t[pos[3]], psi_t[pos[2]], psi_t[pos[4]] ]
        else # flag = 1 => faulty permutation
            vect = [ psi_t[pos[1]], psi_t[pos[2]], psi_t[pos[3]], psi_t[pos[4]] ]
        end

        for m in 1:4
            psi[pos[m]] = vec[m]
            psi_t[pos[m]] = vect[m]
        end # End of m loop

    end # End of i loop
        
    return psi, psi_t

end

# Faulty permutations on entire layer. p is probability of occurrence of faulty permutation
function random_faulty_perm1D(psi0, psi0_t, L, p = 0)

    perm0 = RandomPermutation(L).data
    
    psi = copy(psi0)
    psi_t = copy(psi0_t)

    perms = brick_sort1D(L, perm0)[2] # decompostion of perm0 in 1D perms
    N = length(perms)
    
    for i in 1:N
        flag = 0
        
        if rand(Float64) < p # flag = 1 for faulty
            flag = 1
        end
        
        psi, psi_t = perm_2qubits(perms[i][1], perms[i][2], psi, psi_t, L, flag)
        
    end
    
    return psi, psi_t
    
end

# Applies 2-qb gate between q1 and q2. Matrix U is applied to psi0 and matrix exp(-i*alpha*H)*U to psi0_t
function interaction_2qubit(q1, q2, psi0, psi0_t, L, U, H, alpha = 0)

    if q1 == q2
        error("Cannot have interaction between the same qubit")
    end
    
    nums::Vector{Int} = zeros(L)
    pos::Vector{Int} = zeros(4)
    
    psi = copy(psi0)
    psi_t = copy(psi0_t)
    
    
    for i in 0:2^(L-2) - 1
        Bi = digits(i, base=2, pad=L-2)|>reverse #-> this gives correct binary but not needed
        
        # j loop generates all 4 numbers and puts them in nums vector in the position of the qubits being altered
        for j in 0:3
            Bj = digits(j, base=2, pad=2)|>reverse
            
            nums[q1] = Bj[1]
            nums[q2] = Bj[2]
            
            # k loop fills nums vector
            p::UInt = 1
            
            for k in 1:L
                if k != q1 && k != q2
                    nums[k] = Bi[p]
                    p += 1
                end
            end

            pos[j+1] = BinToDec(nums) + 1
            
        end # End of j loop

        vec = U*[ psi[pos[1]], psi[pos[2]], psi[pos[3]], psi[pos[4]] ]

        vect = exp(-im*alpha*H)*U*[ psi_t[pos[1]], psi_t[pos[2]], psi_t[pos[3]], psi_t[pos[4]] ]


        for m in 1:4
            psi[pos[m]] = vec[m]
            psi_t[pos[m]] = vect[m]
        end # End of m loop

    end # End of i loop
        
    return psi, psi_t

end

function heavy_output(state, state_faulty, L)
    probs::Vector{Float64} = conj(state).*state
    probs_faulty::Vector{Float64} = conj(state_faulty).*state_faulty
    
    med = median(probs)

    # heavy output
    hD = 0
    #count = 0

    for i in 1:2^L
        if probs[i] > med
            hD += probs_faulty[i]
            #count += 1
        end
    end

    #println(count)
    
    return hD
end

# Random Matrix Ensembles

function rand_CUE(N)::Matrix{ComplexF64}
    A = randn(ComplexF64, N, N)
    
    Q, R = qr(A)

    lambda = diagm([R[i,i]/abs(R[i,i]); for i in 1:N])

    return Q*lambda
end


function rand_GUE(n)::Matrix{ComplexF64}
    M = Array{ComplexF64}(undef, n, n)
    for i in 1:n
        M[i, i] = randn(Float64)
        for j in i+1:n
            M[i, j] = randn(ComplexF64)
            M[j, i] = conj(M[i, j])
        end
    end
    
    return M
end
