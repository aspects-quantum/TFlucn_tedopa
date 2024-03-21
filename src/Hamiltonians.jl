###########################################################################################
## this script contains the following functions: 
## "init_Gates", "unit_gates", "unit_dates", "count_evol", "count_devol",
## "total_Hamiltonian", "chain_Hamiltonian"
###########################################################################################


##############################################################################################################
## this function outputs a list of Trotterized gates for imaginary time evolution to get initial thermal state
##############################################################################################################
function exp_xHB(ω_C, ab1, x, s_total)

    #s_total = [(n == S_pos) ? Index(2, "S=3/2") : Index(boson_dim, "Qudit") for n = 1:tot_chain]

    tot_chain = length(s_total)
    S_pos_real = Int(tot_chain/2 + 1)

    gates = ITensor[]
    i = 1
    for j in S_pos_real+1:tot_chain-1
        s1 = s_total[j]
        s2 = s_total[j+1]

        ω_n = ω_C * ab1[i, 1]
        t_n = ω_C * sqrt(ab1[i+1, 2])
        hj = ω_n * op("N", s1) * op("Id", s2) +
            t_n * op("Adag", s1) * op("A", s2) +
            t_n * op("A", s1) * op("Adag", s2)
        
        Gj = exp((x / 2)* hj)
        i=i+1
        push!(gates, Gj)
    end
    append!(gates, reverse(gates))
#= 
    i = 1
    for j in S_pos:-1:2
        s1 = s_total[j]
        s2 = s_total[j-1]

        if i == 1
            #hj = op("zeros_", s_total[j])
            Gj = op("Id",s1)
        else
            ω_n = ω_C * ab2[i-1, 1]
            t_n = ω_C * sqrt(ab2[i, 2])
            hj = -ω_n * op("N", s1) * op("Id", s2) +
                 -t_n * op("Adag", s1) * op("A", s2) +
                 -t_n * op("A", s1) * op("Adag", s2)
            
            Gj = exp((x / 2)* hj)
        end
        i=i+1
        push!(gates, Gj)
    end
    return append!(gates, reverse(gates)) =#
    return gates

end



##########################################################################
function unit_real(ω_0,Ω,c_0,ω_C, ab1, tau, s_total)

    tot_chain = length(s_total)
    S_pos_real = Int(tot_chain/2 + 1)

    gates = ITensor[]
    i=1
    for j in S_pos_real:tot_chain-1
        s1 = s_total[j]
        s2 = s_total[j+1]

        if j == S_pos_real
            hj = ω_0 * op("Sz", s1) * op("Id", s2)  +
                 Ω * op("Sx", s1) * op("Id", s2)   +
                 c_0 * op("Sz", s1) * op("A", s2) +
                 c_0 * op("Sz", s1) * op("Adag", s2)
        else
            ω_n = ω_C * ab1[i-1, 1]
            t_n = ω_C * sqrt(ab1[i, 2])
            hj = ω_n * op("N", s1) * op("Id", s2)  +
                 t_n * op("Adag", s1) * op("A", s2) +
                 t_n * op("A", s1) * op("Adag", s2)
        end

        Gj = exp(-im * (tau / 2) * hj)
        push!(gates, Gj)
        i = i+1
    end
    # Include gates in reverse order too
    # (N,N-1),(N-1,N-2),...
    append!(gates, reverse(gates))

end

#= 

##########################################################################
function unit_gates(N_chain,ω_0,Ω,c_0,ω_C, ab1, ab2, tau, s_total)

    tot_chain = 2*N_chain+1
    S_pos = N_chain+1

    gates = ITensor[]
    i=1
    for j in S_pos:tot_chain-1
        s1 = s_total[j]
        s2 = s_total[j+1]

        if j == S_pos
            hj = ω_0 * op("Sz1", s1) * op("Id", s2)  +
                 Ω * op("Sx1", s1) * op("Id", s2)   +
                 c_0 * op("Sz1", s1) * op("A", s2) +
                 c_0 * op("Sz1", s1) * op("Adag", s2)
        else
            ω_n = ω_C * ab1[i-1, 1]
            t_n = ω_C * sqrt(ab1[i, 2])
            hj = ω_n * op("N", s1) * op("Id", s2)  +
                 t_n * op("Adag", s1) * op("A", s2) +
                 t_n * op("A", s1) * op("Adag", s2)
        end

        Gj = exp(-im * (tau / 2) * hj)
        push!(gates, Gj)
        i = i+1
    end
    # Include gates in reverse order too
    # (N,N-1),(N-1,N-2),...
    append!(gates, reverse(gates))

    i=1
    for j in S_pos:-1:2
        s1 = s_total[j]
        s2 = s_total[j-1]

        if j == S_pos
            hj = (-ω_0) * op("Sz2", s1) * op("Id", s2)  +
                 (-Ω) * op("Sx2", s1) * op("Id", s2)  +
                 c_0 * op("Sz2", s1) * op("A", s2) +
                 c_0 * op("Sz2", s1) * op("Adag", s2)
        else
            ω_n = ω_C * ab2[i-1, 1]
            t_n = ω_C * sqrt(ab2[i, 2])
            hj = -ω_n * op("N", s1) * op("Id", s2)  +
                 (-t_n) * op("Adag", s1) * op("A", s2) +
                 (-t_n) * op("A", s1) * op("Adag", s2)
        end

        Gj = exp(-im * (tau / 2) * hj)
        push!(gates, Gj)
        i = i+1
    end
    # Include gates in reverse order too
    # (N,N-1),(N-1,N-2),...
    return append!(gates, reverse(gates))
end
 =#


#= Parameters #########

w_fn(t)  ## weight function 
supp  ## support of the weight function

N_chain 
N_coeff = N_chain+1;
Nquad          ## Number of quadrature points
ab = recur_coeff(w_fn, supp, N_coeff,Nquad)    ## recurrence coefficients

ω_C      ## bath cutoff
α        ## coupling strength
ω_0      ## spin splitting
Ω        ## independent model if Ω = 0
boson_dim     ## number of level considered for bosonic chain sites
=#
