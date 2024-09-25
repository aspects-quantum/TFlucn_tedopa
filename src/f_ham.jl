
function exp_xHB_comp(ab1, ab2, x, s_total)

    #s_total = [(n == S_pos_t) | (n == S_pos_r) ? Index(2, "S=1/2") : Index(boson_dim, "Qudit") for n = 1:tot_chain]
    #s_tot_comp = s_total[S_pos_r:tot_chain]
    N_ch = Int(length(s_total)/2)

    gates = ITensor[]
    i = 1
    for j in 2:N_ch-1
        s1 = s_total[N_ch+j]
        s2 = s_total[N_ch+j+1]
        s1p = s_total[N_ch+1-j]
        s2p = s_total[N_ch-j]

        ω_n1 = ab1[i, 1]
        t_n1 = sqrt(ab1[i+1, 2])
        ω_n2 = ab2[i, 1]
        t_n2 = sqrt(ab2[i+1, 2])
        hj = ω_n1 * op("N", s1) * op("Id", s2) +
            t_n1 * op("Adag", s1) * op("A", s2) +
            t_n1 * op("A", s1) * op("Adag", s2)
            - ω_n2 * op("N", s1p) * op("Id", s2p) 
            - t_n2 * op("Adag", s1p) * op("A", s2p) 
            - t_n2 * op("A", s1p) * op("Adag", s2p)
        
        Gj = exp((x / 2)*hj)
        i=i+1
        push!(gates, Gj)
    end
    append!(gates, reverse(gates))
    return gates
end



##########################################################################
function tot_gate(ω_0, Ω, c_01, c_02, ab1, ab2, tau, s_total)

    #s_total = [(n == S_pos_t) | (n == S_pos_r) ? Index(2, "S=1/2") : Index(boson_dim, "Qudit") for n = 1:tot_chain]
    #s_tot_comp = s_total[S_pos_r:tot_chain]
    N_ch = Int(length(s_total)/2)

    gates = ITensor[]
    i=1
    for j in 1:N_ch-1
        s1 = s_total[N_ch+j]
        s2 = s_total[N_ch+j+1]
        s1p = s_total[N_ch+1-j]
        s2p = s_total[N_ch-j]

        if j == 1
            hj = ω_0 * op("Sz", s1) * op("Id", s2) +
                 Ω * op("Sx", s1) * op("Id", s2)   +
                 c_01 * op("Sz", s1) * op("A", s2) +
                 c_01 * op("Sz", s1) * op("Adag", s2)
                 #= -ω_0 * op("Sz", s1p) * op("Id", s2p)  
                 -Ω * op("Sx", s1p) * op("Id", s2p)    =#
                 c_02 * op("Sz", s1) * op("A", s2p) 
                 c_02 * op("Sz", s1) * op("Adag", s2p)
        else
            ω_n1 = ab1[i-1, 1]
            t_n1 = sqrt(ab1[i, 2])
            ω_n2 = ab2[i-1, 1]
            t_n2 = sqrt(ab2[i, 2])
            hj = ω_n1 * op("N", s1) * op("Id", s2)  +
                 t_n1 * op("Adag", s1) * op("A", s2) +
                 t_n1 * op("A", s1) * op("Adag", s2)
                -ω_n2 * op("N", s1p) * op("Id", s2p) 
                -t_n2 * op("Adag", s1p) * op("A", s2p) 
                -t_n2 * op("A", s1p) * op("Adag", s2p)
        end

        Gj = exp(-im * (tau / 2) * hj)
        push!(gates, Gj)
        i = i+1
    end
    
    return append!(gates, reverse(gates))
end

##########################################################################
function tot_gate_old(ω_0,Ω,c_01,c_02,ab1, ab2, tau, s_total)

    #s_total = [(n == S_pos_t) | (n == S_pos_r) ? Index(2, "S=1/2") : Index(boson_dim, "Qudit") for n = 1:tot_chain]
    tot_chain = length(s_total)
    S_pos_t = Int(tot_chain/2)
    S_pos_r = S_pos_t+1

    gates = ITensor[]
    i=1
    for j in S_pos_r:tot_chain-1
        s1 = s_total[j]
        s2 = s_total[j+1]

        if j == S_pos_r
            hj = ω_0 * op("Sz", s1) * op("Id", s2)  +
                 Ω * op("Sx", s1) * op("Id", s2)   +
                 c_01 * op("Sz", s1) * op("A", s2) +
                 c_01 * op("Sz", s1) * op("Adag", s2)
        else
            ω_n = ab1[i-1, 1]
            t_n = sqrt(ab1[i, 2])
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
    for j in S_pos_t:-1:2
        s1 = s_total[j]
        s2 = s_total[j-1]

        if j == S_pos_t
            hj = -ω_0 * op("Sz", s1) * op("Id", s2)  +
                  - Ω * op("Sx", s1) * op("Id", s2)   +
                 c_02 * op("Sz", s1) * op("A", s2) +
                 c_02 * op("Sz", s1) * op("Adag", s2)
        else
            ω_n = ab2[i-1, 1]
            t_n = sqrt(ab2[i, 2])
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


##########################################################################


function exp_xHB_comp_old(ab1, x, s_tot_comp)

    #s_total = [(n == S_pos_t) | (n == S_pos_r) ? Index(2, "S=1/2") : Index(boson_dim, "Qudit") for n = 1:tot_chain]
    #s_tot_comp = s_total[S_pos_r:tot_chain]
    N_chain = length(s_tot_comp)

    gates = ITensor[]
    i = 1
    for j in 2:N_chain-1
        s1 = s_tot_comp[j]
        s2 = s_tot_comp[j+1]

        ω_n = ab1[i, 1]
        t_n = sqrt(ab1[i+1, 2])
        hj = ω_n * op("N", s1) * op("Id", s2) +
            t_n * op("Adag", s1) * op("A", s2) +
            t_n * op("A", s1) * op("Adag", s2)
        
        Gj = exp((x / 2)*hj)
        i=i+1
        push!(gates, Gj)
    end
    append!(gates, reverse(gates))

    return gates
end

