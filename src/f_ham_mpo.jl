
function exp_xHB_comp(ab1, ab2, x, s_total)

    #s_total = [(n == S_pos)? Index(2, "S=1/2") : Index(boson_dim, "Qudit") for n = 1:tot_chain]
    
    tot_chain = length(s_total)
    S_pos = Int((tot_chain+1)/2)

    gates = ITensor[]
    i = 1
    for j in S_pos+1:tot_chain-1
        s1 = s_total[j]
        s2 = s_total[j+1]
        s1p = s_total[tot_chain-j+1]
        s2p = s_total[tot_chain-j]

        ω_n1 = ab1[i, 1]
        t_n1 = sqrt(ab1[i+1, 2])
        ω_n2 = ab2[i, 1]
        t_n2 = sqrt(ab2[i+1, 2])
        hj = ω_n1 * op("N", s1) * op("Id", s2) +
            t_n1 * op("Adag", s1) * op("A", s2) +
            t_n1 * op("A", s1) * op("Adag", s2)
            -ω_n2 * op("N", s1p) * op("Id", s2p) 
            -t_n2 * op("Adag", s1p) * op("A", s2p) 
            -t_n2 * op("A", s1p) * op("Adag", s2p)
        Gj = exp((x / 2)*hj)
        i=i+1
        push!(gates, Gj)
    end
    append!(gates, reverse(gates))
    return gates
end



##########################################################################
function tot_gate(ω_0, Ω, c_01, c_02, ab1, ab2, tau, s_total)

    #s_total = [(n == S_pos)? Index(2, "S=1/2") : Index(boson_dim, "Qudit") for n = 1:tot_chain]
    
    tot_chain = length(s_total)
    S_pos = Int((tot_chain+1)/2)

    gates = ITensor[]
    i=1
    for j in S_pos:tot_chain-1
        s1 = s_total[j]
        s2 = s_total[j+1]
        s1p = s_total[tot_chain-j+1]
        s2p = s_total[tot_chain-j]

        if j == S_pos
            hj = ω_0 * op("Sz", s1) * op("Id", s2) +
                 Ω * op("Sx", s1) * op("Id", s2)   +
                 c_01 * op("Sz", s1) * op("A", s2) +
                 c_01 * op("Sz", s1) * op("Adag", s2) 
            Gj = exp(-im * (tau / 2) * hj)
            push!(gates, Gj)
            hj = c_02 * op("Sz", s1) * op("A", s_total[j-1]) +
                 c_02 * op("Sz", s1) * op("Adag", s_total[j-1])
            Gj = exp(-im * (tau / 2) * hj)
            push!(gates, Gj)
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
            Gj = exp(-im * (tau / 2) * hj)
            push!(gates, Gj)
        end
        i = i+1
        
    end
    return append!(gates, reverse(gates))
end
