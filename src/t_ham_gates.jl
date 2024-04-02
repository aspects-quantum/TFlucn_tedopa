
function exp_xHB(ab1, x, s_total)

    #s_total = [(n == S_pos) ? Index(2, "S=1/2") : Index(boson_dim, "Qudit") for n = 1:tot_chain]


    tot_chain = length(s_total)
    S_pos = Int((tot_chain+1)/2)

    gates = ITensor[]
    i = 1
    for j in S_pos+1:tot_chain-1
        s1 = s_total[j]
        s2 = s_total[j+1]

        ω_n = ab1[i, 1]
        t_n = sqrt(ab1[i+1, 2])
        hj = ω_n * op("N", s1) * op("Id", s2) +
            t_n * op("Adag", s1) * op("A", s2) +
            t_n * op("A", s1) * op("Adag", s2)
        
        Gj = exp((x / 2)* hj)
        i=i+1
        push!(gates, Gj)
    end
    append!(gates, reverse(gates))

    return gates

end



##########################################################################
function unit_gates(ω_0,Ω,c_01,c_02,ab1, ab2, tau, s_total)

    #s_total = [(n == S_pos) ? Index(2, "S=1/2") : Index(boson_dim, "Qudit") for n = 1:tot_chain]


    tot_chain = length(s_total)
    S_pos = Int((tot_chain+1)/2)


    gates = ITensor[]
    i=1
    for j in S_pos:tot_chain-1
        s1 = s_total[j]
        s2 = s_total[j+1]

        if j == S_pos
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
    for j in S_pos:-1:2
        s1 = s_total[j]
        s2 = s_total[j-1]

        if j == S_pos
            hj = c_02 * op("Sz", s1) * op("A", s2) +
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


