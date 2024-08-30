
function exp_xHB_tdvp(ab1, x, s_tot_comp)

    #s_total = [(n == S_pos_t) | (n == S_pos_r) ? Index(2, "S=1/2") : Index(boson_dim, "Qudit") for n = 1:tot_chain]
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
        
        Gj = exp((x / 2)* hj)
        i=i+1
        push!(gates, Gj)
    end
    append!(gates, reverse(gates))

    return gates
end
##########################################################################
function HBB_tf(ab1, ab2, s_total)

    #s_total = [(n == S_pos_t) | (n == S_pos_r) ? Index(2, "S=1/2") : Index(boson_dim, "Qudit") for n = 1:tot_chain]
    tot_chain = length(s_total)
    S_pos_t = Int(tot_chain/2)
    S_pos_r = S_pos_t+1

    hbb = OpSum()
    i=1
    for j in S_pos_r+1:tot_chain-1
        ω_n = ab1[i, 1]
        t_n = sqrt(ab1[i+1, 2])
        hbb += ω_n, "N", j
        hbb += t_n, "Adag", j, "A", j+1
        hbb += t_n, "A", j, "Adag", j+1
        i = i+1
    end

    i=1
    for j in S_pos_t-1:-1:2
        s1 = s_total[j]
        s2 = s_total[j-1]
        ω_n = ab2[i, 1]
        t_n = sqrt(ab2[i+1, 2])
        hbb += -ω_n, "N", j
        hbb += -t_n, "Adag", j, "A", j-1
        hbb += -t_n, "A", j, "Adag", j-1
        i = i+1
    end
    
    HBB = MPO(hbb, s_total)

    return HBB
end


##########################################################################
function HTOT_tf(ω_0, Ω, c_01, c_02, ab1, ab2, s_total)

    #s_total = [(n == S_pos_t) | (n == S_pos_r) ? Index(2, "S=1/2") : Index(boson_dim, "Qudit") for n = 1:tot_chain]
    tot_chain = length(s_total)
    S_pos_t = Int(tot_chain/2)
    S_pos_r = S_pos_t+1

    htot = OpSum()
    i=1
    for j in S_pos_r:tot_chain-1

        if j == S_pos_r
            htot += ω_0, "Sz", j
            htot +=   Ω, "Sx", j
            htot +=c_01, "Sz", j, "A", j+1
            htot +=c_01,"Sz", j, "Adag", j+1
        else
            ω_n = ab1[i-1, 1]
            t_n = sqrt(ab1[i, 2])
            htot += ω_n, "N", j
            htot += t_n, "Adag", j, "A", j+1
            htot += t_n, "A", j, "Adag", j+1
        end
        i = i+1
    end

    i=1
    for j in S_pos_t:-1:2
        s1 = s_total[j]
        s2 = s_total[j-1]

        if j == S_pos_t
            htot +=-ω_0, "Sz", j
            htot +=  -Ω, "Sx", j
            htot +=c_02, "Sz", j, "A", j-1
            htot +=c_02,"Sz", j, "Adag", j-1
        else
            ω_n = ab2[i-1, 1]
            t_n = sqrt(ab2[i, 2])
            htot += -ω_n, "N", j
            htot += -t_n, "Adag", j, "A", j-1
            htot += -t_n, "A", j, "Adag", j-1
        end
        i = i+1
    end
    
    HTOT = MPO(htot, s_total)

    return HTOT
end




#= 

function HB_tf(ab1, s_total)

    #s_total = [(n == S_pos_t) | (n == S_pos_r) ? Index(2, "S=1/2") : Index(boson_dim, "Qudit") for n = 1:tot_chain]
    tot_chain = length(s_total)
    S_pos_t = Int(tot_chain/2)
    S_pos_r = S_pos_t+1

    hB = OpSum()
    i = 1
    for j in S_pos_r+1:tot_chain-1
        ω_n = ab1[i, 1]
        t_n = sqrt(ab1[i+1, 2])

        hB += ω_n, "N", j
        hB += t_n, "Adag", j, "A", j+1
        hB += t_n, "A", j, "Adag", j+1
        i=i+1
    end
    
    HB = MPO(hB, s_total)
    return HB
end
 =#
