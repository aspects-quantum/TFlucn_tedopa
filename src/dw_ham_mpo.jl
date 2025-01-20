
##########################################################################

function HB(ab1, ab2, s_total)

    #s_total = [(n == S_pos) ? Index(2, "S=1/2") : Index(boson_dim, "Qudit") for n = 1:tot_chain]

    tot_chain = length(s_total)
    S_pos = Int((tot_chain+1)/2)

    ω_n_REAL = ab1[1:S_pos-1, 1]
    ω_n_TILD = ab2[1:S_pos-1, 1]
    ω_n_total = append!(reverse(ω_n_TILD),[0],ω_n_REAL)
    t_n_REAL = sqrt.(ab1[2:S_pos, 2])
    t_n_TILD = sqrt.(ab2[2:S_pos, 2])
    t_n_total = append!(reverse(t_n_TILD),[0],t_n_REAL)

    H = OpSum()
    
    for j in 2:tot_chain-1
        
        if j < S_pos
            ω_n = ω_n_total[j]
            t_n = t_n_total[j]

            H .+=  (-ω_n, "N", j, "Id", j-1)
            H .+=  (-t_n, "Adag", j, "A", j-1)
            H .+=  (-t_n, "A", j, "Adag", j-1)
        elseif j == S_pos
            continue

        else
            ω_n = ω_n_total[j]
            t_n = t_n_total[j]

            H .+=  (ω_n, "N", j, "Id", j+1)
            H .+=  (t_n, "Adag", j, "A", j+1)
            H .+=  (t_n, "A", j, "Adag", j+1)
        end
    end
    
    return MPO(H, s_total)
end



##########################################################################

function exp_xHB(ab1, ab2, x, s_total)

    #s_total = [(n == S_pos) ? Index(2, "S=1/2") : Index(boson_dim, "Qudit") for n = 1:tot_chain]

    tot_chain = length(s_total)
    S_pos = Int((tot_chain+1)/2)

    ω_n_REAL = ab1[1:S_pos-1, 1]
    ω_n_TILD = ab2[1:S_pos-1, 1]
    ω_n_total = append!(reverse(ω_n_TILD),[0],ω_n_REAL)
    t_n_REAL = sqrt.(ab1[2:S_pos, 2])
    t_n_TILD = sqrt.(ab2[2:S_pos, 2])
    t_n_total = append!(reverse(t_n_TILD),[0],t_n_REAL)

    gates = ITensor[]
    
    for j in 2:tot_chain-1
        
        if j < S_pos
            s1 = s_total[j]
            s2 = s_total[j-1]
    
            ω_n = ω_n_total[j]
            t_n = t_n_total[j]
            hj = (-ω_n) * op("N", s1) * op("Id", s2) +
                (-t_n) * op("Adag", s1) * op("A", s2) +
                (-t_n) * op("A", s1) * op("Adag", s2)
            Gj = exp((x / 2)* hj)
            push!(gates, Gj)

        elseif j == S_pos
            continue

        else
            s1 = s_total[j]
            s2 = s_total[j+1]
    
            ω_n = ω_n_total[j]
            t_n = t_n_total[j]
            hj = ω_n * op("N", s1) * op("Id", s2) +
                t_n * op("Adag", s1) * op("A", s2) +
                t_n * op("A", s1) * op("Adag", s2)
            Gj = exp((x / 2)* hj)
            push!(gates, Gj)
        end
    end
    
    return append!(gates, reverse(gates))
end


##########################################################################
function unit_gates(ω_0,Ω,c_01,c_02,ab1, ab2, tau, s_total)
    tot_chain = length(s_total)
    S_pos = Int((tot_chain+1)/2)

    ω_n_REAL = ab1[1:S_pos-1, 1]
    ω_n_TILD = ab2[1:S_pos-1, 1]
    ω_n_total = append!(reverse(ω_n_TILD),[0],ω_n_REAL)
    t_n_REAL = sqrt.(ab1[2:S_pos, 2])
    t_n_TILD = sqrt.(ab2[2:S_pos, 2])
    t_n_total = append!(reverse(t_n_TILD),[0],t_n_REAL)

    gates = ITensor[]
    
    for j in 2:tot_chain-1
        
        if j < S_pos
            s1 = s_total[j]
            s2 = s_total[j-1]
    
            ω_n = ω_n_total[j]
            t_n = t_n_total[j]
            hj = (-ω_n) * op("N", s1) * op("Id", s2) +
                 (-t_n) * op("Adag", s1) * op("A", s2) +
                 (-t_n) * op("A", s1) * op("Adag", s2)
            Gj = exp(-im * (tau / 2) * hj)
            push!(gates, Gj)

        elseif j == S_pos
            hj = ω_0 * op("Sz", s_total[j]) * op("Id", s_total[j+1])  +
                 Ω * op("Sx", s_total[j]) * op("Id", s_total[j+1])   +
                 c_01 * op("Sx", s_total[j]) * op("A", s_total[j+1]) +
                 c_01 * op("Sx", s_total[j]) * op("Adag", s_total[j+1])
            Gj = exp(-im * (tau / 2) * hj)
            push!(gates, Gj)
            hj = c_02 * op("Sx", s_total[j]) * op("A", s_total[j-1]) +
                 c_02 * op("Sx", s_total[j]) * op("Adag", s_total[j-1])
            Gj = exp(-im * (tau / 2) * hj)
            push!(gates, Gj)

        else
            s1 = s_total[j]
            s2 = s_total[j+1]
    
            ω_n = ω_n_total[j]
            t_n = t_n_total[j]
            hj = ω_n * op("N", s1) * op("Id", s2) +
                 t_n * op("Adag", s1) * op("A", s2) +
                 t_n * op("A", s1) * op("Adag", s2)
            Gj = exp(-im * (tau / 2) * hj)
            push!(gates, Gj)
        end
    end
    
    return append!(gates, reverse(gates))
   
end

