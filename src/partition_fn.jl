
function partition_fn(N_chain, ω_C, ab1, s_total, β)

    S_pos = N_chain+1

    exp_β_HB =  exp_xHB(N_chain, ω_C, ab1, β, s_total)

    for i = 1:S_pos
        exp_β_HB = push!(exp_β_HB, op("Id",s_total[i]))
    end

    EXP_β_HB = MPO(exp_β_HB)

    return tr(EXP_β_HB)

end