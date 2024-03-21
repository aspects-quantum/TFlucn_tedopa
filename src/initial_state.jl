###########################################################################################
######## Outputs an initial state of (spin+chain)
###########################################################################################


function init_state(s_total)
    
    tot_chain = length(s_total)
    S_pos_tilde = Int(tot_chain/2)
    S_pos_real = S_pos_tilde + 1
    
    state = [(n == S_pos_real) ? "ρ" : (n < S_pos_real) ? "Id" : "0" for n = 1:tot_chain]
    RHO_0 = MPO(s_total,state)

    return RHO_0
end

#lvac = vacuum #left_vac(s_total)
    #= RH_tot = ITensor(1.)

    for i=1:S_pos_tilde

        RH = op(state[i],s_total[i]) * op(state[tot_chain+1-i],s_total[tot_chain+1-i])
        x = MPS(RH * delta(s_total[i]',s_total[tot_chain+1-i]'), [s_total[i],s_total[tot_chain+1-i]])
        RH_tot = inner(x 
    end =#

    #= PSI_0 = apply(RHO_0, left_vac(s_total); cutoff, maxdim) =#

#= 
state = [(n == S_pos) ? "ρ_vec" : "0" for n = 1:tot_chain]
    PSI = MPS(s_total,state)  ## (spin+chain) initial MPS

 =#

