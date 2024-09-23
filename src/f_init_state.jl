
function init_state(s_total)
    
    tot_chain = length(s_total)
    S_pos = Int((tot_chain+1)/2)
     
    state = [(n == S_pos) ? "ρ" : "0" for n = 1:tot_chain]
    RHO_0 = MPO(s_total,state)

    return RHO_0
end