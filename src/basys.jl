
function basys(num,N_chain, s_total)
    
    tot_chain = 2*N_chain+1
    S_pos = N_chain+1

    if num == 11
        state = [(n == S_pos) ? "11" : "0" for n = 1:tot_chain]
        PSI = MPS(s_total,state)  ## (spin+chain) initial MPS
    elseif num == 10
        state = [(n == S_pos) ? "10" : "0" for n = 1:tot_chain]
        PSI = MPS(s_total,state)  ## (spin+chain) initial MPS
    elseif num == 01
        state = [(n == S_pos) ? "01" : "0" for n = 1:tot_chain]
        PSI = MPS(s_total,state)  ## (spin+chain) initial MPS
    else
        num == 00
        state = [(n == S_pos) ? "00" : "0" for n = 1:tot_chain]
        PSI = MPS(s_total,state)  ## (spin+chain) initial MPS
    end

    return normalize!(PSI)
end
