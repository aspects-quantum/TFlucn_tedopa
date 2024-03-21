###########################################################################################
######### Outputs characteristic function during the time evolution
###########################################################################################

## this script uses the function "init_state" from "\src\initial_state.jl"


function spin_boson_char_fn(N_chain, ω_0, Ω, ω_C, c_0, ab1, ab2, s_total, tau, u, cutoff, nt)

    tot_chain = 2*N_chain+1
    S_pos = N_chain+1

    count_gates = exp_xHB(N_chain, ω_C, ab1, -im*u, s_total)
    evol = unit_gates(N_chain, ω_0, Ω, c_0, ω_C, ab1, ab2, tau, s_total)    

    ################################################################################################
    ψ = init_state(N_chain, s_total)  ## initial state of (spin+chain)
    ################################################################################################
    cg_ψ0 = []

    s11 = [(n == S_pos) ? "11" : "0" for n = 1:tot_chain]
    bas11 = MPS(s_total,s11)
    push!(cg_ψ0, apply(count_gates, bas11; cutoff))
    s00 = [(n == S_pos) ? "00" : "0" for n = 1:tot_chain]
    bas00 = MPS(s_total,s00)
    push!(cg_ψ0, apply(count_gates, bas00; cutoff))


    
    char_fn=Vector{ComplexF64}()

    cg_ψ = apply(count_gates, ψ; cutoff)
    U_cg_ψ = cg_ψ

    chi = 0
    for s = 1:2
        chi += inner(cg_ψ0[s],cg_ψ)
    end
    push!(char_fn,chi)


    for t in 1:nt
    
        U_cg_ψ = apply(evol, U_cg_ψ; cutoff)
          
        chi = 0
        for s = 1:2
            chi += inner(cg_ψ0[s],U_cg_ψ)
        end
        push!(char_fn,chi)
  
    end

    return char_fn
end


