###########################################################################################
######### Outputs characteristic function during the time evolution
###########################################################################################

function SB_char_1(ω_0, Ω, ω_C, c_0, ab1, s_total, tau, nt, u, cutoff, maxdim)

    count_gates = exp_xHB(ω_C, ab1, -im*u, s_total)
    evol = unit_real(ω_0,Ω,c_0,ω_C, ab1, tau, s_total)   

    ################################################################################################
    ψ = init_state(s_total, cutoff, maxdim)  ## initial state of (spin+chain)
    ι = left_vac(s_total) ## left vacuum state of (spin+chain)
    ################################################################################################

    char_fn=Vector{ComplexF64}()

    cg_ψ = apply(count_gates, ψ; cutoff, maxdim)
    U_cg_ψ = cg_ψ
    U_I = ι
    cg_U_I = apply(count_gates, U_I; cutoff, maxdim)

    chi = inner(cg_U_I,U_cg_ψ)
    push!(char_fn,chi)

    for t in 1:nt
    
        U_cg_ψ = apply(evol, U_cg_ψ; cutoff, maxdim)
        U_I = apply(evol, U_I; cutoff, maxdim)
        cg_U_I = apply(count_gates, U_I; cutoff, maxdim)
          
        chi = inner(cg_U_I,U_cg_ψ)
       
        push!(char_fn,chi)
  
    end

    return char_fn
end


