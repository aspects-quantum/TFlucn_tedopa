###########################################################################################
######### Outputs characteristic function during the time evolution
###########################################################################################

function charfn_tf(ω_0, Ω, c_01, c_02, ab1, ab2, s_total, β, Nbeta, tau, nt, u, cutoff, maxdim)
    
    count_gates = exp_xHB_tf(ab1, -im*u, s_total);
    evol = unit_gates_tf(ω_0, Ω, c_01, c_02, ab1, ab2, tau, s_total)

    ################################################################################################
    ρ_vec = thermal_MPS(ab1, s_total, β, Nbeta, cutoff, maxdim)  ## initial state of (spin+chain);
    ################################################################################################
    #round(Int64,(n1_bsn_dim-2.6)*(i-1)/(N_chain-1))
    char_fn=Vector{ComplexF64}()

    cg_ρ = apply(count_gates, ρ_vec; cutoff, maxdim)
    U_cg_ρ = cg_ρ
    U_ρ = ρ_vec
    cg_U_ρ = cg_ρ

    chi = inner(cg_U_ρ, U_cg_ρ)
    println(chi[1])
    push!(char_fn,chi[1])

    for t in 1:nt
        U_cg_ρ = apply(evol, U_cg_ρ; cutoff, maxdim)
        U_ρ = apply(evol, U_ρ; cutoff, maxdim)
        cg_U_ρ = apply(count_gates, U_ρ; cutoff, maxdim)

        chi = inner(cg_U_ρ, U_cg_ρ)
        println(chi[1])
        push!(char_fn,chi[1])
    end

    return char_fn
end


