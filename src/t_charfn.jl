###########################################################################################
######### Outputs characteristic function during the time evolution
###########################################################################################

function charfn(ω_0, Ω, c_01, c_02, ab1, ab2, s_total, tau, nt, u, cutoff, maxdim)

    count_dates = exp_xHB(ab1, im*u, s_total);
    evol = unit_gates(ω_0, Ω, c_01, c_02, ab1, ab2, tau, s_total)

    ################################################################################################
    ρ = init_state(s_total)  ## initial state of (spin+chain);
    ################################################################################################

    char_fn=Vector{ComplexF64}()

    U_ρ_Ud = ρ
    cd_U_ρ_Ud = apply(count_dates, U_ρ_Ud; cutoff, maxdim)

    chi = tr(cd_U_ρ_Ud)
        
    push!(char_fn,chi[1])


    for t in 1:nt
    
        U_ρ_Ud = apply(evol, U_ρ_Ud; cutoff, maxdim, apply_dag=true)
        cd_U_ρ_Ud = apply(count_dates, U_ρ_Ud; cutoff, maxdim)

        chi = tr(cd_U_ρ_Ud)
        
        push!(char_fn,chi[1])

    end

    return char_fn
end


