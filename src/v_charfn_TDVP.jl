###########################################################################################
######### Outputs characteristic function during the time evolution
###########################################################################################

function charfn_TDVP(ω_0, Ω, c_01, c_02, ab1, ab2, s_total, β, Nbeta, tau, nt, u, cutoff, maxdim)
    
    bath_ham = HB_tf(ab1, s_total)
    tot_ham = HTOT_tf(ω_0,Ω,c_01,c_02,ab1, ab2, s_total)

    ################################################################################################
    ρ_vec = thermal_MPS(ab1, s_total, β, Nbeta, cutoff, maxdim)  ## initial state of (spin+chain);
    ################################################################################################

    char_fn=Vector{ComplexF64}()

    cg_ρ = tdvp(bath_ham, -u*1.0im, ρ_vec; nsteps=5, updater_backend="exponentiate", cutoff=1e-9, maxdim, normalize=true)
    U_cg_ρ = cg_ρ
    U_ρ = ρ_vec
    cg_U_ρ = cg_ρ

    chi = inner(cg_U_ρ, U_cg_ρ)
    println(chi[1])
    push!(char_fn,chi[1])

    for t in 1:nt
        U_cg_ρ = tdvp(tot_ham, -1.0im*tau, U_cg_ρ; nsteps=10, updater_backend="exponentiate", cutoff=1e-9, maxdim, normalize=true)
        U_ρ = tdvp(tot_ham, -1.0im*tau, U_ρ; nsteps=10, updater_backend="exponentiate", cutoff=1e-9, maxdim, normalize=true)
        cg_U_ρ = tdvp(bath_ham, -1.0im*u, U_ρ; nsteps=5, updater_backend="exponentiate", cutoff=1e-9, maxdim, normalize=true)

        chi = inner(cg_U_ρ, U_cg_ρ)
        println(chi[1])
        push!(char_fn,chi[1])
    end

    return char_fn
end


