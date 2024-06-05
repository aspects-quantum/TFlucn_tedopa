

function expectn(ω_0, Ω, c_01,c_02, ab1, ab2, s_total, tau, cutoff, maxdim, nt)

    tot_chain = length(s_total)
    S_pos = Int((tot_chain+1)/2)

    evol = unit_gates(ω_0,Ω,c_01,c_02,ab1, ab2, tau, s_total)

    ################################################################################################
    ρ = init_state(s_total)  ## initial state of (spin+chain)
    ################################################################################################

    obs_list = Vector{Any}()

    i=0;
    for t in 0:1:nt
        obs = tr(apply(op("OBSERVABLE", s_total[S_pos]), ρ; cutoff, maxdim))
        obs_list = push!(obs_list, obs)

        ρ = apply(evol, ρ; apply_dag = true, cutoff, maxdim)
        i=i+1
    end

    return obs_list


end


