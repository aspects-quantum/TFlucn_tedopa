###########################################################################################
######### Outputs characteristic function during the time evolution
###########################################################################################

function SB_char_2(ω_0, Ω, ω_C, c_0, ab1, s_total, tau, nt, u, cutoff, maxdim)

    tot_chain = length(s_total)
    S_pos_tilde = Int(tot_chain/2)

    count_gates = exp_xHB(ω_C, ab1, -im*u, s_total);
    #count_gates = apply(count_g,MPO(s_total,"Id");cutoff,maxdim)
    evol = unit_real(ω_0,Ω,c_0,ω_C, ab1, tau, s_total)   ;
    #evol = apply(evo,MPO(s_total,"Id");cutoff,maxdim)

    ################################################################################################
    ψ_ = init_state(s_total)  ## initial state of (spin+chain);
    #ι = left_vac(s_total) ## left vacuum state of (spin+chain)
    ################################################################################################

    char_fn=Vector{ComplexF64}()
 
    ψ = MPS(S_pos_tilde);
    for i = 1:S_pos_tilde
        ψ[i] = ψ_[i] 
    end
    ψ = swapprime(ψ, 1, 0)
#= 
    for i = 1:S_pos_tilde
        it = ITensor(1.)
        ψ_[i] *= delta(s_total[i],s_total[tot_chain+1-i])
        it = ψ_[i]*ψ_[tot_chain+1-i]
        ψ[i] = it
        indit = inds(it)
        lind = Int(length(indit)/2)
        ψ[i], ψ[tot_chain+1-i] = qr(ψ[i], inds(indit[1:lind]))
    end
    ψ = replaceprime(ψ, 1=>0)
 =#
    cg_ψ = apply(count_gates, ψ; cutoff, maxdim) 
    U_cg_ψ = cg_ψ
    #U_cg_ψ = MPS(S_pos_tilde)
    #= for i = 1:S_pos_tilde
        U_cg_ψ[i] *= delta(s_total[i],s_total[tot_chain+1-i])
        #U_cg_ψ[i] = U_cg_ψ_[i] * U_cg_ψ_[tot_chain+1-i]
    end =#
    #it = ITensor(1.)
    #= U_cg_ψ = MPS(tot_chain)
    for i = 1:S_pos_tilde
        U_cg_ψ_[i] *= delta(s_total[i],s_total[tot_chain+1-i])
        U_cg_ψ[i] = U_cg_ψ_[i]*U_cg_ψ_[tot_chain+1-i]
        U_cg_ψ[i], U_cg_ψ[tot_chain+1-i] = qr(U_cg_ψ[i], (s_total[i]'))
    end
    U_cg_ψ = replaceprime(MPS(it,s_total'), 1=>0) =#


    U_I_ = MPO(s_total, "Id")
    U_I = MPS(S_pos_tilde)
    for i = 1:S_pos_tilde
        U_I[i] =  U_I_[i] 
    end
    U_I = swapprime(U_I, 1, 0)
#= 
    U_I = MPS(tot_chain)
    for i = 1:S_pos_tilde
        it = ITensor(1.)
        U_I_[i] *= delta(s_total[i],s_total[tot_chain+1-i])
        it = U_I_[i]*U_I_[tot_chain+1-i]
        U_I[i] = it
        indit = inds(it)
        lind = Int(length(indit)/2)
        U_I[i], U_I[tot_chain+1-i] = qr(U_I[i], inds(indit[1:lind]))
    end
    U_I = replaceprime(U_I, 1=>0) =#
    #U_I = MPS(S_pos_tilde)
    #= for i = 1:S_pos_tilde
        U_I[i] *= delta(s_total[i],s_total[tot_chain+1-i])
        #U_I[i] = U_I_[i] * U_I_[tot_chain+1-i]
    end =#
    #= it = ITensor(1.)
    for i = 1:tot_chain
        if i<=S_pos_tilde
        U_I_[i] *= delta(s_total[i],s_total[tot_chain+1-i])
        end
        it *= U_I_[i] 
    end
    U_I = replaceprime(MPS(it,s_total'), 1=>0) =#
    cg_U_I = apply(count_gates, U_I; cutoff, maxdim)
#= 
    sU_cg_ψ = MPS(S_pos_tilde)
    scg_U_I = MPS(S_pos_tilde)

    for i = 1:S_pos_tilde
        sU_cg_ψ[i] = U_cg_ψ[i] * U_cg_ψ[tot_chain+1-i]
        scg_U_I[i] = cg_U_I[i] * cg_U_I[tot_chain+1-i]
    end  =#


    chi = inner(cg_U_I,U_cg_ψ)
    push!(char_fn,chi)

    for t in 1:nt
    
        U_cg_ψ = apply(evol, U_cg_ψ; cutoff, maxdim)
        U_I = apply(evol, U_I; cutoff, maxdim)
        cg_U_I = apply(count_gates, U_I; cutoff, maxdim)
          #= 
        sU_cg_ψ = MPS(S_pos_tilde)
        scg_U_I = MPS(S_pos_tilde)
    
        for i = 1:S_pos_tilde
            sU_cg_ψ[i] = U_cg_ψ[i] * U_cg_ψ[tot_chain+1-i]
            scg_U_I[i] = cg_U_I[i] * cg_U_I[tot_chain+1-i]
        end 
     =#
    
        chi = inner(cg_U_I,U_cg_ψ)
       
        push!(char_fn,chi)
  
    end

    return char_fn
end


