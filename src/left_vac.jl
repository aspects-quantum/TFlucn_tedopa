
function left_vac(s_total)

    N = Int(length(s_total)/2)
vacuum = []

    vac1 = delta(s_total[1],s_total[2*N])
    vac = MPS(vac1, [s_total[1],s_total[2*N]]) 

    push!(vacuum,vac)
    
    for i=2:N
        vaci = delta(s_total[i],s_total[2*N+1-i])
        vac = MPS(vaci, [s_total[i],s_total[2*N+1-i]])
        #= vac = apply(vac, vaci; cutoff,maxdim) =#
        push!(vacuum,vac)
    end

   #=  vacuum = MPS(vac,s_total;cutoff,maxdim) =#

    return vacuum
end


#= thermal_vac = basys(11, N_chain, s_total) + basys(00, N_chain, s_total)

    exp_HB_2 =  exp_xHB(N_chain, ω_C, ab1, β/2, s_total)

    vac = apply(exp_HB_2, thermal_vac; cutoff) * sqrt(partition_fn(N_chain, ω_C, ab1, s_total, β)) =#