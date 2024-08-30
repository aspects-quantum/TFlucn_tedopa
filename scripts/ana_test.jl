
using DrWatson
@quickactivate :dissiPT
## https://juliadynamics.github.io/DrWatson.jl/stable/real_world/#Making-your-project-a-usable-module-1

#@disable_warn_order
ITensors.disable_warn_order()



N_chain = 40            ## Number of chain sites for single chain-transformed environment
tot_chain = 2*N_chain


s_total = siteinds("S=1/2", tot_chain)
s_tot_comp = s_total[1:N_chain]




I0 = MPO(s_tot_comp, "Id")
ρ_vec = convert(MPS, I0)
for i = 1:Int(length(s_tot_comp))
    ρ_vec[i] = ρ_vec[i]*delta(s_tot_comp[i]', s_total[N_chain+i])
end
normalize!(ρ_vec)





