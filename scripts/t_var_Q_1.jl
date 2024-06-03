using DrWatson
@quickactivate :flucn_tedopa
## https://juliadynamics.github.io/DrWatson.jl/stable/real_world/#Making-your-project-a-usable-module-1

#@disable_warn_order
ITensors.disable_warn_order()

#################            ########################################################################
################# Parameters ########################################################################
#################            ########################################################################



ρ = [0 0;0 1]
ρ = ρ/tr(ρ)

ITensors.op(::OpName"ρ", ::SiteType"S=1/2") = ρ


## time evolution parameters
cutoff = 1E-8
maxdim = 1000
tau = 10^-2             ## time step duration
nt = 3
ttotal = nt*tau             ## TOTAL TIME evolution

N_chain = 2            ## Number of chain sites for single chain-transformed environment
tot_chain = 2*N_chain+1

boson_dim = 8

# Make an array of 'site' INDICES for the (spin+chain)
S_pos = N_chain+1

s_total = [(n == S_pos) ? Index(2, "S=1/2") : Index(boson_dim, "Qudit") for n = 1:tot_chain]


zer0 = zeros(boson_dim,boson_dim)
zer0[1,1] = 1
ITensors.op(::OpName"0", ::SiteType"Qudit", d::Int) = zer0


ω_C = 5                 ## bath cutoff
ω_0 = 0                 ## spin splitting
Ω = 1                   ## independent model if Ω = 0


###

u = 0.005 # counting field parameter


## 

#= p=plot()
 =#
t_list = collect(0:1:nt)*tau
T_list = [0.1]
α = .01

#cat(a, b) = reshape(append!(vec(a), vec(b)), size(a)[1:end-1]..., :)
support_cutoff = 6*10^2
supp = (0, support_cutoff)                 ## support of the weight function
Nquad = 10^8          ## Number of quadrature points
N_coeff = N_chain + 1



for T = T_list
    β = 1/T 
    n(k) = 1/(exp(β*k) - 1)
    w_fn1(k) = (2*α*k*exp(-k/ω_C))*(1 + n(k))           ## weight function for real space bath
    w_fn2(k) = (2*α*k*exp(-k/ω_C))*n(k)           ## weight function for tilde space bath
    η01 = quadgk(w_fn1, 0, support_cutoff)
    c_01 = sqrt(η01[1])
    η02 = quadgk(w_fn2, 0, support_cutoff)
    c_02 = sqrt(η02[1])
    if N_chain >= 90
        ab1[1:90,1:2] = recur_coeff(w_fn1, supp, 90, Nquad)    ## recurrence coefficients for for real space bath
        ab2[1:90,1:2] = recur_coeff(w_fn2, supp, 90, Nquad)    ## recurrence coefficients for for tilde space bath
        
        a1_100 = ab1[90,1]; b1_100 = ab1[90,2]; a2_100 = ab2[90,1]; b2_100 = ab2[90,2]
        [ab1[i,1] = a1_100  for i = 91:N_coeff]; [ab1[i,2] = b1_100  for i = 91:N_coeff]; 
        [ab2[i,1] = a2_100  for i = 91:N_coeff]; [ab2[i,2] = b2_100  for i = 91:N_coeff]; 
    else
        ab1[1:N_coeff,1:2] = recur_coeff(w_fn1, supp, N_coeff, Nquad)    ## recurrence coefficients for for real space bath
        ab2[1:N_coeff,1:2] = recur_coeff(w_fn2, supp, N_coeff, Nquad)    ## recurrence coefficients for for tilde space bath
    
    end

    chi_2pu = charfn(ω_0, Ω, c_01, c_02, ab1, ab2, s_total, tau, nt, 2*u, cutoff, maxdim)
    chi_pu = charfn(ω_0, Ω, c_01, c_02, ab1, ab2, s_total, tau, nt, u, cutoff, maxdim)
    #chi_0 = ones(length(chi_pu),1)

    var_Q = -(log.(chi_2pu)-2*log.(chi_pu))/(u^2)

    plot!(t_list,real(var_Q),label= "T = $T, α = 0.1")
    @show var_Q

end



plot!(legend=:topright)
xlabel!("t")
title = string("<<Q>>, N_ch = ", N_chain,", b_dim = ", boson_dim,", u = ",u, ", k_max = ", support_cutoff)
title!(title)
display("image/png", p)


a = "varq_bb.png";
savefig(a)