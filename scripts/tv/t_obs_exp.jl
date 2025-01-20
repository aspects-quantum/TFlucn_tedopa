using DrWatson
@quickactivate :flucn_tedopa
## https://juliadynamics.github.io/DrWatson.jl/stable/real_world/#Making-your-project-a-usable-module-1

#@disable_warn_order
ITensors.disable_warn_order()



ρ = (1/2)*[1 1;1 1]
ρ = ρ/tr(ρ)

ITensors.op(::OpName"ρ", ::SiteType"S=1/2") = ρ

observable = [0 1;1 0]
ITensors.op(::OpName"OBSERVABLE", ::SiteType"S=1/2") = observable


## time evolution parameters
cutoff = 1E-8
maxdim = 15
tau = 10^-2             ## time step duration
nt = 500
ttotal = nt*tau             ## TOTAL TIME evolution


N_chain = 15              ## Number of chain sites
tot_chain = 2*N_chain+1
S_pos = N_chain+1

boson_dim = 4           ## Dimension of all chain sites

ω_C = 500                 ## bath cutoff
α = .05                  ## coupling strength
ω_0 = 0                 ## spin splitting
Ω = 0                   ## independent model if Ω = 0

T = .2;                  ## temperature of bath
β = 1/T 
###

N_coeff = N_chain + 1
Nquad = 5000 * N_coeff            ## Number of quadrature points
#(ω_C^2) * (2 * α) * .
#= w_fn1(t) = (t) * exp(-t)/(1-exp(-β*t*ω_C))           ## weight function for real space bath
w_fn2(t) = (t) * exp(-t)/(exp(β*t*ω_C)-1)           ## weight function for tilde space bath =#
n(k) = 1/(exp(β*k) - 1)
w_fn1(k) = (2*α*k*exp(-k/ω_C))*(1 + n(k))           ## weight function for real space bath
w_fn2(k) = (2*α*k*exp(-k/ω_C))*n(k)           ## weight function for tilde space bath
η01 = quadgk(w_fn1, 0, 10000, rtol=1e-4)
c_01 = sqrt(η01[1])
η02 = quadgk(w_fn2, 0, 10000, rtol=1e-4)
c_02 = sqrt(η02[1])
supp = (0, 10000)                 ## support of the weight function
ab1 = recur_coeff(w_fn1, supp, N_coeff, Nquad)    ## recurrence coefficients for for real space bath
ab2 = recur_coeff(w_fn2, supp, N_coeff, Nquad)    ## recurrence coefficients for for tilde space bath

############# Make an array of 'site' INDICES for the (spin+chain) #####################################
s_total = [(n == S_pos) ? Index(2, "S=1/2") : Index(boson_dim, "Qudit") for n = 1:tot_chain]
########################################################################################################

zer0 = zeros(boson_dim,boson_dim)
zer0[1,1] = 1
ITensors.op(::OpName"0", ::SiteType"Qudit", d::Int) = zer0

p=plot();

t_list = collect(0:1:nt)*tau
Out_list = expectn(ω_0, Ω, c_01,c_02, ab1, ab2, s_total, tau, cutoff, maxdim, nt)
plot!(t_list,real(Out_list))


xlabel!("t")
title = string("N_chain = ", N_chain,", boson_dim = ", boson_dim)
title!(title)
display("image/png", p)