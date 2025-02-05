using DrWatson
@quickactivate :flucn_tedopa
## https://juliadynamics.github.io/DrWatson.jl/stable/real_world/#Making-your-project-a-usable-module-1

#@disable_warn_order
ITensors.disable_warn_order()

## this script uses the function "spin_boson_time_evol" from "\src\TEBD_1.jl" ######################

#= 
ITensors.op(::OpName"d_spin", ::SiteType"S=1/2") = [0 0; 0 1]
ITensors.op(::OpName"u_spin", ::SiteType"S=1/2") = [1 0; 0 0]
 =#

#################            ########################################################################
################# Parameters ########################################################################
#################            ########################################################################

spin_state = 1 # spin "Dn"
#= spin_state = 2 # spin "Up"
spin_state = 3 # spin "X+"
spin_state = 4 # spin "X-"
 =#


## time evolution parameters
cutoff = 1E-13
tau = 10^-1             ## time step duration
ttotal = 20            ## TOTAL TIME evolution

N_chain = 12            ## Number of chain sites
boson_dim = 3           ## Dimension of all chain sites

ω_C = 5                 ## bath cutoff
α = 1                   ## coupling strength
ω_0 = 0                 ## spin splitting
Ω = 1                   ## independent model if Ω = 0
c_0 = ω_C * sqrt(2 * α / pi)

TU_list = Vector{Float64}()

T_list = (1:2:9)

for T = T_list

#= T = 10;                  ## temperature of bath
  =#
###
β = 1/T 
N_coeff = N_chain + 1
Nquad = 200 * N_coeff            ## Number of quadrature points
w_fn(t) = (t) * exp(-abs(t))/(1-exp(-β*t*ω_C))   ## weight function leading to associated Laguerre polynomials
supp = (-1000, 1000)        ## support of the weight function

ab = recur_coeff(w_fn, supp, N_coeff, Nquad)    ## recurrence coefficients


u = .001


# Make an array of 'site' INDICES for the (spin+chain)
s_total = [(n == 1) ? Index(2, "S=1/2,n=$n") : Index(boson_dim, "Qudit,n=$n") for n = 1:(N_chain+1)]


    chi_pu = spin_boson_char(spin_state,N_chain, ω_0, Ω, ω_C, c_0, ab, s_total, β, tau, u, cutoff, ttotal)
    chi_mu = spin_boson_char(spin_state,N_chain, ω_0, Ω, ω_C, c_0, ab, s_total, β, tau, -u, cutoff, ttotal)
    chi_0 = 1

    mean_Q = (-im)*(chi_pu-chi_mu)/(2*u)

    var_Q = (-1)*(chi_pu+chi_mu-2*chi_0)/(u^2)

    push!(TU_list,real(var_Q)/real(mean_Q^2))

end

plot(T_list,TU_list)

