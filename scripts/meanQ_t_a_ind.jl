using DrWatson
@quickactivate :flucn_tedopa
## https://juliadynamics.github.io/DrWatson.jl/stable/real_world/#Making-your-project-a-usable-module-1

#@disable_warn_order
ITensors.disable_warn_order()

#################            ########################################################################
################# Parameters ########################################################################
#################            ########################################################################


#ITensors.space(::SiteType"S=3/2") = 4
#= 
ITensors.state(::StateName"11", ::SiteType"S=3/2") = [1, 0, 0, 0]
ITensors.state(::StateName"01", ::SiteType"S=3/2") = [0, 1, 0, 0]
ITensors.state(::StateName"10", ::SiteType"S=3/2") = [0, 0, 1, 0]
ITensors.state(::StateName"00", ::SiteType"S=3/2") = [0, 0, 0, 1]
 =#

#= Id2 = [1 0;0 1]
X = [0 1;1 0]
Z = [1 0;0 -1] =#

ρ = [0 0;0 1]
ρ = ρ/tr(ρ)

ITensors.op(::OpName"ρ", ::SiteType"S=1/2") = ρ

#= 
RHO = kron(ρ, Id2)

ρ_vec = vec(transpose(ρ))  ##vectorized rho after thermofield transformation
 =#
#ITensors.state(::StateName"ρ_vec", ::SiteType"S=3/2") = collect(ρ_vec)
#= 
ITensors.op(::OpName"Id", ::SiteType"S=3/2") = kron(Id2, Id2)
ITensors.op(::OpName"Sz1", ::SiteType"S=3/2") = .5*kron(Z, Id2)
ITensors.op(::OpName"Sz2", ::SiteType"S=3/2") = .5*kron(Id2, Z)
ITensors.op(::OpName"Sx1", ::SiteType"S=3/2") = .5*kron(X, Id2)
ITensors.op(::OpName"Sx2", ::SiteType"S=3/2") = .5*kron(Id2, X)
ITensors.op(::OpName"RHO", ::SiteType"S=3/2") = RHO
 =#
## time evolution parameters
cutoff = 1E-8
maxdim = 100
tau = 10^-2             ## time step duration
nt = 10
ttotal = nt*tau             ## TOTAL TIME evolution

N_chain = 10            ## Number of chain sites for single chain-transformed environment
tot_chain = 2*(N_chain+1)
#= n1_bsn_dim = 6;
boson_dim = [n1_bsn_dim-round(Int64,(n1_bsn_dim-1.6)*(i-1)/(N_chain-1)) for i = 1:N_chain]           ## Dimension of chain sites
 =#
boson_dim = 3

# Make an array of 'site' INDICES for the (spin+chain)
S_pos_tilde = N_chain+1
S_pos_real = N_chain+2
#s_total = [(n == S_pos) ? Index(4, "S=3/2") : Index(boson_dim, "Qudit") for n = 1:tot_chain]
s_total = [(n == S_pos_tilde) | (n == S_pos_real) ? Index(2, "S=1/2") : Index(boson_dim, "Qudit") for n = 1:tot_chain]


zer0 = zeros(boson_dim,boson_dim)
zer0[1,1] = 1
ITensors.op(::OpName"0", ::SiteType"Qudit", d::Int) = zer0


ω_C = 5                 ## bath cutoff
ω_0 = 1                 ## spin splitting
Ω = 0                   ## independent model if Ω = 0

T = 5;                  ## temperature of bath
β = 1/T 
###

u = .001 # counting field parameter


## 

p=plot()

t_list = collect(0:1:nt)*tau
#α_list = [0.1,1.5]
α = 1.5

#cat(a, b) = reshape(append!(vec(a), vec(b)), size(a)[1:end-1]..., :)

#for α = α_list

J(t) = 2*α*t*exp(-t/ω_C)
η0 = quadgk(J, 0, ω_C, rtol=1e-4)


    #c_0 = ω_C * sqrt(2 * α / pi)
    c_0 = sqrt(η0[1]/pi)
    N_coeff = N_chain + 1
    Nquad = 5000 * N_coeff            ## Number of quadrature points
    #(ω_C^2) * (2 * α) * .
    w_fn1(t) = (t) * exp(-t)/(1-exp(-β*t*ω_C))           ## weight function for real space bath
    #w_fn2(t) = (t) * exp(-t)/(exp(β*t*ω_C)-1)           ## weight function for tilde space bath
    supp = (0, 10000)                 ## support of the weight function
    ab1 = recur_coeff(w_fn1, supp, N_coeff, Nquad)    ## recurrence coefficients for for real space bath
    #ab2 = recur_coeff(w_fn2, supp, N_coeff, Nquad)    ## recurrence coefficients for for tilde space bath

    #chi_pu = spin_boson_char_fn(N_chain, ω_0, Ω, ω_C, c_0, ab1, ab2, s_total, tau, u, cutoff, nt)
    chi_pu = SB_char_2(ω_0, Ω, ω_C, c_0, ab1, s_total, tau, nt, u, cutoff, maxdim)
    
    mean_Q = imag.(chi_pu)/u

    plot!(t_list,real(mean_Q),label= "α = $α")

#end


plot!(legend=:topright)
xlabel!("t")
title = string("N_chain = ", N_chain,", boson_dim = ", boson_dim,", u = ",u)
title!(title)
display("image/png", p)

#= 
a = "plots\\meanq_t_a_2.png";
savefig(a)
     =#